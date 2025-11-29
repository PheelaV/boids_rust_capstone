# /// script
# requires-python = ">=3.14"
# dependencies = []
# ///

"""
Benchmark runner that records results per git commit.

Usage:
    uv run scripts/bench_runner.py [--output-dir benchmark_data] [--quick] [--features FEATURES]

Options:
    --quick              Run only 2^14 benchmark for quick iteration/testing
    --features FEATURES  Cargo features to enable (e.g., "simd-neon" or "simd-avx2")
    --nightly           Use cargo +nightly (required for SIMD features)

Results are saved to:
  - benchmark_data/results.json  (structured data for all runs)
  - benchmark_data/logs/<commit>_<timestamp>.log  (raw output)
  - benchmark_data/reports/<commit>_<timestamp>.txt  (human-readable report)

Examples:
    # Run default benchmarks
    uv run scripts/bench_runner.py

    # Run with SIMD optimizations (M3 Mac)
    uv run scripts/bench_runner.py --nightly --features simd-neon

    # Run with SIMD optimizations (AMD Zen 3)
    uv run scripts/bench_runner.py --nightly --features simd-avx2

    # Quick test with SIMD
    uv run scripts/bench_runner.py --quick --nightly --features simd-neon
"""

import subprocess
import json
import re
import sys
import os
from datetime import datetime
from pathlib import Path


def get_git_info():
    """Get current git commit hash and branch."""
    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], text=True
    ).strip()
    commit_short = subprocess.check_output(
        ["git", "rev-parse", "--short", "HEAD"], text=True
    ).strip()
    branch = subprocess.check_output(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"], text=True
    ).strip()
    # Check for uncommitted changes
    status = subprocess.check_output(
        ["git", "status", "--porcelain"], text=True
    ).strip()
    dirty = len(status) > 0

    return {
        "commit": commit,
        "commit_short": commit_short,
        "branch": branch,
        "dirty": dirty,
    }


def run_benchmark(quick=False, features=None, nightly=False):
    """Run cargo bench and stream output while capturing it.

    Returns (output, exit_code, error_message).
    """
    if nightly:
        cmd = ["cargo", "+nightly", "bench", "--bench", "flock_scalability"]
    else:
        cmd = ["cargo", "bench", "--bench", "flock_scalability"]

    if features:
        cmd.extend(["--features", features])

    if quick:
        cmd.extend(["--", r"2\^14"])  # Only run 2^14 benchmark (escape ^ for regex)

    print("=" * 60)
    print(f"Running: {' '.join(cmd)}")
    print("=" * 60)
    print()

    captured_output = []
    error_message = None

    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,  # Line buffered
    )

    # Stream output to console while capturing
    for line in process.stdout:
        print(line, end="")  # Print to console
        captured_output.append(line)

    process.wait()
    output = "".join(captured_output)

    # Check for panic/crash
    if process.returncode != 0:
        # Try to extract panic message
        if "panicked at" in output:
            for line in output.split("\n"):
                if "panicked at" in line or "assertion failed" in line:
                    error_message = line.strip()
                    break
        if not error_message:
            error_message = f"Benchmark failed with exit code {process.returncode}"

    return output, process.returncode, error_message


def parse_criterion_output(output):
    """Parse criterion benchmark output into structured data."""
    results = {}

    # Pattern: benchmark name on its own line, then time on next line
    # flock_scalability/boids/2^14
    #                         time:   [3.2824 s 3.2887 s 3.2951 s]
    #                         thrpt:  [1.2729 Melem/s 1.2754 Melem/s 1.2778 Melem/s]

    name_pattern = r"flock_scalability/boids/(2\^\d+)\s*$"
    time_pattern = r"time:\s+\[([0-9.]+)\s+(\w+)\s+([0-9.]+)\s+(\w+)\s+([0-9.]+)\s+(\w+)\]"
    thrpt_pattern = r"thrpt:\s+\[([0-9.]+)\s+([\w/]+)\s+([0-9.]+)\s+([\w/]+)\s+([0-9.]+)\s+([\w/]+)\]"

    lines = output.split("\n")
    i = 0
    while i < len(lines):
        name_match = re.search(name_pattern, lines[i])
        if name_match:
            name = name_match.group(1)  # e.g., "2^14"

            # Look for time on next few lines
            time_data = None
            thrpt_data = None
            for j in range(i + 1, min(i + 10, len(lines))):
                if time_data is None:
                    time_match = re.search(time_pattern, lines[j])
                    if time_match:
                        time_data = {
                            "low": float(time_match.group(1)),
                            "mid": float(time_match.group(3)),
                            "high": float(time_match.group(5)),
                            "unit": time_match.group(2),
                        }
                        continue

                thrpt_match = re.search(thrpt_pattern, lines[j])
                if thrpt_match:
                    thrpt_data = {
                        "low": float(thrpt_match.group(1)),
                        "mid": float(thrpt_match.group(3)),
                        "high": float(thrpt_match.group(5)),
                        "unit": thrpt_match.group(4),
                    }
                    break

            if time_data:
                results[name] = {
                    "time": time_data,
                    "throughput": thrpt_data,
                }
        i += 1

    return results


def infer_iterations(bench_results):
    """Infer iteration count for each benchmark from throughput and time.

    iterations = (throughput * time) / boid_count
    """
    result = {}
    for name, data in bench_results.items():
        if data["throughput"] is None:
            result[name] = None
            continue

        exp = int(name.split("^")[1])
        boid_count = 1 << exp

        # Convert throughput to elements/s
        thrpt = data["throughput"]["mid"]
        unit = data["throughput"]["unit"]
        if "Melem" in unit:
            thrpt *= 1_000_000
        elif "Kelem" in unit:
            thrpt *= 1_000

        time_s = data["time"]["mid"]
        elements = thrpt * time_s
        iterations = round(elements / boid_count)
        result[name] = iterations

    return result


def generate_report(git_info, bench_results, timestamp, error_message=None, features=None, nightly=False):
    """Generate a human-readable report."""
    lines = []
    lines.append("=" * 60)
    lines.append("BOIDS SCALABILITY BENCHMARK REPORT")
    lines.append("=" * 60)
    lines.append("")
    lines.append(f"Timestamp: {timestamp}")
    lines.append(f"Commit:    {git_info['commit_short']} ({git_info['branch']})")
    lines.append(f"Full SHA:  {git_info['commit']}")
    if features:
        lines.append(f"Features:  {features}")
    if nightly:
        lines.append(f"Toolchain: nightly")
    if git_info["dirty"]:
        lines.append("WARNING:   Working directory has uncommitted changes!")
    lines.append(f"Status:    {'ERROR' if error_message else 'SUCCESS'}")
    if error_message:
        lines.append(f"Error:     {error_message}")
    lines.append("")
    lines.append("-" * 60)
    lines.append("RESULTS" + (" (partial)" if error_message else ""))
    lines.append("-" * 60)
    lines.append("")

    if not bench_results:
        lines.append("No benchmark results collected.")
        lines.append("")
        lines.append("=" * 60)
        return "\n".join(lines)

    # Infer iterations for each benchmark
    iterations_map = infer_iterations(bench_results)

    lines.append(f"{'Benchmark':<12} {'Iters':<8} {'Time (s)':<25} {'Throughput':<15}")
    lines.append(f"{'-'*12} {'-'*8} {'-'*25} {'-'*15}")

    for name in sorted(bench_results.keys(), key=lambda x: int(x.split("^")[1])):
        data = bench_results[name]
        time_str = f"{data['time']['mid']:.3f} [{data['time']['low']:.3f}-{data['time']['high']:.3f}] {data['time']['unit']}"
        iters = iterations_map.get(name)
        iters_str = str(iters) if iters else "?"
        if data["throughput"]:
            thrpt_str = f"{data['throughput']['mid']:.2f} {data['throughput']['unit']}"
        else:
            thrpt_str = "N/A"
        lines.append(f"{name:<12} {iters_str:<8} {time_str:<25} {thrpt_str:<15}")

    lines.append("")

    # Scaling analysis - only compare benchmarks with same iteration count
    if len(bench_results) >= 2:
        # Group by iteration count
        by_iters = {}
        for name, iters in iterations_map.items():
            if iters is not None:
                by_iters.setdefault(iters, []).append(name)

        # Only do scaling analysis for groups with 2+ benchmarks
        comparable_groups = {k: v for k, v in by_iters.items() if len(v) >= 2}

        if comparable_groups:
            lines.append("-" * 60)
            lines.append("SCALING ANALYSIS")
            lines.append("-" * 60)
            lines.append("")

            for iters, names in sorted(comparable_groups.items(), reverse=True):
                sorted_names = sorted(names, key=lambda x: int(x.split("^")[1]))
                lines.append(f"[{iters} iterations]")

                for i in range(1, len(sorted_names)):
                    prev_name = sorted_names[i - 1]
                    curr_name = sorted_names[i]
                    prev_exp = int(prev_name.split("^")[1])
                    curr_exp = int(curr_name.split("^")[1])

                    prev_time = bench_results[prev_name]["time"]["mid"]
                    curr_time = bench_results[curr_name]["time"]["mid"]

                    # Boids doubled = 2x, expected O(n) would be 2x, O(n^2) would be 4x
                    boid_ratio = 2 ** (curr_exp - prev_exp)
                    time_ratio = curr_time / prev_time

                    lines.append(f"  {prev_name} -> {curr_name}:")
                    lines.append(f"    Boid count:  {boid_ratio:.1f}x")
                    lines.append(f"    Time ratio:  {time_ratio:.2f}x")
                    lines.append(f"    Efficiency:  {boid_ratio/time_ratio:.2f}x (1.0 = linear scaling)")
                lines.append("")

        # List benchmarks that couldn't be compared
        uncomparable = [n for n, i in iterations_map.items()
                       if i is None or i not in comparable_groups or len(by_iters.get(i, [])) < 2]
        if uncomparable:
            lines.append("Note: No scaling comparison for: " + ", ".join(sorted(uncomparable)))
            lines.append("")

    lines.append("=" * 60)
    return "\n".join(lines)


def load_results(filepath):
    """Load existing results or return empty list."""
    if filepath.exists():
        with open(filepath) as f:
            return json.load(f)
    return []


def save_results(filepath, results):
    """Save results to JSON file."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w") as f:
        json.dump(results, f, indent=2)


def main():
    output_dir = Path("benchmark_data")
    quick_mode = "--quick" in sys.argv
    nightly_mode = "--nightly" in sys.argv
    features = None

    # Parse args
    if "--output-dir" in sys.argv:
        idx = sys.argv.index("--output-dir")
        if idx + 1 < len(sys.argv):
            output_dir = Path(sys.argv[idx + 1])

    if "--features" in sys.argv:
        idx = sys.argv.index("--features")
        if idx + 1 < len(sys.argv):
            features = sys.argv[idx + 1]

    # Create directories
    logs_dir = output_dir / "logs"
    reports_dir = output_dir / "reports"
    logs_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    # Get git info
    git_info = get_git_info()
    timestamp = datetime.now().isoformat()
    file_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    print()
    print(f"Commit: {git_info['commit_short']} ({git_info['branch']})")
    if git_info["dirty"]:
        print("WARNING: Working directory has uncommitted changes")
    print()

    # Run benchmark (streams to console)
    if quick_mode:
        print("QUICK MODE: Running only 2^14 benchmark\n")
    if features:
        print(f"FEATURES: {features}\n")
    if nightly_mode:
        print("USING: cargo +nightly\n")
    output, exit_code, error_message = run_benchmark(
        quick=quick_mode, features=features, nightly=nightly_mode
    )

    # Save raw log
    log_file = logs_dir / f"{git_info['commit_short']}_{file_timestamp}.log"
    with open(log_file, "w") as f:
        f.write(output)
    print(f"\nRaw log saved to: {log_file}")

    # Parse results (may be partial if benchmark crashed)
    bench_results = parse_criterion_output(output)

    # Generate report (even if partial/failed)
    report = generate_report(
        git_info, bench_results, timestamp, error_message,
        features=features, nightly=nightly_mode
    )

    # Save report
    report_file = reports_dir / f"{git_info['commit_short']}_{file_timestamp}.txt"
    with open(report_file, "w") as f:
        f.write(report)

    # Print report
    print()
    print(report)
    print()
    print(f"Report saved to: {report_file}")

    # Create result entry for JSON
    entry = {
        "timestamp": timestamp,
        "git": git_info,
        "benchmarks": bench_results,
        "status": "error" if error_message else "success",
        "features": features,
        "nightly": nightly_mode,
    }
    if error_message:
        entry["error"] = error_message

    # Load existing results and append
    results_file = output_dir / "results.json"
    all_results = load_results(results_file)
    all_results.append(entry)
    save_results(results_file, all_results)

    print(f"Results appended to: {results_file}")
    print(f"Total recorded runs: {len(all_results)}")

    # Exit with appropriate code
    if error_message:
        sys.exit(1)


if __name__ == "__main__":
    main()
