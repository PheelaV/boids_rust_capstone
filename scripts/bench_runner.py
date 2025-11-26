# /// script
# requires-python = ">=3.14"
# dependencies = []
# ///

"""
Benchmark runner that records results per git commit.

Usage:
    uv run scripts/bench_runner.py [--output-dir benchmark_data]

Results are saved to:
  - benchmark_data/results.json  (structured data for all runs)
  - benchmark_data/logs/<commit>_<timestamp>.log  (raw output)
  - benchmark_data/reports/<commit>_<timestamp>.txt  (human-readable report)
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


def run_benchmark():
    """Run cargo bench and stream output while capturing it.

    Returns (output, exit_code, error_message).
    """
    print("=" * 60)
    print("Running: cargo bench --bench flock_scalability")
    print("=" * 60)
    print()

    captured_output = []
    error_message = None

    process = subprocess.Popen(
        ["cargo", "bench", "--bench", "flock_scalability"],
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

    # Pattern: flock_scalability/boids/2^14
    #                          time:   [3.2824 s 3.2887 s 3.2951 s]
    #                          thrpt:  [1.2729 Melem/s 1.2754 Melem/s 1.2778 Melem/s]

    time_pattern = r"flock_scalability/boids/(2\^\d+)\s+time:\s+\[([0-9.]+)\s+(\w+)\s+([0-9.]+)\s+(\w+)\s+([0-9.]+)\s+(\w+)\]"
    thrpt_pattern = r"thrpt:\s+\[([0-9.]+)\s+(\w+)\s+([0-9.]+)\s+(\w+)\s+([0-9.]+)\s+(\w+)\]"

    lines = output.split("\n")
    i = 0
    while i < len(lines):
        time_match = re.search(time_pattern, lines[i])
        if time_match:
            name = time_match.group(1)  # e.g., "2^14"
            time_low = float(time_match.group(2))
            time_unit = time_match.group(3)
            time_mid = float(time_match.group(4))
            time_high = float(time_match.group(6))

            # Look for throughput on next few lines
            thrpt_data = None
            for j in range(i + 1, min(i + 5, len(lines))):
                thrpt_match = re.search(thrpt_pattern, lines[j])
                if thrpt_match:
                    thrpt_data = {
                        "low": float(thrpt_match.group(1)),
                        "mid": float(thrpt_match.group(3)),
                        "high": float(thrpt_match.group(5)),
                        "unit": thrpt_match.group(4),
                    }
                    break

            results[name] = {
                "time": {
                    "low": time_low,
                    "mid": time_mid,
                    "high": time_high,
                    "unit": time_unit,
                },
                "throughput": thrpt_data,
            }
        i += 1

    return results


def generate_report(git_info, bench_results, timestamp, error_message=None):
    """Generate a human-readable report."""
    lines = []
    lines.append("=" * 60)
    lines.append("BOIDS SCALABILITY BENCHMARK REPORT")
    lines.append("=" * 60)
    lines.append("")
    lines.append(f"Timestamp: {timestamp}")
    lines.append(f"Commit:    {git_info['commit_short']} ({git_info['branch']})")
    lines.append(f"Full SHA:  {git_info['commit']}")
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

    lines.append(f"{'Benchmark':<15} {'Time (s)':<25} {'Throughput':<20}")
    lines.append(f"{'-'*15} {'-'*25} {'-'*20}")

    for name in sorted(bench_results.keys(), key=lambda x: int(x.split("^")[1])):
        data = bench_results[name]
        time_str = f"{data['time']['mid']:.3f} [{data['time']['low']:.3f}-{data['time']['high']:.3f}] {data['time']['unit']}"
        if data["throughput"]:
            thrpt_str = f"{data['throughput']['mid']:.2f} {data['throughput']['unit']}"
        else:
            thrpt_str = "N/A"
        lines.append(f"{name:<15} {time_str:<25} {thrpt_str:<20}")

    lines.append("")

    # Scaling analysis
    if len(bench_results) >= 2:
        lines.append("-" * 60)
        lines.append("SCALING ANALYSIS")
        lines.append("-" * 60)
        lines.append("")

        sorted_names = sorted(
            bench_results.keys(), key=lambda x: int(x.split("^")[1])
        )
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

            lines.append(f"{prev_name} -> {curr_name}:")
            lines.append(f"  Boid count:  {boid_ratio:.1f}x")
            lines.append(f"  Time ratio:  {time_ratio:.2f}x")
            lines.append(f"  Efficiency:  {boid_ratio/time_ratio:.2f}x (1.0 = linear scaling)")
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

    # Parse args
    if "--output-dir" in sys.argv:
        idx = sys.argv.index("--output-dir")
        if idx + 1 < len(sys.argv):
            output_dir = Path(sys.argv[idx + 1])

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
    output, exit_code, error_message = run_benchmark()

    # Save raw log
    log_file = logs_dir / f"{git_info['commit_short']}_{file_timestamp}.log"
    with open(log_file, "w") as f:
        f.write(output)
    print(f"\nRaw log saved to: {log_file}")

    # Parse results (may be partial if benchmark crashed)
    bench_results = parse_criterion_output(output)

    # Generate report (even if partial/failed)
    report = generate_report(git_info, bench_results, timestamp, error_message)

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
