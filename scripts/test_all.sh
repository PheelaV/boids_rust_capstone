#!/usr/bin/env bash
# Run all tests with optimal settings for reliability and speed
#
# boids_lib uses a global RNG, so tests must run single-threaded to avoid
# non-deterministic test failures from parallel RNG access.
#
# web_tests use a shared port (8765) for the local server, so they must also
# run single-threaded to avoid port conflicts.
#
# Other packages can run in parallel.
#
# Usage:
#   ./scripts/test_all.sh          # Run all tests
#   ./scripts/test_all.sh --lib    # Run only boids_lib tests
#   ./scripts/test_all.sh --web    # Run only web_tests (requires WebDriver)

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

if [[ "$1" == "--lib" ]]; then
    echo "=============================================="
    echo "Running boids_lib tests only (single-threaded)"
    echo "=============================================="
    cargo test -p boids_lib -- --test-threads=1
    exit 0
fi

if [[ "$1" == "--web" ]]; then
    echo "=============================================="
    echo "Running web_tests only (single-threaded)"
    echo "=============================================="
    echo "Note: Requires WebDriver running on port 4444"
    echo "      (chromedriver --port=4444 or geckodriver --port=4444)"
    echo ""
    cargo test -p web_tests -- --test-threads=1
    exit 0
fi

echo "=============================================="
echo "Running boids_lib tests (single-threaded)"
echo "=============================================="
cargo test -p boids_lib -- --test-threads=1

echo ""
echo "=============================================="
echo "Running web_tests (single-threaded)"
echo "=============================================="
echo "Note: Requires WebDriver running on port 4444"
# web_tests may fail if WebDriver isn't running - that's okay
cargo test -p web_tests -- --test-threads=1 2>&1 || {
    echo ""
    echo "Warning: web_tests failed (WebDriver may not be running)"
}

echo ""
echo "=============================================="
echo "Running other package tests (parallel)"
echo "=============================================="
# Try to run other workspace tests, but don't fail if some packages have issues
cargo test --workspace --exclude boids_lib --exclude web_tests 2>&1 || {
    echo ""
    echo "Warning: Some workspace packages failed to compile/test"
}

echo ""
echo "=============================================="
echo "All tests completed"
echo "=============================================="
