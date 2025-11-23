#!/bin/bash
#
# Build script for boids_wasm WebAssembly module
#
# This script compiles the Rust code to WebAssembly and generates
# JavaScript bindings using wasm-pack.
#

set -e

echo "Building boids_wasm for WebAssembly..."

# Check if wasm-pack is installed
if ! command -v wasm-pack &> /dev/null; then
    echo "Error: wasm-pack is not installed"
    echo "Install it with: cargo install wasm-pack"
    exit 1
fi

# Check if wasm32-unknown-unknown target is installed
if ! rustup target list | grep -q "wasm32-unknown-unknown (installed)"; then
    echo "Installing wasm32-unknown-unknown target..."
    rustup target add wasm32-unknown-unknown
fi

# Build for web target
echo "Running wasm-pack build..."
wasm-pack build --target web --out-dir ../web/pkg

echo "Build complete! Output is in ../web/pkg/"
echo ""
echo "To run the web demo:"
echo "  cd ../web"
echo "  python -m http.server 8080"
echo "  # Then open http://localhost:8080 in your browser"
