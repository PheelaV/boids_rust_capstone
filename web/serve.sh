#!/bin/bash
#
# Build and serve script for boids WebAssembly demo
#
# This script:
# 1. Builds the WebAssembly module using wasm-pack
# 2. Serves the web demo using miniserve (Rust-based web server)
#

set -e

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== Boids WebAssembly Build & Serve ===${NC}\n"

# Check if wasm-pack is installed
if ! command -v wasm-pack &> /dev/null; then
    echo -e "${YELLOW}Warning: wasm-pack is not installed${NC}"
    echo "Install it with: cargo install wasm-pack"
    exit 1
fi

# Check if miniserve is installed
if ! command -v miniserve &> /dev/null; then
    echo -e "${YELLOW}miniserve is not installed. Installing...${NC}"
    cargo install miniserve --quiet
fi

# Build WASM
echo -e "${GREEN}[1/2] Building WebAssembly module...${NC}"
cd "$(dirname "$0")/../boids_wasm"
wasm-pack build --target web --out-dir ../web/pkg

echo -e "${GREEN}Build complete!${NC}\n"

# Serve the demo
echo -e "${GREEN}[2/2] Starting web server...${NC}"
cd ../web

echo -e "${BLUE}Starting miniserve on http://localhost:8080${NC}"
echo -e "${BLUE}Press Ctrl+C to stop the server${NC}\n"

# Start miniserve with nice defaults
miniserve . \
    --port 8080 \
    --index index.html \
    --color-scheme dark \
    --header "Cache-Control: no-cache" \
    --verbose
