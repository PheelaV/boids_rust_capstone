#!/bin/bash
#
# Alternative build and serve script using basic-http-server
#
# This is a simpler option if you prefer basic-http-server over miniserve
#

set -e

GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${BLUE}=== Boids WebAssembly Build & Serve (basic-http-server) ===${NC}\n"

# Check if wasm-pack is installed
if ! command -v wasm-pack &> /dev/null; then
    echo -e "${YELLOW}Warning: wasm-pack is not installed${NC}"
    echo "Install it with: cargo install wasm-pack"
    exit 1
fi

# Check if basic-http-server is installed
if ! command -v basic-http-server &> /dev/null; then
    echo -e "${YELLOW}basic-http-server is not installed. Installing...${NC}"
    cargo install basic-http-server --quiet
fi

# Build WASM
echo -e "${GREEN}[1/2] Building WebAssembly module...${NC}"
cd "$(dirname "$0")/../boids_wasm"
wasm-pack build --target web --out-dir ../web/pkg

echo -e "${GREEN}Build complete!${NC}\n"

# Serve the demo
echo -e "${GREEN}[2/2] Starting web server...${NC}"
cd ../web

echo -e "${BLUE}Starting server on http://127.0.0.1:4000${NC}"
echo -e "${BLUE}Press Ctrl+C to stop the server${NC}\n"

basic-http-server -a 127.0.0.1:4000 .
