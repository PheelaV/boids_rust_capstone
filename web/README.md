# Boids Simulation - Web Demo

This directory contains a web-based demo of the boids simulation running in WebAssembly.

## Prerequisites

1. **Rust toolchain** with `wasm32-unknown-unknown` target:
   ```bash
   rustup target add wasm32-unknown-unknown
   ```

2. **wasm-pack** - Tool for building WebAssembly:
   ```bash
   cargo install wasm-pack
   ```

3. **A local web server** (choose one):
   - Python: `python -m http.server 8080`
   - Node.js: `npx http-server -p 8080`
   - Rust: `cargo install basic-http-server && basic-http-server`

## Building

From the `boids_wasm` directory, run:

```bash
cd ../boids_wasm
wasm-pack build --target web --out-dir ../web/pkg
```

This will:
- Compile the Rust code to WebAssembly
- Generate JavaScript bindings
- Output everything to `web/pkg/`

## Running

1. Build the WebAssembly module (see above)

2. Start a local web server from the `web` directory:
   ```bash
   cd ../web
   python -m http.server 8080
   ```

3. Open your browser to:
   ```
   http://localhost:8080
   ```

## Features

### Interactive Controls

- **Behavior Toggles**: Enable/disable individual flocking rules
  - Separation: Avoid crowding neighbors
  - Cohesion: Steer towards the average position of neighbors
  - Alignment: Steer towards the average heading of neighbors
  - Wander: Add random exploration behavior

- **Parameter Sliders**: Adjust simulation parameters in real-time
  - Separation, Cohesion, and Alignment coefficients
  - Maximum speed

- **Actions**:
  - Reset: Restart simulation with random positions
  - Pause/Resume: Pause and resume the simulation

### Visual Features

- Boids rendered as triangles pointing in their direction of travel
- Color-coded by speed (red = slow, green = fast)
- Real-time statistics display (FPS, boid count, frame number)

## Performance Notes

- The WASM build is optimized for size (opt-level = "s")
- Clustering features are disabled in WASM builds to reduce bundle size
- Expected performance: 200+ boids at 60 FPS on modern hardware

## Troubleshooting

**"Failed to fetch" or CORS errors:**
- Make sure you're using a local web server, not opening `index.html` directly
- Check that wasm-pack output is in the `web/pkg/` directory

**Blank screen:**
- Check browser console for errors
- Verify WebAssembly is supported (most modern browsers)
- Try refreshing the page

**Low FPS:**
- Reduce the number of boids (edit `init_boids` in `main.js`)
- Check if hardware acceleration is enabled in your browser

## Browser Compatibility

Tested and working on:
- Chrome/Edge 90+
- Firefox 89+
- Safari 15+

WebAssembly is supported in all modern browsers. For older browsers, you may need to use a polyfill.

## Development

To rebuild after making changes to the Rust code:

```bash
cd boids_wasm
wasm-pack build --target web --out-dir ../web/pkg
```

Then refresh your browser (may need a hard refresh: Ctrl+Shift+R or Cmd+Shift+R).

## Architecture

The web demo uses:
- **Rust/WASM**: Simulation logic (boids_lib + boids_wasm)
- **JavaScript**: Rendering and UI controls
- **HTML5 Canvas**: 2D graphics rendering
- **wasm-bindgen**: Rust ↔ JavaScript interop

See `../docs/WASM_ARCHITECTURE.md` for detailed architecture documentation.
