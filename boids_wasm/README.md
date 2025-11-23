# boids_wasm

WebAssembly bindings for the boids simulation library, enabling browser-based flocking simulations.

## Overview

This crate provides JavaScript/TypeScript bindings for `boids_lib`, allowing the boids flocking simulation to run in web browsers via WebAssembly. It uses [wasm-bindgen](https://rustwasm.github.io/wasm-bindgen/) to create a seamless interface between Rust and JavaScript.

## Features

- **Browser-compatible**: Runs in any modern web browser with WebAssembly support
- **High performance**: Compiled to optimized WebAssembly for near-native speed
- **Small bundle size**: Optimized for size with clustering features disabled
- **Type-safe**: Generates TypeScript definitions for type-safe JavaScript integration
- **Real-time control**: Adjust simulation parameters and toggle behaviors on the fly

## Building

### Prerequisites

1. Rust toolchain (1.91.1+)
2. wasm32-unknown-unknown target
3. wasm-pack

```bash
# Install wasm32 target
rustup target add wasm32-unknown-unknown

# Install wasm-pack
cargo install wasm-pack
```

### Build for Web

```bash
# Run the build script
chmod +x build.sh
./build.sh

# Or run wasm-pack directly
wasm-pack build --target web --out-dir ../web/pkg
```

This will generate:
- `boids_wasm_bg.wasm` - The WebAssembly binary
- `boids_wasm.js` - JavaScript bindings
- `boids_wasm.d.ts` - TypeScript type definitions

### Build Options

**Target platforms:**
- `--target web` - For vanilla JavaScript/ES modules (default)
- `--target bundler` - For webpack/rollup/parcel
- `--target nodejs` - For Node.js

**Optimization:**
- Release build (default): Optimized for size
- Debug build: `wasm-pack build --dev`

## Usage

### JavaScript/TypeScript

```javascript
import init, { WasmSimulation } from './pkg/boids_wasm.js';

// Initialize the WASM module
await init();

// Create a simulation with custom config
const config = {
    init_boids: 200,
    window_width: 800,
    window_height: 600,
    separation_coefficient: 1.5,
    cohesion_coefficient: 1.0,
    alignment_coefficient: 1.0,
    max_speed: 4.0,
    min_speed: 2.0,
    sensory_distance: 60.0,
    wander_on: true,
    wander_coefficient: 0.5
};

const sim = WasmSimulation.from_config(config);

// Or use default configuration
const sim = new WasmSimulation();

// Update simulation
sim.update();

// Get boid positions for rendering
const boids = sim.get_boids();
boids.forEach(boid => {
    console.log(`Boid ${boid.id}: (${boid.x}, ${boid.y})`);
});

// Adjust parameters
sim.set_separation_coefficient(2.0);
sim.set_max_speed(5.0);

// Toggle behaviors
sim.toggle_separation();
sim.toggle_cohesion();

// Get statistics
const stats = sim.get_stats();
console.log(`Boids: ${stats.boid_count}, Frame: ${stats.frame_count}`);

// Reset simulation
sim.reset();
```

### TypeScript

TypeScript definitions are automatically generated:

```typescript
import init, { WasmSimulation, SimulationConfig } from './pkg/boids_wasm.js';

const config: SimulationConfig = {
    init_boids: 200,
    window_width: 800,
    window_height: 600,
    // ... other config fields
};

const sim = WasmSimulation.from_config(config);
```

## API Reference

### `WasmSimulation`

Main simulation class.

#### Constructors

- `new WasmSimulation()` - Create with default configuration
- `WasmSimulation.from_config(config: JsValue)` - Create from JavaScript config object

#### Methods

**Simulation Control:**
- `update()` - Advance simulation by one step
- `reset()` - Reset simulation with current configuration
- `get_boid_count(): number` - Get number of boids
- `get_frame_count(): number` - Get current frame number

**Data Access:**
- `get_boids(): JsValue` - Get array of boid states (id, x, y, vx, vy)
- `get_stats(): JsValue` - Get simulation statistics

**Boid Management:**
- `add_boid(x: number, y: number, vx: number, vy: number)` - Add single boid
- `add_boids(count: number)` - Add multiple boids
- `remove_boid(): boolean` - Remove last boid
- `remove_boids(count: number): number` - Remove multiple boids
- `double_boids()` - Double the boid count
- `halve_boids()` - Halve the boid count
- `set_boid_count(count: number)` - Set exact boid count

**Parameter Control:**
- `set_separation_coefficient(value: number)` - Set separation strength
- `set_cohesion_coefficient(value: number)` - Set cohesion strength
- `set_alignment_coefficient(value: number)` - Set alignment strength
- `set_max_speed(value: number)` - Set maximum boid speed

**Behavior Toggles:**
- `toggle_separation()` - Enable/disable separation
- `toggle_cohesion()` - Enable/disable cohesion
- `toggle_alignment()` - Enable/disable alignment
- `toggle_wander()` - Enable/disable wander behavior

### `SimulationConfig`

Configuration object for creating simulations.

```typescript
interface SimulationConfig {
    init_boids: number;           // Number of boids (default: 200)
    window_width: number;          // Width of simulation space
    window_height: number;         // Height of simulation space
    separation_coefficient: number; // Separation strength (default: 1.5)
    cohesion_coefficient: number;  // Cohesion strength (default: 1.0)
    alignment_coefficient: number; // Alignment strength (default: 1.0)
    max_speed: number;            // Maximum speed (default: 4.0)
    min_speed: number;            // Minimum speed (default: 2.0)
    sensory_distance: number;     // Perception radius (default: 60.0)
    wander_on: boolean;           // Enable wander (default: true)
    wander_coefficient: number;   // Wander strength (default: 0.5)
}
```

## Testing

### Unit Tests (Rust)

```bash
cargo test
```

### WASM Tests (Browser)

```bash
wasm-pack test --headless --firefox
wasm-pack test --headless --chrome
```

Tests include:
- Simulation creation and initialization
- Update loop functionality
- Parameter setting and behavior toggles
- Data serialization
- Reset functionality
- Custom configuration

## Performance

**Optimization Strategy:**
- Clustering features disabled (reduces bundle size ~30%)
- Size-optimized compilation (`opt-level = "s"`)
- Link-time optimization (LTO)
- getrandom configured for browser entropy

**Expected Performance:**
- 200 boids @ 60 FPS on modern hardware
- ~500KB WASM bundle (gzipped: ~150KB)
- < 1ms per frame serialization overhead

**Profiling:**
```javascript
const start = performance.now();
sim.update();
const updateTime = performance.now() - start;
console.log(`Update took ${updateTime}ms`);
```

## Bundle Size

- With all optimizations: ~500KB (.wasm)
- Gzipped: ~150KB
- JavaScript glue code: ~50KB

**Size optimization checklist:**
- ✅ Clustering disabled (`--no-default-features`)
- ✅ Size-optimized build (`opt-level = "s"`)
- ✅ LTO enabled
- ✅ Minimal panic infrastructure

## Troubleshooting

**Build fails with "getrandom" error:**
- Make sure the "js" feature is enabled for getrandom
- Check that you're building with `--target wasm32-unknown-unknown`

**"wasm validation error" in browser:**
- Try rebuilding with `--dev` flag for better error messages
- Check browser console for detailed error information

**Poor performance:**
- Reduce number of boids in configuration
- Check if simulation is running in release mode
- Profile to identify bottlenecks

**Memory issues:**
- Avoid holding references to large datasets in JavaScript
- Let Rust manage memory where possible
- Use `get_boids()` only when needed for rendering

## Architecture

See [../docs/WASM_ARCHITECTURE.md](../docs/WASM_ARCHITECTURE.md) for detailed architecture documentation.

## Examples

See the [../web](../web) directory for a complete working example with HTML/CSS/JavaScript.

## License

MIT - See LICENSE.MD in repository root
