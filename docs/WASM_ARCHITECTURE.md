# WebAssembly Architecture for Boids Simulation

## Overview

This document outlines the approach for compiling the boids simulation to WebAssembly (WASM) to enable running it in web browsers via JavaScript.

## Current Architecture Analysis

### Existing Components

1. **`boids_lib`** - Core simulation library
   - Pure Rust logic for boid behaviors
   - Spatial hashing for performance
   - Boundary conditions and field-of-vision
   - Uses `glam` for vector math
   - DBSCAN clustering (via `linfa-clustering`)

2. **`boids_app`** - Desktop visualization
   - Uses `nannou` (wgpu-based graphics)
   - Not browser-compatible
   - Will not be compiled to WASM

### Dependencies Compatibility Assessment

**WASM-Compatible:**
- `glam` - Math library (WASM-friendly)
- `serde` / `serde_json` - Serialization (WASM-friendly)
- `itertools` - Iterator utilities (WASM-friendly)

**Requires Adaptation:**
- `rand` - Needs `wasm-bindgen` feature or getrandom with js feature
- `chrono` - May need js feature enabled

**Potentially Problematic:**
- `linfa-clustering` / `linfa-nn` - May not be WASM-compatible
- `csv` - File I/O doesn't translate to browser (need to adapt)
- `flat_spatial` - Needs verification
- `petal-clustering` - Needs verification

## WebAssembly Compilation Strategy

### Architecture

```
┌─────────────────────────────────────────┐
│           Browser (JavaScript)          │
│  ┌───────────────────────────────────┐  │
│  │   Canvas Rendering / WebGL        │  │
│  └───────────────┬───────────────────┘  │
│                  │                       │
│  ┌───────────────▼───────────────────┐  │
│  │   JavaScript Bindings (glue)     │  │
│  └───────────────┬───────────────────┘  │
└──────────────────┼───────────────────────┘
                   │ wasm-bindgen
┌──────────────────▼───────────────────────┐
│         WebAssembly Module               │
│  ┌───────────────────────────────────┐  │
│  │      boids_wasm crate             │  │
│  │  - Simulation wrapper             │  │
│  │  - JS-friendly API                │  │
│  │  - Memory management              │  │
│  └───────────────┬───────────────────┘  │
│                  │                       │
│  ┌───────────────▼───────────────────┐  │
│  │      boids_lib (core)             │  │
│  │  - Boid behaviors                 │  │
│  │  - Spatial hashing                │  │
│  │  - Flock management               │  │
│  └───────────────────────────────────┘  │
└──────────────────────────────────────────┘
```

### Implementation Plan

#### Phase 1: Create WASM Crate

1. **Create `boids_wasm` crate**
   ```toml
   [package]
   name = "boids_wasm"
   version = "0.1.0"

   [lib]
   crate-type = ["cdylib", "rlib"]

   [dependencies]
   boids_lib = { path = "../boids_lib", default-features = false }
   wasm-bindgen = "0.2"
   serde = { version = "1", features = ["derive"] }
   serde-wasm-bindgen = "0.6"
   getrandom = { version = "0.2", features = ["js"] }
   ```

2. **Feature flags in `boids_lib`**
   - Add `default = ["std"]` feature
   - Add `wasm` feature to disable incompatible dependencies
   - Make clustering optional (not critical for basic simulation)

#### Phase 2: Design JavaScript API

**Core API Surface:**

```rust
// Rust side (boids_wasm)
#[wasm_bindgen]
pub struct WasmSimulation {
    flock: Flock,
    options: RunOptions,
}

#[wasm_bindgen]
impl WasmSimulation {
    #[wasm_bindgen(constructor)]
    pub fn new(config: JsValue) -> Result<WasmSimulation, JsValue>;

    pub fn update(&mut self);

    pub fn get_boids(&self) -> JsValue;  // Returns serialized boid positions/velocities

    pub fn set_option(&mut self, key: String, value: f32);

    pub fn add_boid(&mut self, x: f32, y: f32, vx: f32, vy: f32);

    pub fn remove_boid(&mut self, id: usize);

    pub fn reset(&mut self);

    pub fn get_stats(&self) -> JsValue;  // Returns simulation statistics
}
```

**JavaScript side:**
```javascript
import init, { WasmSimulation } from './boids_wasm.js';

// Initialize WASM
await init();

// Create simulation
const sim = new WasmSimulation({
    init_boids: 200,
    separation_coefficient: 1.5,
    cohesion_coefficient: 1.0,
    alignment_coefficient: 1.0,
    // ... more options
});

// Game loop
function animate() {
    sim.update();
    const boids = sim.get_boids();
    renderBoids(boids);  // Canvas rendering
    requestAnimationFrame(animate);
}
```

#### Phase 3: Handle WASM-Specific Concerns

**Random Number Generation:**
- Use `getrandom` with `js` feature for browser entropy
- Seed-based RNG for deterministic tests
- Replace `rand::thread_rng()` with WASM-compatible alternatives

**Memory Management:**
- Keep large data in Rust (avoid copying to JS)
- Use views/slices where possible
- Serialize only what's needed for rendering

**File I/O:**
- Replace CSV file operations with serialization to JavaScript
- Let browser handle downloads via Blob API

**Clustering (Optional):**
- Feature-gate DBSCAN clustering
- Implement simple distance-based clustering if needed
- Or skip clustering for WASM version

#### Phase 4: Build Process

**Tools:**
```bash
cargo install wasm-pack
```

**Build command:**
```bash
cd boids_wasm
wasm-pack build --target web --out-dir ../web/pkg
```

**Output:**
- `boids_wasm_bg.wasm` - Compiled WASM binary
- `boids_wasm.js` - JavaScript glue code
- `boids_wasm.d.ts` - TypeScript definitions

#### Phase 5: Testing Strategy (TDD)

**Unit Tests (Rust):**
```rust
#[cfg(test)]
mod tests {
    use wasm_bindgen_test::*;

    #[wasm_bindgen_test]
    fn test_simulation_creation() {
        // Test WASM simulation can be created
    }

    #[wasm_bindgen_test]
    fn test_boid_update() {
        // Test update produces valid state
    }
}
```

**Integration Tests (Browser):**
- Use `wasm-bindgen-test` with headless browser
- Test JavaScript ↔ Rust boundary
- Verify serialization/deserialization
- Performance benchmarks

## Hypotheses to Test

### Hypothesis 1: Core boids_lib is WASM-compatible
**Test:** Try compiling boids_lib with `--target wasm32-unknown-unknown`
**Expected:** May fail due to dependencies
**Mitigation:** Add feature flags to disable problematic dependencies

### Hypothesis 2: Spatial hashing performs well in WASM
**Test:** Benchmark WASM vs native simulation with 1000+ boids
**Expected:** 70-90% of native performance
**Mitigation:** Optimize hot paths, reduce allocations

### Hypothesis 3: JS ↔ Rust boundary is a bottleneck
**Test:** Profile time spent in serialization vs simulation
**Expected:** Serialization should be < 10% of frame time
**Mitigation:** Use shared memory, minimize data transfer

### Hypothesis 4: Browser RNG is sufficient for simulation
**Test:** Compare determinism and distribution of WASM RNG
**Expected:** Should match native with same seed
**Mitigation:** Use seed-based RNG, avoid `thread_rng()`

## Performance Considerations

**Optimization Strategies:**
1. Keep simulation state in WASM memory (avoid copying)
2. Only serialize visible boid data
3. Use Float32Array for efficient data transfer
4. Batch updates (update multiple frames per render if needed)
5. Use Web Workers for parallel simulation (advanced)

**Target Performance:**
- 1000 boids @ 60 FPS on modern hardware
- < 16ms per frame (simulation + rendering)
- < 1ms for data serialization

## Success Criteria

1. ✅ Compile boids_lib to WASM with all core features
2. ✅ Expose clean JavaScript API
3. ✅ Achieve > 200 boids @ 60 FPS in browser
4. ✅ Tests pass in both native and WASM targets
5. ✅ Example web application demonstrates all features
6. ✅ Documentation for web integration

## Next Steps

1. Add feature flags to `boids_lib/Cargo.toml`
2. Create `boids_wasm` crate skeleton
3. Write failing tests for WASM API (TDD)
4. Implement WASM bindings incrementally
5. Create minimal HTML/JS example
6. Optimize and benchmark
7. Create comprehensive web demo

## References

- [wasm-bindgen Guide](https://rustwasm.github.io/wasm-bindgen/)
- [Rust and WebAssembly Book](https://rustwasm.github.io/docs/book/)
- [wasm-pack](https://rustwasm.github.io/wasm-pack/)
- [web-sys API](https://rustwasm.github.io/wasm-bindgen/api/web_sys/)
