# Nannou WebAssembly Investigation Report

**Date:** 2025-01-23
**Status:** Experimental / In Progress
**Recommendation:** Use existing canvas-based WASM (boids_wasm) for production; nannou WASM for future exploration

## Executive Summary

Nannou CAN run in WebAssembly, but support is **experimental** and requires significant configuration. After thorough investigation and attempted implementation, here are the findings:

### ✅ What Works
- Simple nannou sketches compile to WASM
- WebGL rendering through wgpu works in browsers
- Basic graphics and animation function
- Working examples exist (astro-nannou-starter, nannou-web-template)

### ⚠️ Challenges Encountered
- Complex web-sys feature configuration (100+ GPU features needed)
- Feature names must be exact snake_case (e.g., `gpu_map_mode` not `GpuMapMode`)
- wgpu version compatibility with web-sys features
- Large bundle sizes (~several MB for WASM module)
- Limited documentation for WASM-specific issues

### ❌ Limitations
- No audio support (CPAL WASM is experimental)
- No file I/O operations
- Cannot use nannou_egui (UI library) in WASM
- Performance lower than native (WebGL vs native graphics)

## Attempted Implementation

### What Was Created
- `boids_app_web/` - New crate for nannou WASM version
- Cargo.toml with extensive web-sys features (WebGPU/WebGL)
- lib.rs with simplified boids rendering using nannou
- HTML/CSS interface matching desktop controls
- Build configuration for wasm-pack

### Build Issues
```
error: failed to select a version for `web-sys`.
package `boids_app_web` depends on `web-sys` with feature `GpuMapMode`
but `web-sys` does not have that feature.
package `web-sys` does have feature `gpu_map_mode`
```

**Root Cause:** web-sys requires exact snake_case feature names, and wgpu@0.17.2 (used by nannou@0.19) needs 100+ specific GPU features to be manually enabled.

## Comparison: Canvas vs Nannou WASM

### Current Working Solution: boids_wasm (Canvas 2D)

**Pros:**
- ✅ Works out of the box
- ✅ Small bundle size (~500KB)
- ✅ Simple to build and deploy
- ✅ Compatible with all browsers
- ✅ Fast development iteration
- ✅ Easy to debug

**Cons:**
- ❌ Manual rendering code (no nannou helpers)
- ❌ Different codebase from desktop app
- ❌ Canvas 2D API limitations (no advanced effects)

### Proposed Nannou Solution: boids_app_web

**Pros:**
- ✅ Code closer to desktop app
- ✅ Same rendering pipeline as desktop
- ✅ WebGPU capabilities (when available)
- ✅ Nannou's graphics abstractions

**Cons:**
- ❌ Complex configuration
- ❌ Larger bundle size (estimated 3-5MB)
- ❌ Experimental/unstable
- ❌ Requires WebGPU or WebGL2
- ❌ Difficult to debug
- ❌ Limited browser support for WebGPU

## Working Examples (Reference)

### 1. astro-nannou-starter (January 2025)
- **Repository:** https://github.com/JulianCataldo/astro-nannou-starter
- **Live Demo:** https://juliancataldo.github.io/astro-nannou-starter/
- **Tech Stack:** Astro + Vite + rsw (Rust to WASM)
- **Features:** Multiple sketches, hot-reload, production builds
- **Use Case:** Best reference for modern nannou WASM setup

### 2. nannou-web-template (2022)
- **Repository:** https://github.com/tomoyanonymous/nannou-web-template
- **Tech Stack:** Webpack-based
- **Status:** Older approach, still functional

## Recommendations

### Short Term (Now)
**Use the existing `boids_wasm` crate (Canvas 2D approach)**

Rationale:
- Already working and tested
- Proven stable for browser deployment
- Fast iteration and debugging
- Good enough for most use cases
- Much smaller bundle size

### Medium Term (3-6 months)
**Monitor nannou WASM maturity**

Actions:
- Watch nannou repository for WASM improvements
- Test astro-nannou-starter approach
- Evaluate if WebGPU browser adoption increases
- Consider contributing to nannou WASM development

### Long Term (6+ months)
**Migrate to nannou WASM when stable**

Conditions for migration:
- nannou documents WASM setup clearly
- wgpu/web-sys feature configuration is simplified
- WebGPU support is widespread (>70% browser coverage)
- Bundle sizes are optimized
- Real-world production examples exist

## Alternative Approaches

If nannou WASM doesn't mature:

### 1. Bevy Engine (Recommended Alternative)
- **Status:** Production-ready WASM support
- **Pros:** Excellent WASM docs, many browser examples, active community
- **Cons:** Different API from nannou, game-focused architecture
- **Use Case:** If you need reliable WASM graphics framework
- **Resources:** https://bevy-cheatbook.github.io/platforms/wasm.html

### 2. Direct wgpu (Advanced)
- **Status:** Full WASM support
- **Pros:** Maximum control, no framework overhead
- **Cons:** Much more boilerplate, steeper learning curve
- **Use Case:** Custom rendering pipelines, advanced graphics

### 3. Macroquad (Simplest)
- **Status:** Excellent WASM support
- **Pros:** Extremely simple API, single-command deploy
- **Cons:** Less powerful than nannou/Bevy
- **Use Case:** Simple 2D games, creative sketches

## Technical Details

### Required web-sys Features (Partial List)
nannou + wgpu require enabling 100+ web-sys features, including:

```toml
web-sys = { version = "0.3", features = [
    # Core
    "Document", "Window", "Navigator",
    # Canvas
    "HtmlCanvasElement", "HtmlElement",
    # WebGPU (100+ features, snake_case names)
    "gpu", "gpu_adapter", "gpu_device",
    "gpu_command_encoder", "gpu_render_pass_encoder",
    "gpu_buffer", "gpu_texture", "gpu_sampler",
    # ... 90+ more GPU features
    # WebGL fallback
    "WebGl2RenderingContext", "WebGlRenderingContext",
]}
```

### Build Command
```bash
cd boids_app_web
wasm-pack build --target web --out-dir web/pkg
```

### Typical Bundle Sizes
- **boids_wasm (Canvas):** ~500KB
- **nannou WASM (estimated):** ~3-5MB
- **Bevy WASM:** ~2-4MB

## Current Project State

### What's Implemented and Working
1. ✅ **boids_wasm** - Canvas 2D WASM (fully functional)
   - Location: `/boids_wasm`
   - Demo: `/web/index.html`
   - Status: Production-ready

2. ✅ **Desktop App** - Nannou native (fully functional)
   - Location: `/boids_app`
   - Status: Production-ready

3. ⚠️ **boids_app_web** - Nannou WASM (experimental, incomplete)
   - Location: `/boids_app_web`
   - Status: Build fails on web-sys features
   - Next Steps: Feature name corrections, extensive testing

### File Structure Created
```
boids_app_web/
├── Cargo.toml          # Dependencies with extensive web-sys features
├── src/
│   └── lib.rs          # Simplified nannou app (no egui)
└── web/
    └── index.html      # HTML interface with controls
```

## Lessons Learned

1. **WASM is Moving Fast:** Browser APIs, wgpu, and nannou are all evolving rapidly
2. **Feature Flags Matter:** Exact feature names are critical for web-sys
3. **Bundle Size Important:** WebAssembly modules can get large quickly
4. **Documentation Gaps:** Experimental features often lack clear guides
5. **Fallbacks Work Well:** Simple Canvas 2D is often good enough

## Resources

### Nannou WASM
- Issue #475: https://github.com/nannou-org/nannou/issues/475
- PR #821: https://github.com/nannou-org/nannou/pull/821
- CPAL WASM Example: https://github.com/nannou-org/cpal_wasm_example

### Working Examples
- astro-nannou-starter: https://github.com/JulianCataldo/astro-nannou-starter
- Live demo: https://juliancataldo.github.io/astro-nannou-starter/

### Alternative Frameworks
- Bevy WASM Guide: https://bevy-cheatbook.github.io/platforms/wasm.html
- Bevy Examples: https://bevy.org/examples/
- wgpu Documentation: https://wgpu.rs/

### WebGPU Support
- Chrome/Edge: ✅ April 2023+
- Safari: ✅ June 2025+ (Safari 26)
- Firefox: ✅ July 2025+ (Firefox 141)

## Conclusion

While nannou WASM is technically possible and working examples exist, the current implementation complexity and experimental nature make it **not recommended for production use** at this time.

**Recommended Path Forward:**
1. Continue using `boids_wasm` (Canvas 2D) for web deployments
2. Maintain `boids_app` (nannou native) for desktop
3. Monitor nannou WASM development
4. Revisit nannou WASM in 6-12 months when ecosystem matures

The good news: **You already have a working WASM solution** that provides excellent browser-based boids simulation!
