# Boids Simulation in Rust

A high-performance implementation of Craig Reynolds' Boids flocking algorithm in Rust, featuring real-time visualization and multiple spatial optimization strategies.

## Overview

This project simulates emergent flocking behavior based on three simple rules:
- **Separation**: Avoid crowding neighbors
- **Alignment**: Steer towards average heading of neighbors
- **Cohesion**: Steer towards average position of neighbors

The simulation supports various boundary conditions, field-of-vision constraints, wander behavior, and includes advanced features like trajectory replay and ghost trails.

## Project Structure

This workspace contains two main crates:

### `boids_lib`
Core simulation library implementing:
- Boid entities with position, velocity, and acceleration
- Flocking rules (separation, cohesion, alignment, wander)
- Multiple spatial indexing strategies:
  - Naive O(n²) tracker for small flocks
  - Spatial hashing for efficient neighbor queries
  - Replay tracker for trajectory playback
- Boundary conditions (toroidal, reflective, absorbing, repulsive)
- DBSCAN clustering for flock identification
- Data recording and CSV export

### `boids_app`
Interactive visualization application using [nannou](https://nannou.cc/):
- Real-time rendering with customizable colors
- Interactive controls via egui GUI
- Keyboard shortcuts for parameter tweaking
- Ghost trail visualization
- Debug overlays (grid, distances, labels)
- Video frame export capability

## Features

- **Performance**: Optimized spatial hashing enables simulation of thousands of boids at 60+ FPS
- **Configurability**: Extensive runtime options via config files or command-line arguments
- **Visualization**: Real-time rendering with field-of-vision cones, neighbor connections, and cluster highlighting
- **Analysis**: Built-in data collection and CSV export for behavioral analysis
- **Research-Ready**: Supports parameter sweeps, reproducible runs, and trajectory replay

## Requirements

- Rust 1.91.1+ (2021 edition)
- Graphics drivers supporting nannou/wgpu

## Building

```bash
cargo build --release
```

## Running

### Interactive Visualization
```bash
cargo run --release --bin boids_app
```

### With Configuration
```bash
cargo run --release --bin boids_app -- --config-path path/to/config.toml
```

## Configuration

The simulation accepts configuration via TOML files. Key parameters:

```toml
no_boids = 200
sensory_distance = 100.0
separation_coefficient = 1.5
cohesion_coefficient = 1.0
alignment_coefficient = 1.0
max_speed = 4.0
min_speed = 2.0
init_width = 1400
init_height = 900
```

## Keyboard Controls

- **Space**: Pause/resume simulation
- **R**: Restart with new random positions
- **1-4**: Toggle individual rules (alignment, cohesion, separation, wander)
- **I/D**: Increase/decrease boid count (double/halve)
- **C**: Toggle control panel
- **V**: Toggle field of vision visualization
- **F**: Toggle flock clustering
- **F8-F12**: Toggle debug overlays
- **Left Click**: Select boid for detailed view
- **Mouse Wheel**: Adjust window zoom

## Dependencies

Key dependencies and their versions:

### Core Simulation
- `glam = "0.17"` - Vector math
- `rand = "0.8"` - Random number generation
- `linfa = "0.6"` - Machine learning (DBSCAN clustering)
- `ndarray = "0.15"` - N-dimensional arrays
- `flat_spatial = "0.6"` - Spatial indexing
- `once_cell = "1.20"` - Lazy static initialization

### Visualization
- `nannou = "0.19"` - Creative coding framework
- `nannou_egui = "0.19"` - GUI integration
- `splines = "4.2.0"` - Trajectory interpolation

### I/O & Serialization
- `csv = "1.3"` - Data export
- `serde = "1"` - Serialization
- `chrono = "0.4"` - Timestamps
- `clap = "4.5"` - CLI argument parsing

## Testing

The project includes a comprehensive integration test suite (49 tests) covering core simulation behavior:

- **`simulation_determinism.rs`** - Deterministic behavior, movement, velocity bounds, boundary containment
- **`boid_behaviors.rs`** - Separation, cohesion, alignment, wander, field-of-vision
- **`boundary_conditions.rs`** - Toroidal, reflective, absorbing, repulsive boundaries
- **`spatial_hashing.rs`** - Spatial optimization correctness vs naive implementation
- **`capture_replay.rs`** - Headless data capture and trajectory replay

Run tests:
```bash
# Run all tests (unit + integration)
cargo test test_

# Or run specific test suites
cargo test --lib                    # Unit tests only
cargo test --test boid_behaviors    # Specific integration test
```

**Note**: Due to a workspace configuration quirk, `cargo test` alone only runs unit tests. Use `cargo test test_` to run all tests including integration tests.

Run benchmarks:
```bash
cargo bench
```

## Performance

On an M1 Pro:
- 200 boids: 120+ FPS
- 1000 boids: 60+ FPS
- 32,768 boids: ~6 FPS

Performance scales roughly O(n log n) with spatial hashing enabled.

## References

- Reynolds, C. W. (1987). "Flocks, herds and schools: A distributed behavioral model." *SIGGRAPH '87*
- [Red Blob Games - Implementation of Influence Maps](https://www.redblobgames.com/grids/hexagons/)
- [Boids Pseudocode](http://www.kfish.org/boids/pseudocode.html)

## License

This project was created as a BSc capstone project. See source files for licensing information.

## Changelog

### 2025-11-23 - Modernization Update
- Updated all dependencies to latest compatible versions
- Replaced `lazy_static` with `once_cell` for better compile times
- Fixed egui API changes (`egui::color` → `egui::ecolor`)
- Fixed spline/glam version compatibility
- Added comprehensive integration test suite (49 tests across 5 test files)
- Cleaned up trailing whitespace and applied rustfmt
- Fixed typo: `trangles()` → `triangles()`
- Improved error handling in several modules

### Original (2022)
- Initial implementation with spatial hashing
- DBSCAN-based flock clustering
- Replay and visualization features
