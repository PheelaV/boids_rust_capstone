# Plan: Introducing SI Units to Boids Simulation

## Executive Summary

This plan outlines the transformation of the boids simulation from a dimensionless, frame-dependent system to a physically-grounded simulation with proper SI units. The goal is to enable real-world physical parameters (speeds in m/s or km/h, accelerations in m/s²) while maintaining frame-rate independence.

## Current State Analysis

### Existing Problems

1. **No Unit System**: All values are dimensionless f32 (pixels, frames)
   - Position: pixels (window coordinates)
   - Velocity: pixels/frame
   - Acceleration: pixels/frame²

2. **Frame-Rate Dependency**: Physics behavior changes with FPS
   - No delta time (dt) tracking
   - `baseline_speed` acts as ad-hoc frame multiplier (line 481 in boid.rs)
   - Different FPS = different simulation behavior

3. **No Physical Scale**: No mapping between screen space and real world
   - Window coordinates are arbitrary
   - No concept of "how big is a boid in reality?"

4. **Mixed Dimensional Analysis**: Parameters lack consistent units
   - `sensory_distance: 60.0` (pixels?)
   - `max_speed: 4.1` (pixels/frame?)
   - `alignment_coefficient: 0.02` (dimensionless?)

### Current Update Loop

```rust
// boid.rs:446-487
pub fn update_location(&mut self, run_options: &RunOptions) {
    // 1. Integrate acceleration
    self.velocity += self.acceleration;

    // 2. Clamp to max speed
    self.velocity = self.velocity.limit_length_sq(max_speed_sq, max_speed);

    // 3. Enforce min speed (with NaN protection)
    if self.velocity.length_squared() < min_speed_sq {
        self.velocity = self.velocity.normalize() * min_speed;
    }

    // 4. Update position (FRAME-DEPENDENT!)
    self.position += self.velocity * baseline_speed;  // ← Problem!

    // 5. Reset acceleration
    self.acceleration *= 0.0;
}
```

**Key Issue**: No `dt` multiplication means physics is coupled to frame rate.

---

## Proposed Solution: Three-Layer Unit System

### Layer 1: Physical World (SI Units)
- Distance: meters (m)
- Time: seconds (s)
- Velocity: meters per second (m/s)
- Acceleration: meters per second squared (m/s²)

### Layer 2: Simulation Space (Internal)
- Uses normalized units for computation
- Independent of screen resolution
- Conversion factors from physical units

### Layer 3: Screen Space (Rendering)
- Pixels for display
- Scale from simulation space to window coordinates

---

## Implementation Approach

### Phase 1: Delta Time Integration

**Goal**: Decouple physics from frame rate

#### 1.1 Add Delta Time Tracking

**File**: `boids_lib/src/options.rs`

Add to `RunOptions`:
```rust
pub struct RunOptions {
    // ... existing fields ...

    // Time management
    pub delta_time: f32,           // Current frame dt in seconds
    pub target_fps: f32,            // Target simulation rate (default: 60.0)
    pub fixed_timestep: bool,       // Use fixed dt vs variable dt
}
```

#### 1.2 Extract Delta Time from Nannou

**File**: `boids_app/src/main.rs`

```rust
fn update(app: &App, model: &mut Model, update: Update) {
    // Extract delta time from Nannou's Update struct
    let dt = update.since_last.as_secs_f32();

    // Clamp dt to avoid instability (max 0.1s = 10 FPS minimum)
    let dt = dt.min(0.1);

    // Update run_options with current dt
    model.run_options.delta_time = dt;

    // ... rest of update logic ...
}
```

#### 1.3 Update Physics Integration

**File**: `boids_lib/src/boid.rs`

Replace `update_location`:
```rust
pub fn update_location(&mut self, run_options: &RunOptions) {
    let dt = if run_options.fixed_timestep {
        1.0 / run_options.target_fps
    } else {
        run_options.delta_time
    };

    // Semi-implicit Euler integration (velocity-first)
    self.velocity += self.acceleration * dt;

    // Clamp velocity
    self.velocity = self.velocity.limit_length_sq(
        run_options.max_speed_sq,
        run_options.max_speed
    );

    // Enforce minimum speed
    if self.velocity.length_squared() < run_options.min_speed_sq {
        self.velocity = self.velocity.normalize() * run_options.min_speed;
        // ... NaN protection ...
    }

    // Update position (FRAME-INDEPENDENT!)
    if !run_options.stop_movement {
        self.position += self.velocity * dt;  // ← Fixed!
    }

    // Reset acceleration
    self.acceleration *= 0.0;

    self.boundaries(run_options);
}
```

**Impact**: Remove `baseline_speed` parameter entirely.

---

### Phase 2: Spatial Scale Definition

**Goal**: Define relationship between screen pixels and physical meters

#### 2.1 Add Scale Parameters

**File**: `boids_lib/src/options.rs`

```rust
pub struct RunOptions {
    // ... existing fields ...

    // Spatial scale
    pub meters_per_pixel: f32,     // Conversion: pixels → meters (default: 0.01 = 1px = 1cm)
    pub pixels_per_meter: f32,     // Conversion: meters → pixels (default: 100)
    pub world_width_meters: f32,   // Physical width of simulation (computed)
    pub world_height_meters: f32,  // Physical height of simulation (computed)
}

impl RunOptions {
    pub fn update_spatial_scale(&mut self) {
        self.pixels_per_meter = 1.0 / self.meters_per_pixel;
        self.world_width_meters = self.window.width() as f32 * self.meters_per_pixel;
        self.world_height_meters = self.window.height() as f32 * self.meters_per_pixel;
    }

    // Conversion helpers
    pub fn pixels_to_meters(&self, pixels: f32) -> f32 {
        pixels * self.meters_per_pixel
    }

    pub fn meters_to_pixels(&self, meters: f32) -> f32 {
        meters * self.pixels_per_meter
    }
}
```

#### 2.2 Recommended Default Scale

For a typical 1200×800 window:
- `meters_per_pixel = 0.01` (1 pixel = 1 centimeter)
- Window represents: 12m × 8m space
- Typical boid size: ~5-10 cm (reasonable for birds/fish)

**User Configurable**: Could represent:
- Microscale: 1px = 1mm → window = 1.2m × 0.8m (fish tank)
- Mesoscale: 1px = 10cm → window = 120m × 80m (bird flock)
- Macroscale: 1px = 1m → window = 1200m × 800m (large swarm)

---

### Phase 3: Unit-Typed Parameters

**Goal**: Make all physical parameters explicit with SI units

#### 3.1 Create Unit System Module

**File**: `boids_lib/src/units.rs` (new file)

```rust
/// Physical units for the simulation
/// All values stored in SI base units (meters, seconds)

use std::ops::{Add, Sub, Mul, Div};
use glam::Vec2;

// ============================================================================
// Distance Units
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Meters(pub f32);

impl Meters {
    pub fn from_kilometers(km: f32) -> Self {
        Meters(km * 1000.0)
    }

    pub fn from_centimeters(cm: f32) -> Self {
        Meters(cm * 0.01)
    }

    pub fn as_meters(&self) -> f32 {
        self.0
    }

    pub fn as_kilometers(&self) -> f32 {
        self.0 / 1000.0
    }

    pub fn as_centimeters(&self) -> f32 {
        self.0 * 100.0
    }
}

// ============================================================================
// Velocity Units
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MetersPerSecond(pub f32);

impl MetersPerSecond {
    pub fn from_kmh(kmh: f32) -> Self {
        MetersPerSecond(kmh / 3.6)
    }

    pub fn as_mps(&self) -> f32 {
        self.0
    }

    pub fn as_kmh(&self) -> f32 {
        self.0 * 3.6
    }
}

impl Mul<f32> for MetersPerSecond {
    type Output = MetersPerSecond;
    fn mul(self, rhs: f32) -> Self::Output {
        MetersPerSecond(self.0 * rhs)
    }
}

impl Div<f32> for MetersPerSecond {
    type Output = MetersPerSecond;
    fn div(self, rhs: f32) -> Self::Output {
        MetersPerSecond(self.0 / rhs)
    }
}

// ============================================================================
// Acceleration Units
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MetersPerSecondSquared(pub f32);

impl MetersPerSecondSquared {
    pub fn as_mps2(&self) -> f32 {
        self.0
    }
}

impl Mul<f32> for MetersPerSecondSquared {
    type Output = MetersPerSecondSquared;
    fn mul(self, rhs: f32) -> Self::Output {
        MetersPerSecondSquared(self.0 * rhs)
    }
}

// ============================================================================
// Vector Quantities (2D)
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Position(pub Vec2);  // meters

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Velocity(pub Vec2);  // meters/second

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Acceleration(pub Vec2);  // meters/second²

// ============================================================================
// Unit Conversions
// ============================================================================

pub struct UnitConverter {
    pub meters_per_pixel: f32,
    pub pixels_per_meter: f32,
}

impl UnitConverter {
    pub fn new(meters_per_pixel: f32) -> Self {
        UnitConverter {
            meters_per_pixel,
            pixels_per_meter: 1.0 / meters_per_pixel,
        }
    }

    // Position conversions
    pub fn position_to_screen(&self, pos: Position) -> Vec2 {
        pos.0 * self.pixels_per_meter
    }

    pub fn screen_to_position(&self, screen: Vec2) -> Position {
        Position(screen * self.meters_per_pixel)
    }

    // Velocity conversions
    pub fn velocity_to_pixels_per_second(&self, vel: Velocity) -> Vec2 {
        vel.0 * self.pixels_per_meter
    }
}
```

#### 3.2 Update Boid Structure

**File**: `boids_lib/src/boid.rs`

**Option A: Full Type Safety** (Recommended for new projects)
```rust
use crate::units::{Position, Velocity, Acceleration};

pub struct Boid {
    pub id: usize,
    pub position: Position,      // SI units: meters
    pub velocity: Velocity,      // SI units: m/s
    acceleration: Acceleration,  // SI units: m/s²
}
```

**Option B: Hybrid Approach** (Easier migration)
```rust
// Keep Vec2 internally, but document units
pub struct Boid {
    pub id: usize,
    pub position: Vec2,      // Internal units: meters
    pub velocity: Vec2,      // Internal units: m/s
    acceleration: Vec2,      // Internal units: m/s²
}

// Add conversion methods
impl Boid {
    pub fn position_meters(&self) -> f32 {
        self.position.length()
    }

    pub fn velocity_mps(&self) -> f32 {
        self.velocity.length()
    }

    pub fn velocity_kmh(&self) -> f32 {
        self.velocity.length() * 3.6
    }
}
```

**Recommendation**: Start with **Option B** for easier migration.

---

### Phase 4: Update RunOptions with SI Parameters

**File**: `boids_lib/src/options.rs`

```rust
pub struct RunOptions {
    // ... existing fields ...

    // PHYSICS PARAMETERS (SI Units)
    // All speeds in m/s, all distances in meters

    // Speed limits
    pub min_speed_mps: f32,           // meters/second (e.g., 0.5 m/s)
    pub max_speed_mps: f32,           // meters/second (e.g., 5.0 m/s = 18 km/h)
    pub max_steering_mps2: f32,       // meters/second² (e.g., 2.0 m/s²)

    // Sensory ranges
    pub sensory_distance_m: f32,      // meters (e.g., 1.0 m)
    pub alignment_range_m: f32,       // meters
    pub cohesion_range_m: f32,        // meters
    pub separation_range_m: f32,      // meters

    // Force coefficients (dimensionless multipliers)
    pub alignment_coefficient: f32,
    pub cohesion_coefficient: f32,
    pub separation_coefficient: f32,

    // Wander behavior
    pub wander_strength: f32,         // m/s² max wander acceleration
    pub wander_rate: f32,             // radians/second
    pub wander_radius_m: f32,         // meters
    pub wander_distance_m: f32,       // meters

    // DEPRECATED (to be removed)
    // pub baseline_speed: f32,       // ← Remove this!
}

impl Default for RunOptions {
    fn default() -> Self {
        // Example defaults for 1px = 1cm scale
        // Window 1200×800 = 12m × 8m world

        RunOptions {
            // Time
            delta_time: 1.0 / 60.0,
            target_fps: 60.0,
            fixed_timestep: true,

            // Spatial scale
            meters_per_pixel: 0.01,  // 1px = 1cm
            pixels_per_meter: 100.0,

            // Physics (realistic bird-like movement)
            min_speed_mps: 0.5,      // 0.5 m/s = 1.8 km/h (slow glide)
            max_speed_mps: 5.0,      // 5.0 m/s = 18 km/h (fast flight)
            max_steering_mps2: 2.0,  // 2.0 m/s² max turn acceleration

            // Sensory ranges
            sensory_distance_m: 1.0, // 1 meter = 100 pixels
            alignment_range_m: 1.15,
            cohesion_range_m: 0.95,
            separation_range_m: 0.35,

            // ... rest of defaults ...
        }
    }
}
```

---

### Phase 5: Update Physics Calculations

#### 5.1 Alignment Rule

**File**: `boids_lib/src/boid.rs` (line 116)

**Before**:
```rust
pub fn alignment(&self, neighbours: &[&Boid], run_options: &RunOptions) -> Vec2 {
    let mut sum = Vec2::ZERO;
    for boid in neighbours {
        sum += boid.velocity;  // dimensionless
    }
    if neighbours.len() > 0 {
        sum /= neighbours.len() as f32;
        sum = self.steer(sum, run_options);
    }
    sum
}
```

**After** (with SI units):
```rust
pub fn alignment(&self, neighbours: &[&Boid], run_options: &RunOptions) -> Vec2 {
    let mut sum = Vec2::ZERO;
    for boid in neighbours {
        sum += boid.velocity;  // m/s
    }
    if neighbours.len() > 0 {
        sum /= neighbours.len() as f32;  // average velocity in m/s
        sum = self.steer(sum, run_options);  // returns m/s²
    }
    sum
}
```

**Key Change**: Return value is now `m/s²` (acceleration) instead of dimensionless.

#### 5.2 Cohesion Rule

**File**: `boids_lib/src/boid.rs` (line 144)

**Before**:
```rust
pub fn cohesion(&self, neighbours: &[&Boid], run_options: &RunOptions) -> Vec2 {
    let mut sum = Vec2::ZERO;
    for boid in neighbours {
        sum += boid.position;  // dimensionless
    }
    if neighbours.len() > 0 {
        sum /= neighbours.len() as f32;
        sum -= self.position;
        sum = self.steer(sum, run_options);
    }
    sum
}
```

**After**:
```rust
pub fn cohesion(&self, neighbours: &[&Boid], run_options: &RunOptions) -> Vec2 {
    let mut sum = Vec2::ZERO;
    for boid in neighbours {
        sum += boid.position;  // meters
    }
    if neighbours.len() > 0 {
        sum /= neighbours.len() as f32;     // center of mass in meters
        sum -= self.position;                // direction vector in meters
        sum = self.steer(sum, run_options);  // returns m/s²
    }
    sum
}
```

#### 5.3 Separation Rule

**File**: `boids_lib/src/boid.rs` (line 172)

**Before**:
```rust
pub fn separation(&self, neighbours: &[&Boid], run_options: &RunOptions) -> Vec2 {
    let mut sum = Vec2::ZERO;
    for boid in neighbours {
        let d = distance_dyn(...);
        if d > 0.0 {
            let mut diff = self.position - boid.position;
            diff = diff.normalize() / d;  // inverse distance weighting
            sum += diff;
        }
    }
    // ...
}
```

**After**:
```rust
pub fn separation(&self, neighbours: &[&Boid], run_options: &RunOptions) -> Vec2 {
    let mut sum = Vec2::ZERO;
    for boid in neighbours {
        let d = distance_dyn(...);  // distance in meters
        if d > 0.0 {
            let mut diff = self.position - boid.position;  // meters
            diff = diff.normalize() / d;  // 1/meters (force magnitude)
            sum += diff;
        }
    }
    // ...
    // Return value is m/s² (acceleration)
}
```

#### 5.4 Steering Function

**File**: `boids_lib/src/boid.rs` (line 322)

**Before**:
```rust
pub fn steer(&self, mut desired: Vec2, run_options: &RunOptions) -> Vec2 {
    desired = desired.normalize();
    if run_options.agent_steering {
        desired *= run_options.max_speed;      // dimensionless
        desired -= self.velocity;              // dimensionless
        desired = desired.limit_length_sq(
            run_options.max_steering_sq,
            run_options.max_steering           // dimensionless
        );
    }
    desired
}
```

**After**:
```rust
pub fn steer(&self, mut desired: Vec2, run_options: &RunOptions) -> Vec2 {
    desired = desired.normalize();
    if run_options.agent_steering {
        desired *= run_options.max_speed_mps;  // m/s (target velocity)
        desired -= self.velocity;               // m/s (velocity error)

        // Convert to acceleration (divide by assumed response time)
        // Using dt as response time: acceleration = velocity_change / dt
        let max_accel = run_options.max_steering_mps2;
        desired = desired.limit_length(max_accel);  // m/s²
    }
    desired  // returns m/s²
}
```

**Important**: The steering function now returns **acceleration** (m/s²) instead of velocity change.

---

### Phase 6: Update Rendering

**File**: `boids_app/src/main.rs`

```rust
fn view(app: &App, model: &Model, frame: Frame) {
    let draw = app.draw();

    let converter = UnitConverter::new(model.run_options.meters_per_pixel);

    for boid in &model.flock.boids {
        // Convert from meters to screen pixels
        let screen_pos = converter.position_to_screen(boid.position);

        // Draw boid at screen position
        draw.ellipse()
            .xy(screen_pos)
            .radius(model.run_options.size)
            .color(model.color);
    }

    draw.to_frame(app, &frame).unwrap();
}
```

---

### Phase 7: Configuration Migration

#### 7.1 Update Default Configs

**File**: `boids_app/configs/*.toml`

**Old Format**:
```toml
baseline_speed = 1.0
max_speed = 4.1
min_speed = 0.65
sensory_distance = 60.0
```

**New Format**:
```toml
# Spatial scale
meters_per_pixel = 0.01  # 1px = 1cm, window = 12m × 8m

# Time
fixed_timestep = true
target_fps = 60.0

# Physics (SI units)
min_speed_mps = 0.5      # 0.5 m/s = 1.8 km/h
max_speed_mps = 5.0      # 5.0 m/s = 18 km/h
max_steering_mps2 = 2.0  # 2.0 m/s²
sensory_distance_m = 1.0 # 1 meter
```

#### 7.2 Conversion Tool

Create a migration utility:

**File**: `boids_lib/src/config_migration.rs` (new file)

```rust
/// Converts old dimensionless configs to SI unit configs
pub fn migrate_old_config(old: &OldRunOptions) -> RunOptions {
    // Assume old config used 60 FPS implicit timing
    let old_dt = 1.0 / 60.0;

    // Assume 1 pixel = 1 cm scale
    let meters_per_pixel = 0.01;

    RunOptions {
        // Convert speeds: (pixels/frame) → (m/s)
        // old_speed_px_per_frame * fps = old_speed_px_per_sec
        // old_speed_px_per_sec * m_per_px = new_speed_m_per_sec
        min_speed_mps: old.min_speed / old_dt * meters_per_pixel,
        max_speed_mps: old.max_speed / old_dt * meters_per_pixel,

        // Convert distances: pixels → meters
        sensory_distance_m: old.sensory_distance * meters_per_pixel,

        // ... rest of conversions ...
    }
}
```

---

## Example Scenarios with Real Units

### Scenario 1: Bird Flock (Mesoscale)

```rust
RunOptions {
    // Scale: 1200×800 window = 120m × 80m
    meters_per_pixel: 0.1,  // 1px = 10cm

    // Bird speeds (starling-like)
    min_speed_mps: 5.0,     // 18 km/h (slow glide)
    max_speed_mps: 20.0,    // 72 km/h (fast flight)
    max_steering_mps2: 5.0, // 5 m/s² (agile turning)

    // Interaction ranges
    sensory_distance_m: 10.0,    // 10 meter awareness (100px)
    separation_range_m: 2.0,     // 2 meter personal space

    // Results in realistic flocking behavior
}
```

### Scenario 2: Fish School (Microscale)

```rust
RunOptions {
    // Scale: 1200×800 window = 1.2m × 0.8m (aquarium)
    meters_per_pixel: 0.001,  // 1px = 1mm

    // Fish speeds (small tropical fish)
    min_speed_mps: 0.05,      // 5 cm/s (slow swim)
    max_speed_mps: 0.3,       // 30 cm/s (burst speed)
    max_steering_mps2: 0.5,   // 0.5 m/s² (quick turn)

    // Interaction ranges
    sensory_distance_m: 0.1,  // 10 cm awareness (100px)
    separation_range_m: 0.02, // 2 cm personal space
}
```

### Scenario 3: Insect Swarm (Microscale, Fast)

```rust
RunOptions {
    // Scale: 1200×800 window = 12m × 8m
    meters_per_pixel: 0.01,  // 1px = 1cm

    // Insect speeds (flies/mosquitoes)
    min_speed_mps: 0.5,       // 0.5 m/s
    max_speed_mps: 2.0,       // 2 m/s = 7.2 km/h
    max_steering_mps2: 10.0,  // 10 m/s² (very agile)

    // Interaction ranges
    sensory_distance_m: 0.5,  // 50 cm awareness
    separation_range_m: 0.1,  // 10 cm personal space
}
```

---

## Implementation Roadmap

### Phase 1: Foundation (Week 1)
- [x] Analyze current codebase
- [ ] Add delta time tracking
- [ ] Update `update_location` with dt
- [ ] Remove `baseline_speed` dependency
- [ ] Test: Verify frame-rate independence

### Phase 2: Spatial Scale (Week 1-2)
- [ ] Add `meters_per_pixel` to `RunOptions`
- [ ] Implement conversion helpers
- [ ] Update rendering with conversions
- [ ] Test: Verify visual output unchanged

### Phase 3: Units Module (Week 2)
- [ ] Create `boids_lib/src/units.rs`
- [ ] Implement typed units (Option B: hybrid approach)
- [ ] Add conversion utilities
- [ ] Document unit conventions

### Phase 4: Parameter Migration (Week 2-3)
- [ ] Add SI unit parameters to `RunOptions`
- [ ] Update default values
- [ ] Create config migration tool
- [ ] Deprecate old parameters

### Phase 5: Physics Updates (Week 3)
- [ ] Update `alignment` rule
- [ ] Update `cohesion` rule
- [ ] Update `separation` rule
- [ ] Update `steer` function
- [ ] Update `wander` behavior

### Phase 6: Testing & Validation (Week 3-4)
- [ ] Unit tests for conversions
- [ ] Integration tests for physics
- [ ] Visual regression tests
- [ ] Performance benchmarks

### Phase 7: Documentation (Week 4)
- [ ] Update README with SI units
- [ ] API documentation
- [ ] Migration guide
- [ ] Example configs

### Phase 8: UI/UX (Week 4)
- [ ] Add UI controls for m/s, km/h display
- [ ] Show real-time velocity/acceleration
- [ ] Add scale visualization
- [ ] Help tooltips with units

---

## Testing Strategy

### Unit Tests

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_velocity_conversion() {
        let vel = MetersPerSecond::from_kmh(36.0);  // 36 km/h
        assert_eq!(vel.as_mps(), 10.0);              // = 10 m/s
    }

    #[test]
    fn test_frame_rate_independence() {
        let mut boid1 = create_test_boid();
        let mut boid2 = boid1.clone();

        // Simulate at 60 FPS
        for _ in 0..60 {
            boid1.update_location(&options_60fps);
        }

        // Simulate at 30 FPS (half speed, same duration)
        for _ in 0..30 {
            boid2.update_location(&options_30fps);
        }

        // Positions should be approximately equal
        assert!((boid1.position - boid2.position).length() < 0.01);
    }

    #[test]
    fn test_physical_limits() {
        let mut boid = create_test_boid();
        boid.velocity = Vec2::new(100.0, 0.0);  // 100 m/s (unrealistic)

        boid.update_location(&run_options);

        // Should be clamped to max_speed_mps
        assert!(boid.velocity.length() <= run_options.max_speed_mps);
    }
}
```

### Integration Tests

1. **Conservation of Energy**: Track total kinetic energy over time
2. **Behavior Invariance**: Same flocking patterns at different scales
3. **Boundary Conditions**: Toroidal/reflective boundaries work correctly
4. **Performance**: No significant slowdown with SI units

---

## Migration Path for Existing Users

### Backwards Compatibility

**Option A**: Auto-detect old configs and migrate
```rust
if config.baseline_speed.is_some() {
    // Old format detected, migrate
    let new_config = migrate_old_config(&config);
    warn!("Old config format detected, migrating to SI units");
    new_config
} else {
    config  // Already in new format
}
```

**Option B**: Provide migration CLI tool
```bash
cargo run --bin migrate-config -- old_config.toml > new_config.toml
```

**Option C**: Support both (with deprecation warning)
```rust
pub struct RunOptions {
    #[deprecated(note = "Use min_speed_mps instead")]
    pub min_speed: Option<f32>,

    pub min_speed_mps: f32,
}
```

---

## Performance Considerations

### Computational Cost

**Added Operations per Frame:**
- Delta time extraction: 1 system call
- Unit conversions: ~10-20 multiplications per boid
- No significant change to core algorithms

**Expected Impact**: < 1% performance overhead

### Optimization Opportunities

1. **Precompute conversions**: Store both pixel and meter values
2. **SIMD**: Use `glam` SIMD features for vector ops
3. **Lookup tables**: For common conversions (if needed)

---

## Benefits Summary

### For Users

✅ **Intuitive parameters**: "max speed = 20 m/s" vs "max speed = 4.1"
✅ **Real-world scales**: Model actual birds, fish, insects
✅ **Frame-rate independence**: 60 FPS = 30 FPS = 144 FPS (same physics)
✅ **Scientific applications**: Export data with meaningful units
✅ **Educational value**: Learn physics concepts

### For Developers

✅ **Dimensional analysis**: Catch unit errors at compile time (with typed units)
✅ **Testability**: Verify physical correctness
✅ **Maintainability**: Clear parameter meanings
✅ **Extensibility**: Easy to add new forces/behaviors

---

## Risks & Mitigation

| Risk | Impact | Mitigation |
|------|--------|------------|
| Breaking changes | High | Provide migration tool + deprecation period |
| Performance regression | Medium | Benchmark early, optimize conversions |
| User confusion | Medium | Clear documentation, examples |
| Bug introduction | High | Extensive testing, gradual rollout |

---

## Alternative Approaches Considered

### 1. Use `uom` Crate (Compile-Time Units)

**Pros:**
- Full type safety
- Zero-cost abstractions
- Automatic conversions

**Cons:**
- Heavy dependency
- Complex type signatures
- Steeper learning curve

**Decision**: Rejected for simplicity. Could be added later.

### 2. Keep Current System, Add Documentation

**Pros:**
- No code changes
- No risk

**Cons:**
- Doesn't solve frame-rate dependency
- Still confusing for users

**Decision**: Rejected. Doesn't address core issues.

### 3. Implicit Units (Current + dt)

**Pros:**
- Minimal changes
- Solves frame-rate issue

**Cons:**
- Still dimensionless
- No real-world scale

**Decision**: Partial solution. Adopted as Phase 1.

---

## Conclusion

This plan provides a **structured, incremental approach** to introducing SI units to the boids simulation. The hybrid implementation (keeping `Vec2` but documenting units) balances:

- **Ease of migration**: Minimal code disruption
- **User benefit**: Real-world parameters
- **Correctness**: Frame-rate independence
- **Future-proofing**: Path to full type safety

**Recommended Start**: Implement Phases 1-2 first (delta time + spatial scale), then evaluate user feedback before proceeding to full SI unit parameters.

**Timeline**: ~4 weeks for full implementation, or ~1 week for minimal viable change (Phases 1-2 only).
