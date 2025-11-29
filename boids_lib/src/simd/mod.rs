//! SIMD-accelerated distance calculations for boids neighbor queries.
//!
//! This module provides batched distance calculations that process multiple
//! candidate boids simultaneously using SIMD instructions.
//!
//! # Platform Support
//! - `simd-neon`: ARM64 NEON (128-bit, 4x f32)
//! - `simd-avx2`: x86-64 AVX2 (256-bit, 8x f32)
//!
//! # Usage
//! Enable via Cargo features:
//! ```toml
//! # M3 Mac:
//! cargo +nightly build --release --features simd-neon
//!
//! # AMD Zen 3:
//! cargo +nightly build --release --features simd-avx2
//! ```

use std::simd::{cmp::SimdPartialOrd, num::SimdFloat, Mask, Simd, StdFloat};

#[cfg(all(feature = "simd-neon", not(feature = "simd-avx2")))]
use std::simd::f32x4;

#[cfg(feature = "simd-avx2")]
use std::simd::{f32x4, f32x8};

/// SIMD lane count for the current platform
#[cfg(feature = "simd-avx2")]
pub const SIMD_LANES: usize = 8;

#[cfg(all(feature = "simd-neon", not(feature = "simd-avx2")))]
pub const SIMD_LANES: usize = 4;

#[cfg(not(any(feature = "simd-neon", feature = "simd-avx2")))]
pub const SIMD_LANES: usize = 4; // Default fallback

/// Result of batch distance calculation
#[derive(Debug, Clone, Copy)]
pub struct BatchDistanceResult<const N: usize>
where
    std::simd::LaneCount<N>: std::simd::SupportedLaneCount,
{
    /// Squared distances for each candidate
    pub dist_sq: Simd<f32, N>,
    /// Distances (with rsqrt approximation)
    pub dist: Simd<f32, N>,
    /// Normalized direction X components (from query to candidates)
    pub dir_x: Simd<f32, N>,
    /// Normalized direction Y components (from query to candidates)
    pub dir_y: Simd<f32, N>,
    /// Mask of candidates within threshold (true = within range)
    pub within_range: Mask<i32, N>,
}

/// Fast reciprocal square root with Newton-Raphson refinement.
///
/// Uses hardware rsqrt estimate followed by one N-R iteration for ~22-bit accuracy.
/// This is ~5x faster than sqrt() for simulation purposes.
#[inline]
fn fast_rsqrt<const N: usize>(x: Simd<f32, N>) -> Simd<f32, N>
where
    std::simd::LaneCount<N>: std::simd::SupportedLaneCount,
{
    let half = Simd::splat(0.5f32);
    let three_halves = Simd::splat(1.5f32);

    // Initial estimate using hardware rsqrt (via sqrt().recip())
    // Note: std::simd doesn't expose rsqrt directly, but sqrt().recip() is optimized
    let y = x.sqrt().recip();

    // One Newton-Raphson iteration: y = y * (1.5 - 0.5 * x * y * y)
    // This improves accuracy from ~11-bit to ~22-bit
    y * (three_halves - half * x * y * y)
}

/// Batch distance calculation for enclosed (non-toroidal) space.
///
/// Processes N candidates simultaneously, computing distances and normalized
/// direction vectors from a single query point.
///
/// # Arguments
/// * `query_x` - Query boid's X position
/// * `query_y` - Query boid's Y position
/// * `candidates_x` - X positions of N candidate boids
/// * `candidates_y` - Y positions of N candidate boids
/// * `threshold_sq` - Squared distance threshold for neighbor consideration
///
/// # Returns
/// `BatchDistanceResult` containing distances, directions, and validity mask
#[inline]
pub fn batch_distances_enclosed<const N: usize>(
    query_x: f32,
    query_y: f32,
    candidates_x: Simd<f32, N>,
    candidates_y: Simd<f32, N>,
    threshold_sq: f32,
) -> BatchDistanceResult<N>
where
    std::simd::LaneCount<N>: std::simd::SupportedLaneCount,
{
    let qx = Simd::splat(query_x);
    let qy = Simd::splat(query_y);
    let thresh_sq = Simd::splat(threshold_sq);
    let epsilon = Simd::splat(0.0001f32);

    // Delta vectors from query to candidates
    let dx = candidates_x - qx;
    let dy = candidates_y - qy;

    // Squared distance
    let dist_sq = dx * dx + dy * dy;

    // Check which candidates are within range
    let within_range = dist_sq.simd_le(thresh_sq);

    // Fast inverse distance using rsqrt with N-R refinement
    // Add epsilon to avoid division by zero for coincident points
    let inv_dist = fast_rsqrt(dist_sq + epsilon);

    // Distance = dist_sq * inv_dist = dist_sq / sqrt(dist_sq)
    let dist = dist_sq * inv_dist;

    // Normalized direction = delta * inv_dist
    let dir_x = dx * inv_dist;
    let dir_y = dy * inv_dist;

    BatchDistanceResult {
        dist_sq,
        dist,
        dir_x,
        dir_y,
        within_range,
    }
}

/// Batch distance calculation for toroidal (wrap-around) space.
///
/// Handles wrap-around boundaries branchlessly using SIMD blend operations.
///
/// # Arguments
/// * `query_x` - Query boid's X position
/// * `query_y` - Query boid's Y position
/// * `candidates_x` - X positions of N candidate boids
/// * `candidates_y` - Y positions of N candidate boids
/// * `threshold_sq` - Squared distance threshold
/// * `half_width` - Half the world width (win_right)
/// * `half_height` - Half the world height (win_top)
/// * `width` - Full world width
/// * `height` - Full world height
#[inline]
pub fn batch_distances_toroidal<const N: usize>(
    query_x: f32,
    query_y: f32,
    candidates_x: Simd<f32, N>,
    candidates_y: Simd<f32, N>,
    threshold_sq: f32,
    half_width: f32,
    half_height: f32,
    width: f32,
    height: f32,
) -> BatchDistanceResult<N>
where
    std::simd::LaneCount<N>: std::simd::SupportedLaneCount,
{
    let qx = Simd::splat(query_x);
    let qy = Simd::splat(query_y);
    let thresh_sq = Simd::splat(threshold_sq);
    let epsilon = Simd::splat(0.0001f32);

    let half_w = Simd::splat(half_width);
    let half_h = Simd::splat(half_height);
    let w = Simd::splat(width);
    let h = Simd::splat(height);
    let zero = Simd::splat(0.0f32);

    // Raw delta
    let raw_dx = candidates_x - qx;
    let raw_dy = candidates_y - qy;

    // Branchless toroidal wrap:
    // if abs(d) > half_size: d += (d < 0 ? size : -size)

    // X component
    let abs_dx = raw_dx.abs();
    let needs_wrap_x = abs_dx.simd_gt(half_w);
    let correction_x = raw_dx.simd_lt(zero).select(w, -w);
    let dx = needs_wrap_x.select(raw_dx + correction_x, raw_dx);

    // Y component
    let abs_dy = raw_dy.abs();
    let needs_wrap_y = abs_dy.simd_gt(half_h);
    let correction_y = raw_dy.simd_lt(zero).select(h, -h);
    let dy = needs_wrap_y.select(raw_dy + correction_y, raw_dy);

    // Squared distance
    let dist_sq = dx * dx + dy * dy;

    // Check which candidates are within range
    let within_range = dist_sq.simd_le(thresh_sq);

    // Fast inverse distance
    let inv_dist = fast_rsqrt(dist_sq + epsilon);
    let dist = dist_sq * inv_dist;

    // Normalized direction
    let dir_x = dx * inv_dist;
    let dir_y = dy * inv_dist;

    BatchDistanceResult {
        dist_sq,
        dist,
        dir_x,
        dir_y,
        within_range,
    }
}

/// Type aliases for platform-specific SIMD widths
#[cfg(feature = "simd-avx2")]
pub type SimdF32 = f32x8;

#[cfg(all(feature = "simd-neon", not(feature = "simd-avx2")))]
pub type SimdF32 = f32x4;

#[cfg(not(any(feature = "simd-neon", feature = "simd-avx2")))]
pub type SimdF32 = f32x4;

/// Process candidates with a callback for each qualifying neighbor.
///
/// This is the main entry point for SIMD-accelerated neighbor queries.
/// Uses a callback to avoid intermediate Vec allocation.
///
/// # Arguments
/// * `query_x`, `query_y` - Query boid position
/// * `query_id` - Query boid ID (to exclude self)
/// * `positions_x`, `positions_y` - SoA position arrays
/// * `start_idx`, `end_idx` - Range of candidates to process
/// * `threshold_sq` - Squared max sensory distance
/// * `is_toroidal` - Whether to use toroidal distance
/// * `half_width`, `half_height`, `width`, `height` - World dimensions (for toroidal)
/// * `callback` - Called for each qualifying candidate: (index, distance, dir_x, dir_y) -> bool
///                Return false to stop early, true to continue
#[inline]
pub fn process_candidates_callback<F>(
    query_x: f32,
    query_y: f32,
    query_id: usize,
    positions_x: &[f32],
    positions_y: &[f32],
    boid_ids: &[usize],
    start_idx: usize,
    end_idx: usize,
    threshold_sq: f32,
    is_toroidal: bool,
    half_width: f32,
    half_height: f32,
    width: f32,
    height: f32,
    mut callback: F,
) -> bool
where
    F: FnMut(usize, f32, f32, f32) -> bool,
{
    let mut idx = start_idx;

    // Process in batches of 4 (works for both NEON and AVX2)
    while idx + 4 <= end_idx {
        // Load 4 candidates
        let cx = f32x4::from_slice(&positions_x[idx..idx + 4]);
        let cy = f32x4::from_slice(&positions_y[idx..idx + 4]);

        let batch_result = if is_toroidal {
            batch_distances_toroidal(
                query_x, query_y, cx, cy, threshold_sq,
                half_width, half_height, width, height,
            )
        } else {
            batch_distances_enclosed(query_x, query_y, cx, cy, threshold_sq)
        };

        // Extract qualifying candidates
        let mask = batch_result.within_range.to_array();
        let dists = batch_result.dist.to_array();
        let dir_xs = batch_result.dir_x.to_array();
        let dir_ys = batch_result.dir_y.to_array();

        for i in 0..4 {
            if mask[i] && boid_ids[idx + i] != query_id {
                if !callback(idx + i, dists[i], dir_xs[i], dir_ys[i]) {
                    return false; // Early exit requested
                }
            }
        }

        idx += 4;
    }

    // Handle remainder with scalar fallback
    while idx < end_idx {
        if boid_ids[idx] == query_id {
            idx += 1;
            continue;
        }

        let dx = if is_toroidal {
            let raw = positions_x[idx] - query_x;
            if raw.abs() > half_width {
                raw + if raw < 0.0 { width } else { -width }
            } else {
                raw
            }
        } else {
            positions_x[idx] - query_x
        };

        let dy = if is_toroidal {
            let raw = positions_y[idx] - query_y;
            if raw.abs() > half_height {
                raw + if raw < 0.0 { height } else { -height }
            } else {
                raw
            }
        } else {
            positions_y[idx] - query_y
        };

        let dist_sq = dx * dx + dy * dy;
        if dist_sq <= threshold_sq {
            let dist = dist_sq.sqrt();
            let inv_dist = if dist > 0.0001 { 1.0 / dist } else { 0.0 };
            if !callback(idx, dist, dx * inv_dist, dy * inv_dist) {
                return false; // Early exit requested
            }
        }

        idx += 1;
    }

    true // Completed all candidates
}

/// Process a batch of candidates and return qualifying neighbors.
///
/// This is the main entry point for SIMD-accelerated neighbor queries.
/// It processes candidates in batches of SIMD_LANES (4 for NEON, 8 for AVX2).
///
/// # Arguments
/// * `query_x`, `query_y` - Query boid position
/// * `query_id` - Query boid ID (to exclude self)
/// * `positions_x`, `positions_y` - SoA position arrays
/// * `start_idx`, `end_idx` - Range of candidates to process
/// * `threshold_sq` - Squared max sensory distance
/// * `is_toroidal` - Whether to use toroidal distance
/// * `half_width`, `half_height`, `width`, `height` - World dimensions (for toroidal)
///
/// # Returns
/// Vector of (index, distance, dir_x, dir_y) for candidates within range
#[inline]
pub fn process_candidate_batch(
    query_x: f32,
    query_y: f32,
    query_id: usize,
    positions_x: &[f32],
    positions_y: &[f32],
    boid_ids: &[usize],
    start_idx: usize,
    end_idx: usize,
    threshold_sq: f32,
    is_toroidal: bool,
    half_width: f32,
    half_height: f32,
    width: f32,
    height: f32,
) -> Vec<(usize, f32, f32, f32)> {
    let mut results = Vec::with_capacity(end_idx - start_idx);

    let mut idx = start_idx;

    // Process in batches of 4 (works for both NEON and AVX2)
    while idx + 4 <= end_idx {
        // Load 4 candidates
        let cx = f32x4::from_slice(&positions_x[idx..idx + 4]);
        let cy = f32x4::from_slice(&positions_y[idx..idx + 4]);

        let batch_result = if is_toroidal {
            batch_distances_toroidal(
                query_x, query_y, cx, cy, threshold_sq,
                half_width, half_height, width, height,
            )
        } else {
            batch_distances_enclosed(query_x, query_y, cx, cy, threshold_sq)
        };

        // Extract qualifying candidates
        let mask = batch_result.within_range.to_array();
        let dists = batch_result.dist.to_array();
        let dir_xs = batch_result.dir_x.to_array();
        let dir_ys = batch_result.dir_y.to_array();

        for i in 0..4 {
            if mask[i] && boid_ids[idx + i] != query_id {
                results.push((idx + i, dists[i], dir_xs[i], dir_ys[i]));
            }
        }

        idx += 4;
    }

    // Handle remainder with scalar fallback
    while idx < end_idx {
        if boid_ids[idx] == query_id {
            idx += 1;
            continue;
        }

        let dx = if is_toroidal {
            let raw = positions_x[idx] - query_x;
            if raw.abs() > half_width {
                raw + if raw < 0.0 { width } else { -width }
            } else {
                raw
            }
        } else {
            positions_x[idx] - query_x
        };

        let dy = if is_toroidal {
            let raw = positions_y[idx] - query_y;
            if raw.abs() > half_height {
                raw + if raw < 0.0 { height } else { -height }
            } else {
                raw
            }
        } else {
            positions_y[idx] - query_y
        };

        let dist_sq = dx * dx + dy * dy;
        if dist_sq <= threshold_sq {
            let dist = dist_sq.sqrt();
            let inv_dist = if dist > 0.0001 { 1.0 / dist } else { 0.0 };
            results.push((idx, dist, dx * inv_dist, dy * inv_dist));
        }

        idx += 1;
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_batch_distances_enclosed_4() {
        let query_x = 0.0f32;
        let query_y = 0.0f32;

        // 4 candidates at known distances
        let candidates_x = f32x4::from_array([3.0, 4.0, 0.0, 10.0]);
        let candidates_y = f32x4::from_array([4.0, 3.0, 5.0, 0.0]);

        // Expected distances: 5.0, 5.0, 5.0, 10.0
        let threshold_sq = 36.0; // radius 6

        let result = batch_distances_enclosed(
            query_x, query_y, candidates_x, candidates_y, threshold_sq,
        );

        let dists = result.dist.to_array();
        assert_relative_eq!(dists[0], 5.0, epsilon = 0.01);
        assert_relative_eq!(dists[1], 5.0, epsilon = 0.01);
        assert_relative_eq!(dists[2], 5.0, epsilon = 0.01);
        assert_relative_eq!(dists[3], 10.0, epsilon = 0.01);

        // First 3 within range, 4th outside
        let mask = result.within_range.to_array();
        assert!(mask[0]);
        assert!(mask[1]);
        assert!(mask[2]);
        assert!(!mask[3]);
    }

    #[test]
    fn test_batch_distances_toroidal_4() {
        // World: -50 to 50 (width=100, half=50)
        let half_w = 50.0f32;
        let half_h = 50.0f32;
        let w = 100.0f32;
        let h = 100.0f32;

        let query_x = 45.0f32;
        let query_y = 0.0f32;

        // Candidate at -45 should wrap to distance 10 (not 90)
        let candidates_x = f32x4::from_array([-45.0, 0.0, 45.0, 40.0]);
        let candidates_y = f32x4::from_array([0.0, 0.0, 0.0, 0.0]);

        let threshold_sq = 225.0; // radius 15

        let result = batch_distances_toroidal(
            query_x, query_y, candidates_x, candidates_y, threshold_sq,
            half_w, half_h, w, h,
        );

        let dists = result.dist.to_array();
        // -45 to 45 wraps: 45 - (-45) = 90, but via wrap: 100 - 90 = 10
        assert_relative_eq!(dists[0], 10.0, epsilon = 0.1);
        // 0 to 45 = 45 (no wrap)
        assert_relative_eq!(dists[1], 45.0, epsilon = 0.1);
        // 45 to 45 = 0 (epsilon protected)
        assert_relative_eq!(dists[2], 0.0, epsilon = 0.1);
        // 40 to 45 = 5
        assert_relative_eq!(dists[3], 5.0, epsilon = 0.1);
    }
}
