
use crate::{
    boid::Boid, options::{WindowSize, RunOptions, Distance},
    options::Boundary
};

pub fn distance_dyn_boid(b1: &Boid, b2: &Boid, run_options: &RunOptions) -> f32 {
    match run_options.distance {
        Distance::EucEnclosed => simple_distance_boid(b1, b2, &run_options.window),
        Distance::EucThoroidal => thoroidal_distance_boid(b1, b2, &run_options.window),
    }
}

pub fn distance_dyn(x1: f32, x2: f32, y1: f32, y2: f32, run_options: &RunOptions) -> f32 {
    match run_options.distance {
        Distance::EucEnclosed => simple_distance(x1, x2, y1, y2, &run_options.window),
        Distance::EucThoroidal => thoroidal_distance(x1, x2, y1, y2, &run_options.window),
    }
}

#[allow(dead_code)]
pub fn simple_distance_boid(b1: &Boid, b2: &Boid, _window_size: &WindowSize) -> f32 {
    simple_distance_sq_boid(b1, b2, _window_size).sqrt()
}

#[allow(dead_code)]
pub fn simple_distance_sq_boid(b1: &Boid, b2: &Boid, _window_size: &WindowSize) -> f32 {
    simple_distance_sq(b1.position.x, b2.position.x, b1.position.y, b2.position.y, _window_size)
}

#[allow(dead_code)]
pub fn simple_distance(x1: f32, x2: f32, y1: f32, y2: f32, _window_size: &WindowSize) -> f32 {
    simple_distance_sq(x1, x2, y1, y2, _window_size).sqrt()
}

pub fn simple_distance_sq(x1: f32, x2: f32, y1: f32, y2: f32, _window_size: &WindowSize) -> f32 {
    (x1 - x2).powi(2) + (y1 - y2).powi(2)
}

pub fn thoroidal_distance_boid(b1: &Boid, b2: &Boid, window_size: &WindowSize) -> f32 {
    thoroidal_distance_sq_boid(b1, b2, window_size).sqrt()
}

pub fn thoroidal_distance_sq_boid(b1: &Boid, b2: &Boid, window_size: &WindowSize) -> f32 {
    thoroidal_distance_sq(b1.position.x, b2.position.x, b1.position.y, b2.position.y, window_size)
    // thoroidal_distance(b1.pos_x, b2.pos_x, b1.pos_y, b2.pos_y)
}

#[allow(dead_code)]
pub fn thoroidal_distance(x1: f32, x2: f32, y1: f32, y2: f32, window_size: &WindowSize) -> f32 {
    thoroidal_distance_sq(x1, x2, y1, y2, window_size).sqrt()
}

// Thoroidal distance as per https://blog.demofox.org/2017/10/01/calculating-the-distance-between-points-in-wrap-around-toroidal-space/
/// Get distance in between points in a thoroidal space "continuous", that "wraps around"
pub fn thoroidal_distance_sq(x1: f32, x2: f32, y1: f32, y2: f32, window_size: &WindowSize) -> f32 {
    let mut dx = (x1 - x2).abs();
    let mut dy = (y1 - y2).abs();

    // coordinates are centered w.r.t. (0, 0)
    // that means that the highest x and y values are exactly half of width and height respectively

    // e.g., if on dimension x we got an absolute difference, dx = |x1 - x2| larger than |x|/2
    // we went the wrong direction and should consider the other side to get the minimal distance

    if dx > window_size.win_right {
        // dx = 1.0 - dx;
        dx = window_size.win_w - dx;
    }

    if dy > window_size.win_top {
        // dy = 1.0 - dy;
        dy = window_size.win_h - dy;
    }

    return dx.powi(2) + dy.powi(2);
}

// TODO: test

#[cfg(test)]
mod tests {
    use approx::{assert_relative_eq};
    use nannou::prelude::{vec2, Vec2};

    use crate::options::{WindowSize};

    use super::thoroidal_distance;

    macro_rules! assert_eqf32 {
        ($x:expr, $y:expr) => {
            assert_relative_eq!($x, $y, epsilon = 1e-3f32)
        };
    }

    struct ThoroidialTestData {
        a: Vec2,
        b: Vec2,
        c: Vec2,
        d: Vec2,
    }

    impl ThoroidialTestData {
        fn new() -> Self {
            ThoroidialTestData {
                 a: vec2(3., 2.),
                 b: vec2(0., 0.),
                 c: vec2(-3., 3.),
                 d: vec2(-3., -3.5)
            }
        }
    }

    fn _thoroidal_distance(v1: Vec2, v2: Vec2) -> f32 {
        let window_size = WindowSize{
            win_left: -4.,
            win_right: 4.,
            win_top: 4.,
            win_bottom: -4.,
            win_h: 8.,
            win_w: 8.,
        };

        thoroidal_distance(v1.x, v2.x, v1.y, v2.y, &window_size)
    }

    #[test]
    fn thoroidal_distance_a_b() {
        let data = ThoroidialTestData::new();
        let result = _thoroidal_distance(data.a, data.b);
        assert_eqf32!(result, 3.606);
    }

    #[test]
    fn thoroidal_distance_a_c() {
        let data = ThoroidialTestData::new();
        let result = _thoroidal_distance(data.a, data.c);
        assert_eqf32!(result, 2.236);
    }

    #[test]
    fn thoroidal_distance_a_d() {
        let data = ThoroidialTestData::new();
        let result = _thoroidal_distance(data.a, data.d);
        assert_eqf32!(result, 3.202);
    }

    #[test]
    fn thoroidal_distance_c_d() {
        let data = ThoroidialTestData::new();
        let result = _thoroidal_distance(data.c, data.d);
        assert_eqf32!(result, 1.5);
    }
}