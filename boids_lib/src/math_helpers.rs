
use glam::Vec2;

use crate::{
    boid::Boid, options::{WindowSize, RunOptions, Distance}
};

// pub fn position_dyn_boid(b1: &Boid, b2: &Boid, run_options: &RunOptions) -> Vec2 {
//     match run_options.distance {
//         Distance::EucToroidal => tor_vec(b1.position, b2.position, &run_options.window),
//         Distance::EucEnclosed => ,
//     }
// }

pub fn distance_dyn_boid(b1: &Boid, b2: &Boid, run_options: &RunOptions) -> f32 {
    match run_options.distance {
        Distance::EucEnclosed => simple_distance_boid(b1, b2, &run_options.window),
        Distance::EucToroidal => thoroidal_distance_boid(b1, b2, &run_options.window),
    }
}

pub fn distance_dyn(x1: f32, x2: f32, y1: f32, y2: f32, run_options: &RunOptions) -> f32 {
    match run_options.distance {
        Distance::EucEnclosed => simple_distance(x1, x2, y1, y2, &run_options.window),
        Distance::EucToroidal => toroidal_distance(x1, x2, y1, y2, &run_options.window),
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
    toroidal_distance_sq_boid(b1, b2, window_size).sqrt()
}

pub fn toroidal_distance_sq_boid(b1: &Boid, b2: &Boid, window_size: &WindowSize) -> f32 {
    toroidal_distance2_sq(b1.position.x, b2.position.x, b1.position.y, b2.position.y, window_size)
    // thoroidal_distance(b1.pos_x, b2.pos_x, b1.pos_y, b2.pos_y)
}

#[allow(dead_code)]
pub fn toroidal_distance(x1: f32, x2: f32, y1: f32, y2: f32, window_size: &WindowSize) -> f32 {
    toroidal_distance2_sq(x1, x2, y1, y2, window_size).sqrt()
}

// Thoroidal distance as per https://blog.demofox.org/2017/10/01/calculating-the-distance-between-points-in-wrap-around-toroidal-space/
/// Get distance in between points in a thoroidal space "continuous", that "wraps around"
pub fn toroidal_distance_sq(x1: f32, x2: f32, y1: f32, y2: f32, window_size: &WindowSize) -> f32 {
    let mut dx = (x1 - x2).abs();
    let mut dy = (y1 - y2).abs();

    // coordinates are centered w.r.t. (0, 0)
    // that means that the highest x and y values are exactly half of width and height respectively

    // e.g., if on dimension x we got an absolute difference, dx = |x1 - x2| larger than |x|/2
    // we went the wrong direction and should consider the other side to get the minimal distance

    if dx > window_size.win_right as f32 {
        // dx = 1.0 - dx;
        dx = window_size.win_w as f32 - dx;
    }

    if dy > window_size.win_top as f32 {
        // dy = 1.0 - dy;
        dy = window_size.win_h as f32 - dy;
    }

    return dx.powi(2) + dy.powi(2);
}

pub fn toroidal_distance2_sq(x1: f32, x2: f32, y1: f32, y2: f32, window_size: &WindowSize) -> f32 {
    // let px_1 = x1 + window_size.win_left as f32;
    // let px_2 = x2 +window_size.win_right as f32;
    // let dx = 
    // (px_1 - px_2).abs()
    // .min((px_1 - (px_2 + window_size.win_w as f32)).abs())
    // .min((px_1 - (px_2 - window_size.win_w as f32)).abs());

    // let py_1 = y1 + window_size.win_top as f32;
    // let py_2 = y2 +window_size.win_top as f32;
    // let dy = 
    // (py_1 - py_2).abs()
    // .min((py_1 - (py_2 + window_size.win_h as f32)).abs())
    // .min((py_1 - (py_2 - window_size.win_h as f32)).abs());

    let (dx, dy) = tor_vec_p(x1, x2, y1, y2, window_size);

    return dx.powi(2) + dy.powi(2);
}

pub fn tor_vec(v1: Vec2, v2: Vec2, window_size: &WindowSize) -> Vec2 {
    // let v = tor_vec_p(v1.x, v2.x, v1.y, v2.y, window_size);
    // vec2(v)
    Vec2::from(tor_vec_p(v1.x, v2.x, v1.y, v2.y, window_size))
}

/// this producess a vector along the shortest path from p1(x1,y1) to p2 (x2, t2)
pub fn tor_vec_p(x1: f32, x2: f32, y1: f32, y2: f32, window_size: &WindowSize) -> (f32, f32) {
    // components pointing towards P2
    let mut dx_p2 = x2 - x1;
    let mut dy_p2 = y2 - y1;

    if dx_p2.abs() > window_size.win_right as f32 {
        dx_p2 = if dx_p2 < 0. { 
            window_size.win_w 
        } else { 
            -window_size.win_w 
        } as f32 + dx_p2;
    }

    if dy_p2.abs() > window_size.win_top as f32 {
        dy_p2 = if dy_p2 < 0. { 
            window_size.win_h
        } else { 
            -window_size.win_h 
        } as f32 + dy_p2;
    }

    (dx_p2, dy_p2)
}



// TODO: test


// pub fn rand_circle_point(){

// }

#[cfg(test)]
mod tests {
    use approx::{assert_relative_eq};
    use glam::Vec2;
    // use nannou::prelude::{vec2, Vec2};

    use crate::options::{WindowSize};

    use super::{toroidal_distance, tor_vec_p};

    macro_rules! assert_eqf32 {
        ($x:expr, $y:expr) => {
            assert_relative_eq!($x, $y, epsilon = 1e-3f32)
        };
    }

    #[test]
    fn test_tor_vec_p() {
        let window_size = WindowSize{
            win_left: -35,
            win_right: 35,
            win_top: 35,
            win_bottom: -35,
            win_h: 70,
            win_w: 70,
        };

        let x1 = 30.;
        let y1 = -30.;
        let x2 = -30.;
        let y2 = 30.;

        let res = tor_vec_p(x1, x2, y1, y2, &window_size);

        assert_eqf32!(res.0, 10.);
        assert_eqf32!(res.1, -10.);

        let res2 = tor_vec_p(x2, x1, y2, y1, &window_size);

        assert_eqf32!(res2.0, -10.);
        assert_eqf32!(res2.1, 10.);
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
                 a: Vec2::new(3., 2.),
                 b: Vec2::new(0., 0.),
                 c: Vec2::new(-3., 3.),
                 d: Vec2::new(-3., -3.5)
            }
        }
    }

    fn _thoroidal_distance(v1: Vec2, v2: Vec2) -> f32 {
        let window_size = WindowSize{
            win_left: -4,
            win_right: 4,
            win_top: 4,
            win_bottom: -4,
            win_h: 8,
            win_w: 8,
        };

        toroidal_distance(v1.x, v2.x, v1.y, v2.y, &window_size)
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