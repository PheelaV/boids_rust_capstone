use glam::Vec2;



use crate::{
    boid::Boid,
    options::{Distance, RunOptions, WindowSize},
};

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
    simple_distance_sq(
        b1.position.x,
        b2.position.x,
        b1.position.y,
        b2.position.y,
        _window_size,
    )
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
    toroidal_distance_sq(
        b1.position.x,
        b2.position.x,
        b1.position.y,
        b2.position.y,
        window_size,
    )
}

#[allow(dead_code)]
pub fn toroidal_distance(x1: f32, x2: f32, y1: f32, y2: f32, window_size: &WindowSize) -> f32 {
    toroidal_distance_sq(x1, x2, y1, y2, window_size).sqrt()
}

/// Get distance in between points in a thoroidal space "continuous", that "wraps around"
// pub fn toroidal_distance_sq(x1: f32, x2: f32, y1: f32, y2: f32, window_size: &WindowSize) -> f32 {
//     let mut dx = (x1 - x2).abs();
//     let mut dy = (y1 - y2).abs();

//     if dx > window_size.win_right as f32 {
//         // dx = 1.0 - dx;
//         dx = window_size.win_w as f32 - dx;
//     }

//     if dy > window_size.win_top as f32 {
//         // dy = 1.0 - dy;
//         dy = window_size.win_h as f32 - dy;
//     }
//     return dx.powi(2) + dy.powi(2);
// }

pub fn toroidal_distance_sq(x1: f32, x2: f32, y1: f32, y2: f32, window_size: &WindowSize) -> f32 {
    let (dx, dy) = tor_vec_p(x1, x2, y1, y2, window_size);

    return dx.powi(2) + dy.powi(2);
}

pub fn tor_vec(v1: Vec2, v2: Vec2, window_size: &WindowSize) -> Vec2 {
    Vec2::from(tor_vec_p(v1.x, v2.x, v1.y, v2.y, window_size))
}

/// Producess a vector along the shortest path from p1(x1,y1) to p2 (x2, t2) in a toroidal space.
/// Originally inspired by toroidal distance as per [source](https://blog.demofox.org/2017/10/01/calculating-the-distance-between-points-in-wrap-around-toroidal-space/),
/// with a modification to preserve directionality (sign)
#[inline]
pub fn tor_vec_p(x1: f32, x2: f32, y1: f32, y2: f32, window_size: &WindowSize) -> (f32, f32) {
    // The idea big idea is:
    // coordinates are centered w.r.t. (0, 0)
    // e.g., if on dimension x we got an absolute difference, dx = |x1 - x2| larger than |x|/2
    // we went the wrong direction and should consider the other side to get the minimal distance
    // this additionally preserves the sign of the shortest path

    let dx_p2 = tor_vec_pc(x1, x2, window_size.win_w as f32, window_size.win_right as f32);
    let dy_p2 = tor_vec_pc(y1, y2, window_size.win_h as f32, window_size.win_top as f32);

    (dx_p2, dy_p2)
}

/// pc - point's component
#[inline]
pub fn tor_vec_pc(x1: f32, x2: f32, size: f32, max: f32) -> f32 {
    // component pointing towards x2
    let d_p2 = x2 - x1;
    // if we are crossing more than half the space, the other way it shall be
    if d_p2.abs() > max {
        d_p2 + if d_p2 < 0. { size } else { -size } as f32
    } else {
        d_p2
    }
}

// this does not work as I would need to provide a custom implementation of the nn datastructure (could be flock.rs, but
// it would neet a large separation of conerns, i.e. it is just doing way too much atm)
// use linfa_nn::distance::Distance as lDistance;
// #[derive(Debug, Clone, PartialEq, Eq)]
// pub struct TorL2Dist{
//     pub window: WindowSize
// }
// impl TorL2Dist {
//     pub fn new(window_size: WindowSize) -> Self {
//         TorL2Dist { window: window_size }
//     }
// }
// impl<F: Float> lDistance<F> for TorL2Dist {
//     fn distance<D: ndarray::Dimension>(&self, a: ndarray::ArrayView<F, D>, b: ndarray::ArrayView<F, D>) -> F {
//         let mut result = F::zero();

//         // additional checks in debug
//         if cfg!(debug_assertions){
//             let a_shape = a.shape();
//             let b_shape = b.shape();

//             if a_shape != b_shape {
//                 panic!("shapes must equal")
//             }

//             if a_shape[0] != 2 {
//                 panic!("only two dimensions are supported")
//             }
//         }
//         // const first_axis: Axis = Axis(0);

//         // for (aa, bb) in a.axis_iter(first_axis).zip(b.axis_iter(first_axis)) {

//         // }
//         // let s = *a.get::<NdIndex<Ix1>>(0.into()).unwrap();

//         // let t =  tor_vec_p(
//         //     *(a.get),
//         //     *b.get(0.into()).unwrap(),
//         //     *a.get(1.into()).unwrap(),
//         //     *b.get(1.into()).unwrap(),
//         //      &self.window);

//         let mut odd = false;
//         Zip::from(a).and(b).for_each(|ai,  bi| {
//             let dist: F;
//             if !odd {
//                 dist = F::from(tor_vec_pc(ai.clone().to_f32().unwrap(), bi.clone().to_f32().unwrap(), self.window.win_w as f32, self.window.win_right as f32)).unwrap();
//             } else {
//                 dist = F::from(tor_vec_pc(ai.clone().to_f32().unwrap(), bi.clone().to_f32().unwrap(), self.window.win_h as f32, self.window.win_top as f32)).unwrap();
//             }
//             odd = !odd;
//             result += dist * dist;
//         });
//         // Zip::from(a).and(b).apply(function)

//         result.sqrt()
//     }
// }

// because of nannou I have to keep glam version at .17
// as glam is not stable yet, there are breaking changes
// in between decimal versions, somewhere between .17 and .23
pub(crate) trait MyRotate {
    fn rotate(&self, rhs: Vec2) -> Self;
}

impl MyRotate for Vec2 {
    #[inline]
    fn rotate(&self, rhs: Vec2) -> Self {
        Vec2::new(
            self.x * rhs.x - self.y * rhs.y,
            self.y * rhs.x + self.x * rhs.y,
        )
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use glam::Vec2;

    use super::{tor_vec_p, toroidal_distance};
    use crate::options::{get_window_size, WindowSize};

    macro_rules! assert_eqf32 {
        ($x:expr, $y:expr) => {
            assert_relative_eq!($x, $y, epsilon = 1e-3_f32)
        };
    }

    #[test]
    fn test_tor_vec_p() {
        let window_size = get_window_size(70, 70);

        let x1 = 30_f32;
        let y1 = -30_f32;
        let x2 = -30_f32;
        let y2 = 30_f32;

        let res = tor_vec_p(x1, x2, y1, y2, &window_size);

        assert_eqf32!(res.0, 10.);
        assert_eqf32!(res.1, -10.);

        let res2 = tor_vec_p(x2, x1, y2, y1, &window_size);

        assert_eqf32!(res2.0, -10.);
        assert_eqf32!(res2.1, 10.);
    }

    struct ToroidialTestData {
        a: Vec2,
        b: Vec2,
        c: Vec2,
        d: Vec2,
    }

    impl ToroidialTestData {
        fn new() -> Self {
            ToroidialTestData {
                a: Vec2::new(3., 2.),
                b: Vec2::new(0., 0.),
                c: Vec2::new(-3., 3.),
                d: Vec2::new(-3., -3.5),
            }
        }
    }

    fn _toroidal_distance(v1: Vec2, v2: Vec2) -> f32 {
        let window_size = WindowSize {
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
        let data = ToroidialTestData::new();
        let result = _toroidal_distance(data.a, data.b);
        assert_eqf32!(result, 3.606);
    }

    #[test]
    fn thoroidal_distance_a_c() {
        let data = ToroidialTestData::new();
        let result = _toroidal_distance(data.a, data.c);
        assert_eqf32!(result, 2.236);
    }

    #[test]
    fn thoroidal_distance_a_d() {
        let data = ToroidialTestData::new();
        let result = _toroidal_distance(data.a, data.d);
        assert_eqf32!(result, 3.202);
    }

    #[test]
    fn thoroidal_distance_c_d() {
        let data = ToroidialTestData::new();
        let result = _toroidal_distance(data.c, data.d);
        assert_eqf32!(result, 1.5);
    }
}
