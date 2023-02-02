use boids_lib::{
    birdwatcher::Birdwatcher,
    flock::Flock,
    options::{
        RunOptions, SaveOptions, self, Boundary
    },
};
use extendr_api::prelude::*;

use geo::{ConvexHull, MultiPoint, Point, Area};

#[derive(Debug, IntoDataFrameRow)]
struct TestData {
    id: usize,
    x: f32,
    y: f32,
    cluster_id: usize,
    n_neighbours: usize,
    time: u64
}

#[extendr]
fn flock(no_iter: u32, init_boids: u32, save_locations_path: String) -> () {

    let mut run_options: RunOptions = Default::default();

    run_options.init_boids = init_boids as usize;

    run_options.save_options = get_save_options(save_locations_path);
    flock_base(no_iter, run_options);
}

#[extendr]
/// executes flocking and returns a dataframe with the location data
fn flock_return(no_iter: u32, init_boids: u32, save_locations_path: String, sample_rate: u64, init_width: f32, init_height: f32) -> Robj {

    let mut run_options: RunOptions = Default::default();

    run_options.init_boids = init_boids as usize;
    run_options.save_options = get_save_options(save_locations_path);
    run_options.sample_rate = sample_rate;

    run_options.window = options::get_window_size(init_width, init_height);
    
    flock_base(no_iter, run_options)
}

#[extendr]
/// executes flocking and returns a dataframe with the location data
fn flock_detailed(
    no_iter: u32,
    init_boids: u32,
    save_locations_path: String,
    sample_rate: u64,
    init_width: f32,
    init_height: f32,
    sensory_distance: f32,
    allignment_coef: f32,
    cohesion_coef: f32,
    separation_coef: f32,
    allignment_trs_coef:f32,
    cohesion_trs_coef:f32,
    separation_trs_coef: f32,
    min_speed: f32,
    max_speed: f32,
    max_steering: f32,
    dbscan_clustering: bool,
    boundary_config: Option<&str>) -> Robj {

    let mut run_options: RunOptions = Default::default();

    run_options.init_boids = init_boids as usize;
    run_options.save_options = get_save_options(save_locations_path);
    run_options.sample_rate = sample_rate;

    run_options.window = options::get_window_size(init_width, init_height);

    run_options.sensory_distance = sensory_distance;
    run_options.allignment_coefficient = allignment_coef;
    run_options.cohesion_coefficient = cohesion_coef;
    run_options.separation_coefficient = separation_coef;
    run_options.allignment_treshold_coefficient = allignment_trs_coef;
    run_options.cohesion_treshold_coefficient = cohesion_trs_coef;
    run_options.separation_treshold_coefficient = separation_trs_coef;
    run_options.min_speed = min_speed;
    run_options.max_speed = max_speed;
    run_options.max_steering = max_steering;
    run_options.dbscan_flock_clustering_on = dbscan_clustering;
   
    if let Some(boundary_config_string) = boundary_config{
        match serde_json::from_str::<Boundary>(boundary_config_string) {
            Ok(boundary) => run_options.boundary = boundary,
            Err(err) => panic!("Error, boundary deserialization failed: {}", err)
        }
    }
    
    // let boundary: Boundary = ;
    
    flock_base(no_iter, run_options)
}

/// private function for "overloading"
fn flock_base(no_iter: u32, run_options: RunOptions) -> Robj {

    let mut flock = Flock::new(&run_options);
    let mut bird_watcher = Birdwatcher::new(run_options.sample_rate);

    (0..no_iter).for_each(|_| {
        flock.update(&run_options);
        bird_watcher.watch(&flock);
    });

    let data = bird_watcher
    .pop_data_save(&run_options.save_options);
    // this might force recompillation? 

    data
    .iter()
    .map(|bd| TestData {
        id: bd.id,
        x: bd.x,
        y: bd.y,
        cluster_id: bd.cluster_id,
        n_neighbours: bd.n_neighbours,
        time: bd.time
    })
    .collect::<Vec<TestData>>()
    .into_dataframe()
    .unwrap()
    .as_robj()
    .to_owned()
        
}
fn get_save_options(save_locations_path: String) -> SaveOptions {
    let path: Option<String> ;
    if save_locations_path.is_empty(){
        path = None;
    } else {
        path = Some(save_locations_path);
    } 

    SaveOptions {save_locations: true, save_locations_path: path, save_locations_timestamp: true}
}

#[extendr]
fn get_convex_hull(x: &[f64], y: &[f64]) -> f64 {
    if x.len() != y.len() { panic!("Inputs of unequal length are not allowed!"); }

    let points = x.iter()
    .zip(y.iter())
    .map(|(x, y)| {
        // Coord { x: *x, y: *y};
        Point::new(*x, *y)
    }).collect::<Vec<Point>>();

    let hull = MultiPoint::new(points).convex_hull();

    hull.unsigned_area()
}

// #[extendr]
// fn force_recompile() -> &'static str{
//     return r#"hello world"#
// }
// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod boidr;
    fn flock;
    fn flock_return;
    fn flock_detailed;
    // fn force_recompile;
    fn get_convex_hull;
}

#[cfg(test)]
mod tests {
    use boids_lib::options::Boundary;

    use super::*;

    #[test]
    fn test_boundary_setting_serializaion() {
        let data = r#"
        {
            "type": "Repulsive",
            "distance": 100,
            "force": 0.05
        }
        "#;
        
        let boundary: Boundary = serde_json::from_str(data).unwrap();

        match boundary {
            // Boundary::Thoroidal => todo!(),
            // Boundary::Absorbing => todo!(),
            // Boundary::Reflective => todo!(),
            Boundary::Repulsive { distance, force } => {
                assert_eq!(100., distance);
                assert_eq!(0.05, force);
            },
            _ => panic!("Different type of boundary inputted than is being tested.")
        }
    }
    #[test]
    fn test_boundary_setting_serializaion2() {
        let data = r#"{"type": "Repulsive", "distance": 100, "force": 0.05}"#;
        
        let boundary: Boundary = serde_json::from_str(data).unwrap();

        match boundary {
            // Boundary::Thoroidal => todo!(),
            // Boundary::Absorbing => todo!(),
            // Boundary::Reflective => todo!(),
            Boundary::Repulsive { distance, force } => {
                assert_eq!(100., distance);
                assert_eq!(0.05, force);
            },
            _ => panic!("Different type of boundary inputted than is being tested.")
        }
    }
}