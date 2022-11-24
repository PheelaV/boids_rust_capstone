use boids_lib::{
    birdwatcher::Birdwatcher,
    flock::Flock,
    options::{
        get_run_options, RunOptions, SaveOptions, self
    },
};
use extendr_api::prelude::*;

#[derive(Debug, IntoDataFrameRow)]
struct TestData {
    id: u32,
    x: f32,
    y: f32,
}

#[extendr]
fn flock(no_iter: u32, init_boids: u32, save_locations_path: String) -> () {

    let mut run_options = get_run_options();

    run_options.init_boids = init_boids;

    run_options.save_options = get_save_options(save_locations_path);
    flock_base(no_iter, run_options);
}

#[extendr]
/// executes flocking and returns a dataframe with the location data
fn flock_return(no_iter: u32, init_boids: u32, save_locations_path: String, sample_rate: u32, init_width: f32, init_height: f32) -> Robj {

    let mut run_options = get_run_options();

    run_options.init_boids = init_boids;
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
    sample_rate: u32,
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
    max_steering: f32) -> Robj {

    let mut run_options = get_run_options();

    run_options.init_boids = init_boids;
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

    data
    .iter()
    .map(|bd| TestData {
        id: bd.id,
        x: bd.x,
        y: bd.y,
    })
    .collect::<Vec<TestData>>()
    .into_dataframe()
    .unwrap()
    .as_robj()
    .to_owned()
        
}
fn get_save_options(save_locations_path: String) -> SaveOptions {
    let mut path : Option<String> = None;
    if save_locations_path.is_empty(){
        path = None;
    } else {
        path = Some(save_locations_path);
    } 

    SaveOptions {save_locations: true, save_locations_path: path, save_locations_timestamp: true}
}
// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod boidr;
    // fn hello_world;
    fn flock;
    fn flock_return;
    fn flock_detailed;
}

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn burn_test() {
//         hello_world2();
//     }
// }
