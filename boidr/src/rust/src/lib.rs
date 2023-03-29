use std::path::Path;
use boids_lib::{
    birdwatcher::{Birdwatcher, BoidData},
    flock::Flock,
    math_helpers::tor_vec_pc,
    options::{self, Boundary, Distance, RunOptions, SaveOptions, NoiseModel},
};
use csv::{Reader, Writer};
use extendr_api::prelude::*;

use geo::{Area, ConvexHull, LineString, MultiPoint, Point, Polygon};
use serde::{Deserialize, Serialize};
use voronoice::{BoundingBox, VoronoiBuilder};

#[derive(Debug, IntoDataFrameRow)]
struct TestData {
    id: usize,
    x: f32,
    y: f32,
    cluster_id: usize,
    n_neighbours: usize,
    time: u64,
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
fn flock_return(
    no_iter: u32,
    init_boids: u32,
    save_locations_path: String,
    sample_rate: f32,
    init_width: u32,
    init_height: u32,
) -> Robj {
    let mut run_options: RunOptions = Default::default();

    run_options.init_boids = init_boids as usize;
    run_options.save_options = get_save_options(save_locations_path);
    run_options.sample_rate = sample_rate as u16; // this intentional is because of R shenanigans

    run_options.window = options::get_window_size(init_width, init_height);

    flock_base(no_iter, run_options)
}

#[extendr]
/// executes flocking and returns a dataframe with the location data
fn flock_detailed_no_return(
    no_iter: u32,
    init_boids: u32,
    save_locations_path: String,
    sample_rate: f32,
    init_width: u32,
    init_height: u32,
    sensory_distance: f32,
    alignment_coef: f32,
    cohesion_coef: f32,
    separation_coef: f32,
    alignment_trs_coef: f32,
    cohesion_trs_coef: f32,
    separation_trs_coef: f32,
    min_speed: f32,
    agent_steering: bool,
    max_speed: f32,
    max_steering: f32,
    dbscan_clustering: bool,
    boundary_config: Option<&str>,
    distance_config: Option<&str>,
    field_of_vision: f32,
    rules_impl: bool,
    wander_on: bool,
    wander_coef: f32,
    wander_rate: f32,
    wander_radius: f32,
    wander_distance: f32,
    baseline_speed: f32,
    wander_random: bool,
) -> () {
    let mut run_options: RunOptions = Default::default();

    run_options.init_boids = init_boids as usize;
    run_options.save_options = get_save_options(save_locations_path);
    run_options.sample_rate = sample_rate as u16; //  

    run_options.window = options::get_window_size(init_width, init_height);

    run_options.sensory_distance = sensory_distance;
    run_options.alignment_coefficient = alignment_coef;
    run_options.cohesion_coefficient = cohesion_coef;
    run_options.separation_coefficient = separation_coef;
    run_options.alignment_treshold_coefficient = alignment_trs_coef;
    run_options.cohesion_treshold_coefficient = cohesion_trs_coef;
    run_options.separation_treshold_coefficient = separation_trs_coef;
    run_options.min_speed = min_speed;
    run_options.max_speed = max_speed;
    run_options.agent_steering = agent_steering;
    run_options.max_steering = max_steering;
    run_options.dbscan_flock_clustering_on = dbscan_clustering;
    run_options.field_of_vision_deg = field_of_vision;
    run_options.wander_on = wander_on;
    run_options.rules_impl = rules_impl;
    run_options.wander_rate = wander_rate;
    run_options.wander_radius = wander_radius;
    run_options.wander_coefficient = wander_coef;
    run_options.wander_distance = wander_distance;
    run_options.baseline_speed = baseline_speed;
    run_options.noise_model = match wander_random {
        true => NoiseModel::Vicsek,
        false => NoiseModel::Reynolds,
    };
    // attempts to retreive the boundary, if none is set, uses reflective as default
    if let Some(boundary_config_string) = boundary_config {
        match serde_json::from_str::<Boundary>(boundary_config_string) {
            Ok(boundary) => run_options.boundary = boundary,
            Err(err) => panic!("Error, boundary deserialization failed: {}", err),
        }
    } else {
        run_options.boundary = Boundary::Reflective
    }

    // attempts to retreive the distance, if none is set uses:
    // - euclidean toroidial if space is toroidial
    // - euclidean
    if let Some(distance_config_str) = distance_config {
        match serde_json::from_str::<Distance>(distance_config_str) {
            Ok(distance) => run_options.distance = distance,
            Err(err) => panic!("Error, distance deserialization failed: {}", err),
        }
    } else {
        run_options.distance = match run_options.boundary {
            Boundary::Toroidal => Distance::EucToroidal,
            _ => Distance::EucEnclosed,
        }
    }

    flock_base_no_return(no_iter, run_options)
}

#[extendr]
/// executes flocking and returns a dataframe with the location data
fn flock_detailed(
    no_iter: u32,
    init_boids: u32,
    save_locations_path: String,
    sample_rate: f32,
    init_width: u32,
    init_height: u32,
    sensory_distance: f32,
    alignment_coef: f32,
    cohesion_coef: f32,
    separation_coef: f32,
    alignment_trs_coef: f32,
    cohesion_trs_coef: f32,
    separation_trs_coef: f32,
    min_speed: f32,
    agent_steering: bool,
    max_speed: f32,
    max_steering: f32,
    dbscan_clustering: bool,
    boundary_config: Option<&str>,
    distance_config: Option<&str>,
    field_of_vision: f32,
    rules_impl: bool,
    wander_on: bool,
    wander_coef: f32,
    wander_rate: f32,
    wander_radius: f32,
    wander_distance: f32,
    baseline_speed: f32,
    wander_random: bool,
) -> Robj {
    let mut run_options: RunOptions = Default::default();

    run_options.init_boids = init_boids as usize;
    run_options.save_options = get_save_options(save_locations_path);
    run_options.sample_rate = sample_rate as u16; //  

    run_options.window = options::get_window_size(init_width, init_height);

    run_options.sensory_distance = sensory_distance;
    run_options.alignment_coefficient = alignment_coef;
    run_options.cohesion_coefficient = cohesion_coef;
    run_options.separation_coefficient = separation_coef;
    run_options.alignment_treshold_coefficient = alignment_trs_coef;
    run_options.cohesion_treshold_coefficient = cohesion_trs_coef;
    run_options.separation_treshold_coefficient = separation_trs_coef;
    run_options.min_speed = min_speed;
    run_options.max_speed = max_speed;
    run_options.agent_steering = agent_steering;
    run_options.max_steering = max_steering;
    run_options.dbscan_flock_clustering_on = dbscan_clustering;
    run_options.field_of_vision_deg = field_of_vision;
    run_options.wander_on = wander_on;
    run_options.rules_impl = rules_impl;
    run_options.wander_rate = wander_rate;
    run_options.wander_radius = wander_radius;
    run_options.wander_coefficient = wander_coef;
    run_options.wander_distance = wander_distance;
    run_options.baseline_speed = baseline_speed;
    run_options.noise_model = match wander_random {
        true => NoiseModel::Vicsek,
        false => NoiseModel::Reynolds,
    };
    // attempts to retreive the boundary, if none is set, uses reflective as default
    if let Some(boundary_config_string) = boundary_config {
        match serde_json::from_str::<Boundary>(boundary_config_string) {
            Ok(boundary) => run_options.boundary = boundary,
            Err(err) => panic!("Error, boundary deserialization failed: {}", err),
        }
    } else {
        run_options.boundary = Boundary::Reflective
    }

    // attempts to retreive the distance, if none is set uses:
    // - euclidean toroidial if space is toroidial
    // - euclidean
    if let Some(distance_config_str) = distance_config {
        match serde_json::from_str::<Distance>(distance_config_str) {
            Ok(distance) => run_options.distance = distance,
            Err(err) => panic!("Error, distance deserialization failed: {}", err),
        }
    } else {
        run_options.distance = match run_options.boundary {
            Boundary::Toroidal => Distance::EucToroidal,
            _ => Distance::EucEnclosed,
        }
    }

    flock_base(no_iter, run_options)
}

/// private function for "overloading"
fn flock_base(no_iter: u32, run_options: RunOptions) -> Robj {
    let mut ro = run_options;
    let mut flock = Flock::new(&ro);
    let mut bird_watcher = Birdwatcher::new(ro.sample_rate);

    (0..no_iter).for_each(|_| {
        flock.update(&mut ro);
        bird_watcher.watch(&flock);
    });

    let data = bird_watcher.pop_data_save(&ro.save_options);

    data.iter()
        .map(|bd| TestData {
            id: bd.id,
            x: bd.x,
            y: bd.y,
            cluster_id: bd.cluster_id,
            n_neighbours: bd.n_neighbours,
            time: bd.time,
        })
        .collect::<Vec<TestData>>()
        .into_dataframe()
        .unwrap()
        .as_robj()
        .to_owned()
}

fn flock_base_no_return(no_iter: u32, run_options: RunOptions) -> () {
    let mut ro = run_options;
    let mut flock = Flock::new(&ro);
    let mut bird_watcher = Birdwatcher::new(ro.sample_rate);

    (0..no_iter).for_each(|_| {
        flock.update(&mut ro);
        bird_watcher.watch(&flock);
    });

    _ = bird_watcher.pop_data_save(&ro.save_options);
}
fn get_save_options(save_locations_path: String) -> SaveOptions {
    let path: Option<String>;
    if save_locations_path.is_empty() {
        path = None;
    } else {
        path = Some(save_locations_path);
    }

    SaveOptions {
        save_locations: true,
        save_locations_path: path,
        save_locations_timestamp: true,
    }
}

#[extendr]
fn get_convex_hull(x: &[f64], y: &[f64]) -> f64 {
    if x.len() != y.len() {
        panic!("Inputs of unequal length are not allowed!");
    }

    let points = x
        .iter()
        .zip(y.iter())
        .map(|(x, y)| {
            // Coord { x: *x, y: *y};
            Point::new(*x, *y)
        })
        .collect::<Vec<Point>>();

    let hull = MultiPoint::new(points).convex_hull();

    hull.unsigned_area()
}

#[extendr]
fn get_voronoi_areas(x: &[f64], y: &[f64], width: i32, height: i32) -> Vec<f64> {
    if x.len() != y.len() {
        panic!("Inputs of unequal length are not allowed!");
    }

    let points: Vec<voronoice::Point> = x
        .iter()
        .zip(y.iter())
        .map(|(x, y)| voronoice::Point { x: *x, y: *y })
        .collect();

    let my_voronoi = VoronoiBuilder::default()
        .set_sites(points)
        .set_bounding_box(BoundingBox::new_centered(width as f64, height as f64))
        .set_lloyd_relaxation_iterations(5)
        .build();

    if let Some(v) = my_voronoi {
        v.iter_cells()
            .map(|c| {
                let geom_cell_points = c
                    .iter_vertices()
                    .map(|p| Point::new(p.x, p.y))
                    .collect::<Vec<Point>>();
                Polygon::new(LineString::from(geom_cell_points), vec![]).unsigned_area()
            })
            .collect()
    } else {
        panic!("Somethign went wrong with the voronoid creation!")
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy, Default)]
pub struct BoidDataPrep {
    pub id: usize,
    pub x: f32,
    pub y: f32,
    pub dx: f32,
    pub dy: f32,
    pub cluster_id: usize,
    pub n_neighbours: usize,
    pub time: u64,
}

impl BoidDataPrep {
    fn new(data: BoidData) -> Self {
        Self {
            id: data.id,
            x: data.x,
            y: data.y,
            dx: Default::default(),
            dy: Default::default(),
            cluster_id: data.cluster_id,
            n_neighbours: data.n_neighbours,
            time: data.time,
        }
    }
}

#[extendr]
fn preprocess_file(file_path: &str, init_width: u32, init_height: u32, no_boids: u32) -> String {
    let p = Path::new(file_path);

    if !p.exists() {
        panic!("Path does not exist! Path: {}", p.to_str().unwrap());
    }

    let mut results: Vec<BoidDataPrep> = Vec::new();

    // we need the records sorted by time and id
    // unless we want to get into external sorting on disk, we have to load everything and sort
    // to ensure that, at the moment birdwatcher records everything as is given 
    // (spathash1d does not guarantee order)
    let mut rdr = Reader::from_path(p).expect("Failed to open CSV file!");
    // let mut wtr2 = Writer::from_writer(vec![]);
    // let mut count: usize = 0;

    let win_max_w  = init_width as f32/ 2.;
    let win_max_h  = init_height as f32/ 2.;
    let win_w = init_width as f32;
    let win_h = init_height as f32;

    for result in rdr.deserialize() {
        let record: BoidData = result.expect("single record");
        let updated = BoidDataPrep::new(record);
        results.push(updated);
    }

    results.sort_unstable_by(|a, b| 
        a.time.cmp(&b.time).then(a.id.cmp(&b.id)));

    let time_steps = results.len() as u32 / no_boids;

    let diff = 1;
    for i in 0..no_boids {
        for t in diff..time_steps{
            let b_cur = (i + t * no_boids) as usize;
            let b_pas = (i + (t - diff) * no_boids) as usize;

            results[b_cur].dx = tor_vec_pc(results[b_pas].x, results[b_cur].x, win_w, win_max_w);
            results[b_cur].dy = tor_vec_pc(results[b_pas].y, results[b_cur].y, win_h, win_max_h);
        }
    }
    let location = match p.parent(){
        Some(parent) => parent.to_path_buf(),
        None => Path::new("./").to_path_buf()
    };

    let file_name = vec!["prepro_", p.file_name().unwrap().to_str().unwrap()].join("");
    let mut wtr2 = Writer::from_path(location.join(file_name)).unwrap();
    // let mut wtr = Writer::from_writer(vec![]);
    results.iter().for_each(|d| {
        wtr2.serialize(*d).expect("boid data prep failed");
    });

    results.len().to_string()
    // data
}

#[extendr]
fn force_recompile() -> &'static str {
    return r#"hello world"#; 
}

#[extendr]
fn test_execution() -> usize {

    let no_iter = 10000_u32;
    let mut run_options: RunOptions = Default::default();

    run_options.init_boids = 1024;
    run_options.sample_rate = 32;
    run_options.window = options::get_window_size(1000, 1000);
    run_options.sensory_distance = 45.;
    run_options.alignment_coefficient = 0.26;
    run_options.cohesion_coefficient = 0.9;
    run_options.separation_coefficient = 1.2;
    run_options.alignment_treshold_coefficient = 1.;
    run_options.cohesion_treshold_coefficient = 0.8;
    run_options.separation_treshold_coefficient = 0.3;
    run_options.min_speed = 2.;
    run_options.max_speed = 2.;
    run_options.agent_steering = true;
    run_options.max_steering = 0.7;
    run_options.dbscan_flock_clustering_on = true;
    run_options.field_of_vision_deg = 280.;
    run_options.wander_on = true;
    run_options.rules_impl = true;
    run_options.wander_rate = 1.;
    run_options.wander_radius = 16.0;
    run_options.wander_coefficient = 0.17;
    run_options.wander_distance = 23.0;
    run_options.baseline_speed = 1.0;
    run_options.noise_model = NoiseModel::Reynolds;
    run_options.boundary = Boundary::Toroidal;
    run_options.distance = Distance::EucToroidal;

    let mut flock = Flock::new(&run_options);
    let mut bird_watcher = Birdwatcher::new(run_options.sample_rate);

    (0..no_iter).for_each(|_| {
        flock.update(&mut run_options);
        bird_watcher.watch(&flock);
    });

    let data = bird_watcher.pop_data_save(&run_options.save_options);

    let converted = data.iter()
    .map(|bd| TestData {
        id: bd.id,
        x: bd.x,
        y: bd.y,
        cluster_id: bd.cluster_id,
        n_neighbours: bd.n_neighbours,
        time: bd.time,
    })
    .collect::<Vec<TestData>>();

     converted.len()
}
// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod boidr;
    fn flock;
    fn flock_return;
    fn flock_detailed;
    fn flock_detailed_no_return;
    fn force_recompile;
    fn get_convex_hull;
    fn get_voronoi_areas;
    fn preprocess_file;
    fn test_execution;
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
            // Boundary::Toroidal => todo!(),
            // Boundary::Absorbing => todo!(),
            // Boundary::Reflective => todo!(),
            Boundary::Repulsive { distance, force } => {
                assert_eq!(100., distance);
                assert_eq!(0.05, force);
            }
            _ => panic!("Different type of boundary inputted than is being tested."),
        }
    }
    #[test]
    fn test_boundary_setting_serializaion2() {
        let data = r#"{"type": "Repulsive", "distance": 100, "force": 0.05}"#;

        let boundary: Boundary = serde_json::from_str(data).unwrap();

        match boundary {
            // Boundary::Toroidal => todo!(),
            // Boundary::Absorbing => todo!(),
            // Boundary::Reflective => todo!(),
            Boundary::Repulsive { distance, force } => {
                assert_eq!(100., distance);
                assert_eq!(0.05, force);
            }
            _ => panic!("Different type of boundary inputted than is being tested."),
        }
    }
    #[test]
    fn test_execution() {

        let no_iter = 10000_u32;
        let mut run_options: RunOptions = Default::default();

        run_options.init_boids = 1024;
        run_options.sample_rate = 32;
        run_options.window = options::get_window_size(1000, 1000);
        run_options.sensory_distance = 45.;
        run_options.alignment_coefficient = 0.26;
        run_options.cohesion_coefficient = 0.9;
        run_options.separation_coefficient = 1.2;
        run_options.alignment_treshold_coefficient = 1.;
        run_options.cohesion_treshold_coefficient = 0.8;
        run_options.separation_treshold_coefficient = 0.3;
        run_options.min_speed = 2.;
        run_options.max_speed = 2.;
        run_options.agent_steering = true;
        run_options.max_steering = 0.7;
        run_options.dbscan_flock_clustering_on = true;
        run_options.field_of_vision_deg = 280.;
        run_options.wander_on = true;
        run_options.rules_impl = true;
        run_options.wander_rate = 1.0;
        run_options.wander_radius = 16.0;
        run_options.wander_coefficient = 0.2;
        run_options.wander_distance = 23.0;
        run_options.baseline_speed = 1.0;
        run_options.noise_model = NoiseModel::Reynolds;
        run_options.boundary = Boundary::Toroidal;
        run_options.distance = Distance::EucToroidal;

        let mut flock = Flock::new(&run_options);
        let mut bird_watcher = Birdwatcher::new(run_options.sample_rate);
    
        (0..no_iter).for_each(|_| {
            flock.update(&mut run_options);
            bird_watcher.watch(&flock);
        });
    
        let data = bird_watcher.pop_data_save(&run_options.save_options);
    
        let converted = data.iter()
        .map(|bd| TestData {
            id: bd.id,
            x: bd.x,
            y: bd.y,
            cluster_id: bd.cluster_id,
            n_neighbours: bd.n_neighbours,
            time: bd.time,
        })
        .collect::<Vec<TestData>>();

        println!("{}", converted.len()) 
        
    }
}
