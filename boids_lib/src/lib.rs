use birdwatcher::{Birdwatcher, BoidData};
use flock::Flock;
use options::RunOptions;

pub mod boid;
pub mod flock;

pub mod math_helpers;
pub mod options;
pub mod birdwatcher;


pub fn flock_base(no_iter: u64, run_options: RunOptions) -> Vec<BoidData> {

    let mut flock = Flock::new(&run_options);
    let mut bird_watcher = Birdwatcher::new(run_options.sample_rate);

    (0..no_iter).for_each(|_| {
        flock.update(&run_options);
        bird_watcher.watch(&flock);
    });

    let data = bird_watcher
    .pop_data_save(&run_options.save_options);

    data
}