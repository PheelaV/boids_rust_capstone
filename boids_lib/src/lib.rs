use birdwatcher::{Birdwatcher, BoidData};
use flock::Flock;
use options::RunOptions;

pub mod boid;
pub mod flock;

pub mod birdwatcher;
pub mod math_helpers;
pub mod options;

pub fn flock_base(no_iter: u64, run_options: RunOptions) -> Vec<BoidData> {
    let mut ro = run_options;
    let mut flock = Flock::new(&ro);
    let mut bird_watcher = Birdwatcher::new(ro.sample_rate);

    (0..no_iter).for_each(|_| {
        flock.update(&mut ro);
        bird_watcher.watch(&flock);
    });

    let data = bird_watcher.pop_data_save(&ro.save_options);

    data
}
