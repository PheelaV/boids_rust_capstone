use crate::options::RunOptions;

use super::tracker::Tracker;

pub struct ReplayTracker {

}

impl ReplayTracker {
    pub fn from_path(replay_path: &String, run_options: &RunOptions) -> Self {
        todo!();
    }
}

impl Tracker for ReplayTracker {
    fn new(entities: &[crate::boid::Boid], run_options: &RunOptions) -> Self
    where
        Self: Sized {
        todo!()
    }

    fn update(&mut self, run_options: &RunOptions) {
        todo!()
    }

    fn insert_single(&mut self, _: crate::boid::Boid, _: &RunOptions) {
        
    }

    fn restart(&mut self, _: &[crate::boid::Boid]) {
        
    }

    fn delete_last(&mut self) -> Option<crate::boid::Boid> {
        todo!()
    }

    fn delete_multiple(&mut self, _: &[usize], _: &RunOptions) {
        todo!()
    }

    fn insert_multiple(&mut self, _: &[crate::boid::Boid], _: &RunOptions) {
        todo!()
    }

    fn view(&self) -> (&Vec<crate::boid::Boid>, &Vec<crate::boid::BoidMetadata>) {
        todo!()
    }

    fn view2<'a>(&'a self) -> Box<dyn Iterator<Item = (&'a crate::boid::Boid, &'a crate::boid::BoidMetadata)> + 'a> {
        todo!()
    }

    fn get_neighbours<'a>(&'a self, _: &crate::boid::Boid, _: &RunOptions) -> Vec<&'a crate::boid::Boid> {
        todo!()
    }
}