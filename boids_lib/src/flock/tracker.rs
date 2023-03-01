use crate::{boid::{Boid, BoidMetadata}, options::RunOptions};

// a tracker takes ownership of a flock of boids
// a tracker will follow an implementation consisting of choosing a datastructure to
// hold spatial information about boids in and have traits for querying 2D, later 3D
// information
pub trait Tracker {
    fn new(entities: &[Boid], run_options: &RunOptions) -> Self
    where
        Self: Sized;
    fn update(&mut self, run_options: &RunOptions);
    fn insert_single(&mut self, entity: Boid, _: &RunOptions);
    fn restart(&mut self, entities: &[Boid]);
    fn delete_last(&mut self) -> Option<Boid>;
    fn delete_multiple(&mut self, ids_delete: &[usize], run_options: &RunOptions);
    fn insert_multiple(&mut self, entities: &[Boid], _: &RunOptions);
    fn view(&self) -> (&Vec<Boid>, &Vec<BoidMetadata>);
    fn view2<'a>(&'a self) -> Box<dyn Iterator<Item = (&'a Boid, &'a BoidMetadata)> + 'a>;
    fn get_neighbours<'a>(&'a self, boid: &Boid, run_options: &RunOptions) -> Vec<&'a Boid>;
}