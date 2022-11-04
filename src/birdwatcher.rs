use std::mem;

use serde::Serialize;

use crate::flock::Flock;

// todo: it would be cool to have bird watcher store a reference to the flock, both would have to have the same lifetime anotations
// the problem is that the nannou Model would also have to be annotated 

// so right now, this is more of a bird data acummulator than a birdwatcher
#[derive(Serialize, Debug)]
pub struct BoidData {
    pub id: u32,
    pub x: f32,
    pub y: f32,
}

pub struct Birdwatcher{
    locations: Vec<BoidData>,
    // flock: &'a Flock,
}

impl Birdwatcher {
    pub fn new() -> Self {

        Birdwatcher
        { 
            locations: Vec::new(), 
            // flock
        }
    }

    /// Triggers data collection
    pub fn watch(&mut self, flock: &Flock) -> () {
        let mut current_locations : Vec<BoidData> = 
            flock.boids.iter()
            .map(|b| {
                BoidData{id: b.id, x: b.position.x, y:b.position.y}
            })
            .collect();
        
        self.locations.append(&mut current_locations);
    }

    pub fn pop_data(&mut self) -> Vec<BoidData> {
        mem::take(&mut self.locations)
    }
}