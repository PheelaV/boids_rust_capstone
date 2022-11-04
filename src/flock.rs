use nannou::prelude::*;
use nannou::Draw;
use rand::Rng;

use crate::boid::Boid;
use crate::math_helpers::distance_dyn_boid;
use crate::options::InitiationStrategy;
use crate::options::Options;
use crate::options::RunOptions;

fn get_boids(no_boids: u32, options: &Options) -> Vec<Boid> {

    (0..no_boids).map(|id|{
        get_boid(options, id)
    })
    .collect()
}

fn get_boid(options: &Options, id: u32) -> Boid {
    let mut rng = rand::thread_rng();

    match options.initiation_strat {
        InitiationStrategy::CircleCenterOut => todo!(),
        InitiationStrategy::CircleCircumferenceIn => {
            let r = options.init_height as f32 / 2.0;

            let init_vel: f32 = rng.gen::<f32>() / 3. + 0.5;
            let init_pos: f32 = rng.gen::<f32>() / 3. + 2. / 3.;
            let deg = rng.gen::<f32>() * 360.0;
    
            let x = r * deg.cos();
            let y = r * deg.sin();
    
            let mut init_vec = vec2(-x, -y);
            init_vec = init_vec.clamp_length_max(0.5 * init_vel);
    
            Boid::new(x * init_pos, y * init_pos, init_vec, id)
        },
        InitiationStrategy::TwoWalls => todo!(),
        InitiationStrategy::RectangleIn => todo!(),
        InitiationStrategy::RandomIn => todo!(),
        InitiationStrategy::RandomRandom => {
            // x in [-width/2, width/2)
            let x = (rng.gen::<f32>() - 0.5) * options.init_width as f32;
            // y in [-height/2, height/2)
            let y = (rng.gen::<f32>() - 0.5) * options.init_height as f32;
            // x in [-max_speed, max_speed)
            let x_vel = (rng.gen::<f32>() * 2. - 1.) * options.run_options.max_speed;
            // y in [-max_speed, max_speed)
            let y_vel = (rng.gen::<f32>() * 2. - 1.) * options.run_options.max_speed;

            let mut init_vec = vec2(x_vel, y_vel);

            init_vec = init_vec.clamp_length(options.run_options.min_speed, options.run_options.max_speed);

            Boid::new(x, y, init_vec, id)
        },
    }
}

pub struct Flock {
    pub boids: Vec<Boid>,
}

impl Flock {
    pub fn new(no_boids: u32, options: &Options) -> Self {
        let boids = get_boids(no_boids, &options);

        Flock { 
            boids,
        }
    }

    pub fn update(&mut self, options: &Options) {
        let clone = self.boids.to_owned();

        let max_sensory_distance = 
        options.run_options.allignment_treshold_distance
        .max(options.run_options.cohesion_treshold_distance
            .max(options.run_options.separation_treshold_distance));

        for b_cur in self.boids.iter_mut() {
            let mut neighbours: Vec<&Boid> = Vec::new();

            for b_other in clone.iter() {
                if b_other.id == b_cur.id {
                    continue;
                }

                let distance = distance_dyn_boid(b_cur, b_other, &options.run_options);
                if distance < max_sensory_distance {
                    neighbours.push(b_other);
                }
            }

            b_cur.run_rules(neighbours, &options.run_options);
            b_cur.update_location(&options.run_options)
        }
    }

    pub fn insert(&mut self, options: &Options){
        let b = get_boid(options, self.boids.len() as u32);
        self.insert_boid(b);
    }

    pub fn insert_boid(&mut self, b: Boid) {
        self.boids.push(b);
    }

    pub fn delete_last(&mut self) -> Option<u32> {
        let b_delete = self.boids.get(self.boids.len() - 1);

        match b_delete {
            Some(x) => Some(self.delete_boid(x.id)),
            None => None,
        }
    }

    /// Deletes a boid with given id and reindexes boids
    pub fn delete_boid(&mut self, id_delete: u32) -> u32{
        self.boids = self.boids
        .iter_mut()
        .enumerate()
        .filter(|(_, b)|{ 
            b.id != id_delete
        })
        .map(|(index, b)| {
            b.id = index as u32;
            *b
        })
        .collect();

        id_delete
    }

    pub fn draw(&self, draw: &Draw, color: &Hsv, run_options: &RunOptions) {

        for b in self.boids.iter() {
            b.display(&draw, color, run_options);
        }
    }
}
