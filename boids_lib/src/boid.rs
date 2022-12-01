

// use nannou::prelude::{Vec2, vec2, PI};

use std::{f32::consts::PI};

use glam::Vec2;

use crate::{
    math_helpers::distance_dyn_boid,
    options::{
        RunOptions, 
        Boundary
    },
};


#[derive(Debug, Clone)]
pub struct BoidMetadata {
    pub id: usize,
    pub cluster_id: usize,
    pub clicked_neighbour_id: usize,
    pub n_neighbours: usize
}

impl BoidMetadata {
    pub fn new(boid: &Boid) -> Self {
        let mut default: BoidMetadata = Default::default();
        default.id = boid.id;
        default
    }
}

impl Default for BoidMetadata {
    fn default() -> Self {
        Self { id: std::usize::MAX, clicked_neighbour_id: std::usize::MAX, cluster_id: 0, n_neighbours: 0 }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Boid {
    // sequential id starting from 0
    pub id: usize,
    pub position: Vec2,
    pub velocity: Vec2,
    acceleration: Vec2,
    // pub n_neighbours: i32,
    // pub cluster_id: usize
}

impl Boid {
    /// Creates a new [`Boid`].
    pub fn new(x: f32, y: f32, velocity: Vec2, id: usize) -> Self {
        let acceleration = Vec2::ZERO;
        let position = Vec2::new(x, y);

        Boid {
            id,
            position,
            velocity,
            acceleration,
            // n_neighbours: 0,
            // cluster_id: 0,
        }
    }

    pub fn run_rules(& self, nearest_boids: &Vec<&Boid>, run_options: &RunOptions) -> Vec2 {
        let mut sum = Vec2::ZERO;
        let filtered = if !run_options.field_of_vision_on 
            { nearest_boids.to_owned() }  // todo: watchout for the copy here
        else 
            { self.filter_sight2(&nearest_boids, run_options)};

        // self.n_neighbours = filtered.len() as i32;

        if run_options.separation_on {
            sum += self.separation(&filtered, run_options);
        }

        if run_options.cohesion_on {
            sum += self.cohesion(&filtered, run_options);
        }

        if run_options.alignment_on {
            sum += self.alignment(&filtered, run_options);
        }

        sum
    }

    pub fn filter_sight<'a>(&self, others: &[&'a Boid], run_options: &RunOptions) -> Vec<&'a Boid> {
        let res: Vec<&Boid> = others.iter()
        .filter(|b_other| {
            if self.id == b_other.id {
                return false; 
            }
            
            let vec_to_other = b_other.position - self.position;
            
            if vec_to_other.length() < 0.01 {
                return false;
            }
            
            let atan2_self = self.velocity.y.atan2(self.velocity.x);
            let atan2_other = vec_to_other.y.atan2(vec_to_other.x);
            let mut atan2_diff = atan2_other - atan2_self;
            
            // normalize
            atan2_diff = 
            if atan2_diff > PI { 
                atan2_diff - 2.*PI
            } else if atan2_diff <= -1.* PI{
                atan2_diff + 2.*PI
            } else {atan2_diff};

            run_options.field_of_vision_on && atan2_diff.abs() < run_options.field_of_vision_half_rad  
          })
          .map(|b| {*b})
        .collect();

        res
    }

    pub fn filter_sight2<'a>(&self, others: &Vec<&'a Boid>, run_options: &RunOptions) -> Vec<&'a Boid> {
        let res: Vec<&Boid> = others.iter()
        .filter(|b_other| {
            if self.id == b_other.id {
                return false; 
            }
            
            let vec_to_other = b_other.position - self.position;
            
            if vec_to_other.length() < 0.01 {
                return false;
            }
            
            let atan2_self = self.velocity.y.atan2(self.velocity.x);
            let atan2_other = vec_to_other.y.atan2(vec_to_other.x);
            let mut atan2_diff = atan2_other - atan2_self;
            
            // normalize
            atan2_diff = 
            if atan2_diff > PI { 
                atan2_diff - 2.*PI
            } else if atan2_diff <= -1.* PI{
                atan2_diff + 2.*PI
            } else {atan2_diff};

            run_options.field_of_vision_on && atan2_diff.abs() < run_options.field_of_vision_half_rad  
          })
          .map(|b| {*b})
        .collect();

        res
    }

    pub fn separation(&self, others: &Vec<&Boid>, run_options: &RunOptions) -> Vec2 {
        let mut res = Vec2::new(0.0, 0.0);
        let mut count = 0;
        if others.len() != 0 {
            for other in others {
                let distance = distance_dyn_boid(self, other, &run_options);
                if distance < run_options.separation_treshold_distance {
                    count += 1;
                    if !run_options.separation_impl_mode {
                        res += (self.position - other.position).normalize() / distance.powi(2)
                    } else {
                        res += (self.position - other.position).normalize() / distance
                    }
                }
            }

            // if run_options.separation_impl_mode {
            //     res /= count as f32;
            // }

            res *= run_options.separation_coefficient;
        }

        if count > 0 {
            res.clamp_length_max(run_options.max_steering)
        } else {
            Vec2::ZERO
        }
    }

    pub fn cohesion(&self, others: &Vec<&Boid>, run_options: &RunOptions) -> Vec2 {
        let mut center = Vec2::new(0.0, 0.0);
        let mut count = 0;

        for other in others {
            let distance = distance_dyn_boid(self, other, &run_options);
            if distance > 0. && distance < run_options.cohesion_treshold_distance {
                center += other.position;
                // center += vec2(other.pos_x, other.pos_y);
                count += 1;
            }
        }

        if count > 0 {
            center /= count as f32;
            self.seek(center, run_options) * run_options.cohesion_coefficient
        } else {
            Vec2::ZERO
        }
    }

    pub fn seek(&self, target: Vec2, run_options: &RunOptions) -> Vec2 {
        let mut desired = target - self.position;        

        if desired.length() == 0. {
            Vec2::new(0., 0.)
        } else {
            if run_options.cohesion_impl_mode {
                desired = desired.normalize();
                desired *= run_options.max_speed;
            }
            desired = desired.clamp_length_max(run_options.max_steering);
            desired
        }
    }

    pub fn alignment(&self, others: &Vec<&Boid>, run_options: &RunOptions) -> Vec2 {
        let mut avg = Vec2::ZERO;
        let mut count = 0.;

        for other in others {
            let distance = distance_dyn_boid(self, other, &run_options);
            if distance < run_options.allignment_treshold_distance {
                avg += other.velocity;
                count += 1.;
            }
        }

        if count > 0. {
            avg /= count;
            avg = avg.normalize();
            avg *= run_options.max_speed;
            avg = (avg - self.velocity) * run_options.allignment_coefficient;
            avg.clamp_length_max(run_options.max_steering)
        } else {
            Vec2::new(0.0, 0.0)
        }
    }

    // pub fn avoid(&self, x: f32, y:f32, run_options: &RunOptions) {
    //     let distance = distance_dyn(self.position.x, x, self.position.y, y, &run_options);

    // 

    // Actually shifts the individual's location
    pub fn update_location(&mut self, run_options: &RunOptions) {
        self.velocity += self.acceleration;
        self.velocity = self.velocity.clamp_length_max(run_options.max_speed);

        if self.velocity.length() < run_options.min_speed {
            self.velocity = self.velocity.normalize() * run_options.min_speed;
        }
        
        if !run_options.stop_movement {
            self.position += self.velocity * run_options.baseline_speed;
        }

        self.acceleration *= 0.0;

        self.boundaries(run_options)
    }

    pub fn apply_force(&mut self, force: Vec2) {
        self.acceleration += force;
    }

    fn boundaries(&mut self, run_options: &RunOptions) {
        match run_options.boundary {
            Boundary::Thoroidal => {
                // wrap around the barrier
                if self.position.x < run_options.window.win_left {
                    self.position.x = run_options.window.win_right;
                } else if self.position.x > run_options.window.win_right {
                    self.position.x = run_options.window.win_left;
                }

                if self.position.y > run_options.window.win_top {
                    self.position.y = run_options.window.win_bottom;
                } else if self.position.y < run_options.window.win_bottom {
                    self.position.y = run_options.window.win_top;
                }
            }
            Boundary::Absorbing => todo!(),
            Boundary::Reflective => {
                // flip velocity if it is going beyond the edge
                if (self.position.x < run_options.window.win_left && self.velocity.x < 0.) || (self.position.x > run_options.window.win_right && self.velocity.x > 0.) {
                    self.velocity.x = -self.velocity.x;
                }

                if (self.position.y > run_options.window.win_top && self.velocity.y > 0.) || (self.position.y < run_options.window.win_bottom && self.velocity.y < 0.) {
                    self.velocity.y = -self.velocity.y;
                }
            },
            Boundary::Repulsive { distance, force } => {
                // add accelleration as a vector pointing away from the barrier
                if self.position.x < run_options.window.win_left + distance {
                    self.acceleration.x += force;
                } else if self.position.x > run_options.window.win_right - distance {
                    self.acceleration.x -= force;
                }

                if self.position.y > run_options.window.win_top - distance {
                    self.acceleration.y -= force;
                } else if self.position.y < run_options.window.win_bottom + distance {
                    self.acceleration.y += force;
                }
            }
        };
    }
}
