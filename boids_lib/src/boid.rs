use std::{f32::consts::E, fmt::Debug};

use glam::f32::Vec2;
use rand::Rng;

use crate::{
    math_helpers::{distance_dyn_boid, tor_vec, MyRotate},
    options::{Boundary, Distance, RunOptions},
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoidType {
    Mob,
    Disruptor,
    Ghost,
}

#[derive(Debug, Clone, Copy)]
pub struct BoidMetadata {
    pub id: usize,
    pub cluster_id: usize,
    pub clicked_neighbour_id: usize,
    pub n_neighbours: usize,
    pub accelleration_update: Vec2,
    pub boid_type: BoidType,
    pub wander_direction: f32,
    pub wander_next: f32,
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
        Self {
            id: std::usize::MAX,
            clicked_neighbour_id: std::usize::MAX,
            cluster_id: 0,
            n_neighbours: 0,
            accelleration_update: Default::default(),
            boid_type: BoidType::Mob,
            wander_direction: 0.,
            wander_next: 0.,
        }
    }
}

#[derive(Clone, Copy, Default)]
pub struct Boid {
    // sequential id starting from 0
    pub id: usize,
    pub position: Vec2,
    pub velocity: Vec2,
    acceleration: Vec2,
    // pub wander_direction: f32
    // pub n_neighbours: i32,
    // pub cluster_id: usize
}

impl Debug for Boid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Boid")
            .field("id", &self.id)
            .field("position", &self.position)
            .finish()
    }
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
            // wander_direction: 0.,
            // n_neighbours: 0,
            // cluster_id: 0,
        }
    }

    // TODO: action selection will be important as the number of rules grows
    pub fn run_rules(
        &self,
        nearest_boids: &Vec<&Boid>,
        metadata: &Vec<BoidMetadata>,
        run_options: &RunOptions,
    ) -> Vec2 {
        let mut sum = Vec2::ZERO;

        // if run_options.rules_impl {
        let filtered = if !run_options.field_of_vision_on {
            nearest_boids.to_owned()
        } else {
            self.filter_sight2(&nearest_boids, run_options)
        };

        if run_options.separation_on {
            // let sep = self.separation(&filtered, run_options);
            // if self.id == run_options.clicked_boid_id {
            //     println!("separation: {:?}, len: {:?}", sep, sep.length());
            // }
            // sum += sep
            sum += self.separation(&filtered, run_options);
            Boid::check_for_unruly_rule(sum, "separation");
        }

        if run_options.cohesion_on && (!run_options.rules_impl || sum == Vec2::ZERO) {
            // let cohes = self.cohesion(&filtered, run_options);
            // if self.id == run_options.clicked_boid_id {
            //     println!("cohesion: {:?}, len: {:?}", cohes, cohes.length());
            // }
            // sum += cohes;
            sum += self.cohesion(&filtered, run_options);
            Boid::check_for_unruly_rule(sum, "cohesion");
        }

        if run_options.alignment_on {
            let align = self.allignment(&filtered, run_options);
            if self.id == run_options.clicked_boid_id {
                println!("align: {:?}, len: {:?}", align, align.length());
            }
            sum += align;
            // sum += self.allignment(&filtered, run_options);
            Boid::check_for_unruly_rule(sum, "allignment");
        }

        if run_options.wander_on {
            let wander = self.wander(metadata, run_options);
            if self.id == run_options.clicked_boid_id {
                println!("wander: {:?}, len: {:?}", wander, wander.length());
            }
            sum += wander;
            // sum += self.wander(metadata, run_options);
            Boid::check_for_unruly_rule(sum, "wander");
        }

        // this is mostly for manual manipulation when testing, is not used in experiments atm
        if run_options.seek_target_on {
            match run_options.seek_location {
                Some(target) => sum = self.seek(target, run_options),
                None => (),
            }
        }
        // } else {
        // if run_options.separation_on {
        //     sum += self.separation(&self.filter_sight3(&nearest_boids, run_options, run_options.separation_fov_half_cos), run_options);
        // }

        // if run_options.cohesion_on {
        //     sum += self.cohesion(&self.filter_sight3(&nearest_boids, run_options, run_options.cohesion_fov_half_cos), run_options);
        // }

        // if run_options.alignment_on {
        //     sum += self.alignment(&self.filter_sight3(&nearest_boids, run_options, run_options.allignment_fov_half_cos), run_options);
        // }
        // let filtered = if !run_options.field_of_vision_on {
        //     nearest_boids.to_owned()
        // } else {
        //     self.filter_sight2(&nearest_boids, run_options)
        // };

        // if run_options.separation_on {
        //     sum += self.separation(&filtered, run_options);
        // }

        // if run_options.cohesion_on {
        //     sum += self.cohesion(&filtered, run_options);
        // }

        // if run_options.alignment_on {
        //     sum += self.alignment(&filtered, run_options);
        // }
        // }

        // if run_options.wander_on {
        //     sum += self.wander(&filtered, metadata, run_options);
        // }

        // if run_options.cohesion_on {
        //     sum += self.cohesion(&filtered, run_options);
        // }

        // if run_options.alignment_on {
        //     sum += self.alignment(&filtered, run_options);
        // }

        // match metadata[self.id].boid_type {
        //     BoidType::Mob => (),
        //     BoidType::Disruptor => {
        //         let disrupt = self.disrupt(&filtered, metadata,  run_options);
        //         sum += disrupt;
        //     },
        // }

        // dbg!("{:?}", sum);
        // sum = sum.
        // if self.id == run_options.clicked_boid_id {
        //     println!("sum: {:?}, len: {:?}", sum, sum.length());
        // }
        sum
    }

    pub fn filter_sight2<'a>(
        &self,
        others: &Vec<&'a Boid>,
        run_options: &RunOptions,
    ) -> Vec<&'a Boid> {
        let res: Vec<&Boid> = others
            .iter()
            .filter(|b_other| {
                if self.id == b_other.id {
                    return false;
                }

                let vec_to_other = match run_options.distance {
                    Distance::EucToroidal => {
                        tor_vec(self.position, b_other.position, &run_options.window)
                    }
                    Distance::EucEnclosed => b_other.position - self.position,
                };

                if vec_to_other.length() < 0.01 {
                    return false;
                }

                let vel_norm = self.velocity.normalize();
                let vec_to_other_norm = vec_to_other.normalize();

                // let rad_to_other = vel_norm.dot(vec_to_other_norm).acos();
                // this calculates v•u = |v||u|cos(ß), which is cos(ß) because of v and u being unit vectors
                vel_norm.dot(vec_to_other_norm) > run_options.field_of_vision_cos
                // if self.id == run_options.clicked_boid_id {
                //     println!("-----------");
                //     println!("dot: {:?}", vel_norm.dot(vec_to_other_norm));
                //     println!(
                //         "{:?}",
                //         vel_norm.dot(vec_to_other_norm) > run_options.field_of_vision_cos
                //     );
                //     println!(
                //         "fov cos: {:?}",
                //         run_options.field_of_vision_cos
                //     );
                //     println!("dot.acos: {:?}", rad_to_other);
                //     println!("{:?}", rad_to_other < run_options.field_of_vision_half_rad);
                //     println!("fov half rad: {:?}", run_options.field_of_vision_half_rad);
                // }
                // if self.id == run_options.clicked_boid_id {
                //     println!("cos: {:?}", vel_norm.dot(vec_to_other_norm));
                //     println!("acos: {:?}", rad_to_other)
                // }

                // rad_to_other < run_options.field_of_vision_half_rad
                // let atan2_self = self.velocity.y.atan2(self.velocity.x);
                // let atan2_other = vec_to_other.y.atan2(vec_to_other.x);
                // let mut atan2_diff = atan2_other - atan2_self;

                // // normalize
                // atan2_diff = if atan2_diff > PI {
                //     atan2_diff - 2. * PI
                // } else if atan2_diff <= -1. * PI {
                //     atan2_diff + 2. * PI
                // } else {
                //     atan2_diff
                // };

                // run_options.field_of_vision_on
                //     && atan2_diff.abs() < run_options.field_of_vision_half_rad
            })
            .map(|b| *b)
            .collect();

        res
    }

    pub fn filter_sight3<'a>(
        &self,
        others: &Vec<&'a Boid>,
        run_options: &RunOptions,
        half_cos_treshold: f32,
    ) -> Vec<&'a Boid> {
        let res: Vec<&Boid> = others
            .iter()
            .filter(|b_other| {
                if self.id == b_other.id {
                    return false;
                }

                let vec_to_other = match run_options.distance {
                    Distance::EucToroidal => {
                        tor_vec(self.position, b_other.position, &run_options.window)
                    }
                    Distance::EucEnclosed => b_other.position - self.position,
                };

                if vec_to_other.length() < 0.01 {
                    return false;
                }

                let vel_norm = self.velocity.normalize();
                let vec_to_other_norm = vec_to_other.normalize();

                // let rad_to_other = vel_norm.dot(vec_to_other_norm).acos();
                // this calculates v•u = |v||u|cos(ß), which is cos(ß) because of v and u being unit vectors
                vel_norm.dot(vec_to_other_norm) > half_cos_treshold
            })
            .map(|b| *b)
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

                    let value = match run_options.distance {
                        Distance::EucToroidal => {
                            tor_vec(other.position, self.position, &run_options.window)
                        }
                        Distance::EucEnclosed => self.position - other.position,
                    };
                    
                    // res += value;
                    
                    if !run_options.separation_impl_mode {
                        // value /= distance;
                        // value /= distance;
                        // res += value * Self::sep_scale(distance, run_options)
                        res += value.normalize() / distance;
                    } else {
                        res += value.normalize() * Self::scale_sigmoid(distance, run_options);
                    }

                    // if !run_options.separation_impl_mode {
                    //     res += ((self.position - other.position) / distance).normalize()
                    // } else {
                    //     res += (self.position - other.position).normalize() / distance
                    // }
                    // if !run_options.separation_impl_mode {
                    //     res += (self.position - other.position).normalize() / distance.powi(2)
                    // } else {
                    //     res += (self.position - other.position).normalize() / distance
                    // }
                }
            }

            // if run_options.separation_impl_mode {
            //     res /= count as f32;
            // }

            // res *= run_options.separation_coefficient;
        }

        if count > 0 {
            res /= count as f32;
            res = res.normalize();
            // res.normalize() * run_options.separation_coefficient;
            res *= run_options.max_speed;
            res -= self.velocity;
            res.clamp_length_max(run_options.max_steering);
            res * run_options.separation_coefficient
        } else {
            Vec2::ZERO
        }
    }

    fn scale_fall_off(distance: f32, run_options: &RunOptions) -> f32 {
        let dist_ration = distance / run_options.separation_treshold_distance;

        // 0 < beta <= 1
        // 0 < alpha < 1
        // produces a curve that assymptotically touches y = 1, then falls off to intercept x = 1 at 1 - beta

        // i.e. soft transition from y = 1 with a sharp turn as we approach x = 1

        const ALPHA: f32 = 4.;
        const BETA: f32 = 1.;

        1. - BETA * E.powf(ALPHA * dist_ration - ALPHA)
    }

    fn scale_sigmoid(distance: f32, run_options: &RunOptions) -> f32 {
        let dist_ration = distance / run_options.separation_treshold_distance;
        const ALPHA: f32 = 15.;
        const BETA: f32 = 0.6;

        1. - (1. / ( 1. + E.powf(-ALPHA * (dist_ration - BETA))))
    }

    pub fn cohesion(&self, others: &Vec<&Boid>, run_options: &RunOptions) -> Vec2 {
        let mut center = Vec2::new(0.0, 0.0);
        let mut count = 0;

        for other in others {
            let distance = distance_dyn_boid(self, other, &run_options);
            if distance > 0. && distance < run_options.cohesion_treshold_distance {
                center += match run_options.distance {
                    Distance::EucToroidal => {
                        tor_vec(self.position, other.position, &run_options.window)
                    }
                    Distance::EucEnclosed => other.position,
                };
                count += 1;
            }
        }

        if count > 0 {
            center /= count as f32;
            // if run_options.distance == Distance::EucToroidal {
            //     center += self.position;
            // }
            match run_options.distance {
                Distance::EucToroidal => {
                    self.steer(center, run_options) * run_options.cohesion_coefficient
                }
                Distance::EucEnclosed => {
                    self.steer(center - self.position, run_options)
                        * run_options.cohesion_coefficient
                }
            }
            // self.steer(center, run_options) * run_options.cohesion_coefficient
        } else {
            Vec2::ZERO
        }
    }

    pub fn seek(&self, target: Vec2, run_options: &RunOptions) -> Vec2 {
        let desired = match run_options.distance {
            Distance::EucToroidal => tor_vec(self.position, target, &run_options.window),
            Distance::EucEnclosed => target - self.position,
        };

        self.steer(desired, run_options)
    }

    pub fn steer(&self, mut desired: Vec2, run_options: &RunOptions) -> Vec2 {
        if desired.length() == 0. {
            Vec2::new(0., 0.)
        } else {
            desired = desired.normalize();
            if !run_options.cohesion_impl_mode { 
                // this is used on cohesion, Reynolds does not do this (but has other mechanics in place later in his pipeline)
                desired *= run_options.max_speed;
                desired -= self.velocity;
                desired = desired.clamp_length_max(run_options.max_steering);
            }
            desired
        }
    }

    pub fn allignment(&self, others: &Vec<&Boid>, run_options: &RunOptions) -> Vec2 {
        let mut avg_vel = Vec2::ZERO;
        let mut count = 0.;

        for other in others {
            let distance = distance_dyn_boid(self, other, &run_options);
            if distance < run_options.allignment_treshold_distance {
                avg_vel += other.velocity;
                count += 1.;
            }
        }

        if count > 0. {
            avg_vel /= count;
            // if !run_options.allignment_impl_mode {
            //     avg_vel = avg_vel.normalize();
            //     avg_vel *= run_options.max_speed;

            // }
            // //  else {
            // // }
            // avg_vel = (avg_vel - self.velocity).normalize_or_zero();
            // avg_vel.clamp_length_max(run_options.max_steering);
            // avg_vel * run_options.allignment_coefficient
            avg_vel
        } else {
            Vec2::new(0.0, 0.0)
        }
    }

    pub fn disrupt(
        &self,
        others: &Vec<&Boid>,
        metadata: &Vec<BoidMetadata>,
        run_options: &RunOptions,
    ) -> Vec2 {
        let mut nearest: Option<&&Boid> = None;
        let mut nearest_distance = f32::MAX;

        for other in others {
            match metadata[other.id].boid_type {
                BoidType::Disruptor => break,
                _ => (),
            }

            let distance = distance_dyn_boid(self, other, &run_options);
            if distance < nearest_distance {
                nearest_distance = distance;
                nearest = Some(other);
            }
        }

        match nearest {
            Some(target) => (target.position - self.position) * 10.,
            None => Vec2::new(0., 0.),
        }
    }

    pub fn wander(&self, metadata: &Vec<BoidMetadata>, run_options: &RunOptions) -> Vec2 {
        match run_options.noise_model {
            crate::options::NoiseModel::Reynolds => {
                // the current velocity vector normalized
                let heading = self.velocity.normalize();
                // gives the center of the circle driving the locomotion
                // let loco_center = self.position +  heading * (run_options.wander_radius * (2_f32).sqrt());
                let loco_center =
                    self.position + heading * (run_options.wander_distance * (2_f32).sqrt());

                // vector pointing at the point on circumference
                let wander_point = heading.rotate(Vec2::new(
                    run_options.wander_radius * metadata[self.id].wander_direction.cos(),
                    run_options.wander_radius * metadata[self.id].wander_direction.sin(),
                ));

                // places the point onto the locomotion circle with respect to agent's location
                let wander_f = loco_center + wander_point;

                // steer towards the point on the locomotion circle
                self.seek(wander_f, run_options) * run_options.wander_coefficient
            }
            crate::options::NoiseModel::Vicsek => {
                let wander_vector = Vec2::new(
                    metadata[self.id].wander_direction.cos(),
                    metadata[self.id].wander_direction.sin(),
                );
                let heading = self.velocity.normalize().rotate(
                    wander_vector.clamp_length(run_options.max_speed, run_options.max_speed),
                );

                self.seek(self.position + heading, run_options) * run_options.wander_coefficient
            }
        }
    }

    // pub fn avoid(&self, x: f32, y:f32, run_options: &RunOptions) {
    //     let distance = distance_dyn(self.position.x, x, self.position.y, y, &run_options);

    //

    // Actually shifts the individual's location
    pub fn update_location(&mut self, run_options: &RunOptions) {
        // let vel_prior = self.velocity;
        self.velocity += self.acceleration;

        self.velocity = self.velocity.clamp_length_max(run_options.max_speed);

        // this is problematic
        if self.velocity.length() < run_options.min_speed {
            self.velocity = self.velocity.normalize() * run_options.min_speed;

            // This was extremely hard to spot and debug.
            // In this very unlikely case where agents have meat head on with equal
            // but opposite velocities and accelerations, the normalization above
            // produces NaN vector.

            // Presuming we have no notion of communication and social interaction:
            // Real life equivalent would be if you meet someone exactly head on, at the same speed,
            // and you slow down (at the same rate), to avoid collision, but then you both have to make
            // a decision, which way to go to avoid standing in place or bumbing into each other again
            // because it does not matter whether or not you CAN move(accelerate) what matters is to
            // start moving somewhere (e.g. make a random choice).
            // This is the resolution, arbitrary in this case:
            if self.velocity.is_nan() {
                let mut rng = rand::thread_rng();
                let x_vel = (rng.gen::<f32>() * 2. - 1.) * run_options.max_speed;
                let y_vel = (rng.gen::<f32>() * 2. - 1.) * run_options.max_speed;
                self.velocity.x = x_vel;
                self.velocity.y = y_vel;
            }
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

    fn check_for_unruly_rule(v: Vec2, rule_name: &str) {
        if v.x.is_infinite() || v.x.is_nan() || v.x.is_infinite() || v.x.is_nan() {
            panic!("the faulty rule is: {}, the faulty bits are: x fin: {}, x nan: {}, y fin: {}, y nan; {}", rule_name, v.x.is_infinite(), v.x.is_nan(), v.x.is_infinite(), v.x.is_nan())
        }
    }

    fn boundaries(&mut self, run_options: &RunOptions) {
        match run_options.boundary {
            Boundary::Toroidal => {
                if self.position.x < run_options.window.win_left as f32 {
                    self.position.x = run_options.window.win_right as f32 + self.position.x
                        - run_options.window.win_left as f32;
                } else if self.position.x > run_options.window.win_right as f32 {
                    self.position.x = run_options.window.win_left as f32 + self.position.x
                        - run_options.window.win_right as f32;
                }

                if self.position.y > run_options.window.win_top as f32 {
                    self.position.y = run_options.window.win_bottom as f32 + self.position.y
                        - run_options.window.win_top as f32;
                } else if self.position.y < run_options.window.win_bottom as f32 {
                    self.position.y = run_options.window.win_top as f32 + self.position.y
                        - run_options.window.win_bottom as f32;
                }
            }
            Boundary::Absorbing => {
                if self.position.x < run_options.window.win_left as f32 {
                    self.position.x = run_options.window.win_left as f32;
                    self.velocity.x = 0.;
                } else if self.position.x > run_options.window.win_right as f32 {
                    self.position.x = run_options.window.win_right as f32;
                    self.velocity.x = 0.;
                }

                if self.position.y > run_options.window.win_top as f32 {
                    self.position.y = run_options.window.win_top as f32;
                    self.velocity.y = 0.;
                } else if self.position.y < run_options.window.win_bottom as f32 {
                    self.position.y = run_options.window.win_bottom as f32;
                    self.velocity.y = 0.;
                }
            }
            Boundary::Reflective => {
                // flip velocity if it is going beyond the edge
                if (self.position.x < run_options.window.win_left as f32 && self.velocity.x < 0.)
                    || (self.position.x > run_options.window.win_right as f32
                        && self.velocity.x > 0.)
                {
                    self.velocity.x = -self.velocity.x;
                }

                if (self.position.y > run_options.window.win_top as f32 && self.velocity.y > 0.)
                    || (self.position.y < run_options.window.win_bottom as f32
                        && self.velocity.y < 0.)
                {
                    self.velocity.y = -self.velocity.y;
                }
            }
            Boundary::Repulsive { distance, force } => {
                // add accelleration as a vector pointing away from the barrier
                if self.position.x < run_options.window.win_left as f32 + distance {
                    self.acceleration.x += force;
                } else if self.position.x > run_options.window.win_right as f32 - distance {
                    self.acceleration.x -= force;
                }

                if self.position.y > run_options.window.win_top as f32 - distance {
                    self.acceleration.y -= force;
                } else if self.position.y < run_options.window.win_bottom as f32 + distance {
                    self.acceleration.y += force;
                }
            }
            Boundary::RepulsiveCircle { radius: _ } => todo!(),
        };
    }
}
