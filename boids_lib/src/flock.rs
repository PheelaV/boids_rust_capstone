use glam::Vec2;
use linfa::DatasetBase;
use linfa::dataset::Labels;
use linfa::traits::Transformer;
use linfa_clustering::Dbscan;
use ndarray::Array2;
// use nannou::prelude::vec2;
use rand::Rng;

use crate::boid::Boid;
use crate::math_helpers::distance_dyn_boid;
use crate::options::InitiationStrategy;
use crate::options::RunOptions;

fn get_boids(run_options: &RunOptions) -> Vec<Boid> {

    (0..run_options.init_boids).map(|id|{
        get_boid(run_options, id)
    })
    .collect()
}

fn get_boid(run_options: &RunOptions, id: u32) -> Boid {
    let mut rng = rand::thread_rng();

    match run_options.initiation_strat {
        InitiationStrategy::CircleCenterOut => todo!(),
        InitiationStrategy::CircleCircumferenceIn => {
            let r = run_options.window.win_h / 2.0;

            let init_vel: f32 = rng.gen::<f32>() / 3. + 0.5;
            let init_pos: f32 = rng.gen::<f32>() / 3. + 2. / 3.;
            let deg = rng.gen::<f32>() * 360.0;
    
            let x = r * deg.cos();
            let y = r * deg.sin();
    
            let mut init_vec = Vec2::new(-x, -y);
            init_vec = init_vec.clamp_length_max(0.5 * init_vel);
    
            Boid::new(x * init_pos, y * init_pos, init_vec, id)
        },
        InitiationStrategy::TwoWalls => todo!(),
        InitiationStrategy::RectangleIn => todo!(),
        InitiationStrategy::RandomIn => todo!(),
        InitiationStrategy::RandomRandom => {
            // x in [-width/2, width/2)
            let x = (rng.gen::<f32>() - 0.5) * run_options.window.win_w as f32;
            // y in [-height/2, height/2)
            let y = (rng.gen::<f32>() - 0.5) * run_options.window.win_h as f32;
            // x in [-max_speed, max_speed)
            let x_vel = (rng.gen::<f32>() * 2. - 1.) * run_options.max_speed;
            // y in [-max_speed, max_speed)
            let y_vel = (rng.gen::<f32>() * 2. - 1.) * run_options.max_speed;

            let mut init_vec = Vec2::new(x_vel, y_vel);

            init_vec = init_vec.clamp_length(run_options.min_speed, run_options.max_speed);

            Boid::new(x, y, init_vec, id)
        },
    }
}

pub struct Flock {
    pub boids: Vec<Boid>,
}

impl Flock {
    pub fn new(run_options: &RunOptions) -> Self {
        let boids = get_boids(&run_options);

        Flock { 
            boids,
        }
    }

    pub fn update(&mut self, run_options: &RunOptions) {
        let mut clone = self.boids.to_owned();

        let test_data: Array2<f32> = clone.iter()
        .map(|row| [row.position.x, row.position.y])
        .collect::<Vec<_>>()
        .into();
        // let dataset: DatasetBase<_, _> = generate::blobs(20, &expected_centroids, &mut rng).into();
        let test_dataset: DatasetBase<_, _> = test_data.into();
        // let expected_centroids = array![[10., 10.], [1., 12.], [20., 30.], [-20., 30.],];        
        let res =  Dbscan::params(3)
        .tolerance(run_options.sensory_distance)
        .transform(test_dataset)
        .unwrap();

        self.boids = res.targets.iter().zip(clone.iter_mut()).map(|pair| {
            match pair.0 {
                Some(cluster) => {
                    pair.1.cluster_id = *cluster + 1;
                    pair.1.to_owned()
                },
                None => {
                    pair.1.cluster_id = 0;
                    pair.1.to_owned()
                },
            }
        }).collect::<Vec<Boid>>();

        // println!("clusters: {count:2.}", count = res.label_count().len());
        let label_count = res.label_count().remove(0);

        println!();
        println!("Result: ");
        for (label, count) in label_count {
            match label {
                None => println!(" - {} noise points", count),
                Some(i) => println!(" - {} points in cluster {}", count, i),
            }
        }

        let max_sensory_distance = 
        run_options.allignment_treshold_distance
        .max(run_options.cohesion_treshold_distance
            .max(run_options.separation_treshold_distance));

        for b_cur in self.boids.iter_mut() {
            let mut neighbours: Vec<&Boid> = Vec::new();

            for b_other in clone.iter() {
                if b_other.id == b_cur.id {
                    continue;
                }

                let distance = distance_dyn_boid(b_cur, b_other, &run_options);
                if distance < max_sensory_distance {
                    neighbours.push(b_other);
                }
            }

            b_cur.run_rules(neighbours, &run_options);
            b_cur.update_location(&run_options)
        }
    }

    pub fn insert(&mut self, run_options: &RunOptions){
        let b = get_boid(run_options, self.boids.len() as u32);
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

    pub fn restart(&mut self, run_options: &RunOptions) -> (){
        self.boids = get_boids(&run_options);
    }
}
