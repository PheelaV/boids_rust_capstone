use std::{collections::HashSet, f32::consts::PI};

use glam::Vec2;
use rand::Rng;

use crate::{
    boid::{Boid, BoidMetadata},
    math_helpers::distance_dyn_boid,
    options::RunOptions,
};

use super::{get_flock_ids, tracker::Tracker, MY_RNG};

/// A naive implementation of boids tracking, which uses an O(N^2) algorithm for
/// finding boid's neighbours.
pub struct NaiveTracker {
    boids: Vec<Boid>,
    metadata: Vec<BoidMetadata>,
}


// TODO: naive tracker is behind at the moment in terms of noise and in terms of the acceleration update
impl NaiveTracker {
    // pub(self) fn get_neighbours<'a>(
    //     &'a self,
    //     boid: &Boid,
    //     run_options: &RunOptions,
    //     neighbours: &mut Vec<&'a Boid>
    // ) {
    //     for b_other in self.boids.iter() {
    //         if b_other.id == boid.id {
    //             continue;
    //         }

    //         let distance = distance_dyn_boid(&boid, b_other, &run_options);
    //         if distance < run_options.max_sensory_distance {
    //             neighbours.push(b_other);
    //         }
    //     }
    // }

    pub fn get_neighbours_naive<'a>(
        boid: &Boid,
        all_boids: &'a Vec<Boid>,
        run_options: &RunOptions,
        neighbours: &mut Vec<&'a Boid>,
    ) {
        for b_other in all_boids.iter() {
            if b_other.id == boid.id {
                continue;
            }

            let distance = distance_dyn_boid(&boid, b_other, &run_options);
            if distance < run_options.max_sensory_distance {
                neighbours.push(b_other);
            }
        }
    }
}

impl Tracker for NaiveTracker {
    fn new(entities: &[Boid], run_options: &RunOptions) -> Self {
        NaiveTracker {
            boids: entities.to_vec(),
            metadata: vec![Default::default(); run_options.init_boids],
        }
    }

    fn get_neighbours<'a>(&'a self, boid: &Boid, run_options: &RunOptions) -> Vec<&'a Boid> {
        let mut res = Vec::<&'a Boid>::new();

        NaiveTracker::get_neighbours_naive(boid, &self.boids, run_options, &mut res);

        res
    }

    fn update(&mut self, run_options: &RunOptions) {
        let mut accelleration: Vec<Vec2> = Vec::with_capacity(self.boids.len());
        let mut metadata: Vec<BoidMetadata> = vec![Default::default(); self.boids.len()];
        let mut neighbours: Vec<&Boid> = Vec::new();
        // let mut clicked_neighbours: Vec<&Boid> =
        //     Vec::with_capacity(run_options.neighbours_cosidered);

        // calculation loop
        for i_cur in 0..self.boids.len() {
            neighbours.clear();

            let b_current = &self.boids[i_cur];

            metadata[i_cur].id = b_current.id;
                        // this should really be in one of the boid rules,
            // but then there goes borrow checker and boids having
            // to mutate their own state for which the current pipeline
            // simply is not built for, so here is a workaround
            if run_options.wander_on {
                let wander_next = MY_RNG
                    .lock()
                    .unwrap()
                    .gen_range(-(run_options.wander_rate/2.)..(run_options.wander_rate/2.));
                metadata[i_cur].wander_direction = match run_options.noise_model {
                    crate::options::NoiseModel::Vicsek => wander_next,
                    crate::options::NoiseModel::Reynolds => (metadata[i_cur].wander_direction + wander_next) % (2. * PI),
                }
                    
            }

            NaiveTracker::get_neighbours_naive(b_current, &self.boids, run_options, &mut neighbours);
            accelleration.push(self.boids[i_cur].run_rules(&neighbours, &metadata, &run_options));
            // for cn in &clicked_neighbours {
            //     metadata[cn.id].clicked_neighbour_id = run_options.clicked_boid_id;
            // }
            if metadata[i_cur].id == run_options.clicked_boid_id {
                for cn in neighbours.iter() {
                    metadata[cn.id].clicked_neighbour_id = metadata[i_cur].id;
               }
            }
        }
        // update loop
        for i_cur in 0..self.boids.len() {
            // update metadata
            self.metadata[i_cur].id = i_cur;
            self.metadata[i_cur].cluster_id = Default::default();
            self.metadata[i_cur].clicked_neighbour_id = metadata[i_cur].clicked_neighbour_id;
            self.metadata[i_cur].n_neighbours = metadata[i_cur].n_neighbours;

            // progress
            self.boids[i_cur].apply_force(accelleration[i_cur]);
            self.boids[i_cur].update_location(&run_options)
        }

        // also try https://docs.rs/cogset/latest/cogset/struct.Optics.html
        if run_options.dbscan_flock_clustering_on {
            let flock_ids = get_flock_ids(self, &self.boids, run_options);

            self.boids
                .iter()
                .enumerate()
                .for_each(|(index, boid)| self.metadata[boid.id].cluster_id = flock_ids[index])
        }
    }

    fn insert_single(&mut self, entity: Boid, _: &RunOptions) {
        let mut metadata: BoidMetadata = Default::default();
        metadata.id = entity.id;
        self.metadata.push(metadata);

        self.boids.push(entity);
    }

    fn restart(&mut self, entities: &[Boid]) {
        self.boids = entities.to_vec();
    }

    fn delete_last(&mut self) -> Option<Boid> {
        if self.boids.len() == 0 {
            None
        } else {
            let index = self.boids.len() - 1;
            self.metadata.swap_remove(index);
            Some(self.boids.swap_remove(index))
        }
    }

    fn delete_multiple(&mut self, ids_delete: &[usize], _: &RunOptions) {
        let ids_delete_set: HashSet<usize> = HashSet::from_iter(ids_delete.iter().cloned());

        // removing values from the 1D table
        self.boids.retain(|boid| !ids_delete_set.contains(&boid.id));

        self.metadata.resize(self.boids.len(), Default::default())
    }

    fn insert_multiple(&mut self, entities: &[Boid], _: &RunOptions) {
        self.boids.extend_from_slice(entities);

        self.metadata.extend(entities.iter().map(|e| {
            let mut metadata: BoidMetadata = Default::default();
            metadata.id = e.id;
            metadata
        }))
    }

    fn view(&self) -> (&Vec<Boid>, &Vec<BoidMetadata>) {
        (&&self.boids, &self.metadata)
    }

    // return a map iterator
    // https://stackoverflow.com/questions/31904842/return-a-map-iterator-which-is-using-a-closure-in-rust
    fn view2<'a>(&'a self) -> Box<dyn Iterator<Item = (&'a Boid, &'a BoidMetadata)> + 'a> {
        // self.table.iter().map
        Box::new(self.boids.iter().map(|e| (e, &self.metadata[e.id])))
    }

    fn signal(&mut self, signal: super::tracker::TrackerSignal) {
        () // noop
    }
}
