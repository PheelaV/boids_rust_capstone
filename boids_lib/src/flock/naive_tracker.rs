use std::{any::Any, collections::HashSet, f32::consts::PI};

use glam::Vec2;
use rand::Rng;

use crate::{
    boid::{Boid, BoidMetadata, NeighborData},
    math_helpers::distance_and_direction_dyn_boid,
    options::RunOptions,
};

use super::{
    get_flock_ids,
    tracker::{Tracker, TrackerSignal},
    MY_RNG,
};

/// A naive implementation of boids tracking, which uses an O(N^2) algorithm for
/// finding boid's neighbours.
pub struct NaiveTracker {
    entities: Vec<Boid>,
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

    /// Returns neighbors with cached distance and direction data.
    /// Also performs FOV filtering if enabled.
    pub fn get_neighbours_naive_with_data<'a>(
        boid: &Boid,
        all_boids: &'a Vec<Boid>,
        run_options: &RunOptions,
        neighbours: &mut Vec<NeighborData<'a>>,
    ) {
        // Pre-compute velocity direction for FOV check (if enabled)
        let vel_norm = if run_options.field_of_vision_on {
            Some(boid.velocity.normalize())
        } else {
            None
        };

        for b_other in all_boids.iter() {
            if b_other.id == boid.id {
                continue;
            }

            let (distance, direction) = distance_and_direction_dyn_boid(boid, b_other, run_options);

            if distance >= run_options.max_sensory_distance {
                continue;
            }

            // FOV check using already-computed direction
            if let Some(ref vel_dir) = vel_norm {
                if vel_dir.dot(direction) <= run_options.field_of_vision_cos {
                    continue; // Outside field of view
                }
            }

            neighbours.push(NeighborData {
                boid: b_other,
                distance,
                direction,
            });

            if run_options.neighbours_cosidered != 0
                && neighbours.len() >= run_options.neighbours_cosidered
            {
                return; // Return early once limit is reached
            }
        }
    }
}

impl Tracker for NaiveTracker {
    fn new(entities: &[Boid], run_options: &RunOptions) -> Self {
        NaiveTracker {
            entities: entities.to_vec(),
            metadata: vec![Default::default(); run_options.init_boids],
        }
    }

    fn get_neighbours<'a>(&'a self, boid: &Boid, run_options: &RunOptions) -> Vec<NeighborData<'a>> {
        let mut res = Vec::<NeighborData<'a>>::new();

        NaiveTracker::get_neighbours_naive_with_data(boid, &self.entities, run_options, &mut res);

        res
    }

    fn update(&mut self, run_options: &RunOptions) {
        let mut accelleration: Vec<Vec2> = Vec::with_capacity(self.entities.len());
        let mut metadata: Vec<BoidMetadata> = vec![Default::default(); self.entities.len()];
        let mut neighbours: Vec<NeighborData> = Vec::new();
        // let mut clicked_neighbours: Vec<&Boid> =
        //     Vec::with_capacity(run_options.neighbours_cosidered);

        // calculation loop
        for i_cur in 0..self.entities.len() {
            neighbours.clear();

            let b_current = &self.entities[i_cur];

            metadata[i_cur].id = b_current.id;
            // this should really be in one of the boid rules,
            // but then there goes borrow checker and boids having
            // to mutate their own state for which the current pipeline
            // simply is not built for, so here is a workaround
            // the gist is: "Reynolds" noise is stateful and needs to write metadata
            if run_options.wander_on {
                let wander_next = MY_RNG
                    .lock()
                    .unwrap()
                    .gen_range(-(run_options.wander_rate)..(run_options.wander_rate))
                    * PI;
                metadata[i_cur].wander_direction = match run_options.noise_model {
                    crate::options::NoiseModel::Vicsek => wander_next,
                    crate::options::NoiseModel::Reynolds => {
                        (self.metadata[i_cur].wander_direction + wander_next) % (2. * PI)
                    }
                }
            }

            NaiveTracker::get_neighbours_naive_with_data(
                b_current,
                &self.entities,
                run_options,
                &mut neighbours,
            );
            // FOV filtering is now done inside get_neighbours_naive_with_data
            accelleration.push(self.entities[i_cur].run_rules(
                &neighbours,
                &metadata,
                &run_options,
            ));
            // for cn in &clicked_neighbours {
            //     metadata[cn.id].clicked_neighbour_id = run_options.clicked_boid_id;
            // }
            if metadata[i_cur].id == run_options.clicked_boid_id {
                for cn in neighbours.iter() {
                    metadata[cn.boid.id].clicked_neighbour_id = metadata[i_cur].id;
                }
            }
        }
        // update loop
        for i_cur in 0..self.entities.len() {
            // update metadata
            self.metadata[i_cur].id = i_cur;
            self.metadata[i_cur].cluster_id = Default::default();
            self.metadata[i_cur].clicked_neighbour_id = metadata[i_cur].clicked_neighbour_id;
            self.metadata[i_cur].n_neighbours = metadata[i_cur].n_neighbours;
            self.metadata[i_cur].wander_direction = metadata[i_cur].wander_direction;

            // progress
            self.entities[i_cur].apply_force(accelleration[i_cur]);
            self.entities[i_cur].update_location(&run_options)
        }

        // also try https://docs.rs/cogset/latest/cogset/struct.Optics.html
        if run_options.dbscan_flock_clustering_on {
            let flock_ids = get_flock_ids(self, &self.entities, run_options);

            self.entities
                .iter()
                .enumerate()
                .for_each(|(index, boid)| self.metadata[boid.id].cluster_id = flock_ids[index])
        }
    }

    fn insert_single(&mut self, entity: Boid, _: &RunOptions) {
        let mut metadata: BoidMetadata = Default::default();
        metadata.id = entity.id;
        self.metadata.push(metadata);

        self.entities.push(entity);
    }

    fn restart(&mut self, entities: &[Boid]) {
        self.entities = entities.to_vec();
    }

    fn delete_last(&mut self, _run_options: &RunOptions) -> Option<Boid> {
        if self.entities.len() == 0 {
            None
        } else {
            let index = self.entities.len() - 1;
            self.metadata.swap_remove(index);
            Some(self.entities.swap_remove(index))
        }
    }

    fn delete_multiple(&mut self, ids_delete: &[usize], _: &RunOptions) {
        let ids_delete_set: HashSet<usize> = HashSet::from_iter(ids_delete.iter().cloned());

        // removing values from the 1D table
        self.entities
            .retain(|boid| !ids_delete_set.contains(&boid.id));

        self.metadata
            .resize(self.entities.len(), Default::default())
    }

    fn insert_multiple(&mut self, entities: &[Boid], _: &RunOptions) {
        self.entities.extend_from_slice(entities);

        self.metadata.extend(entities.iter().map(|e| {
            let mut metadata: BoidMetadata = Default::default();
            metadata.id = e.id;
            metadata
        }))
    }

    // return a map iterator
    // https://stackoverflow.com/questions/31904842/return-a-map-iterator-which-is-using-a-closure-in-rust
    fn view<'a>(&'a self) -> Box<dyn Iterator<Item = (&'a Boid, &'a BoidMetadata)> + 'a> {
        // self.table.iter().map
        Box::new(self.entities.iter().map(|e| (e, &self.metadata[e.id])))
    }

    fn signal(&mut self, _: TrackerSignal) {
        // noop
    }

    fn get_no_entities(&self) -> usize {
        self.entities.len()
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}
