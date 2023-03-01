use std::collections::HashMap;
use std::collections::HashSet;
use std::sync::Mutex;

use glam::Vec2;
use itertools::Itertools;
use linfa::traits::Transformer;
use linfa::DatasetBase;
use linfa_clustering::Dbscan;
use linfa_nn::distance::L2Dist;
use linfa_nn::CommonNearestNeighbour;
// use linfa_clustering::Optics;
use ndarray::Array2;
use petal_clustering::{
    //  Dbscan as pDbscan,
    Fit as pFit,
    Optics as pOptics,
};
use petal_neighbors::distance::Euclidean;
use rand::Rng;

use lazy_static::lazy_static;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro128Plus;
lazy_static! {
    static ref MY_RNG: Mutex<Xoshiro128Plus> = Mutex::new(Xoshiro128Plus::from_entropy());
}

use crate::boid::Boid;
use crate::boid::BoidMetadata;
use crate::math_helpers::distance_dyn_boid;
use crate::options::Boundary;
use crate::options::InitiationStrategy;
use crate::options::RunOptions;

use self::spathash_tracker::SpatHash1D;
use self::tracker::Tracker;

pub mod naive_tracker;
pub mod spathash_tracker;
pub mod tracker;

// pub mod tracker;

// pub use tracker::*;
// type BoidTracker = SpatHash1D;
// type BoidTracker = BoidTracker;

// type FlockTracker = BoidTracker;
type FlockTracker = SpatHash1D;
pub struct Flock {
    pub tracker: FlockTracker,
}

impl Flock {
    // const GREY_SCALE: &str =
    //     " .\'`^\",:;Il!i><~+_-?][}{1)(|\\/tfjrxnuvczXYUJCLQ0OZmwqpdbkhao*#MW&8%B@$";
    pub fn new(run_options: &RunOptions) -> Self {
        let boids = get_boids(&run_options);

        Flock {
            tracker: FlockTracker::new(&boids, run_options),
            // tracker: BoidTracker::new(&boids, run_options),
        }
    }

    pub fn get_random() {
        todo!()
    }

    pub fn view(&self) -> (&Vec<Boid>, &Vec<BoidMetadata>) {
        self.tracker.view()
    }

    pub fn view2<'a>(&'a self) -> Box<dyn Iterator<Item = (&'a Boid, &'a BoidMetadata)> + 'a> {
        self.tracker.view2()
    }

    pub fn update(&mut self, run_options: &RunOptions) {
        self.tracker.update(run_options);
    }

    pub fn insert_single(&mut self, run_options: &RunOptions) {
        let b = get_boid(run_options, self.tracker.view().0.len());
        self.tracker.insert_single(b, run_options);
    }

    pub fn insert_multiple(&mut self, run_options: &RunOptions) {
        let no_boids = self.tracker.view().0.len();
        let boid_difference = run_options.init_boids - no_boids;
        let boids = get_boids_n(boid_difference, no_boids, run_options);
        self.tracker.insert_multiple(&boids, run_options);
    }

    pub fn delete_last(&mut self) -> Option<Boid> {
        self.tracker.delete_last()
    }

    pub fn delete_boid(&mut self, id_delete: usize, run_options: &RunOptions) {
        self.tracker.delete_multiple(&[id_delete], run_options)
    }

    pub fn delete_multiple(&mut self, ids_delete: &[usize], run_options: &RunOptions) {
        self.tracker.delete_multiple(ids_delete, run_options)
    }

    pub fn restart(&mut self, run_options: &RunOptions) -> () {
        let boids = get_boids(run_options);
        self.tracker.restart(&boids);
    }
}

fn get_boids(run_options: &RunOptions) -> Vec<Boid> {
    // todo: this was an arbitrary set up, remove
    // let mut boid: Vec<Boid> = Vec::new();
    // boid.push(
    //     Boid::new(
    //         run_options.window.win_left as f32 + 1.,
    //         0.,
    //         vec2(run_options.max_speed, 0.),
    //         0
    //     )
    // );
    // boid.push(
    //     Boid::new(
    //         run_options.window.win_right as f32 - 1.,
    //         0.,
    //         vec2(-run_options.max_speed, 0.),
    //         1
    //     )
    // );
    // return boid;

    get_boids_n(run_options.init_boids, 0, run_options)
}

fn get_boids_n(n: usize, start_id: usize, run_options: &RunOptions) -> Vec<Boid> {
    (start_id..(start_id + n))
        .map(|id| get_boid(run_options, id))
        .collect()
}

fn get_boid(run_options: &RunOptions, id: usize) -> Boid {
    let mut rng = rand::thread_rng();

    match run_options.initiation_strat {
        InitiationStrategy::CircleCenterOut => todo!(),
        InitiationStrategy::CircleCircumferenceIn => {
            let r = run_options.window.win_h as f32 / 2.0;

            let init_vel: f32 = rng.gen::<f32>() / 3. + 0.5;
            let init_pos: f32 = rng.gen::<f32>() / 3. + 2. / 3.;
            let deg = rng.gen::<f32>() * 360.0;

            let x = r * deg.cos();
            let y = r * deg.sin();

            let mut init_vec = Vec2::new(-x, -y);
            init_vec = init_vec.clamp_length_max(0.5 * init_vel);

            Boid::new(x * init_pos, y * init_pos, init_vec, id)
        }
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
        }
    }
}

// fn get_flock_ids_gmm(entities: &[Boid], run_options: &RunOptions) -> Vec<usize> {
//     let test_data: Array2<f32> = entities
//         .iter()
//         .map(|row| [row.position.x, row.position.y])
//         .collect::<Vec<_>>()
//         .into();

//     let test_dataset: DatasetBase<_, _> = test_data.into();

//     // let gmm = GaussianMixtureModel::params(n_clusters)
//     //         .n_runs(10)
//     //         .tolerance(1e-4)
//     //         .with_rng(rng)
//     //         .fit(&dataset).expect("GMM fitting");
//     todo!()
// }

fn get_flock_ids_optics(entities: &[Boid], run_options: &RunOptions) -> Vec<usize> {
    let test_data: Array2<f32> = entities
        .iter()
        .map(|row| [row.position.x, row.position.y])
        .collect::<Vec<_>>()
        .into();

    // let test_dataset: DatasetBase<_, _> = test_data.into();

    let (clusters, noise) =
        pOptics::new(run_options.sensory_distance, 3, Euclidean::default()).fit(&test_data);

    // for map in clusters.iter() {
    //     println!("-----------------");
    //     println!("label: {:?}, len: {:?}", map.0, map.1.len());
    //     for point in map.1 {
    //         println!("    point: {:?}", point)
    //     }
    // }

    // println!("-----------------");
    // println!("label: {:?}, len: {:?}", -1, noise.len());
    // for noise_point in noise.iter() {
    //     println!("    point: {:?}", noise_point)
    // }

    // let res_map: HashMap<usize, usize> = HashMap::new();

    let mut output = vec![0 as usize; entities.len()];

    for (label, points) in clusters {
        for p in points {
            output[p] = label + 1;
        }
    }

    for noise_point in noise {
        output[noise_point] = 0;
    }

    return output;

    // let res = Optics::params(3)
    //     .tolerance(run_options.sensory_distance)
    //     .transform(test_data.view())
    //     .unwrap();

    // println!();
    // println!("OPTICS Result: ");
    // for sample in res.iter() {
    //     println!("{:?}", sample.index());
    // }
    // println!();
}

fn get_flock_ids(tracker: &dyn Tracker, entities: &[Boid], run_options: &RunOptions) -> Vec<usize> {
    // if run_options.dbscan_flock_clustering_on {
    let test_data: Array2<f32> = entities
        .iter()
        // .sorted_by(|b1, b2| b1.id.cmp(&b2.id))
        .map(|row| [row.position.x, row.position.y])
        .collect::<Vec<_>>()
        .into();
    let test_dataset: DatasetBase<_, _> = test_data.into();
    // let linfa_dist_metric = TorL2Dist::new(run_options.window.clone());

    // let res = Dbscan::params_with(3, linfa_dist_metric, CommonNearestNeighbour::KdTree)
    let res = Dbscan::params_with(3, L2Dist, CommonNearestNeighbour::BallTree)
        .tolerance(run_options.sensory_distance)
        .transform(test_dataset)
        .unwrap();

    let res_labels = res
        .targets
        .iter()
        .map(|clust_res| match clust_res {
            Some(cluster) => *cluster + 1,
            None => 0,
        })
        .collect::<Vec<usize>>();
    // res_labels

    return match run_options.boundary {
        Boundary::Toroidal => {
            let output_labels = join_adjacent_flocks(tracker, entities, res_labels, run_options);
            output_labels
        }
        _ => res_labels,
    };
}

/// Helps to join adjacent flocks in toroidal space, using flock labels from an euclidean space
fn join_adjacent_flocks(
    tracker: &dyn Tracker,
    entities: &[Boid],
    labels: Vec<usize>,
    run_options: &RunOptions,
) -> Vec<usize> {
    type ClusterId = usize;
    type BoidId = usize;
    // prepare result, to be modified
    let mut resul_labels = labels.clone();
    // keeping track of which flocks are to be replaced with which ids
    // key = old flock_id, value = new flock_id
    let mut replaced_map: HashMap<ClusterId, HashSet<ClusterId>> = HashMap::new();
    // keeping track of noise entities, as these could cause a chain of joins
    // entity id followed by a set of all the flocks it intersects
    let mut noise_entities: HashMap<BoidId, (&Boid, HashSet<ClusterId>)> = HashMap::new();
    // we will need this to join all members of a flock if part of them participated in boundary joins
    // maps flock onto usize of entities as they are coming in
    let mut inverse_map: HashMap<ClusterId, Vec<BoidId>> = HashMap::new();
    for (index, label) in labels.iter().enumerate() {
        inverse_map
            .entry(*label)
            .or_insert_with(|| Vec::with_capacity(1))
            .push(index)
    }

    // we are essentially trying to create a unidirected graph of all flocks that should be merged,
    // then traverse it and make those merges
    // get boundary entities

    let boids_cluster_lookup = entities
        .iter()
        .zip(resul_labels.iter())
        .map(|tuple| (tuple.0.id, tuple))
        .collect::<HashMap<BoidId, (&Boid, &ClusterId)>>();

    let boundary_boids = boids_cluster_lookup
        .iter()
        .filter(|(_, (b, _))| {
            run_options.window.win_right as f32 - b.position.x.abs() < run_options.sensory_distance
                || run_options.window.win_top as f32 - b.position.y.abs()
                    < run_options.sensory_distance
        })
        .map(|tuple| (*tuple.0, *tuple.1))
        .collect::<HashMap<BoidId, (&Boid, &ClusterId)>>();


    'outer: for (_, (b1, c1)) in boundary_boids.iter() {
        'inner: for b2 in tracker.get_neighbours(b1, run_options) {
            if run_options.clicked_boid_id == b1.id || run_options.clicked_boid_id == b2.id {
                println!("jackpot")
            }
            let c2 = if boundary_boids.contains_key(&b2.id) {
                boundary_boids.get(&b2.id).unwrap().1
            } else {
                continue 'inner;
            };
            // check whether we are have identitcal entities or clusters (could be any flock, even 0 and 0 = noise)
            if b1.id == b2.id || **c1 == *c2 {
                continue 'inner;
            }
            // check whether they are within sensory distance
            if distance_dyn_boid(b1, b2, run_options) > run_options.max_sensory_distance {
                continue 'inner;
            }

            // always use higher cluster_id as replacer, otherwise preserve fetch order,
            // later merge starting with lowest old_cluster_id to ensure a chain merge, e.g.
            // we have cluster_ids: 1, 2, 3 and 4; 1 -> 2 -> 3 -> 4, now we have all as memberers of cluster 4
            let c1_is_replacer = *c1 > c2;
            let (replacer, replaced) = if c1_is_replacer { (*c1, c2) } else { (c2, *c1) };
            // check whether it has already been used
            if replaced_map.contains_key(replaced) {
                let ref mut key = replaced_map.get_mut(&replaced).unwrap();

                if key.contains(replacer) {
                    // already logged, continue loop
                    continue 'inner;
                    //if !c1_is_replacer { continue 'outer} else { continue 'inner }
                } else {
                    // we found a new boundary relationship
                    // replacer = A, replaced = B and key = C
                    // we have B <- C and B <- A, thus B <- A|B <- B|A
                    key.insert(*replacer);
                }
            }
            // if one is noise, add to tracked nosie entities and deal with them later
            // we actually fake them as if they were part of a cluster already as it makes
            // updating everything possible in one sweep
            // todo: duplicated branches, we can reduce this
            if *replacer == 0 {
                noise_entities
                    .entry(if c1_is_replacer { b1.id } else { b2.id })
                    .and_modify(|(_, clusters)| {
                        clusters.insert(*replaced);
                    })
                    .or_insert_with(|| {
                        (
                            if c1_is_replacer { b1 } else { b2 },
                            HashSet::from([*replaced]),
                        )
                    });
                continue 'inner;
                //if c1_is_replacer { continue 'outer} else { continue 'inner }
            } else if *replaced == 0 {
                noise_entities
                    .entry(if !c1_is_replacer { b1.id } else { b2.id })
                    .and_modify(|(_, clusters)| {
                        clusters.insert(*replacer);
                    })
                    .or_insert_with(|| {
                        (
                            if !c1_is_replacer { b1 } else { b2 },
                            HashSet::from([*replacer]),
                        )
                    });
                continue 'inner;
                //if !c1_is_replacer { continue 'outer} else { continue 'inner }
            }

            replaced_map
                .entry(*replaced)
                .or_insert(HashSet::new())
                .insert(*replacer);
        }
    }

    // Now we need to deal with noise entities, the tricky part is multiple noise entities
    // can create a chain that connects up to a multiple of existing clusters. That chain is limited
    // by min_group clustering parameter.
    // In a special case where there is noise, that does not connect with any cluster, but
    // it does connect with noise that neighbours clusters, we have to again keep track of that.
    // One thing to keep in mind is that here we are just dealing with noise and entities
    // that exist on the boundary (|sensory_distance| away from border), everything else should be picked
    // up by the root clustering algorithm this funciton complements. An especially simplifying assumtion
    // is that we are using |min_group| = 3, i.e. when there is a cluster entity and a neighbouring noise
    // the only other noise can be at maximum 1 distance away from that neighbouring noise.
    // F1E1 | <-> NN1 <-> NN2, here NN1 is still a core point, NN2 being edge
    // F1E1 | <-> NN! <-> NN2 <-> NN3 is not possible because under min_group = 3 it would become
    // F1E1 | <-> F2E1 <-> F2E2 <-> F2E3, which this function crunches all into the same cluster F2
    let mut assigned_noise_entities: HashMap<BoidId, ClusterId> = HashMap::new();
    for (nn1_id, (nn1_entity, nn1_clusters)) in noise_entities.iter() {
        // asign the NN1 to a cluster
        let chosen_cluster = nn1_clusters.iter().max().unwrap();
        assigned_noise_entities.insert(*nn1_id, *chosen_cluster);
        // find all NN2 candidates and do like wise
        tracker
            .get_neighbours(nn1_entity, run_options)
            .into_iter()
            .filter(|b| *boids_cluster_lookup[&b.id].1 == 0)
            .for_each(|b| {
                assigned_noise_entities.insert(b.id, *chosen_cluster);
            });
    }

    // now that the noise points have been assignmned, we need to make sure that they are included in merging
    // three options:
    // 1) noise point(s) are added onto existing clusters
    // 2) noise point(s) are added onto existing clusters and cause two clusters to merge where
    //      2a) at least one has participated in current replacement
    //      2b) none have participated in current replacement
    // we can check whether one or the other are participating by symmetry if A and B are merging,
    // then there must be that either B <- A or A <- B (the later in our case as higher values get preference)
    // if none are merging, we'll create a new merging record lower <- higher
    // we can check which clusters are merging by checking which NN1 have more than one neighbour cluster
    for (entity_id, chosen_cluster_id) in assigned_noise_entities.iter() {
        // get the original NN1
        let noise_entry = noise_entities.get(&entity_id);
        // nn1_clusters
        //     .iter()
        //     .filter(|c| **c != *chosen_cluster_id && !replaced_map.contains_key(c))
        //     .collect::<Vec<&usize>>();

        // we do not even have to check the other way around as we know that chosen cluster must nescessarily
        // be greater (in ID number) than the alternative, hence it was chosen earlier
        // thus we either ensure the merging relationship is there or we create it

        // if it has been included in the initial scan for noise entities, it will produce Some
        // otherwise it is a border point and we just assign it to the correct cluster later on
        if let Some((_, nn1_clusters)) = noise_entry {
            // it has been included as part of the initial scan for noise entities
            for nn1_clust in nn1_clusters.iter().filter(|c| **c != *chosen_cluster_id) {
                if replaced_map.contains_key(nn1_clust) {
                    replaced_map
                        .entry(*nn1_clust)
                        .and_modify(|e| {
                            e.insert(*chosen_cluster_id);
                        })
                        .or_insert(HashSet::from([*chosen_cluster_id]));
                }
            }
        }
        //  else {
        //     // it has not been included as part of the initial scan, it is a border point
        // }
    }

    // finally, regarding noise, we get to spook their cluster membership so that they are mapped
    // correctly by what follows afterwards
    // let maybe_noise_entities = inverse_map.get(&0);
    // if let Some(noise_entities) = maybe_noise_entities{
    //     for noise_entity in noise_entities {
    //         if !assigned_noise_entities.contains_key(noise_entity) { continue; }
    //         let assigned_cluster = assigned_noise_entities.get(noise_entity).unwrap();
    //         inverse_map.entry(*assigned_cluster).and_modify(|e| e.push(*noise_entity));
    //     }
    // }

    for (entity_id, assigned_cluster) in assigned_noise_entities.iter() {
        // WARNING: to future self: we won't bother with removing them from the inverse_map[0], this might shoot me in the foot
        // inverse_map.entry(0).and_modify(|e| { e.remove(*entity_id); });
        inverse_map
            .entry(*assigned_cluster)
            .and_modify(|e| e.push(*entity_id));

        // label the noise entity in the result, this might get shifted later if it's cluster gets merged onto
        resul_labels[*entity_id] = *assigned_cluster;
    }

    // using the replace_map we traverse the unidirected graph (in case there is a cycle we simply break it)
    // this gives us all the available clusters to merge, we simply pick the largest label to be the final and voala
    // dbg!(&replaced_map);
    while replaced_map.len() > 0 {
        let mut to_be_replaced: Vec<usize> = Vec::new();
        let mut current_replacer: usize;

        fn traverse_replaced(
            replacer: &usize,
            replaced_map: &mut HashMap<usize, HashSet<usize>>,
        ) -> HashSet<usize> {
            let values = replaced_map.remove(replacer);
            let mut ret = match values {
                Some(replacers) => replacers
                    .iter()
                    .fold(HashSet::<usize>::new(), |mut acc, x| {
                        if !acc.contains(x) {
                            acc.extend(traverse_replaced(x, replaced_map).iter());
                        }
                        acc
                    }),
                None => HashSet::new(),
            };

            ret.insert(*replacer);
            ret
        }

        for (replaced, replacers) in replaced_map.clone().iter().sorted_by_key(|x| x.0) {
            if replaced_map.contains_key(replaced) {
                replaced_map.remove(replaced);
            } else {
                continue;
            }
            to_be_replaced.push(*replaced);

            for replacer in replacers {
                let mut wtf = traverse_replaced(replacer, &mut replaced_map)
                    .into_iter()
                    .collect::<Vec<usize>>();
                to_be_replaced.append(&mut wtf);
            }

            current_replacer = match to_be_replaced.iter().max() {
                Some(max_replacer) => *max_replacer,
                None => panic!("Invalid state!"),
            };

            // println!("replaced: {:?}", &to_be_replaced);
            // println!("replacer: {:?}", &current_replacer);

            for replaced in to_be_replaced.iter().filter(|r| **r != current_replacer) {
                let old_flock_members = inverse_map.get(replaced).expect("should always contain values as no single flock should be merged more than once, exception being noise");
                for old_member in old_flock_members {
                    resul_labels[*old_member] = current_replacer;
                }
            }
            to_be_replaced.clear();
        }
    }

    return resul_labels;
}

#[cfg(test)]
mod tests {

    use std::collections::HashSet;
    use std::f32::consts::PI;

    use glam::Vec2;
    use rand::{Rng, SeedableRng};
    use rand_xoshiro::{self, Xoshiro128StarStar};

    use crate::flock::spathash_tracker::{SpatHash1D, SpatialHashingTableSettings};
    // use crate::flock::{BoidTracker, SpatialHashingTableSettings};
    use crate::math_helpers::distance_dyn;
    use crate::{
        boid::Boid,
        // flock::SpatHash1D,
        options::{self, RunOptions, WindowSize},
    };

    use super::naive_tracker::BoidTracker;
    use super::Tracker;

    #[test]
    fn should_produce_1d_index_positive_point() {
        let p_x = 300.;
        let p_y = 300.;

        let min_x = -400.;
        let max_x = 400.;
        let min_y = -400.;
        let max_y = 400.;
        let cs_x = 200.;
        let cs_y = 200.;
        let table_height = 4.;

        let index = SpatHash1D::get_table_index(
            p_x,
            p_y,
            min_x,
            max_x,
            min_y,
            max_y,
            cs_x,
            cs_y,
            table_height,
        );
        assert_eq!(index, 15);
    }

    #[test]
    fn should_produce_1d_index_positive_extreme() {
        let p_x = 400.;
        let p_y = 400.;

        let min_x = -400.;
        let max_x = 400.;
        let min_y = -400.;
        let max_y = 400.;
        let cs_x = 200.;
        let cs_y = 200.;
        let table_height = 4.;

        let index = SpatHash1D::get_table_index(
            p_x,
            p_y,
            min_x,
            max_x - 1.,
            min_y,
            max_y - 1.,
            cs_x,
            cs_y,
            table_height,
        );
        assert_eq!(index, 15);
    }

    #[test]
    fn should_produce_1d_index_negative() {
        let p_x = -400.;
        let p_y = -400.;

        let min_x = -400.;
        let max_x = 400.;
        let min_y = -400.;
        let max_y = 400.;
        let cs_x = 200.;
        let cs_y = 200.;
        let table_height = 4.;

        let index = SpatHash1D::get_table_index(
            p_x,
            p_y,
            min_x,
            max_x - 1.,
            min_y,
            max_y - 1.,
            cs_x,
            cs_y,
            table_height,
        );
        assert_eq!(index, 0);
    }

    #[test]
    fn test_get_neighbour_single() {
        let clicked = Boid::new(-111.31945, 31.2412548, Vec2::new(0., 0.), 73);
        let neighbour = Boid::new(-80.3626938, 104.281799, Vec2::new(0., 0.), 67);

        let boids = vec![clicked.to_owned(), neighbour.to_owned()];

        let mut run_options = RunOptions::default();

        // the only set at the thime of occurence
        run_options.sensory_distance = 100.;
        run_options.update_sensory_distances();
        run_options.allignment_treshold_coefficient = 1.14999998;
        run_options.cohesion_treshold_coefficient = 0.949999988;
        run_options.separation_treshold_coefficient = 0.349999994;
        run_options.window = options::get_window_size(1400, 900);

        // let flock = Flock::new(&run_options);

        let mut neighbours: Vec<&Boid> = Vec::with_capacity(2);
        let mut tracker = SpatHash1D::new(&boids, &run_options);

        tracker.update_table(&run_options);
        SpatHash1D::get_neighbours(&tracker, &clicked, 17, &run_options, &mut neighbours);

        assert_eq!(neighbours.len(), 1);
    }

    #[test]
    fn should_get_neighbour_single() {
        let clicked = Boid::new(-111.31945, 31.2412548, Vec2::new(0., 0.), 73);
        let neighbour = Boid::new(-80.3626938, 104.281799, Vec2::new(0., 0.), 67);

        let boids = vec![clicked.to_owned(), neighbour.to_owned()];

        let mut run_options = RunOptions::default();

        // the only set at the time of occurence
        run_options.sensory_distance = 100.;
        run_options.update_sensory_distances();
        run_options.allignment_treshold_coefficient = 1.14999998;
        run_options.cohesion_treshold_coefficient = 0.949999988;
        run_options.separation_treshold_coefficient = 0.349999994;
        run_options.window = options::get_window_size(1400, 900);

        // let flock = Flock::new(&run_options);

        let mut neighbours: Vec<&Boid> = Vec::with_capacity(2);
        let mut tracker = SpatHash1D::new(&boids, &run_options);

        tracker.update_table(&run_options);
        SpatHash1D::get_neighbours(&tracker, &clicked, 17, &run_options, &mut neighbours);

        assert_eq!(neighbours.len(), 1);
    }

    /// Get a random point within a distance of a boid
    fn get_boid_points(
        centroid: &Boid,
        max_distance: f32,
        win: &WindowSize,
        n: usize,
        seed: u64,
    ) -> Vec<Vec2> {
        let mut rng = Xoshiro128StarStar::seed_from_u64(seed);

        let mut res: Vec<Vec2> = Vec::with_capacity(n);
        while res.len() < n {
            // https://programming.guide/random-point-within-circle.html

            // random point on a circle with r = max_distance
            let a = rng.gen::<f32>() * 2. * PI;
            let r = max_distance * rng.gen::<f32>().sqrt();

            let x = r * a.cos() + centroid.position.x;
            let y = r * a.sin() + centroid.position.y;

            // if x and y are within boundary, add them to the result
            if !(x < win.win_left as f32
                || x > win.win_right as f32
                || y < win.win_bottom as f32
                || y > win.win_top as f32)
            {
                res.push(Vec2::new(x, y));
            }
        }

        res
    }

    fn get_neighbourhood_setup(
        run_options: &RunOptions,
        sensory_distance: f32,
        no_neighbours: usize,
    ) -> (Vec<Boid>, Vec<Vec<Boid>>, Vec<Boid>) {
        let left = run_options.window.win_left as f32 + sensory_distance - 1.;
        let right = run_options.window.win_right as f32 - sensory_distance + 1.;
        let bottom = run_options.window.win_bottom as f32 + sensory_distance - 1.;
        let top = run_options.window.win_top as f32 - sensory_distance + 1.;

        // set up the boids
        // create a few clicked boids to get the neighbours of
        let clicked = vec![
            // center
            Boid::new(0., 0., Default::default(), 0),
            // left bottom edge
            Boid::new(left, bottom, Default::default(), 1),
            //right top edge
            Boid::new(right, top, Default::default(), 2),
            //left middle
            Boid::new(left, 0., Default::default(), 3),
            // bottom middle
            Boid::new(0., bottom, Default::default(), 4),
        ];

        // for each clicked boid, create a few randomly distributed flockmates around, within sensory distance
        let mut id: usize = 4;
        let mut clicked_neighbours: Vec<Vec<Boid>> = vec![Default::default(); clicked.len()];

        clicked.iter().enumerate().for_each(|(index, b)| {
            clicked_neighbours[index] = get_boid_points(
                b,
                run_options.max_sensory_distance * 0.9,
                &run_options.window,
                no_neighbours,
                11551155,
            )
            .iter()
            .map(|p| {
                id += 1;
                Boid::new(p.x, p.y, Default::default(), id)
            })
            .collect();
        });

        // flatten the neighbours, mix in the clicked boids and create the tracker instance

        let boids: Vec<Boid> = clicked_neighbours
            .to_owned()
            .into_iter()
            .flatten()
            .chain(clicked.to_owned().into_iter())
            .collect();

        (clicked, clicked_neighbours, boids)
    }

    /// integrates all of the flocking mechanics to retrieve neighbours and tests correct neighbours are retrieved
    /// systems integrated - spatial subdivison, determining cell index of a given entity, checking it's distance to a centroid
    // #[ignore]
    #[test]
    fn should_get_neighbours_sweep_spathhash() {
        let mut run_options = RunOptions::default();

        // set up the boid space
        let sensory_distance = 100.;
        let no_neighbours: usize = 51;

        // 5 from set up, no_neighbours for each set up boid
        run_options.init_boids = 5 + 5 * no_neighbours;
        run_options.sensory_distance = sensory_distance;
        run_options.allignment_treshold_coefficient = 1.14999998;
        run_options.cohesion_treshold_coefficient = 0.949999988;
        run_options.separation_treshold_coefficient = 0.349999994;
        run_options.update_sensory_distances();
        run_options.window = options::get_window_size(1400, 900);
        run_options.neighbours_cosidered = no_neighbours;

        let (clicked, clicked_neighbours, boids) =
            get_neighbourhood_setup(&run_options, sensory_distance, no_neighbours);

        let mut tracker = SpatHash1D::new(&boids, &run_options);
        tracker.update_table(&run_options);

        fn get_index(
            entity: &Boid,
            run_options: &RunOptions,
            settings: &SpatialHashingTableSettings,
        ) -> usize {
            SpatHash1D::get_table_index(
                entity.position.x,
                entity.position.y,
                run_options.window.win_left as f32,
                (run_options.window.win_right - 1) as f32,
                run_options.window.win_bottom as f32,
                (run_options.window.win_top - 1) as f32,
                settings.x_cell_res,
                settings.y_cell_res,
                settings.x_cell_count as f32,
            )
        }

        // check all flockmates have been fetched
        for i in 0..clicked.len() {
            let clicked_boid = &clicked[i];

            let clicked_boid_index = get_index(clicked_boid, &run_options, &tracker.settings);

            let mut neighbours: Vec<&Boid> = Vec::with_capacity(no_neighbours);

            // neighbours.clear();
            tracker.get_neighbours(
                clicked_boid,
                clicked_boid_index,
                &run_options,
                &mut neighbours,
            );

            // let retrieved_neighbour_id_set: HashSet<usize> =
            //     neighbours.iter().map(|nb| nb.id).collect();

            // let neighbour_difference: Vec<&Boid> = clicked_neighbours[i]
            //     // .to_owned()
            //     .iter()
            //     .filter(|nb| !retrieved_neighbour_id_set.contains(&nb.id))
            //     .collect();

            // let neighbour_difference_indexes: Vec<usize> = neighbour_difference
            //     .iter()
            //     .map(|b| get_index(*b, &run_options, &tracker.settings))
            //     .collect();

            // let neighbour_difference_dist: Vec<f32> = neighbour_difference
            //     .iter()
            //     .map(|b| distance_dyn(b.position.x, clicked_boid.position.x, b.position.y,clicked_boid.position.y, &run_options))
            //     .collect();

            assert_eq!(clicked_neighbours[i].len(), neighbours.len()) //, "centroid id: {}, ND: {:#?}, NDI: {:#?}, NDD: {:#?}", clicked_boid.id, neighbour_difference, neighbour_difference_indexes, neighbour_difference_dist)
        }
    }

    #[test]
    fn should_get_neighbours_sweep_boidtracker() {
        let mut run_options = RunOptions::default();

        // set up the boid space
        let sensory_distance = 100.;
        let no_neighbours: usize = 51;

        run_options.sensory_distance = sensory_distance;
        run_options.update_sensory_distances();
        run_options.allignment_treshold_coefficient = 1.14999998;
        run_options.cohesion_treshold_coefficient = 0.949999988;
        run_options.separation_treshold_coefficient = 0.349999994;
        run_options.window = options::get_window_size(1400, 900);
        run_options.neighbours_cosidered = no_neighbours;

        let (clicked, clicked_neighbours, boids) =
            get_neighbourhood_setup(&run_options, sensory_distance, no_neighbours);

        // let tracker = BoidTracker::new(&boids, &run_options);
        // check all flockmates have been fetched

        let mut neighbours: Vec<&Boid> = Vec::with_capacity(no_neighbours);

        for i in 0..clicked.len() {
            neighbours.clear();

            let clicked_boid = &clicked[i];

            BoidTracker::get_neighbours_naive(clicked_boid, &boids, &run_options, &mut neighbours);

            let retrieved_neighbour_id_set: HashSet<usize> =
                neighbours.iter().map(|nb| nb.id).collect();

            let neighbour_difference: Vec<&Boid> = clicked_neighbours[i]
                // .to_owned()
                .iter()
                .filter(|nb| !retrieved_neighbour_id_set.contains(&nb.id))
                .collect();

            let neighbour_difference_dist: Vec<f32> = neighbour_difference
                .iter()
                .map(|b| {
                    distance_dyn(
                        b.position.x,
                        clicked_boid.position.x,
                        b.position.y,
                        clicked_boid.position.y,
                        &run_options,
                    )
                })
                .collect();

            assert_eq!(
                clicked_neighbours[i].len(),
                neighbours.len(),
                "centroid id: {}, ND: {:#?}, NDD: {:#?}",
                clicked_boid.id,
                neighbour_difference,
                neighbour_difference_dist
            )
        }
    }

    // // there is a problem, which makes the SpatHash1D 'mostly' work, at times it does not fetch all of the neighbours, and is extremely difficult to reproduce
    // // here I set up a full boids environment, where all of the entities will be controlled and moved by SpatHash1D, and the double checked by copying over
    // // all of the boids to the naive BoidsTracker to see which entities and under what conditions get omitted
    // #[test]
    // fn should_retrieve_the_same_entities() {
    //     let mut run_options = RunOptions::default();

    //     // set up the boid space
    //     let sensory_distance = 100.;
    //     let no_neighbours: usize = 51;

    //     // 5 from set up, no_neighbours for each set up boid
    //     run_options.init_boids = 5 + 5 * no_neighbours;
    //     run_options.sensory_distance = sensory_distance;
    //     run_options.update_sensory_distances();
    //     run_options.allignment_treshold_coefficient = 1.14999998;
    //     run_options.cohesion_treshold_coefficient = 0.949999988;
    //     run_options.separation_treshold_coefficient = 0.349999994;
    //     run_options.window = options::get_window_size(1400., 900.);
    //     run_options.neighbours_cosidered = no_neighbours;

    //     let (
    //         clicked,
    //         clicked_neighbours,
    //         boids
    //     ) = get_neighbourhood_setup(&run_options, sensory_distance, no_neighbours);

    //     let mut tracker = SpatHash1D::new(&boids, &run_options);
    //     tracker.update_table(&run_options);

    //     fn get_index(entity: &Boid, run_options: &RunOptions, settings: &SpatialHashingTableSettings) -> usize {
    //          SpatHash1D::get_table_index(
    //             entity.position.x,
    //             entity.position.y,
    //             run_options.window.win_left,
    //             run_options.window.win_right - 1.,
    //             run_options.window.win_bottom,
    //             run_options.window.win_top - 1.,
    //             settings.x_cell_res,
    //             settings.y_cell_res,
    //             settings.x_cell_count,
    //         )
    //     }

    //     // check all flockmates have been fetched
    //     for i in 0..clicked.len() {
    //         let clicked_boid = &clicked[i];

    //         let clicked_boid_index = get_index(clicked_boid, &run_options, &tracker.settings);

    //         let mut clicked_boid_neighbours: Vec<&Boid> = Vec::with_capacity(no_neighbours);

    //         clicked_boid_neighbours.clear();
    //         tracker.get_neighbours(
    //             clicked_boid,
    //             clicked_boid_index,
    //             &run_options,
    //             &mut clicked_boid_neighbours,
    //         );

    //         let retrieved_neighbour_id_set: HashSet<usize> = clicked_boid_neighbours.iter().map(|nb| nb.id).collect();

    //         let neighbour_difference: Vec<&Boid> = clicked_neighbours[i]
    //             // .to_owned()
    //             .iter()
    //             .filter(|nb| !retrieved_neighbour_id_set.contains(&nb.id))
    //             .collect();

    //         let neighbour_difference_indexes: Vec<usize> = neighbour_difference
    //             .iter()
    //             .map(|b| get_index(*b, &run_options, &tracker.settings))
    //             .collect();

    //         let neighbour_difference_dist: Vec<f32> = neighbour_difference
    //             .iter()
    //             .map(|b| distance_dyn(b.position.x, clicked_boid.position.x, b.position.y,clicked_boid.position.y, &run_options))
    //             .collect();

    //         assert_eq!(clicked_neighbours[i].len(), clicked_boid_neighbours.len(), "centroid id: {}, ND: {:#?}, NDI: {:#?}, NDD: {:#?}", clicked_boid.id, neighbour_difference, neighbour_difference_indexes, neighbour_difference_dist)
    //     }
    // }
}
