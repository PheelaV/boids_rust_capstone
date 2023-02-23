use std::collections::HashSet;
use std::f32::consts::PI;
use std::iter;
use std::sync::Mutex;

use glam::Vec2;
use linfa::traits::Transformer;
use linfa::DatasetBase;
use linfa_clustering::Dbscan;
use ndarray::Array2;
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
use crate::options::InitiationStrategy;
use crate::options::RunOptions;

// pub mod tracker;

// pub use tracker::*;

// a tracker takes ownership of a flock of boids
// a tracker will follow an implementation consisting of choosing a datastructure to
// hold spatial information about boids in and have traits for querying 2D, later 3D
// information
pub trait Tracker {
    fn new(entities: &[Boid], run_options: &RunOptions) -> Self;
    fn update(&mut self, run_options: &RunOptions);
    fn insert_single(&mut self, entity: Boid, _: &RunOptions);
    fn restart(&mut self, entities: &[Boid]);
    fn delete_last(&mut self) -> Option<Boid>;
    fn delete_multiple(&mut self, ids_delete: &[usize], run_options: &RunOptions);
    fn insert_multiple(&mut self, entities: &[Boid], _: &RunOptions);
    fn view(&self) -> (&Vec<Boid>, &Vec<BoidMetadata>);
    fn view2<'a>(&'a self) -> Box<dyn Iterator<Item = (&'a Boid, &'a BoidMetadata)> + 'a>;
}

// todo: maybe there is more performance if we stop can discriminate entities that are already in the correct cell
/// Implementation follows [paper](https://www.researchgate.net/publication/277870601_A_Hash_Table_Construction_Algorithm_for_Spatial_Hashing_Based_on_Linear_Memory)
/// It has been modified into a form of unstable insert
pub struct SpatHash1D {
    // for each cell, there is a pivot with ranges into table
    // index of pivots corresponds to cell representation
    pivots: Vec<SpatHashPiv>,
    // 1D representation of cells
    // indexes have no special meaning
    table: Vec<Boid>,
    // everything that a Boid does not nescesarily need for its own functioning
    // indexes are effectively boid ids
    metadata: Vec<BoidMetadata>,
    // for each boid, records which cell it's at
    // indexes are effectively boid ids
    index: Vec<usize>,
    // contains information about the table, used for updating and reconstruction
    settings: SpatialHashingTableSettings,
}

#[derive(Clone)]
struct SpatHashPiv {
    usg: usize,
    init: Option<usize>,
    fin: Option<usize>,
}

impl Default for SpatHashPiv {
    fn default() -> Self {
        SpatHashPiv {
            usg: 0,
            init: None,
            fin: None,
        }
    }
}

impl Tracker for SpatHash1D {
    fn new(entities: &[Boid], run_options: &RunOptions) -> Self {
        let settings = SpatHash1D::get_tracker_settings(run_options);
        let metadata: Vec<BoidMetadata> =
            entities.iter().map(|e| BoidMetadata::new(e)).collect();
        // metadata[0].boid_type = BoidType::Disruptor;

        SpatHash1D {
            // initialize vector with both capacity and values prefilled to simplify code in the update_table
            pivots: vec![Default::default(); settings.cell_count],
            table: entities.to_vec(),
            metadata: metadata,
            // initialize vector with both capacity and values prefilled to simplify code in the update_table
            index: vec![0; run_options.init_boids as usize],
            settings,
        }
    }

    fn update(&mut self, run_options: &RunOptions) {
        let new_settings = SpatHash1D::get_tracker_settings(run_options);
        let cell_difference = self.settings.cell_count as i64 - new_settings.cell_count as i64;
        if cell_difference != 0 {
            self.pivots
                .resize(new_settings.cell_count, Default::default());
            self.settings = new_settings;
            self.update_table(run_options);
        }

        self.update_location(run_options);
        self.update_table(run_options);

        if run_options.dbscan_flock_clustering_on {
            let flock_ids = get_flock_ids(&self.table, run_options);

            self.table
                .iter()
                .enumerate()
                .for_each(|(index, boid)| self.metadata[boid.id].cluster_id = flock_ids[index])
        }
    }

    fn insert_single(&mut self, entity: Boid, run_options: &RunOptions) {
        self.insert(entity, run_options);

        // make sure the table updates its structure
        self.update_table(run_options);
    }

    fn restart(&mut self, entities: &[Boid]) {
        self.table = entities.to_vec();
    }

    fn delete_last(&mut self) -> Option<Boid> {
        let last_id = self.index.len() - 1;
        match self.table.iter().position(|e| e.id == last_id) {
            Some(table_index) => Some(self.table.swap_remove(table_index)),
            None => None,
        }
    }

    /// deletes a certain boid from the table and its auxilary arrays,
    fn delete_multiple(&mut self, ids_delete: &[usize], run_options: &RunOptions) {
        let ids_delete_set: HashSet<usize> = HashSet::from_iter(ids_delete.iter().cloned());

        // removing values from the 1D table
        self.table.retain(|boid| !ids_delete_set.contains(&boid.id));

        // resize auxilary arrays
        self.metadata.resize(self.table.len(), Default::default());
        self.index.resize(self.table.len(), Default::default());

        // update data structure (index, pivots)
        self.update_table(run_options);
    }

    fn insert_multiple(&mut self, entities: &[Boid], run_options: &RunOptions) {
        entities
            .into_iter()
            .for_each(|e| self.insert(*e, run_options));

        // make sure the table updates its structure
        self.update_table(run_options);
    }

    fn view<'a>(&self) -> (&Vec<Boid>, &Vec<BoidMetadata>) {
        (&self.table, &self.metadata)
    }

    // https://stackoverflow.com/questions/31904842/return-a-map-iterator-which-is-using-a-closure-in-rust
    fn view2<'a>(&'a self) -> Box<dyn Iterator<Item = (&'a Boid, &'a BoidMetadata)> + 'a> {
        // self.table.iter().map
        Box::new(self.table.iter().map(|e| (e, &self.metadata[e.id])))
    }
}

impl SpatHash1D {
    // get the cell vectors pointing towards the neighbourhood cell
    // in eucledian space
    const HOME: [i32; 2] = [0, 0];
    const S: [i32; 2] = [0, -1];
    const SW: [i32; 2] = [-1, -1];
    const SE: [i32; 2] = [1, -1];
    const W: [i32; 2] = [-1, 0];
    const E: [i32; 2] = [1, 0];
    const N: [i32; 2] = [0, 1];
    const NE: [i32; 2] = [1, 1];
    const NW: [i32; 2] = [-1, 1];
    
    #[rustfmt::skip]
    const LOOKUP: [[[i32; 2]; 9]; 16] = [
        // 0 = middle
        [Self::NW, Self::N, Self::NE, Self::W, Self::HOME, Self::E, Self::SW, Self::S, Self::SE,],
        // 1 = top
        [Self::HOME, Self::HOME, Self::HOME, Self::E, Self::HOME, Self::W, Self::SW, Self::S, Self::SE],
        // 2 = bottom
        [Self::NW, Self::N, Self::NE, Self::W, Self::HOME, Self::E, Self::HOME, Self::HOME, Self::HOME],
        // placeholder
        [[0; 2]; 9],
        // 4 = right
        [Self::NW, Self::N, Self::HOME, Self::W, Self::HOME, Self::HOME, Self::SW, Self::S, Self::HOME],
        // 5 = right | top
        [Self::HOME, Self::HOME, Self::HOME, Self::W, Self::HOME, Self::HOME, Self::SW, Self::S, Self::HOME],
        // 6 = right | bottom
        [Self::NW, Self::N, Self::HOME, Self::W, Self::HOME, Self::HOME, Self::HOME, Self::HOME, Self::HOME],
        // placeholder
        [[0; 2]; 9],
        // 8 = left
        [Self::HOME, Self::N, Self::NE, Self::HOME, Self::HOME, Self::E, Self::HOME, Self::S, Self::SE],
        // 9 = left | top
        [Self::HOME, Self::HOME, Self::HOME, Self::HOME, Self::HOME, Self::E, Self::HOME, Self::S, Self::SE],
        // placeholder
        [[0; 2]; 9],
        // 10 = left | botom
        [Self::HOME, Self::N, Self::NE, Self::HOME, Self::HOME, Self::E, Self::HOME, Self::HOME, Self::HOME],
        // placeholder
        [[0; 2]; 9],
        // placeholder
        [[0; 2]; 9],
        // placeholder
        [[0; 2]; 9],
        // placeholder
        [[0; 2]; 9],
    ];

    // https://stackoverflow.com/questions/31904842/return-a-map-iterator-which-is-using-a-closure-in-rust
    // fn view2<'a>(&'a self) -> Box<dyn Iterator<Item= (&'a Boid, &'a BoidMetadata)> + 'a> {
    //     // self.table.iter().map
    //     Box::new(self.table.iter().map(|e| (e, &self.metadata[e.id])))
    // }

    fn insert(&mut self, entity: Boid, run_options: &RunOptions) {
        let index = SpatHash1D::get_table_index(
            entity.position.x,
            entity.position.y,
            run_options.window.win_left,
            run_options.window.win_right - 1.,
            run_options.window.win_bottom,
            run_options.window.win_top - 1.,
            self.settings.x_cell_res,
            self.settings.y_cell_res,
            self.settings.x_cell_count,
        );

        // expand table and index vectors
        self.index.push(index);
        self.metadata.push(BoidMetadata::new(&entity));
        self.table.push(entity);
    }

    pub(self) fn update_location(&mut self, run_options: &RunOptions) {
        // reset metadata
        let metadata_default: BoidMetadata = Default::default();
        for b in self.table.iter() {
            self.metadata[b.id].id = b.id;
            self.metadata[b.id].clicked_neighbour_id = metadata_default.clicked_neighbour_id;
            self.metadata[b.id].cluster_id = metadata_default.cluster_id;
            self.metadata[b.id].n_neighbours = metadata_default.n_neighbours;
            self.metadata[b.id].accelleration_update = metadata_default.accelleration_update;
        }

        let mut neighbours: Vec<&Boid> = Vec::with_capacity(32);
        let mut accellerations: Vec<Vec2> = vec![Default::default(); run_options.init_boids];
        let mut metadata: Vec<BoidMetadata> = vec![Default::default(); run_options.init_boids];

        let mut clicked_neighbours: Vec<usize> =
            Vec::with_capacity(run_options.neighbours_cosidered);
        // todo: look at ways to reduce the vector allocations above

        for e in 0..self.table.len() {
            neighbours.clear();

            let b: &Boid = &self.table[e];
            self.get_neighbours(b, self.index[b.id], run_options, &mut neighbours);

            // this should really be in one of the boid rules,
            // but then there goes borrow checker and boids having
            // to mutate their own state for which the current pipeline
            // simply is not buit for, so here is a workaround
            if run_options.wander_on {
                let wander_next = MY_RNG
                    .lock()
                    .unwrap()
                    .gen_range(-run_options.wander_rate..run_options.wander_rate);
                metadata[self.table[e].id].wander_direction =
                    (self.metadata[self.table[e].id].wander_direction + wander_next) % (2. * PI);
            }

            metadata[self.table[e].id].n_neighbours = neighbours.len();

            // clicked neighbour
            if self.table[e].id == run_options.clicked_boid_id {
                for n in &neighbours {
                    clicked_neighbours.push(n.id)
                }
            }

            // this blows up, because the vector is allocated with appropriate capacity, but is empty
            // accellerations[e] = self.table[e].run_rules3(&neighbours, run_options);
            // I will attempt to patch it with push to get rid of the overhead in pre-creating the Vec2s
            // and rely on implicit order, which might cause issues, but should not as the for loop goes
            //sequentually through the entities
            let accelleration = self.table[e].run_rules(&neighbours, &self.metadata, run_options);
            metadata[self.table[e].id].accelleration_update = accelleration;
            accellerations[e] = accelleration;
        }

        clicked_neighbours
            .iter()
            .for_each(|cn_id| metadata[*cn_id].clicked_neighbour_id = run_options.clicked_boid_id);

        // propagate metadata update
        for id in 0..self.metadata.len() {
            self.metadata[id].clicked_neighbour_id = metadata[id].clicked_neighbour_id;
            self.metadata[id].n_neighbours = metadata[id].n_neighbours;
            self.metadata[id].accelleration_update = metadata[id].accelleration_update;
            self.metadata[id].wander_direction = metadata[id].wander_direction
        }

        // apply the forces
        for e in 0..self.table.len() {
            self.table[e].apply_force(accellerations[e]);
            self.table[e].update_location(&run_options)
        }
        accellerations.clear();
    }

    /// assumes all the auxilary arrays are of correct size
    pub(self) fn update_table(&mut self, run_options: &RunOptions) {
        // reset pivot metadata as we are in the next iterration, e.g., we do
        // not want pivot.usg to linger along from previous iterration
        self.pivots.iter_mut().for_each(|pivot| {
            pivot.usg = Default::default();
            pivot.init = Default::default();
            pivot.fin = Default::default();
        });

        for e in 0..self.table.len() {
            let index = SpatHash1D::get_table_index(
                self.table[e].position.x,
                self.table[e].position.y,
                run_options.window.win_left,
                run_options.window.win_right - 1.,
                run_options.window.win_bottom,
                run_options.window.win_top - 1.,
                self.settings.x_cell_res,
                self.settings.y_cell_res,
                self.settings.x_cell_count,
            );

            // if we are dealing with the same entity, do not tick usage
            self.index[self.table[e].id] = index;
            self.pivots[index].usg += 1;
        }

        let mut accum: usize = 0;
        for p in 0..self.pivots.len() {
            if self.pivots[p].usg > 0 {
                self.pivots[p].init = Some(accum);
                accum += self.pivots[p].usg;
                self.pivots[p].fin = Some(accum);
            } else {
                self.pivots[p].init = None;
                self.pivots[p].fin = None;
            }
        }

        let mut pivots = self.pivots.to_owned();

        let mut sorted = vec![false; self.table.len()];
        let mut destination: usize;

        let mut e = 0;
        while e < self.table.len() {
            // here is where the algorithm differs slightly from the original
            // as we are reusing the same table as the set of entities, we
            // have to make a swap this collapses into a form of unstable
            // insertion sort...

            let id = self.table[e].id;

            destination = usize::wrapping_sub(pivots[self.index[id]].fin.unwrap(), 1);

            // skip if the destination has alreade been sorted, current table[e] is in right place
            // not checking for 0 but for MAX (18446744073709551615) as we take advantage of wraparound
            if destination == std::usize::MAX || sorted[destination] || e == destination {
                e += 1;
                continue;
            }
            self.table.swap(e, destination);

            // update the pivot now that one entity of the cell has been mapped to the correct place in the 1D table
            // each unique value of fin during the iterration effectively represents a unique location for the next
            // entity in that cell
            // we consume it during swap and then decrease to point at the location of cell's next entity
            // until it points at the first(swapped, last in terms of the 1D table) entity of the previous cell
            pivots[self.index[id]].fin = Some(pivots[self.index[id]].fin.unwrap() - 1);

            // set destination as sorted
            sorted[destination] = true;
        }
    }

    /// Returns "hashed" value representing an index for spatial subdivision, handles a zero centered coordinate system
    ///
    /// ## Arguments
    ///
    /// - `p_x` - x coordinate of a point of interest
    /// - `p_y` - y coordinate of a point of interest
    /// - `min_x` - minimum x axis value
    /// - `min_y` - minimum y axis value
    /// - `cs_x` - cell size along the x axis
    /// - `cs_y` - cell size along the y axis
    /// - `table_width` - the table width, how many cells make up a row
    ///
    /// ## Returns
    /// an index as _usize_, that takes values from 0..(x_cells * y_cells)
    ///
    /// **watch out**, in case of a symmetric/centered coordinate system (like nannou)
    /// the recommendation is to subtract 1 from the x/y max values to reduce the
    /// range of returned indexes by one as the centering will cause the dimensions
    /// to have uneven counts, returning index that would be out of bounds for the
    /// grid
    pub fn get_table_index(
        p_x: f32,
        p_y: f32,
        min_x: f32,
        max_x: f32,
        min_y: f32,
        max_y: f32,
        cs_x: f32,
        cs_y: f32,
        table_width: f32,
    ) -> usize {
        (((p_x.clamp(min_x, max_x) - min_x) / cs_x).floor()
            + (((p_y.clamp(min_y, max_y) - min_y) / cs_y).floor() * table_width)) as usize
    }
    // pub fn get_table_index2(
    //     p_x: f32,
    //     p_y: f32,
    //     min_x: f32,
    //     max_x: f32,
    //     min_y: f32,
    //     max_y: f32,
    //     cs_x: f32,
    //     cs_y: f32,
    //     table_width: f32,
    // ) -> usize {
    //     (((p_x.clamp(min_x, max_x) - min_x) / cs_x).ceil()
    //         + (((p_y.clamp(min_y, max_y) - min_y) / cs_y).ceil() * table_width)) as usize
    // }

    // pub fn lattice_query_to_grid_cell(
    //     lattice_cell: i32,
    //     min_x: f32,
    //     max_x: f32,
    //     min_y: f32,
    //     max_y: f32
    // ) -> usize {
    //     let lc = lattice_cell as f32;

    //     todo!();
    // }

    pub fn get_tracker_settings(run_options: &RunOptions) -> SpatialHashingTableSettings {
        let x_range = (run_options.window.win_right - run_options.window.win_left).ceil();
        let y_range = (run_options.window.win_top - run_options.window.win_bottom).ceil();

        // let x_cell_count = (x_range as f32 / (100.)).ceil();

        // here we multiply the sensory distance by two, because it is a radius and the method
        // for fetching neighbours is implemented with the assumption that the sensory radius
        // fits into one cell of the grid
        // value any higher than that will an unnescessary ammount of boids
        // any lower, and it won't fetch all boids that should have been considered

        let mut x_cell_count: f32;
        let mut y_cell_count: f32;

        // the problems I am hitting here are related to the gauss circle problem
        // which consists of counting the number of lattice points N(r) inside the boundary of
        // a circle of radius r with center at the origin
        // there is an exact solution: https://mathworld.wolfram.com/GausssCircleProblem.html
        // - considers the floored radius
        // tests all discrete points within positive quartal
        if x_range.max(y_range) / 2. < run_options.max_sensory_distance {
            x_cell_count = 1.;
            y_cell_count = 1.;
        } else {
            x_cell_count = (x_range as f32 / (run_options.max_sensory_distance * 2.)).ceil();
            y_cell_count = (y_range as f32 / (run_options.max_sensory_distance * 2.)).ceil();
        }

        if y_cell_count < 2. || x_cell_count < 2. {
            y_cell_count = 1.;
            x_cell_count = 1.;
        }

        let cell_count = (x_cell_count * y_cell_count) as usize;

        let x_cell_res = x_range / x_cell_count;
        let y_cell_res = y_range / y_cell_count;

        SpatialHashingTableSettings {
            x_range,
            y_range,
            x_cell_count,
            y_cell_count,
            cell_count,
            x_cell_res,
            y_cell_res,
        }
    }

    pub(self) fn get_neighbours<'a>(
        &'a self,
        boid: &Boid,
        cell_index: usize,
        run_options: &RunOptions,
        neighbours: &mut Vec<&'a Boid>,
    ) {
        // if the table grid is less than 4x4
        if self.settings.x_cell_count <= 3. && self.settings.y_cell_count <= 3. {
            // fall back to naive method
            return BoidTracker::get_neighbours_naive(boid, &self.table, run_options, neighbours);
        }

        // only add it once
        let is_left = cell_index % self.settings.x_cell_count as usize == 0;
        let is_right = (cell_index % self.settings.x_cell_count as usize) + 1
            == self.settings.x_cell_count as usize;
        let is_bottom = cell_index / self.settings.x_cell_count as usize == 0;
        let is_top =
            cell_index as f32 >= (self.settings.y_cell_count - 1f32) * self.settings.x_cell_count;

        // self.settings
        // if there is any piece of code in this code base that deserves an internet trophy it is this
        // and LOOKUP - unreadable, unreasonable, non-idiomatic, hella fun and blazingly fast?
        let lookup_index = {
            (is_left as u8) << 3 | (is_right as u8) << 2 | (is_bottom as u8) << 1 | is_top as u8
        } as usize;
        // depending on the match, iterate the vectors, presuming, serves as a
        // lookup for a potential (1..3)*(1..3) area on a n*n grid
        for cell in Self::LOOKUP[lookup_index]
            // depending on the match, iterate the vectors
            .iter()
            // filter out the out-of-bound vectors, the home vector is used as a placeholder
            .filter(|l| **l != Self::HOME)
            // always consider the boid's current 'home' vector
            .chain(iter::once(&Self::HOME))
            // convert vectors to the table's 1D indexes pointing to cells
            .map(|l| (cell_index as i32 + l[0] + l[1] * self.settings.x_cell_count as i32) as usize)
        {
            if self.pivots[cell].usg == 0 {
                continue;
            }
            for index in self.pivots[cell].init.unwrap()..self.pivots[cell].fin.unwrap() {
                if self.table[index].id != boid.id
                    && distance_dyn_boid(boid, &self.table[index], run_options)
                        <= run_options.max_sensory_distance
                {
                    neighbours.push(&self.table[index]);
                    if run_options.neighbours_cosidered != 0
                        && neighbours.len() == run_options.neighbours_cosidered
                    {
                        break;
                    }
                }
            }
        }
    }
}
/// A naive implementation of boids tracking, which uses an O(N^2) algorithm for
/// finding boid's neighbours.
pub struct BoidTracker {
    boids: Vec<Boid>,
    metadata: Vec<BoidMetadata>,
}

impl BoidTracker {
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

impl Tracker for BoidTracker {
    fn new(entities: &[Boid], run_options: &RunOptions) -> Self {
        BoidTracker {
            boids: entities.to_vec(),
            metadata: vec![Default::default(); run_options.init_boids],
        }
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

            BoidTracker::get_neighbours_naive(b_current, &self.boids, run_options, &mut neighbours);
            accelleration.push(self.boids[i_cur].run_rules(&neighbours, &metadata, &run_options));
            // for cn in &clicked_neighbours {
            //     metadata[cn.id].clicked_neighbour_id = run_options.clicked_boid_id;
            // }
            if metadata[i_cur].id == run_options.clicked_boid_id {
                for cn in neighbours.iter() {
                    metadata[cn.id].clicked_neighbour_id = metadata[i_cur].id
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

        if run_options.dbscan_flock_clustering_on {
            let flock_ids = get_flock_ids(&self.boids, run_options);

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

    // https://stackoverflow.com/questions/31904842/return-a-map-iterator-which-is-using-a-closure-in-rust
    fn view2<'a>(&'a self) -> Box<dyn Iterator<Item = (&'a Boid, &'a BoidMetadata)> + 'a> {
        // self.table.iter().map
        Box::new(self.boids.iter().map(|e| (e, &self.metadata[e.id])))
    }
}

/// Uses a spatial hashing space division method, where all cells of the underlying
/// table are stored in a 1D array, with the individual cell's being allocated
/// dynamically, reused or updated as needed.
///
/// Uses an Θ(n) algorithm for finding the boid's neighbours.
/// ref: Optimization_of_large-scale_real-time_simulations
pub struct SpatialHashingTableSettings {
    pub x_range: f32,
    pub y_range: f32,
    pub x_cell_count: f32,
    pub y_cell_count: f32,
    pub cell_count: usize,
    pub x_cell_res: f32,
    pub y_cell_res: f32,
    // run_options: &'a RunOptions
}
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
            let r = run_options.window.win_h / 2.0;

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

fn get_flock_ids(entities: &[Boid], run_options: &RunOptions) -> Vec<usize> {
    // if run_options.dbscan_flock_clustering_on {
    let test_data: Array2<f32> = entities
        .iter()
        .map(|row| [row.position.x, row.position.y])
        .collect::<Vec<_>>()
        .into();
    let test_dataset: DatasetBase<_, _> = test_data.into();
    let res = Dbscan::params(3)
        .tolerance(run_options.sensory_distance)
        .transform(test_dataset)
        .unwrap();

    let t = res
        .targets
        .iter()
        .map(|clust_res| match clust_res {
            Some(cluster) => *cluster + 1,
            None => 0,
        })
        .collect::<Vec<usize>>();
    t
    // let label_count = res.label_count().remove(0);

    // println!();
    // println!("Result: ");
    // for (label, count) in label_count {
    //     match label {
    //         None => println!(" - {} noise points", count),
    //         Some(i) => println!(" - {} points in cluster {}", count, i),
    //     }
    // }
    // };
}

#[cfg(test)]
mod tests {

    use std::collections::HashSet;
    use std::f32::consts::PI;

    use glam::Vec2;
    use rand::{Rng, SeedableRng};
    use rand_xoshiro::{self, Xoshiro128StarStar};

    use crate::flock::{BoidTracker, SpatialHashingTableSettings};
    use crate::math_helpers::distance_dyn;
    use crate::{
        boid::Boid,
        flock::SpatHash1D,
        options::{self, RunOptions, WindowSize},
    };

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
        run_options.window = options::get_window_size(1400., 900.);

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
        run_options.window = options::get_window_size(1400., 900.);

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
            if !(x < win.win_left || x > win.win_right || y < win.win_bottom || y > win.win_top) {
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
        let left = run_options.window.win_left + sensory_distance - 1.;
        let right = run_options.window.win_right - sensory_distance + 1.;
        let bottom = run_options.window.win_bottom + sensory_distance - 1.;
        let top = run_options.window.win_top - sensory_distance + 1.;

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
        run_options.window = options::get_window_size(1400., 900.);
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
                run_options.window.win_left,
                run_options.window.win_right - 1.,
                run_options.window.win_bottom,
                run_options.window.win_top - 1.,
                settings.x_cell_res,
                settings.y_cell_res,
                settings.x_cell_count,
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

            let retrieved_neighbour_id_set: HashSet<usize> =
                neighbours.iter().map(|nb| nb.id).collect();

            let neighbour_difference: Vec<&Boid> = clicked_neighbours[i]
                // .to_owned()
                .iter()
                .filter(|nb| !retrieved_neighbour_id_set.contains(&nb.id))
                .collect();

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
        run_options.window = options::get_window_size(1400., 900.);
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
