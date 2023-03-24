use std::{collections::{HashSet, HashMap}, f32::consts::PI, iter, any::Any};

use glam::Vec2;
use itertools::Itertools;
use rand::Rng;

use crate::{options::{RunOptions, Distance}, boid::{BoidMetadata, Boid}, math_helpers::distance_dyn_boid};

use super::{Tracker, MY_RNG, get_flock_ids, naive_tracker::NaiveTracker, tracker::TrackerSignal};

/// Uses a spatial hashing space division method, where all cells of the underlying
/// table are stored in a 1D array, with the individual cell's being allocated
/// dynamically, reused or updated as needed.
///
/// Uses an Θ(n) algorithm for finding the boid's neighbours.
/// ref: Optimization_of_large-scale_real-time_simulations
pub struct SpatialHashingTableSettings {
    /// environment x range
    pub x_range: usize,
    /// environment y range
    pub y_range: usize,
    /// grid x range
    pub x_cell_count: usize,
    /// grid y range
    pub y_cell_count: usize,
    /// n cells of grid
    pub cell_count: usize,
    /// environment x units per grid x cell
    pub x_cell_res: f32,
    /// environment y units per grid y cell
    pub y_cell_res: f32,
    // run_options: &'a RunOptions
}

// todo: maybe there is more performance if we stop can discriminate entities that are already in the correct cell
/// Implementation follows [paper](https://www.researchgate.net/publication/277870601_A_Hash_Table_Construction_Algorithm_for_Spatial_Hashing_Based_on_Linear_Memory)
/// It has been modified into a form of an unstable insertion sort
pub struct SpatHash1D {
    // for each cell, there is a pivot with ranges into table
    // index of pivots corresponds to cell representation
    pivots: Vec<SpatHashPiv>,
    /// 1D representation of cells.
    /// Indexes have no special meaning, as a matter of fact they will be dependent
    /// on spatial information and are not stable.
    table: Vec<Boid>,
    /// This works in conjunction with table, provides references for boids in table.
    /// Indexes are agennt IDs.
    view: Vec<usize>,
    /// Everything that a Boid does not nescesarily need for its own functioning.
    /// Indexes are agent IDs.
    metadata: Vec<BoidMetadata>,
    /// For each boid, records its cell index.
    /// Indexes are agent IDs.
    index: Vec<usize>,
    /// Contains information about the table, used for updating and reconstruction of spatial
    /// index.
    pub(crate)settings: SpatialHashingTableSettings,
    /// Querying metadata for neighbourhood, gives lookup vectors depending on the current
    /// environment and table dimensions
    query_metadata: Option<[[[i32; 2]; 9]; 16]>,

    distance: Distance,
    // clicked_neighbours: Option<Vec::<usize>>
}

#[derive(Clone, Debug)]
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
        let metadata: Vec<BoidMetadata> = entities.iter().map(|e| BoidMetadata::new(e)).collect();
        // metadata[0].boid_type = BoidType::Disruptor;

        SpatHash1D {
            // initialize vector with both capacity and values prefilled to simplify code in the update_table
            pivots: vec![Default::default(); settings.cell_count],
            table: entities.to_vec(),
            metadata,
            // initialize vector with both capacity and values prefilled to simplify code in the update_table
            index: vec![0; run_options.init_boids as usize],
            settings,
            query_metadata: None,
            distance: run_options.distance,
            view: entities.iter().map(|b| b.id).collect_vec(),
        }
    }

    fn update(&mut self, run_options: &RunOptions) {
        // check if table grid needs updating (expanding or shrinking)
        let new_settings = SpatHash1D::get_tracker_settings(run_options);
        let table_resized = self.settings.cell_count as i64 - new_settings.cell_count as i64 != 0;
        if table_resized {
            self.pivots
                .resize(new_settings.cell_count, Default::default());
            self.settings = new_settings;
            self.update_table(run_options);
        }

        // check if we need to regenerate query metadata for toroidial lookup
        // of the table grid, that is when:
        //   distance has changed
        //   or table was resized
        //   or we begun with toroidial distance
        if self.distance != run_options.distance
            || table_resized
            || (run_options.distance == Distance::EucToroidal && self.query_metadata == None)
        {
            self.distance = run_options.distance;
            self.query_metadata = if run_options.distance == Distance::EucEnclosed {
                None
            } else {
                Some(self.get_query_metadata())
            };
        }

        self.update_location(run_options);
        self.update_table(run_options);

        if run_options.dbscan_flock_clustering_on {
            // make a copy of the boid data as order very much matters in terms of density based clustering
            // and we do not want to disturb the order in table
            let data: Vec<Boid> = self
                .table
                .iter()
                .map(|b| b.to_owned())
                .sorted_by(|b1, b2| b1.id.cmp(&b2.id))
                .collect();

            let mut clust_res_map: HashMap<usize, usize> = HashMap::new();

            let flock_ids: Vec<usize> =get_flock_ids(self, &data, run_options);

            for (b, c) in data.iter().zip(flock_ids.iter()) {
                clust_res_map.insert(b.id, *c);
            }

            self.table.iter().for_each(|boid| {
                self.metadata[boid.id].cluster_id = *clust_res_map.get(&boid.id).unwrap_or(&0)
            })
        }
    }

    fn insert_single(&mut self, entity: Boid, run_options: &RunOptions) {
        self.insert(entity, run_options);

        // make sure the table updates its structure
        self.update_table(run_options);
    }

    fn restart(&mut self, entities: &[Boid]) {
        self.table = entities.to_vec();
        self.view = entities.iter().map(|e| e.id).collect_vec();
    }

    fn delete_last(&mut self) -> Option<Boid> {
        let last_id = self.index.len() - 1;
        match self.table.iter().position(|e| e.id == last_id) {
            Some(table_index) => {
                self.view.swap_remove(last_id);
                Some(self.table.swap_remove(table_index))
            },
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

    // https://stackoverflow.com/questions/31904842/return-a-map-iterator-which-is-using-a-closure-in-rust
    fn view<'a>(&'a self) -> Box<dyn Iterator<Item = (&'a Boid, &'a BoidMetadata)> + 'a> {
        Box::new(self.metadata.iter().map(|m| (&self.table[self.view[m.id]], m)))
    }

    fn get_neighbours<'a>(&'a self, boid: &Boid, run_options: &RunOptions) -> Vec<&'a Boid> {
        let cell_index = SpatHash1D::get_table_index(
            boid.position.x,
            boid.position.y,
            run_options.window.win_left as f32,
            (run_options.window.win_right - 1) as f32,
            run_options.window.win_bottom as f32,
            (run_options.window.win_top - 1) as f32,
            self.settings.x_cell_res,
            self.settings.y_cell_res,
            self.settings.x_cell_count as f32,
        );
        let mut result = Vec::<&Boid>::new();

        self.get_neighbours(boid, cell_index, run_options, &mut result);

        result
    }

    fn signal(&mut self, _: TrackerSignal) {
        () // noop
    }

    fn get_no_entities(&self) -> usize {
       self.table.len() 
    }

    fn as_any(&self) -> &dyn Any {
        self
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
        // placeholder
        [[0; 2]; 9],
    ];

    // set of helper functions for toroidially wrapped environment space

    /// if current cell is left (W), returns vector pointing towards it's W neighbour
    /// according to toroidal wrapping
    ///
    /// 'vo' - vertical offset in (-1..=1)
    fn gt_l(&self, vo: i32) -> [i32; 2] {
        [self.settings.x_cell_count as i32 - 1, vo]
    }

    /// if current cell is right, returns the cell that would have been on the right
    /// according to toroidal wrappin
    ///
    /// 'vo' - vertical offset in (-1..=1)g
    fn gt_r(&self, vo: i32) -> [i32; 2] {
        [-(self.settings.x_cell_count as i32) + 1, vo]
    }

    /// if current cell is topdners, returns the cell that would have been on the top
    /// according to toroidal wrappin
    ///
    /// 'ho' - horizontal offset in (-1..=1)
    fn gt_t(&self, ho: i32) -> [i32; 2] {
        [ho, -(self.settings.y_cell_count as i32) + 1]
    }

    /// if current cell is bottom, returns the cell that would have been on the bottom
    /// according to toroidal wrapping
    ///
    /// 'ho' - horizontal offset in (-1..=1)
    fn gt_b(&self, ho: i32) -> [i32; 2] {
        [ho, (self.settings.y_cell_count as i32) - 1]
    }

    /// if current cell is left top, returns the cell that would have been on the
    /// following the  left top diagonal
    fn gt_lt(&self) -> [i32; 2] {
        [self.gt_l(0)[0], self.gt_t(0)[1]]
    }

    /// if current cell is left bottom, returns the cell that would have been on the
    /// following the  left bottom diagonal
    fn gt_lb(&self) -> [i32; 2] {
        [self.gt_l(0)[0], self.gt_b(0)[1]]
    }

    /// if current cell is right top, returns the cell that would have been on the
    /// following the  roght top diagonal
    fn gt_rt(&self) -> [i32; 2] {
        [self.gt_r(0)[0], self.gt_t(0)[1]]
    }

    /// if current cell is right bottom, returns the cell that would have been on the
    /// following the  right bottom diagonal
    fn gt_rb(&self) -> [i32; 2] {
        [self.gt_r(0)[0], self.gt_b(0)[1]]
    }

    // https://stackoverflow.com/questions/31904842/return-a-map-iterator-which-is-using-a-closure-in-rust
    // fn view2<'a>(&'a self) -> Box<dyn Iterator<Item= (&'a Boid, &'a BoidMetadata)> + 'a> {
    //     // self.table.iter().map
    //     Box::new(self.table.iter().map(|e| (e, &self.metadata[e.id])))
    // }

    #[rustfmt::skip]
    fn get_query_metadata(&self) -> [[[i32; 2]; 9]; 16] {
        [
            // 0 = middle
            [Self::NW, Self::N, Self::NE, Self::W, Self::HOME, Self::E, Self::SW, Self::S, Self::SE],
            // 1 = top
            [self.gt_t(-1), self.gt_t(0), self.gt_t(1), Self::E, Self::HOME, Self::W, Self::SW, Self::S, Self::SE],
            // 2 = bottom
            [Self::NW, Self::N, Self::NE, Self::W, Self::HOME, Self::E, self.gt_b(-1), self.gt_b(0), self.gt_b(1)],
            // placeholder
            [[0; 2]; 9],
            // 4 = right
            [Self::NW, Self::N, self.gt_r(1), Self::W, Self::HOME, self.gt_r(0), Self::SW, Self::S, self.gt_r(-1)],
            // 5 = right | top
            [self.gt_t(-1), self.gt_t(0), self.gt_rt(), Self::W, Self::HOME, self.gt_r(0), Self::SW, Self::S, self.gt_r(-1)],
            // 6 = right | bottom
            [Self::NW, Self::N, self.gt_r(1), Self::W, Self::HOME, self.gt_r(0), self.gt_b(-1), self.gt_b(0), self.gt_rb()],
            // placeholder
            [[0; 2]; 9],
            // 8 = left
            [self.gt_l(1), Self::N, Self::NE, self.gt_l(0), Self::HOME, Self::E, self.gt_l(-1), Self::S, Self::SE],
            // 9 = left | top
            [self.gt_lt(), self.gt_t(0), self.gt_t(1), self.gt_l(0), Self::HOME, Self::E, self.gt_l(-1), Self::S, Self::SE],
            // 10 = left | botom
            [self.gt_l(1), Self::N, Self::NE, self.gt_l(0), Self::HOME, Self::E, self.gt_lb(), self.gt_b(0), self.gt_b(1)],
            // placeholder
            [[0; 2]; 9],
            // placeholder
            [[0; 2]; 9],
            // placeholder
            [[0; 2]; 9],
            // placeholder
            [[0; 2]; 9],
            // placeholder
            [[0; 2]; 9],
        ]
    }

    fn insert(&mut self, entity: Boid, run_options: &RunOptions) {
        let index = SpatHash1D::get_table_index(
            entity.position.x,
            entity.position.y,
            run_options.window.win_left as f32,
            (run_options.window.win_right - 1) as f32,
            run_options.window.win_bottom as f32,
            (run_options.window.win_top - 1) as f32,
            self.settings.x_cell_res,
            self.settings.y_cell_res,
            self.settings.x_cell_count as f32,
        );

        // expand table and index vectors
        self.index.push(index);
        self.metadata.push(BoidMetadata::new(&entity));
        self.view.push(self.table.len());
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
            // simply is not built for, so here is a workaround
            // the gist is: "Reynolds" noise is stateful and needs to write metadata 
            if run_options.wander_on {
                let wander_next = MY_RNG
                    .lock()
                    .unwrap()
                    .gen_range(-(run_options.wander_rate)..(run_options.wander_rate)) * PI;
                metadata[self.table[e].id].wander_direction = match run_options.noise_model {
                    crate::options::NoiseModel::Vicsek => wander_next,
                    crate::options::NoiseModel::Reynolds => (self.metadata[self.table[e].id].wander_direction + wander_next) % (2. * PI),
                }
            }

            metadata[self.table[e].id].n_neighbours = neighbours.len();

            // clicked neighbour
            if self.table[e].id == run_options.clicked_boid_id {
                for n in &neighbours {
                    clicked_neighbours.push(n.id)
                }
            }

            // todo: this is 
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
            self.metadata[id].wander_direction = metadata[id].wander_direction;
        }

        // apply the forces
        for e in 0..self.table.len() {
            self.table[e].apply_force(accellerations[e]);
            self.table[e].update_location(&run_options)
        }
        accellerations.clear();
    }

    /// updates table grid, assumes all the auxilary arrays are of correct size
    pub(crate) fn update_table(&mut self, run_options: &RunOptions) {
        if self.distance != run_options.distance {}
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
                run_options.window.win_left as f32,
                (run_options.window.win_right - 1) as f32,
                run_options.window.win_bottom as f32,
                (run_options.window.win_top - 1) as f32,
                self.settings.x_cell_res,
                self.settings.y_cell_res,
                self.settings.x_cell_count as f32,
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

        // Now we finally sort through the existing agents using newly created pivots and indeces
        while e < self.table.len() {

            // it has already been taken care of, skip
            if sorted[e] {
                e += 1;
                continue;
            }

            let id = self.table[e].id;

            // destination is the last range of a given pivot
            // one pivot represents a cell
            // when we ocupy a cell's range with one agent or find an agent already in place
            // we decrease the cell's range -> since it is a copy, we are keeping track
            // of which agen'ts need yet to
            destination = match pivots[self.index[id]].fin {
                Some(d) => {
                    if d == 0 {
                        panic!("We are trying to sort through too many units!");
                    } else {
                        d - 1
                    }
                },
                None => {
                    panic!("Fin was none, invalid state!"); 
                },
            };

            // if boid is not already sorted
            if e != destination
            {
                // swap the agent into place
                self.table.swap(e, destination);
                // // make a record of where he is
                // self.view[self.table[destination].id] = destination;
                // mark as done
                sorted[destination] = true;
            } else {
                // check if this index has not been put in place
                // todo: this might not be required
                // if !sorted[e] {
                //     self.view[self.table[e].id] = e;
                // }
                e += 1;
            }
            
            pivots[self.index[id]].fin = Some(destination);
        }

        for e in 0..self.table.len() {
            self.view[self.table[e].id] = e;
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
        let x_range = (run_options.window.win_right - run_options.window.win_left) as usize;
        let y_range = (run_options.window.win_top - run_options.window.win_bottom) as usize;

        // let x_cell_count = (x_range as f32 / (100.)).ceil();

        // here we multiply the sensory distance by two, because it is a radius and the method
        // for fetching neighbours is implemented with the assumption that the sensory radius
        // fits into one cell of the grid
        // value any higher than that will an unnescessary ammount of boids
        // any lower, and it won't fetch all boids that should have been considered

        let mut x_cell_count: usize;
        let mut y_cell_count: usize;

        // the problems I am hitting here are related to the gauss circle problem
        // which consists of counting the number of lattice points N(r) inside the boundary of
        // a circle of radius r with center at the origin
        // there is an exact solution: https://mathworld.wolfram.com/GausssCircleProblem.html
        // - considers the floored radius
        // tests all discrete points within positive quartal
        if x_range.max(y_range) as f32 / 2. < run_options.max_sensory_distance {
            x_cell_count = 1;
            y_cell_count = 1;
        } else {
            x_cell_count =
                x_range as usize / (run_options.max_sensory_distance.ceil() as usize);
            y_cell_count =
                y_range as usize / (run_options.max_sensory_distance.ceil() as usize);
            // x_cell_count =
            //     x_range as usize / (run_options.max_sensory_distance.ceil() as usize * 2);
            // y_cell_count =
            //     y_range as usize / (run_options.max_sensory_distance.ceil() as usize * 2);
        }

        if y_cell_count < 2 || x_cell_count < 2 {
            y_cell_count = 1;
            x_cell_count = 1;
        }

        let cell_count = (x_cell_count * y_cell_count) as usize;

        let x_cell_res = x_range as f32 / x_cell_count as f32;
        let y_cell_res = y_range as f32 / y_cell_count as f32;

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

    pub(crate) fn get_neighbours<'a>(
        &'a self,
        boid: &Boid,
        cell_index: usize,
        run_options: &RunOptions,
        neighbours: &mut Vec<&'a Boid>,
    ) {
        // if the table grid is less than 4x4
        if self.settings.x_cell_count <= 3 && self.settings.y_cell_count <= 3 {
            // fall back to naive method
            return NaiveTracker::get_neighbours_naive(boid, &self.table, run_options, neighbours);
        }

        // only add it once
        let is_left = cell_index % self.settings.x_cell_count as usize == 0;
        let is_right = (cell_index % self.settings.x_cell_count as usize) + 1
            == self.settings.x_cell_count as usize;
        let is_bottom = cell_index / self.settings.x_cell_count as usize == 0;
        let is_top = cell_index >= (self.settings.y_cell_count - 1) * self.settings.x_cell_count;

        // self.settings
        // if there is any piece of code in this code base that deserves an internet trophy it is this
        // and LOOKUP - unreadable, unreasonable, non-idiomatic, hella fun and blazingly fast?
        let lookup_index = {
            (is_left as u8) << 3 | (is_right as u8) << 2 | (is_bottom as u8) << 1 | is_top as u8
        } as usize;
        // depending on the match, iterate the vectors, presuming, serves as a
        // lookup for a potential (1..3)*(1..3) area on a n*n grid

        let m = match (self.query_metadata, self.distance) {
            (Some(metadata), Distance::EucToroidal) => metadata[lookup_index],
            (None, Distance::EucEnclosed) => Self::LOOKUP[lookup_index],
            (_, _) => {
                panic!("it should have never come to this")
            }
        };

        // #[cfg(debug)]
        // if boid.id == run_options.clicked_boid_id {
        //     println !("lookup vectors: {:?}", m);
        // }

        for cell in m
            // depending on the match, iterate the vectors
            .iter()
            // filter out the out-of-bound vectors, the home vector is used as a placeholder
            .filter(|l| **l != Self::HOME)
            // always consider the boid's current 'home' vector
            .chain(iter::once(&Self::HOME))
            // convert vectors to the table's 1D indexes pointing to cells
            .map(|l| (cell_index as i32 + l[0] + l[1] * self.settings.x_cell_count as i32) as usize)
        {
            // if run_options.clicked_boid_id == boid.id {
            //     println!("{:?}", cell);
            // }
            // if cell > 60 {
            //     panic!();
            // }
            if self.pivots[cell].usg == 0 {
                continue;
            }
            for index in self.pivots[cell].init.unwrap()..self.pivots[cell].fin.unwrap() {
                // #[cfg(debug)]
                // if run_options.clicked_boid_id == boid.id && self.table[index].id != boid.id {
                //     let a = distance_dyn_boid(boid, &self.table[index], run_options);
                //     println!("from:{:?}, to:{:?}, distance:{:?}", boid.id, self.table[index].id, a);
                // }
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

        // #[cfg(debug)]
        // if run_options.clicked_boid_id == boid.id {
        //     println!("neighbours: {:?}", neighbours);
        // }
    }
}