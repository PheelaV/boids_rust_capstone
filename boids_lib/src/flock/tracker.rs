// use std::{collections::HashSet, iter};

// use glam::Vec2;

// use crate::{options::RunOptions, boid::{Boid, BoidMetadata}, math_helpers::distance_dyn_boid};

// use super::get_flock_ids;

// // a tracker takes ownership of a flock of boids
// // a tracker will follow an implementation consisting of choosing a datastructure to
// // hold spatial information about boids in and have traits for querying 2D, later 3D
// // information
// pub trait Tracker {
//     fn new(entities: &[Boid], run_options: &RunOptions) -> Self;
//     fn update(&mut self, run_options: &RunOptions);
//     fn insert_single(&mut self, entity: Boid, _: &RunOptions);
//     fn restart(&mut self, entities: &[Boid]);
//     fn delete_last(&mut self) -> Option<Boid>;
//     fn delete_multiple(&mut self, ids_delete: &[usize], run_options: &RunOptions);
//     fn insert_multiple(&mut self, entities: &[Boid], _: &RunOptions);
//     fn view(&self) -> (&Vec<Boid>, &Vec<BoidMetadata>);
//     fn view2<'a>(&'a self) -> Box<dyn Iterator<Item = (&'a Boid, &'a BoidMetadata)> + 'a>;
// }

// // todo: maybe there is more performance if we stop can discriminate entities that are already in the correct cell
// /// Implementation follows [paper](https://www.researchgate.net/publication/277870601_A_Hash_Table_Construction_Algorithm_for_Spatial_Hashing_Based_on_Linear_Memory)
// /// It has been modified into a form of unstable insert
// pub struct SpatHash1D {
//     // for each cell, there is a pivot with ranges into table
//     // index of pivots corresponds to cell representation
//     pivots: Vec<SpatHashPiv>,
//     // 1D representation of cells
//     // indexes have no special meaning
//     table: Vec<Boid>,
//     // everything that a Boid does not nescesarily need for its own functioning
//     metadata: Vec<BoidMetadata>,
//     // for each boid, records which cell it's at
//     // indexes are effectively boid ids
//     index: Vec<usize>,
//     // contains information about the table, used for updating and reconstruction
//     settings: SpatialHashingTableSettings,
// }

// #[derive(Clone)]
// struct SpatHashPiv {
//     usg: usize,
//     init: Option<usize>,
//     fin: Option<usize>,
// }

// impl Default for SpatHashPiv {
//     fn default() -> Self {
//         SpatHashPiv {
//             usg: 0,
//             init: None,
//             fin: None,
//         }
//     }
// }

// impl Tracker for SpatHash1D {
//     fn new(entities: &[Boid], run_options: &RunOptions) -> Self {
//         let settings = SpatHash1D::get_tracker_settings(run_options);

//         SpatHash1D {
//             // pivots: Vec::with_capacity(settings.cell_count),
//             // initialize vector with both capacity and values prefilled to simplify code in the update_table
//             pivots: vec![Default::default(); settings.cell_count],
//             table: entities.to_vec(),
//             metadata: entities.iter().map(|e| BoidMetadata::new(e)).collect(),
//             // index: Vec::with_capacity(run_options.init_boids as usize),
//             // initialize vector with both capacity and values prefilled to simplify code in the update_table
//             index: vec![0; run_options.init_boids as usize],
//             settings,
//             // accellerations: Vec::with_capacity(run_options.init_boids as usize),
//         }
//     }

//     /// deletes a certain boid from the table and its auxilary arrays,
//     fn delete_multiple(&mut self, ids_delete: &[usize], run_options: &RunOptions) {
//         let ids_delete_set: HashSet<usize> = HashSet::from_iter(ids_delete.iter().cloned());

//         // removing values from the 1D table
//         self.table.retain(|boid| !ids_delete_set.contains(&boid.id));

//         // resize auxilary arrays
//         self.metadata.resize(self.table.len(), Default::default());
//         self.index.resize(self.table.len(), Default::default());

//         // update data structure (index, pivots)
//         self.update_table(run_options);
//     }

//     fn update(&mut self, run_options: &RunOptions) {
//         let new_settings = SpatHash1D::get_tracker_settings(run_options);
//         let cell_difference = self.settings.cell_count as i64 - new_settings.cell_count as i64;
//         if cell_difference != 0 {
//             self.pivots
//                 .resize(new_settings.cell_count, Default::default());
//             self.settings = new_settings;
//             self.update_table(run_options);
//         }
//         // if self.settings.x_cell_res != new_settings.x_cell_res || self.settings.y_cell_res != new_settings.y_cell_res {

//         //     // self.settings = new_settings;
//         //     // pivots: vec![Default::default(); settings.cell_count],
//         //     // self.update_table(run_options);
//         // } else {
//         //     self.settings = new_settings;

//         // }

//         // if run_options.clicked_boid_id != std::usize::MAX {
//         //     let clicked_neighbours = self.get_neighbours(, cell_index, run_options, neighbours)
//         // }

//         self.update_loation(run_options);
//         self.update_table(run_options);

//         if run_options.dbscan_flock_clustering_on {
//             let flock_ids = get_flock_ids(&self.table, run_options);
        
//             self.table
//             .iter()
//             .enumerate()
//             .for_each(|(index, boid)| {
//                 self.metadata[boid.id].cluster_id = flock_ids[index]
//             })
//         }

//     }

//     fn insert_single(&mut self, entity: Boid, run_options: &RunOptions) {
//         self.insert(entity, run_options);

//         // make sure the table updates its structure
//         self.update_table(run_options);
//     }

//     fn restart(&mut self, entities: &[Boid]) {
//         self.table = entities.to_vec();
//     }

//     fn delete_last(&mut self) -> Option<Boid> {
//         let last_id = self.index.len() - 1;
//         match self.table.iter().position(|e| e.id == last_id) {
//             Some(table_index) => Some(self.table.swap_remove(table_index)),
//             None => None,
//         }
//     }

//     fn insert_multiple(&mut self, entities: &[Boid], run_options: &RunOptions) {
//         entities
//             .into_iter()
//             .for_each(|e| self.insert(*e, run_options));

//         // make sure the table updates its structure
//         self.update_table(run_options);
//     }

//     fn view<'a>(&self) -> (&Vec<Boid>, &Vec<BoidMetadata>) {
//         (&self.table, &self.metadata)
//     }

//     // https://stackoverflow.com/questions/31904842/return-a-map-iterator-which-is-using-a-closure-in-rust
//     fn view2<'a>(&'a self) -> Box<dyn Iterator<Item = (&'a Boid, &'a BoidMetadata)> + 'a> {
//         // self.table.iter().map
//         Box::new(self.table.iter().map(|e| (e, &self.metadata[e.id])))
//     }
// }

// // get the cell vectors pointing towards the neighbourhood cell
// // in eucledian space
// const HOME: [i64; 2] = [0, 0];
// const S: [i64; 2] = [0, -1];
// const SW: [i64; 2] = [-1, -1];
// const SE: [i64; 2] = [1, -1];
// const W: [i64; 2] = [-1, 0];
// const E: [i64; 2] = [1, 0];
// const N: [i64; 2] = [0, 1];
// const NE: [i64; 2] = [1, 1];
// const NW: [i64; 2] = [-1, 1];

// impl SpatHash1D {
//     // https://stackoverflow.com/questions/31904842/return-a-map-iterator-which-is-using-a-closure-in-rust
//     // fn view2<'a>(&'a self) -> Box<dyn Iterator<Item= (&'a Boid, &'a BoidMetadata)> + 'a> {
//     //     // self.table.iter().map
//     //     Box::new(self.table.iter().map(|e| (e, &self.metadata[e.id])))
//     // }

//     fn insert(&mut self, entity: Boid, run_options: &RunOptions) {
//         let index = SpatHash1D::get_table_index(
//             entity.position.x,
//             entity.position.y,
//             run_options.window.win_left,
//             run_options.window.win_right - 1.,
//             run_options.window.win_bottom,
//             run_options.window.win_top - 1.,
//             self.settings.x_cell_res,
//             self.settings.y_cell_res,
//             self.settings.x_cell_count,
//         );

//         // expand table and index vectors
//         self.index.push(index);
//         self.metadata.push(BoidMetadata::new(&entity));
//         self.table.push(entity);
//     }

//     pub(self) fn update_loation(&mut self, run_options: &RunOptions) {
//         // reset metadata
//         let metadata_default: BoidMetadata = Default::default();
//         for e in 0..self.table.len() {
//             self.metadata[e].id = e;
//             self.metadata[e].clicked_neighbour_id = metadata_default.clicked_neighbour_id;
//             self.metadata[e].cluster_id = metadata_default.cluster_id;
//             self.metadata[e].n_neighbours = metadata_default.n_neighbours;
//         }
//         // todo: temp print
//         // println!("{:#?}", self.index);
//         // println!("--------");
//         // println!("tracker entity count: {}", self.table.len());

//         let mut neighbours: Vec<&Boid> = Vec::with_capacity(32);
//         let mut accellerations: Vec<Vec2> = Vec::with_capacity(run_options.init_boids);
//         let mut metadata: Vec<BoidMetadata> = vec![Default::default(); run_options.init_boids];

//         let mut clicked_neighbours: Vec<usize> =
//             Vec::with_capacity(run_options.neighbours_cosidered);
//         // todo: this is likely where the performance degrades because of the creation of the vectors above
//         // get neighbours and final accelleration per boid
//         for e in 0..self.table.len() {
//             neighbours.clear();

//             let b: &Boid = &self.table[e];
//             self.get_neighbours(b, self.index[b.id], run_options, &mut neighbours);

//             // neighbours length
//             metadata[self.table[e].id].n_neighbours = neighbours.len();

//             // clicked neighbour
//             if self.table[e].id == run_options.clicked_boid_id {
//                 for n in &neighbours {
//                     // println!("neighbour id: {}", n.id);
//                     clicked_neighbours.push(n.id)
//                 }
//             }
//             // todo: temp print
//             // println!("tracker neighbours: {}", neighbours.len());

//             // this blows up, because the vector is allocated with appropriate capacity, but is empty
//             // accellerations[e] = self.table[e].run_rules3(&neighbours, run_options);
//             // I will attempt to patch it with push to get rid of the overhead in pre-creating the Vec2s
//             // and rely on implicit order, which might cause issues, but should not as the for loop goes
//             //sequentually through the entities
//             accellerations.push(self.table[e].run_rules3(&neighbours, run_options));
//         }

//         clicked_neighbours.iter().for_each(|cn_id| {
//             // println!("cn_id: {}",cn_id);
//             metadata[*cn_id].clicked_neighbour_id = run_options.clicked_boid_id
//         });

//         // update metadata
//         for e in 0..self.index.len() {
//             self.metadata[e].clicked_neighbour_id = metadata[e].clicked_neighbour_id;
//             self.metadata[e].n_neighbours = metadata[e].n_neighbours;
//         }

//         // apply the forces
//         for e in 0..self.table.len() {
//             self.table[e].apply_force(accellerations[e]);
//             self.table[e].update_location(&run_options)
//         }

//         // println!("{:#?}", self.index);
//         // println!("{:#?}", self.metadata);
//         // todo: temp print
//         // println!("tracker first pos: {}", self.table[0].position);
//         // println!("tracker first acc: {}", accellerations.first().unwrap_or(&Vec2::new(0., 0.)));
//         // println!("tracker accs length: {}", accellerations.len());
//         // println!("########");
//         accellerations.clear();
//         // self.update_table(&self.table, run_options)
//     }

//     /// assumes all the auxilary arrays are of correct size
//     pub(self) fn update_table(&mut self, run_options: &RunOptions) {
//         // reset pivot metadata as we are in the next iterration, e.g., we do
//         // not want pivot.usg to linger along from previous iterration
//         self.pivots.iter_mut().for_each(|pivot| {
//             pivot.usg = Default::default();
//             pivot.init = Default::default();
//             pivot.fin = Default::default();
//         });

//         for e in 0..self.table.len() {
//             let index = SpatHash1D::get_table_index(
//                 self.table[e].position.x,
//                 self.table[e].position.y,
//                 run_options.window.win_left,
//                 run_options.window.win_right - 1.,
//                 run_options.window.win_bottom,
//                 run_options.window.win_top - 1.,
//                 self.settings.x_cell_res,
//                 self.settings.y_cell_res,
//                 self.settings.x_cell_count,
//             );

//             // if we are dealing with the same entity, do not tick usage
//             // self.table[e].id ==
//             self.index[self.table[e].id] = index;
//             self.pivots[index].usg += 1;
//         }

//         let mut accum: usize = 0;
//         for p in 0..self.pivots.len() {
//             if self.pivots[p].usg > 0 {
//                 self.pivots[p].init = Some(accum);
//                 accum += self.pivots[p].usg;
//                 self.pivots[p].fin = Some(accum);
//             } else {
//                 self.pivots[p].init = None;
//                 self.pivots[p].fin = None;
//             }
//         }

//         let mut pivots = self.pivots.to_owned();

//         // if self
//         //     .pivots
//         //     .iter()
//         //     .any(|x| x.usg != 0 && (x.fin == None || x.fin.unwrap() == 0))
//         // {
//         //     print!("something fishy going on here!");
//         // }

//         let mut sorted = vec![false; self.table.len()];
//         let mut destination: usize;

//         let mut e = 0;
//         while e < self.table.len() {
//             // here is where the algorithm differs slightly from the original
//             // as we are reusing the same table as the set of entities, we
//             // have to make a swap this collapses into a form of unstable
//             // insertion sort...

//             // self.table[self.pivots[self.index[e]].fin.unwrap() - 1] = self.table[e];
//             // mem::swap(&mut self.table[self.pivots[self.index[e]].fin.unwrap() - 1], &mut self.table[e]);
//             let id = self.table[e].id;

//             // if id == std::usize::MAX {
//             //     println!("something fishy about id");
//             // }
//             // if self.pivots[self.index[id]].fin.unwrap() == 0 {
//             //     println!("something fishy about usize");
//             // }

//             destination = usize::wrapping_sub(pivots[self.index[id]].fin.unwrap(), 1);

//             // skip if the destination has alreade been sorted, current table[e] is in right place
//             // not checking for 0 but for MAX (18446744073709551615) as we take advantage of wraparound
//             if destination == std::usize::MAX || sorted[destination] || e == destination {
//                 e += 1;
//                 continue;
//             }
//             // self.table[e].;
//             self.table.swap(e, destination);

//             // update the pivot now that one entity of the cell has been mapped to the correct place in the 1D table
//             // each unique value of fin during the iterration effectively represents a unique location for the next
//             // entity in that cell
//             // we consume it during swap and then decrease to point at the location of cell's next entity
//             // until it points at the first(swapped, last in terms of the 1D table) entity of the previous cell
//             pivots[self.index[id]].fin = Some(pivots[self.index[id]].fin.unwrap() - 1);
//             // if self.pivots.iter().any(|x| x.usg != 0 && (x.fin == None || x.fin.unwrap() == 0)) {
//             //     print!("something fishy going on here!");
//             // }

//             // set destination as sorted
//             sorted[destination] = true;
//         }
//     }

//     /// Returns "hashed" value representing an index for spatial subdivision, handles a zero centered coordinate system
//     ///
//     /// ## Arguments
//     ///
//     /// - `p_x` - x coordinate of a point of interest
//     /// - `p_y` - y coordinate of a point of interest
//     /// - `min_x` - minimum x axis value
//     /// - `min_y` - minimum y axis value
//     /// - `cs_x` - cell size along the x axis
//     /// - `cs_y` - cell size along the y axis
//     /// - `table_width` - the table width, how many cells make up a row
//     ///
//     /// ## Returns
//     /// an index as _usize_, that takes values from 0..(x_cells * y_cells)
//     ///
//     /// **watch out**, in case of a symmetric/centered coordinate system (like nannou)
//     /// the recommendation is to subtract 1 from the x/y max values to reduce the
//     /// range of returned indexes by one as the centering will cause the dimensions
//     /// to have uneven counts, returning index that would be out of bounds for the
//     /// grid
//     pub fn get_table_index(
//         p_x: f32,
//         p_y: f32,
//         min_x: f32,
//         max_x: f32,
//         min_y: f32,
//         max_y: f32,
//         cs_x: f32,
//         cs_y: f32,
//         table_width: f32,
//     ) -> usize {
//         (((p_x.clamp(min_x, max_x) - min_x) / cs_x).floor()
//             + (((p_y.clamp(min_y, max_y) - min_y) / cs_y).floor() * table_width)) as usize
//     }
//     pub fn get_table_index2(
//         p_x: f32,
//         p_y: f32,
//         min_x: f32,
//         max_x: f32,
//         min_y: f32,
//         max_y: f32,
//         cs_x: f32,
//         cs_y: f32,
//         table_width: f32,
//     ) -> usize {
//         (((p_x.clamp(min_x, max_x) - min_x) / cs_x).ceil()
//             + (((p_y.clamp(min_y, max_y) - min_y) / cs_y).ceil() * table_width)) as usize
//     }

//     pub fn get_tracker_settings(run_options: &RunOptions) -> SpatialHashingTableSettings {
//         let x_range = (run_options.window.win_right - run_options.window.win_left).ceil();
//         let y_range = (run_options.window.win_top - run_options.window.win_bottom).ceil();

//         // let x_cell_count = (x_range as f32 / (100.)).ceil();

//         // here we multiply the sensory distance by two, because it is a radius and the method
//         // for fetching neighbours is implemented with the assumption that the sensory radius
//         // fits into one cell of the grid
//         // value any higher than that will an unnescessary ammount of boids
//         // any lower, and it won't fetch all boids that should have been considered

//         let x_cell_count: f32;
//         let y_cell_count: f32;

//         // the problems I am hitting here are related to the gauss circle problem
//         // which consists of counting the number of lattice points N(r) inside the boundary of
//         // a circle of radius r with center at the origin
//         // there is an exact solution: https://mathworld.wolfram.com/GausssCircleProblem.html
//         // - considers the floored radius
//         // tests all discrete points within positive quartal
//         if x_range.max(y_range) / 2. < run_options.max_sensory_distance {
//             x_cell_count = 1.;
//             y_cell_count = 1.;
//         } else {
//             x_cell_count = (x_range as f32 / (run_options.max_sensory_distance * 2.)).ceil();
//             y_cell_count = (y_range as f32 / (run_options.max_sensory_distance * 2.)).ceil();
//         }
//         // let x_cell_count = (x_range as f32 / (run_options.max_sensory_distance * 2.)).ceil();
//         // let y_cell_count = x_cell_count;
//         // let x_cell_count = 4.;
//         // let y_cell_count = (y_range as f32 / (100.)).ceil();
//         // let y_cell_count = (y_range as f32 / (run_options.max_sensory_distance * 2.)).ceil();
//         // let y_cell_count = 4.;

//         let cell_count = (x_cell_count * y_cell_count) as usize;

//         let x_cell_res = x_range / x_cell_count;
//         let y_cell_res = y_range / y_cell_count;

//         SpatialHashingTableSettings {
//             x_range,
//             y_range,
//             x_cell_count,
//             y_cell_count,
//             cell_count,
//             x_cell_res,
//             y_cell_res,
//         }
//     }

//     pub(self) fn get_neighbours<'a>(
//         &'a self,
//         boid: &Boid,
//         cell_index: usize,
//         run_options: &RunOptions,
//         neighbours: &mut Vec<&'a Boid>,
//     ) {
//         // only add it once
//         let is_left = cell_index % self.settings.x_cell_count as usize == 0;
//         let is_right = (cell_index % self.settings.x_cell_count as usize) + 1
//             == self.settings.x_cell_count as usize;
//         let is_bottom = cell_index / self.settings.x_cell_count as usize == 0;
//         let is_top =
//             cell_index as f32 >= (self.settings.y_cell_count - 1f32) * self.settings.x_cell_count;
//         if boid.id == run_options.clicked_boid_id {
//             println!("cell_index: {}", cell_index);
//             println!(
//                 "x_y: {x:04.1}_{y:04.1}",
//                 x = boid.position.x,
//                 y = boid.position.y
//             );
//         }
//         // let matched_cells = match (is_left, is_right, is_bottom, is_top) {
//         //     (false, false, false, false) => [NW, N, NE, W, HOME, E, SW, S, SE],
//         //     (true, false, true, false) => [HOME, N, NE, HOME, HOME, E, HOME, HOME, HOME],
//         //     (true, false, false, true) => [HOME, HOME, HOME, HOME, HOME, E, HOME, S, S],
//         //     (true, false, false, false) => [HOME, N, NE, HOME, HOME, E, HOME, S, SE],
//         //     (false, true, true, false) => [NW, N, HOME, W, HOME, HOME, HOME, HOME, HOME],
//         //     (false, true, false, true) => [HOME, HOME, HOME, W, HOME, HOME, SW, S, HOME],
//         //     (false, true, false, false) => [NW, N, HOME, W, HOME, HOME, SW, S, HOME],
//         //     (false, false, true, false) => [NW, N, NE, W, HOME, E, HOME, HOME, HOME],
//         //     (false, false, false, true) => [HOME, HOME, HOME, E, HOME, W, SW, S, SE],
//         //     _ => panic!("it should not have come to this :"),
//         // };
//         // matched_cells
//         //     .iter()
//         //     .filter(|l| **l != HOME)
//         //     // always consider the boid's current 'home' vector
//         //     .chain(iter::once(&HOME))
//         //     // convert vectors to the table's 1D indexes pointing to cells
//         //     .map(|lookup| {
//         //         // println!("{}, {}", lookup[0], lookup[1]);
//         //         let cell =
//         //             (cell_index as i64 + lookup[0] + lookup[1] * self.settings.x_cell_count as i64)
//         //                 as usize;
//         //         cell
//         //     })
//         //     .filter(|cell| self.pivots[*cell].usg != 0)
//         //     .for_each(|cell| {
//         //         // println!("range_cell: {}", cell);
//         //         // acum += 1;

//         //         // if self.pivots[cell].usg == 0 {
//         //         //     break;
//         //         // }
//         //         for index in self.pivots[cell].init.unwrap()
//         //             ..(self.pivots[cell].init.unwrap() + self.pivots[cell].usg)
//         //         {
//         //             if self.table[index].id != boid.id
//         //                 && distance_dyn_boid(boid, &self.table[index], run_options)
//         //                     <= run_options.max_sensory_distance
//         //             {
//         //                 neighbours.push(&self.table[index]);
//         //                 if neighbours.len() == run_options.neighbours_cosidered {
//         //                     break;
//         //                 }
//         //             }
//         //         }
//         //     })

//         // depending on the match, iterate the vectors
//         // .iter()
//         for cell in match (is_left, is_right, is_bottom, is_top) {
//             (false, false, false, false) => [NW, N, NE, W, HOME, E, SW, S, SE],
//             (true, false, true, false) => [HOME, N, NE, HOME, HOME, E, HOME, HOME, HOME],
//             (true, false, false, true) => [HOME, HOME, HOME, HOME, HOME, E, HOME, S, S],
//             (true, false, false, false) => [HOME, N, NE, HOME, HOME, E, HOME, S, SE],
//             (false, true, true, false) => [NW, N, HOME, W, HOME, HOME, HOME, HOME, HOME],
//             (false, true, false, true) => [HOME, HOME, HOME, W, HOME, HOME, SW, S, HOME],
//             (false, true, false, false) => [NW, N, HOME, W, HOME, HOME, SW, S, HOME],
//             (false, false, true, false) => [NW, N, NE, W, HOME, E, HOME, HOME, HOME],
//             (false, false, false, true) => [HOME, HOME, HOME, E, HOME, W, SW, S, SE],
//             _ => panic!("it should not have come to this :"),
//         }
//         // depending on the match, iterate the vectors
//         .iter()
//         // filter out the out-of-bound vectors, the home vector is used as a placeholder
//         .filter(|l| **l != HOME)
//         // always consider the boid's current 'home' vector
//         .chain(iter::once(&HOME))
//         // convert vectors to the table's 1D indexes pointing to cells
//         .map(|lookup| {
//             // println!("{}, {}", lookup[0], lookup[1]);
//             (cell_index as i64 + lookup[0] + lookup[1] * self.settings.x_cell_count as i64) as usize
//         }) {
//             // println!("range_cell: {}", cell);
//             // acum += 1;

//             if self.pivots[cell].usg == 0 {
//                 continue;
//             }
//             for index in self.pivots[cell].init.unwrap()
//                 ..(self.pivots[cell].init.unwrap() + self.pivots[cell].usg)
//             {
//                 if self.table[index].id != boid.id
//                     && distance_dyn_boid(boid, &self.table[index], run_options)
//                         <= run_options.max_sensory_distance
//                 {
//                     neighbours.push(&self.table[index]);
//                     if neighbours.len() == run_options.neighbours_cosidered {
//                         break;
//                     }
//                 }
//             }
//         }
//         // todo: temp print
//         // for _ in 0..(9 - acum){
//         //     println!(" - ");
//         // }
//     }
// }
// /// A naive implementation of boids tracking, which uses an O(N^2) algorithm for
// /// finding boid's neighbours.
// pub struct BoidTracker {
//     boids: Vec<Boid>,
//     metadata: Vec<BoidMetadata>,
// }

// impl Tracker for BoidTracker {
//     fn new(entities: &[Boid], run_options: &RunOptions) -> Self {
//         BoidTracker {
//             boids: entities.to_vec(),
//             metadata: vec![Default::default(); run_options.init_boids],
//         }
//     }

//     fn update(&mut self, run_options: &RunOptions) {
//         let mut accelleration: Vec<Vec2> = Vec::with_capacity(self.boids.len());
//         let mut metadata: Vec<BoidMetadata> = vec![Default::default(); self.boids.len()];
//         let mut neighbours: Vec<&Boid> = Vec::new();
//         let mut clicked_neighbours: Vec<&Boid> =
//             Vec::with_capacity(run_options.neighbours_cosidered);

//         // calculation loop
//         for i_cur in 0..self.boids.len() {
//             neighbours.clear();

//             let b_current = &self.boids[i_cur];

//             metadata[i_cur].id = b_current.id;

//             for b_other in self.boids.iter() {
//                 if b_other.id == b_current.id {
//                     continue;
//                 }

//                 let distance = distance_dyn_boid(&b_current, b_other, &run_options);
//                 if distance < run_options.max_sensory_distance {
//                     neighbours.push(b_other);
//                 }

//                 if b_current.id == run_options.clicked_boid_id {
//                     for n in neighbours.iter() {
//                         // metadata == n.id;
//                         metadata[n.id].clicked_neighbour_id = b_current.id;
//                         clicked_neighbours.push(n);
//                     }
//                 }
//             }
//             // todo: temp print
//             // if self.boids[i_cur].id == 0 {
//             //     println!("flock neighbours: {}", neighbours.len());
//             // }
//             accelleration.push(self.boids[i_cur].run_rules3(&neighbours, &run_options));

//             // todo: temp print
//             // println!("flock first pos: {}", self.boids[0].position);

//             for cn in &clicked_neighbours {
//                 metadata[cn.id].clicked_neighbour_id = run_options.clicked_boid_id;
//             }

//             // update metadata
//             for e in 0..self.metadata.len() {
//                 self.metadata[e].id = e;
//                 self.metadata[e].cluster_id = Default::default();
//                 self.metadata[e].clicked_neighbour_id = metadata[e].clicked_neighbour_id;
//                 self.metadata[e].n_neighbours = metadata[e].n_neighbours;
//             }
//         }
//         // todo: temp print
//         // println!("flock first acc: {}", accelleration.first().unwrap_or(&Vec2::new(0., 0.)));
//         // println!("flock accs length: {}", accelleration.len());
//         // update loop
//         for i_cur in 0..self.boids.len() {
//             self.boids[i_cur].apply_force(accelleration[i_cur]);
//             self.boids[i_cur].update_location(&run_options)
//         }

//         if run_options.dbscan_flock_clustering_on {
//             let flock_ids = get_flock_ids(&self.boids, run_options);
        
//             self.boids.iter()
//             .enumerate()
//             .for_each(|(index, boid)| {
//                 self.metadata[boid.id].cluster_id = flock_ids[index]
//             })
//         }
//     }

//     fn insert_single(&mut self, entity: Boid, _: &RunOptions) {
//         let mut metadata: BoidMetadata = Default::default();
//         metadata.id = entity.id;
//         self.metadata.push(metadata);

//         self.boids.push(entity);
//     }

//     fn restart(&mut self, entities: &[Boid]) {
//         self.boids = entities.to_vec();
//     }

//     fn insert_multiple(&mut self, entities: &[Boid], _: &RunOptions) {
//         self.boids.extend_from_slice(entities);

//         self.metadata.extend(entities.iter().map(|e| {
//             let mut metadata: BoidMetadata = Default::default();
//             metadata.id = e.id;
//             metadata
//         }))
//     }

//     fn delete_last(&mut self) -> Option<Boid> {
//         if self.boids.len() == 0 {
//             None
//         } else {
//             let index = self.boids.len() - 1;
//             self.metadata.swap_remove(index);
//             Some(self.boids.swap_remove(index))
//         }
//     }

//     fn view(&self) -> (&Vec<Boid>, &Vec<BoidMetadata>) {
//         (&&self.boids, &self.metadata)
//     }

//     // https://stackoverflow.com/questions/31904842/return-a-map-iterator-which-is-using-a-closure-in-rust
//     fn view2<'a>(&'a self) -> Box<dyn Iterator<Item = (&'a Boid, &'a BoidMetadata)> + 'a> {
//         // self.table.iter().map
//         Box::new(self.boids.iter().map(|e| (e, &self.metadata[e.id])))
//     }

//     fn delete_multiple(&mut self, ids_delete: &[usize], run_options: &RunOptions) {
//         let ids_delete_set: HashSet<usize> = HashSet::from_iter(ids_delete.iter().cloned());

//         // removing values from the 1D table
//         self.boids.retain(|boid| !ids_delete_set.contains(&boid.id));

//         self.metadata.resize(self.boids.len(), Default::default())
//     }
// }
// /// Uses a spatial hashing space division method, where all cells of the underlying
// /// table are stored in a 1D array, with the individual cell's being allocated
// /// dynamically, reused or updated as needed.
// ///
// /// Uses an Θ(n) algorithm for finding the boid's neighbours.
// /// ref: Optimization_of_large-scale_real-time_simulations
// pub struct SpatialHashingTableSettings {
//     pub x_range: f32,
//     pub y_range: f32,
//     pub x_cell_count: f32,
//     pub y_cell_count: f32,
//     pub cell_count: usize,
//     pub x_cell_res: f32,
//     pub y_cell_res: f32,
//     // run_options: &'a RunOptions
// }

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn test_retreive_neighbour(){

//     }
// }