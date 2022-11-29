use std::env::home_dir;
use std::iter;

use glam::Vec2;
use linfa::traits::Transformer;
use linfa::DatasetBase;
use linfa_clustering::Dbscan;
use ndarray::Array2;
use rand::Rng;

use crate::boid::Boid;
use crate::math_helpers::distance_dyn_boid;
use crate::options::InitiationStrategy;
use crate::options::RunOptions;

// a tracker takes ownership of a flock of boids
// a tracker will follow an implementation consisting of choosing a datastructure to
// hold spatial information about boids in and have traits for querying 2D, later 3D
// information

// todo: maybe there is more performance if we stop can discriminate entities that are already in the correct cell
/// Implementation follows [paper](https://www.researchgate.net/publication/277870601_A_Hash_Table_Construction_Algorithm_for_Spatial_Hashing_Based_on_Linear_Memory)
pub struct SpatHash1D {
    pivots: Vec<SpatHashPiv>,
    table: Vec<Boid>,
    index: Vec<usize>,
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

// get the cell vectors pointing towards the neighbourhood cell
// in eucledian space
const HOME: [i64; 2] = [0, 0];
const S: [i64; 2] = [0, -1];
const SW: [i64; 2] = [-1, -1];
const SE: [i64; 2] = [1, -1];
const W: [i64; 2] = [-1, 0];
const E: [i64; 2] = [1, 0];
const N: [i64; 2] = [0, 1];
const NE: [i64; 2] = [1, 1];
const NW: [i64; 2] = [-1, 1];

impl SpatHash1D {
    fn new(entities: &[Boid], run_options: &RunOptions) -> Self {
        let settings = SpatHash1D::get_tracker_settings(run_options);

        SpatHash1D {
            // pivots: Vec::with_capacity(settings.cell_count),
            // initialize vector with both capacity and values prefilled to simplify code in the update_table
            pivots: vec![Default::default(); settings.cell_count],
            table: entities.to_vec(),
            // index: Vec::with_capacity(run_options.init_boids as usize),
            // initialize vector with both capacity and values prefilled to simplify code in the update_table
            index: vec![0; run_options.init_boids as usize],
            settings,
            // accellerations: Vec::with_capacity(run_options.init_boids as usize),
        }
    }

    pub fn insert(&mut self, entity: Boid, run_options: &RunOptions) {
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
        self.table.push(entity);
        // make sure the table updates its structure
        self.update_table(run_options);
    }

    pub fn view(&self) -> &Vec<Boid> {
        &self.table
    }

    pub fn update(&mut self, run_options: &RunOptions) {
        let new_settings = SpatHash1D::get_tracker_settings(run_options);
        let cell_difference = self.settings.cell_count as i64 - new_settings.cell_count as i64;
        if cell_difference != 0 {
            self.pivots.resize(new_settings.cell_count, Default::default());
            self.settings = new_settings;
        }
        // if self.settings.x_cell_res != new_settings.x_cell_res || self.settings.y_cell_res != new_settings.y_cell_res {

        //     // self.settings = new_settings;
        //     // pivots: vec![Default::default(); settings.cell_count],
        //     // self.update_table(run_options);
        // } else {
        //     self.settings = new_settings;

        // }

        self.update_loation(run_options);
        self.update_table(run_options);
    }

    fn update_loation(&mut self, run_options: &RunOptions) {
        // todo: temp print
        // println!("{:#?}", self.index);
        // println!("--------");
        // println!("tracker entity count: {}", self.table.len());

        let mut neighbours: Vec<&Boid> = Vec::with_capacity(32);
        let mut accellerations: Vec<Vec2> = Vec::with_capacity(run_options.init_boids as usize);

        // get neighbours and final accelleration per boid
        for e in 0..self.table.len() {
            neighbours.clear();

            self.get_neighbours(&self.table[e], self.index[e], run_options, &mut neighbours);
            // todo: temp print
            // println!("tracker neighbours: {}", neighbours.len());

            // this blows up, because the vector is allocated with appropriate capacity, but is empty
            // accellerations[e] = self.table[e].run_rules3(&neighbours, run_options);
            // I will attempt to patch it with push to get rid of the overhead in pre-creating the Vec2s
            // and rely on implicit order, which might cause issues, but should not as the for loop goes
            //sequentually through the entities
            accellerations.push(self.table[e].run_rules3(&neighbours, run_options));
        }

        // apply the forces
        for e in 0..self.table.len() {
            self.table[e].apply_force(accellerations[e]);
            self.table[e].update_location(&run_options)
        }

        // todo: temp print
        // println!("tracker first pos: {}", self.table[0].position);
        // println!("tracker first acc: {}", accellerations.first().unwrap_or(&Vec2::new(0., 0.)));
        // println!("tracker accs length: {}", accellerations.len());
        // println!("########");
        accellerations.clear();
        // self.update_table(&self.table, run_options)
    }

    pub fn update_table(&mut self, run_options: &RunOptions) {
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
            // self.table[e].id ==
            self.index[e] = index;
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

        for e in 0..self.table.len() {
            // here is where the algorithm differs slightly from the original
            // as we are reusing the same table as the set of entities, we
            // have to make a swap this collapses into a form of unstable
            // insertion sort...

            // self.table[self.pivots[self.index[e]].fin.unwrap() - 1] = self.table[e];
            // mem::swap(&mut self.table[self.pivots[self.index[e]].fin.unwrap() - 1], &mut self.table[e]);
            self.table
                .swap(self.pivots[self.index[e]].fin.unwrap() - 1, e);

            // update the pivot now that one entity of the cell has been mapped to the correct place in the 1D table
            self.pivots[self.index[e]].fin = Some(self.pivots[self.index[e]].fin.unwrap() - 1)
        }
    }

    pub fn get_neighbours<'a>(
        &'a self,
        boid: &Boid,
        cell_index: usize,
        run_options: &RunOptions,
        neighbours: &mut Vec<&'a Boid>,
    ) {
        // only add it once
        let is_left = cell_index % self.settings.x_cell_count as usize == 0;
        let is_right = (cell_index % self.settings.x_cell_count as usize) + 1
            == self.settings.x_cell_count as usize;
        let is_bottom = cell_index / self.settings.x_cell_count as usize == 0;
        let is_top =
            cell_index as f32 >= (self.settings.y_cell_count - 1f32) * self.settings.x_cell_count;
        // todo: temp print
        // {
        //     println!("no cells: {}", self.settings.cell_count);
        //     println!("no x cells: {}", self.settings.x_cell_count);
        //     println!("no y cells: {}", self.settings.y_cell_count);
        //     println!("cell: {}", cell_index);
        //     println!("left: {}", is_left);
        //     println!("right: {}", is_right);
        //     println!("bottom: {}", is_bottom);
        //     println!("top: {}", is_top);
        // }
        // let mut acum = 0u32;
        for cell in match (is_left, is_right, is_bottom, is_top) {
            (false, false, false, false) => [NW, N, NE, W, HOME, E, SW, S, SE],
            (true, false, true, false) => [HOME, N, NE, HOME, HOME, E, HOME, HOME, HOME],
            (true, false, false, true) => [HOME, HOME, HOME, HOME, HOME, E, HOME, S, S],
            (true, false, false, false) => [HOME, N, NE, HOME, HOME, E, HOME, S, SE],
            (false, true, true, false) => [NW, N, HOME, W, HOME, HOME, HOME, HOME, HOME],
            (false, true, false, true) => [HOME, HOME, HOME, W, HOME, HOME, SW, S, HOME],
            (false, true, false, false) => [NW, N, HOME, W, HOME, HOME, SW, S, HOME],
            // (false, false, true, true) => todo!(),
            (false, false, true, false) => [NW, N, NE, W, HOME, E, HOME, HOME, HOME],
            (false, false, false, true) => [HOME, HOME, HOME, E, HOME, E, SW, S, SE],
            _ => panic!("it should not have come to this :"),
        }
        // depending on the match, iterate the vectors
        .iter()
        // filter out the out-of-bound vectors, the home vector is used as a placeholder
        .filter(|l| **l != HOME)
        // always consider the boid's current 'home' vector
        .chain(iter::once(&HOME))
        // convert vectors to the table's 1D indexes pointing to cells
        .map(|lookup| {
            // println!("{}, {}", lookup[0], lookup[1]);
            (cell_index as i64 + lookup[0] + lookup[1] * self.settings.x_cell_count as i64) as usize
        }) {
            // println!("range_cell: {}", cell);
            // acum += 1;

            if self.pivots[cell].usg == 0 {
                continue;
            }
            for index in self.pivots[cell].init.unwrap()
                ..(self.pivots[cell].init.unwrap() + self.pivots[cell].usg)
            {
                if self.table[index].id != boid.id
                    && distance_dyn_boid(boid, &self.table[index], run_options)
                        <= run_options.max_sensory_distance
                {
                    neighbours.push(&self.table[index]);
                    if neighbours.len() == run_options.neighbours_cosidered {
                        break;
                    }
                }
            }
        }
        // todo: temp print
        // for _ in 0..(9 - acum){
        //     println!(" - ");
        // }
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

    fn get_tracker_settings(run_options: &RunOptions) -> SpatialHashingTableSettings {
        let x_range = (run_options.window.win_right - run_options.window.win_left).ceil();
        let y_range = (run_options.window.win_top - run_options.window.win_bottom).ceil();

        let x_cell_count = (x_range as f32 / (run_options.max_sensory_distance * 2.)).ceil();
        let y_cell_count = x_cell_count;
        // let x_cell_count = 128.;
        // let y_cell_count = (y_range as f32 / (run_options.max_sensory_distance)).ceil();
        // let y_cell_count = 128.;

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
}

// if let Some(fin) = self.pivots[self.index[e]].fin {
//     self.pivots[self.index[e]].fin = Some(fin - 1);
// }
/// A naive implementation of boids tracking, which uses an O(N^2) algorithm for
/// finding boid's neighbours.
struct BoidTracker {
    boids: Vec<Boid>,
}
/// Uses a spatial hashing space division method, where all cells of the underlying
/// table are stored in a 1D array, with the individual cell's being allocated
/// dynamically, reused or updated as needed.
///
/// Uses an Θ(n) algorithm for finding the boid's neighbours.
/// ref: Optimization_of_large-scale_real-time_simulations
struct SpatialHashingTableSettings {
    pub x_range: f32,
    pub y_range: f32,
    pub x_cell_count: f32,
    pub y_cell_count: f32,
    pub cell_count: usize,
    pub x_cell_res: f32,
    pub y_cell_res: f32,
    // run_options: &'a RunOptions
}

pub struct Flock {
    boids: Vec<Boid>,
    // pub boids2: Vec<Boid>,
    // tracker: SpatialHashingTracker15D,
    pub tracker: SpatHash1D,
}

impl Flock {
    // const GREY_SCALE: &str = "$@B%8&WM#*oahkbdpqwmZO0QLCJUYXzcvunxrjft/\\|()1{}[]?-_+~<>i!lI;:,\"^`'. ";
    const GREY_SCALE: &str =
        " .\'`^\",:;Il!i><~+_-?][}{1)(|\\/tfjrxnuvczXYUJCLQ0OZmwqpdbkhao*#MW&8%B@$";
    pub fn new(run_options: &RunOptions) -> Self {
        let boids = get_boids(&run_options);

        Flock {
            // boidss: Cell::new(boids.iter().map(|b| b.to_owned()).collect()),
            tracker: SpatHash1D::new(&boids, run_options),
            // boids
            boids: Vec::new(),
        }
    }

    pub fn view(&self) -> &Vec<Boid> {
        &self.boids
    }

    pub fn update(&mut self, run_options: &RunOptions) {
        self.tracker.update(run_options);
        // test_me();
        // panic!();

        // let table = SpatialHashingTracker15D::get_spatial_subdivision_table(self.boids.to_owned(), run_options);
        // println!("table size: {}", table.len());

        if run_options.dbscan_flock_clustering_on {
            let test_data: Array2<f32> = self
                .boids
                .iter()
                .map(|row| [row.position.x, row.position.y])
                .collect::<Vec<_>>()
                .into();
            let test_dataset: DatasetBase<_, _> = test_data.into();
            let res = Dbscan::params(3)
                .tolerance(run_options.sensory_distance)
                .transform(test_dataset)
                .unwrap();

            self.boids = res
                .targets
                .iter()
                .zip(self.boids.iter_mut())
                .map(|pair| match pair.0 {
                    Some(cluster) => {
                        pair.1.cluster_id = *cluster + 1;
                        pair.1.to_owned()
                    }
                    None => {
                        pair.1.cluster_id = 0;
                        pair.1.to_owned()
                    }
                })
                .collect::<Vec<Boid>>();

            // let label_count = res.label_count().remove(0);

            // println!();
            // println!("Result: ");
            // for (label, count) in label_count {
            //     match label {
            //         None => println!(" - {} noise points", count),
            //         Some(i) => println!(" - {} points in cluster {}", count, i),
            //     }
            // }
        };

        // for i_cur in 0..self.boidss.get_mut().len() {
        //     for i_other in 0..self.boidss.get_mut().len() {
        //         if i_cur == i_other {continue}

        //         self.boidss[i_cur]
        //     }
        // }
        // for b_cur in self.boidss.iter(){
        //     let mut neighbours: Vec<&Boid> = Vec::new();

        //     for b_other in self.boidss.iter(){
        //         if b_other.get().id == b_cur.get().id {
        //             continue;
        //         }

        //         let distance = distance_dyn_boid(b_cur.get(), b_other, &run_options);
        //         if distance < max_sensory_distance {
        //             neighbours.push(b_other);
        //         }
        //     }
        // }
        // let clone = self.boids.to_owned();

        // for b_cur in self.boids.iter_mut() {
        //     let mut neighbours: Vec<&Boid> = Vec::new();

        //     for b_other in clone.iter() {
        //         if b_other.id == b_cur.id {
        //             continue;
        //         }

        //         let distance = distance_dyn_boid(b_cur, b_other, &run_options);
        //         if distance < max_sensory_distance {
        //             neighbours.push(b_other);
        //         }
        //     }

        //     b_cur.run_rules(neighbours, &run_options);
        //     b_cur.update_location(&run_options)
        // }

        // let table = HashMap::new()
        // use flat_spatial::Grid;

        // let mut g: Grid<(), [f32; 2]> = Grid::new(10);
        // let a = g.insert([3.0, 3.0], ());
        // let _b = g.insert([12.0, -8.0], ());

        // let around: Vec<_> = g
        //     .query_around([2.0, 2.0], 5.0)
        //     .map(|(id, _pos)| id)
        //     .collect();

        let mut accelleration: Vec<Vec2> = Vec::with_capacity(self.boids.len());

        // calculation loop
        for i_cur in 0..self.boids.len() {
            let mut neighbours: Vec<&Boid> = Vec::new();

            for b_other in self.boids.iter() {
                if b_other.id == self.boids[i_cur].id {
                    continue;
                }

                let distance = distance_dyn_boid(&self.boids[i_cur], b_other, &run_options);
                if distance < run_options.max_sensory_distance {
                    neighbours.push(b_other);
                }
            }
            // todo: temp print
            // if self.boids[i_cur].id == 0 {
            //     println!("flock neighbours: {}", neighbours.len());
            // }
            accelleration.push(self.boids[i_cur].run_rules2(neighbours, &run_options));

            // todo: temp print
            // println!("flock first pos: {}", self.boids[0].position);
        }
        // todo: temp print
        // println!("flock first acc: {}", accelleration.first().unwrap_or(&Vec2::new(0., 0.)));
        // println!("flock accs length: {}", accelleration.len());
        // update loop
        for i_cur in 0..self.boids.len() {
            self.boids[i_cur].apply_force(accelleration[i_cur]);
            self.boids[i_cur].update_location(&run_options)
        }
    }

    pub fn insert(&mut self, run_options: &RunOptions) {
        let b = get_boid(run_options, self.boids.len() as u32);
        // self.insert_boid(b.to_owned());
        self.tracker.insert(b, run_options)
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
    pub fn delete_boid(&mut self, id_delete: u32) -> u32 {
        self.boids = self
            .boids
            .iter_mut()
            .enumerate()
            .filter(|(_, b)| b.id != id_delete)
            .map(|(index, b)| {
                b.id = index as u32;
                *b
            })
            .collect();

        id_delete
    }

    pub fn restart(&mut self, run_options: &RunOptions) -> () {
        self.boids = get_boids(&run_options);
    }
}

fn get_boids(run_options: &RunOptions) -> Vec<Boid> {
    (0..run_options.init_boids)
        .map(|id| get_boid(run_options, id))
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
// todo create a datastructure to build u preferences of the local boids for searching neighbours
// begin by getting rid of the copying above and put that into the datastructure
// continue on by usinr the same interface for the spatial hashingi

// the structure should be able to take a copy of the data to build its internal guts
// then a function for querying n-nearest neighbours should be made publicly available

// struct NaiveNeighbourFinder {
//     pub boids: Vec<Boid>
// }

// impl NaiveNeighbourFinder {
//     pub fn new(run_options: &RunOptions) -> Self {
//         let boids = get_boids(&run_options);

//         Flock {
//             boids,
//         }
//     }
// }

#[cfg(test)]
mod tests {

    use crate::flock::{Flock, SpatHash1D};

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
}
