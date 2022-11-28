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
struct SpatialHashingTracker15D {
    settings: SpatialHashingTableSettings,
    table: Vec<Vec<Boid>>,
    // pub accellerations: Vec<Vec2>

}

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
impl SpatialHashingTracker15D {
    // The inspiration for mapping 
    pub fn new(run_options: &RunOptions) -> Self {
        // let SpatialHashingTableSettings {
        //     x_range,
        //     y_range,
        //     x_cell_count,
        //     y_cell_count,
        //     cell_count,
        //     x_cell_res,
        //     y_cell_res,
        // } = SpatialHashingTracker15D::get_tracker_settings(run_options);

        let settings = SpatialHashingTracker15D::get_tracker_settings(run_options);
        
        SpatialHashingTracker15D {
            table: vec![Vec::with_capacity(30); settings.cell_count],
            // accellerations: Vec::with_capacity(run_options.init_boids as usize),
            settings,

        }
    }

    pub fn update(&mut self, run_options: &RunOptions){

        // todo tune the capacity so it matches more with the actual value
        let mut neighbours: Vec<&Boid> = Vec::with_capacity(20);

        // self.accellerations.clear();

        let mut accellerations :Vec<Vec2> = Vec::with_capacity(run_options.init_boids as usize);
        
        let mut boid_counter: usize = 0;

        // self.table.iter().enumerate()
        // .for_each(|(i, b)| {
        //     b.iter()
        //     .for_each(|b| {
        //         neighbours.clear();
        //         self.get_neighbours(&b.to_owned(), i, run_options, &mut neighbours);
        //         accellerations[boid_counter] = b.run_rules3(&neighbours, run_options);
        //         boid_counter += 1; 

        //     })
        // });

        
        for cell_i in 0..self.table.len(){
            for boid_i in 0..self.table[cell_i].len(){
                let boid = self.table[cell_i][boid_i];
                neighbours.clear();
                self.get_neighbours(&boid, cell_i, run_options, &mut neighbours);
                // neighbours.append(other)
                accellerations[boid_counter] = boid.run_rules3(&neighbours, run_options);
                boid_counter += 1; 
            }
        }

        boid_counter = 0;
        self.table.iter_mut()
        .enumerate()
        .for_each(|(cell_i, cell)|{
            // apply accelleration
            cell.iter_mut()
            .for_each(|b|{
                b.apply_force(accellerations[boid_counter]);

            });

            // // recalculate cell the boid belongs to
            // cell.iter_mut()
            // .for_each(|b|{
            //     let target_cell = Flock::get_table_index(
            //         b.position.x,
            //         b.position.y,
            //         run_options.window.win_left,
            //         run_options.window.win_right - 1.,
            //         run_options.window.win_bottom,
            //         run_options.window.win_top - 1.,
            //         self.settings.x_cell_res,
            //         self.settings.y_cell_res,
            //         self.settings.x_cell_count,
            //     );

            //     if cell_i != target_cell {
                    
            //     }
            // })
        });


        for cell_i in 0..self.table.len(){
            let mut boid_i = 0;
            while boid_i < self.table[cell_i].len() {
                let boid = self.table[cell_i][boid_i];
                // boid.apply_force(accellerations[boid_counter]);
                // boid_counter += 1; 
                let target_cell = SpatialHashingTracker15D::get_table_index(
                    boid.position.x,
                    boid.position.y,
                    run_options.window.win_left,
                    run_options.window.win_right - 1.,
                    run_options.window.win_bottom,
                    run_options.window.win_top - 1.,
                    self.settings.x_cell_res,
                    self.settings.y_cell_res,
                    self.settings.x_cell_count,
                );

                if cell_i != target_cell {
                    _ = self.table[cell_i].swap_remove(cell_i);
                    self.table[target_cell].push(boid)
                } else {
                    boid_i += 1;
                }
            }
        }

        // self.table
        // .iter()
        // .enumerate()
        // .map(|(i, b)|{
        //     b.iter()
        //     .map(|b|{
        //         // neighbours.clear();
        //         self.get_neighbours(&b.to_owned(), i, run_options, &mut neighbours);
 
        //     })
        // });

        // self.table
        // .iter_mut()
        // .enumerate()
        // // enumerate to apply forces
        // .for_each(|(cell_index, cell_boids)| {
        //     cell_boids
        //     .iter_mut()
        //     .for_each(|boid| {
        //         neighbours.clear();
        //         self.get_neighbours(boid, cell_index, run_options, &mut neighbours);

        //         boid.run_rules(neighbours, &run_options);
        //     })
        // })
        // iterate through all cells, noting the current cell number
        // go through all the boids
        // collect the accelleration
        // 
        // then iterate again, apply acceleration
        // then update the datastructure
    }

    pub fn get_neighbours<'a>(&'a self, boid: &Boid, cell_index: usize, run_options: &RunOptions, neighbours:&mut Vec<&'a Boid>) {

        // // todo: get rid of this, as whoever will be iterating the table will know the current boid's cell
        // let i = Flock::get_table_index(
        //     boid.position.x,
        //     boid.position.y,
        //     run_options.window.win_left,
        //     run_options.window.win_right - 1.,
        //     run_options.window.win_bottom,
        //     run_options.window.win_top - 1.,
        //     self.settings.x_cell_res,
        //     self.settings.y_cell_res,
        //     self.settings.x_cell_count,
        // );
        
        // todo: do something about the ammount of casting here
        // todo: this might cause rounding error issues
        let is_left = cell_index % self.settings.x_cell_count as usize == 0;
        let is_right = (cell_index % 5) + 1 == self.settings.x_cell_count as usize;
        let is_bottom = cell_index / self.settings.x_cell_count as usize == 0;
        let is_top = cell_index as f32 >= (self.settings.y_cell_count - 1f32) * self.settings.x_cell_count;

        const HOME: (i32, i32) = (0, 0);
        // #[rustfmt::skip]
        // mod unformatted{
        //     const N:  (i32, i32) = (0, 1);
        //     const NW: (i32, i32) = (-1, 1);
        //     const W:  (i32, i32) = (-1, 0);
        //     const WS: (i32, i32) = (-1, -1);
        //     const S:  (i32, i32) = (0, -1);
        //     const SE: (i32, i32) = (1, -1);
        //     const E:  (i32, i32) = (1, 0);
        // }

        // use unformatted::*;

        {
            // Depending on the position of a boid
            // (left, right, bottom top) in cells rows/columns,
            //  take a subset of cells for consideration when 
            // determining its neighbours
            // E.g., when boid is on the left most edge in the middle,
            // it will ignore anything left of it as it woul have been
            // out of bounds and only check the other 5 left, right,
            // bottom, top cells rows/columns individually
            // todo: this could be made a bit more efficient if it considered the boid's angle and FOV
            #[rustfmt::skip]
            match (is_left, is_right, is_bottom, is_top) {
                // (true, true, true, true) => todo!(),
                // (true, true, true, false) => todo!(),
                // (true, true, false, true) => todo!(),
                // (true, true, false, false) => todo!(),
                // (true, false, true, true) => todo!(),

                // &[
                //     (-1, 1), (0, 1), (1, 1), 
                //     (-1, 0), (0, 0), (1, 0), 
                //     (-1, -1), (0, -1), (1, -1),
                // ],
                (true, false, true, false) => &[
                    (0 as i32, 0 as i32), (0, 1), (1, 1), 
                    (0, 0), (0, 0), (1, 0), 
                    (0, 0), (0, 0), (0, 0),
                ],
                (true, false, false, true) => &[
                    (0, 0), (0, 0), (0, 0), 
                    (0, 0), (0, 0), (1, 0), 
                    (0, 0), (0, -1), (0, -1),
                ],
                (true, false, false, false) => &[
                    (0, 0), (0, 1), (1, 1), 
                    (0, 0), (0, 0), (1, 0), 
                    (0, 0), (0, -1), (1, -1),
                ],
                (false, true, true, false) => &[
                    (-1, 1), (0, 1), (0, 0), 
                    (-1, 0), (0, 0), (0, 0), 
                    (0, 0), (0, 0), (0, 0),
                ],
                (false, true, false, true) => &[
                    (0, 0), (0, 0), (0, 0), 
                    (-1, 0), (0, 0), (0, 0), 
                    (-1, -1), (0, -1), (0, 0),
                ],
                (false, true, false, false) => &[
                    (-1, 1), (1, 1), (0, 0), 
                    (-1, 0), (0, 0), (0, 0), 
                    (-1, -1), (0, -1), (0, 0),
                ],
                // (false, false, true, true) => todo!(),
                (false, false, true, false) => &[
                    (-1, 1), (0, 1), (1, 1), 
                    (-1, 0), (0, 0), (1, 0), 
                    (0, 0), (0, 0), (0, 0),
                ],
                (false, false, false, true) => &[
                    (0, 0), (0, 0), (0, 0), 
                    (1, 0), (0, 0), (1, 0), 
                    (-1, -1), (0, -1), (1, 0-1),
                ],
                (false, false, false, false) => &[
                    (-1, 1), (0, 1), (1, 1), 
                    (-1, 0), (0, 0), (1, 0), 
                    (-1, -1), (0, -1), (1, -1),
                ],
                _ => panic!("it should not have come to this :")

            }
            // depending on the match, iterate neighbouring cells
            .iter()
            // filter out the out-of-bound cells, the home cells is used as a placeholder
            .filter(|l| **l != (0, 0)) 
            // always consider the boid's current 'home' cell
            .chain(iter::once(&HOME)) // only add it once
            // convert to the table's 1D index
            .map(|lookup| {
                cell_index as i32 + lookup.0 + lookup.1 * self.settings.x_cell_count as i32
            })
            // get the cell vectors as iterators
            .map(|index| -> std::slice::Iter<Boid> {
                self.table[index as usize].iter()
            })
            // flatten into a single iterator
            .flatten()
            // filter boids within max sensory distance
            .filter(|b| distance_dyn_boid(boid, b, run_options) <= run_options.max_sensory_distance)
            // take out the neighbours and push them to the vector of neighbours
            .for_each(|b| neighbours.push(b));
        };
 
    }

    fn get_spatial_subdivision_table(boids: Vec<Boid>, run_options: &RunOptions) -> Vec<Vec<Boid>> {
        // ref: Optimization_of_large-scale_real-time_simulations
        // float min = 0
        // float max = 100
        // int cell size = 25
        // int width = (max-min)/cell size = 4
        // int number of buckets = width
        // int grid_cell = ( floor(x/cell size) ) + ( floor(y/cell size) )*width
        
        // for the time being I am just using these hardcoded values for the division to make debugging easier
        
        let x_range = (run_options.window.win_right - run_options.window.win_left).ceil();
        let y_range = (run_options.window.win_top - run_options.window.win_bottom).ceil();
        
        // let width =  (x_range / x_cell_res) as usize;
        let x_cell_count = (x_range as f32 / (run_options.max_sensory_distance / 2.)).ceil();
        // let height = (y_range / y_cell_res) as usize;
        let y_cell_count = (y_range as f32 / (run_options.max_sensory_distance / 2.)).ceil();
        
        let cell_count = (x_cell_count * y_cell_count) as usize;
        
        let mut table: Vec<Vec<Boid>> = vec![Vec::with_capacity(30); cell_count];
        
        let x_cell_res = x_range / x_cell_count;
        let y_cell_res = y_range / y_cell_count;
        
        // conversion_factor is 1/size
        // int grid_cell =  x*conversion_factor +  y*conversion_factor*width
        
        // let conversion_factor_x = 1. / x_cell_res;
        // let conversion_factor_y = 1. / y_cell_res;
        for b in boids.iter() {
            let table_cell_i = SpatialHashingTracker15D::get_table_index(
                b.position.x,
                b.position.y,
                run_options.window.win_left,
                run_options.window.win_right - 1.,
                run_options.window.win_bottom,
                run_options.window.win_top - 1.,
                x_cell_res,
                y_cell_res,
                x_cell_count,
            );
            table[table_cell_i].push(*b)
        }
        
        // fun print
        // for y in (0..(y_cell_count as usize)).rev() {
        //     // 1..=16
        //     for x in 0..(x_cell_count as usize) {
        //         // 0..16
        //         let len = table[x + y * x_cell_count as usize].len();
                
        //         // print!("{:>3}", len);
        //         // println!("len: {:>3}", Flock::GREY_SCALE.len());
        //         print!(
        //             "{}",
        //             Flock::GREY_SCALE
        //             .chars()
        //             .nth(len % 70)
        //             .expect("out of character range!")
        //         );
        //     }
        //     println!();
        // }
        table
    }
    
    /// Returns "hashed" value representing an index for spatial subdivision, handles a zero centered coordinate system
    ///
    /// ## Arguments
    ///
    /// `p_x` - x coordinate of a point of interest
    /// `p_y` - y coordinate of a point of interest
    /// `min_x` - minimum x axis value
    /// `min_y` - minimum y axis value
    /// `cs_x` - cell size along the x axis
    /// `cs_y` - cell size along the y axis
    /// `table_width` - the table width, how many cells make up a row
    /// 
    /// ## Returns
    /// an index as _usize_, that takes values from 0..(x_cells * y_cells)
    ///
    /// watch out, in case of a symmetric/centered coordinate system (like nannou)
    /// the recommendation is to subtract 1 from the x/y max values to reduce the range of returned indexes
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

    pub fn get_tracker_settings(run_options: &RunOptions) -> SpatialHashingTableSettings {
        let x_range = (run_options.window.win_right - run_options.window.win_left).ceil();
        let y_range = (run_options.window.win_top - run_options.window.win_bottom).ceil();

        let x_cell_count = (x_range as f32 / (run_options.max_sensory_distance / 2.)).ceil();
        let y_cell_count = (y_range as f32 / (run_options.max_sensory_distance / 2.)).ceil();

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


pub struct Flock {
    pub boids: Vec<Boid>,
    tracker: SpatialHashingTracker15D,
}

impl Flock {
    // const GREY_SCALE: &str = "$@B%8&WM#*oahkbdpqwmZO0QLCJUYXzcvunxrjft/\\|()1{}[]?-_+~<>i!lI;:,\"^`'. ";
    const GREY_SCALE: &str =
    " .\'`^\",:;Il!i><~+_-?][}{1)(|\\/tfjrxnuvczXYUJCLQ0OZmwqpdbkhao*#MW&8%B@$";
    pub fn new(run_options: &RunOptions) -> Self {
        let boids = get_boids(&run_options);
        
        Flock {
            // boidss: Cell::new(boids.iter().map(|b| b.to_owned()).collect()),
            boids,
            tracker: SpatialHashingTracker15D::new(run_options)
        }
    }
    
    
    pub fn update(&mut self, run_options: &RunOptions) {
        // test_me();
        // panic!();
        let table = SpatialHashingTracker15D::get_spatial_subdivision_table(self.boids.to_owned(), run_options);
        
        println!("table size: {}", table.len());
        
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
            
            accelleration.push(self.boids[i_cur].run_rules2(neighbours, &run_options));
        }
        
        // update loop
        for i_cur in 0..self.boids.len() {
            self.boids[i_cur].apply_force(accelleration[i_cur]);
            self.boids[i_cur].update_location(&run_options)
        }
    }
    
    pub fn insert(&mut self, run_options: &RunOptions) {
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

    use crate::flock::{Flock, SpatialHashingTracker15D};

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

        let index = SpatialHashingTracker15D::get_table_index(
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

        let index = SpatialHashingTracker15D::get_table_index(
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

        let index = SpatialHashingTracker15D::get_table_index(
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
