use std::{any::Any, path::Path};

use super::tracker::{Tracker, TrackerSignal};
use crate::{
    birdwatcher::BoidData,
    boid::{self, Boid, BoidMetadata, BoidType},
    math_helpers::tor_vec_p,
    options::{get_window_size, RunOptions},
};
use csv::Reader;
use glam::{vec2, Mat2};
use itertools::Itertools;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum ReplayTrackerError {
    // #[error("data store disconnected")]
    // Disconnect(#[from] io::Error),
    // #[error("the data for key `{0}` is not available")]
    // Redaction(String),
    // #[error("invalid header (expected {expected:?}, found {found:?})")]
    // InvalidHeader {
    //     expected: String,
    //     found: String,
    // },
    #[error("replay tracker could not load the next buffer, msg: `{0}`")]
    BufferLoadError(String),
    #[error("data has non-uniform or uneven structures, number of boids is expected to stay constant for all timepoints, msg: `{0}`")]
    DataBlockError(String),
}

#[derive(Clone, Debug, PartialEq, Copy)]
pub enum ReplayDirection {
    Forward,
    Backward,
}

pub struct ReplayTracker {
    pub locations: Vec<BoidData>,
    pub buff_boid: Vec<Boid>,
    pub buff_meta: Vec<BoidMetadata>,
    pub buffer_size: usize,
    pub no_boids: usize,
    // ghosts make up trails behind individuals -> helps signify movement at stop-motion
    pub no_ghosts: i32,
    pub run_options: RunOptions,
    // points at the current position before the about-to-be-played data block
    buff_playhead: usize,
    // makes sure we do not hit the ends so that we can calculate headings via differecnces of positions
    // also used to make sure we have enough ghosts to display
    // todo: there is a bug where if playhead reaches bumber end, then backwards tick is issued, it plays a frame ahead and then goes backward,
    // the inverse is true for going backward then ticking forward
    // todo: bumper logic assumed reverse direction will literally reverse the velocities, which we decided is not how we want it to work
    // so bumber should not be stopping playhead at +bumper on the left and len-bumper on the right (providing bumped space on either side)
    // instead it should only do so on the left, and leave the right to play to the end as ghosts will come from the "future" when going backward,
    // instead of from the past when going forward
    buff_bumper: usize,
    replay_direction: ReplayDirection,
    // rescale: f32,
    // for rescaling, in case we have too small or too large of a window for the current replay
    rescale_mat: Mat2,
    ghosts_on: bool,
    paused: bool,
}

impl ReplayTracker {
    pub fn from_path(
        replay_path: &String,
        run_options: &RunOptions,
    ) -> Result<Self, ReplayTrackerError> {
        let path = Path::new(replay_path);

        if !path.exists() {
            panic!("can't find replay file")
        }

        // load all the data into memory
        let mut rdr = Reader::from_path(path).unwrap();
        let locations = rdr
            .deserialize::<BoidData>()
            .map(|i| i.expect("boid data record"))
            .sorted_by(|a, b| {
                if a.time != b.time {
                    a.time.cmp(&b.time)
                } else {
                    a.id.cmp(&b.id)
                }
            })
            .collect::<Vec<BoidData>>();

        let mut ro = run_options.to_owned();

        // find the whole max x and y values of the data at hand
        let (mut max_x, mut max_y) = locations
            .iter()
            .fold::<(f32, f32), _>((0., 0.), |acc: (f32, f32), cur| {
                (cur.x.max(acc.0), cur.y.max(acc.1))
            });
        max_x = max_x.ceil();
        max_y = max_y.ceil();

        // let max_x = 450.;
        // let max_y = 450.;

        // get ratios
        // these will be between (0, inf]
        let rdx = run_options.window.win_right as f32 / max_x;
        let rdy = run_options.window.win_top as f32 / max_y;

        // in case (0, 1) we are scaling down, we have to use the smaller of the two to make sure
        // either side fits

        // in case (1, inf] we are scaling up, again we use the smaller one to make sure we do not
        // over stretch the larger dimension

        // = 1 we do not care
        let rescale = rdx.min(rdy);
        // ˜λI 2x2 rescaling matrix
        let rescale_mat = rescale * Mat2::IDENTITY;

        ro.window = get_window_size((max_x * 2.) as u32, (max_y * 2.) as u32);
        // scan the items until one repetition is found, assume the number of
        // items skipped is the number of boids per iteration
        let no_boids = locations
            .iter()
            .enumerate()
            .find_or_first(|(index, d)| index != &0 && d.id == locations[0].id)
            .expect("")
            .0;

        // overwrite defaults
        ro.init_boids = no_boids;

        let no_ghosts: i32 = 4;
        let ghosts_on = true;

        let (buff_bumper, buffer_size) = match ghosts_on {
            true => (
                2 + no_ghosts.abs() as usize,
                no_boids + no_boids * no_ghosts.abs() as usize,
            ),
            false => (1 as usize, no_boids),
        };

        // let buff_bumper: usize = 1;

        let mut ret = ReplayTracker {
            locations,
            buff_boid: vec![Boid::default(); buffer_size],
            buff_meta: vec![BoidMetadata::default(); buffer_size],
            buffer_size,
            buff_playhead: ReplayTracker::get_start_playhead(buff_bumper, ro.init_boids),
            run_options: ro,
            buff_bumper,
            replay_direction: ReplayDirection::Forward,
            rescale_mat,
            no_ghosts,
            ghosts_on,
            paused: false,
            no_boids,
        };

        let res = ret.load_buffer();

        if let Err(err) = res {
            return Err(err);
        }

        if ret.buff_boid.len() % ret.run_options.init_boids != 0 {
            return Err(ReplayTrackerError::DataBlockError(
                "failed to initialise".to_owned(),
            ));
        }

        return Ok(ret);
    }

    fn get_start_playhead(buff_bumper: usize, init_boids: usize) -> usize {
        buff_bumper * init_boids
    }

    fn is_right(&self) -> bool {
        self.buff_playhead
            == self.locations.len() - self.run_options.init_boids * (self.buff_bumper)
    }

    fn is_left(&self) -> bool {
        self.buff_playhead == self.run_options.init_boids * (self.buff_bumper)
    }

    fn playhead_at_end(&self) -> bool {
        match self.replay_direction {
            ReplayDirection::Forward => self.is_right(),
            ReplayDirection::Backward => self.is_left(),
        }
    }

    fn playhead_at_start(&self) -> bool {
        match self.replay_direction {
            ReplayDirection::Forward => self.is_left(),
            ReplayDirection::Backward => self.is_right(),
        }
    }
    fn load_buffer(&mut self) -> Result<bool, ReplayTrackerError> {

        let RunOptions {
            init_boids, window, ..
        } = self.run_options;

        let finished = self.playhead_at_end();
        let starting = self.playhead_at_start();

        // does the next step have enough data points?
        if init_boids
            > match self.replay_direction {
                ReplayDirection::Backward => self.buff_playhead + 1,
                ReplayDirection::Forward => self.locations.len() - self.buff_playhead,
            }
        {
            return Err(ReplayTrackerError::DataBlockError(
                "next step does not have enough data points".to_owned(),
            ));
        }

        let time = self.locations[self.buff_playhead].time;

        let to_skip = match self.ghosts_on {
            // true => self.buff_playhead - self.no_ghosts as usize * init_boids,
            true => match self.replay_direction {
                // ReplayDirection::Forward => self.buff_playhead - (self.no_ghosts + 1) as usize * init_boids,
                // ReplayDirection::Backward => self.buff_playhead - (self.no_ghosts + 2) as usize * init_boids,
                ReplayDirection::Forward => {
                    (init_boids * (time as usize - 1)) - (self.no_ghosts as usize * init_boids)
                }
                ReplayDirection::Backward => {
                    (init_boids * (time as usize - 2)) - (self.no_ghosts as usize * init_boids)
                }
            },
            false => match self.replay_direction {
                ReplayDirection::Forward => init_boids * (time as usize - 1),
                ReplayDirection::Backward => init_boids * (time as usize - 2),
            },
        };

        for (data, index) in self
            .locations
            .iter()
            .skip(to_skip)
            .take(self.buffer_size)
            // .enumerate()
            .zip(
                // this generates a sequence of 0..init_boids repeating along buffer_size
                // the purpose is to check whether ids are
                //  1. consecutive, and
                //  2. whether time points contain the same number of agents
                (0..init_boids).into_iter().cycle().take(self.buffer_size),
            )
        {
            if index != data.id {
                return Err(ReplayTrackerError::DataBlockError(
                    "id missalignment".to_owned(),
                ));
            }

            let time_diff_abs = ((data.time as i128 - time as i128) as i32).abs();

            let check_time = match self.replay_direction {
                ReplayDirection::Forward => time,
                ReplayDirection::Backward => time - 1,
            };

            if self.ghosts_on
                && time_diff_abs
                    > match self.replay_direction {
                        ReplayDirection::Forward => self.no_ghosts.abs(),
                        ReplayDirection::Backward => self.no_ghosts.abs() + 1,
                    }
                || !self.ghosts_on && data.time != check_time
            {
                return Err(ReplayTrackerError::DataBlockError(
                    "time missalignment".to_owned(),
                ));
            }

            let buff_index = data.id
                + match self.ghosts_on {
                    true => {
                        init_boids
                            * match self.replay_direction {
                                ReplayDirection::Forward => time_diff_abs,
                                ReplayDirection::Backward => time_diff_abs - 1,
                            } as usize
                    }
                    false => 0,
                };


            let boid = self.buff_boid.get_mut(buff_index).unwrap();
            let meta = self.buff_meta.get_mut(buff_index).unwrap();

            boid.id = data.id;
            boid.position = self.rescale_mat * vec2(data.x, data.y);
            boid.velocity = vec2(2., 2.);

            // depending on replay direction and on whether or not we are in a starting position,
            // set the velocity by taking a lag(-1) of the difference of positions
            boid.velocity = self.rescale_mat
                * match (&self.replay_direction, starting) {
                    // (ReplayDirection::Forward, _) => {
                    (_, _) => {
                        let last_index =
                            (data.time as usize - 1) * init_boids - init_boids + data.id;
                        // println!("las_index: {:?}", last_index);
                        let last_location = self.locations[last_index];
                        // println!("last_location: {:?}", last_location);
                        let vec =
                            tor_vec_p(last_location.x, data.x, last_location.y, data.y, &window);
                        vec2(vec.0, vec.1)
                    } // (ReplayDirection::Backward, _) => {
                      //     let a = match self.ghosts_on {
                      //         // true => time_diff_abs as usize,
                      //         true => 1,
                      //         false => 0,
                      //     };
                      //     let last_index = (data.time as usize - a as usize) * init_boids + data.id;
                      //     println!("las_index: {:?}", last_index);
                      //     let last_location = self.locations[last_index];
                      //     let vec =
                      //         tor_vec_p(last_location.x, data.x, last_location.y, data.y, &window);
                      //     vec2(vec.0, vec.1)
                      // }
                };
            // println!("position: {:?}", boid.position);
            // println!("velocity: {:?}", boid.velocity);
            meta.id = data.id;
            meta.n_neighbours = data.n_neighbours;
            meta.cluster_id = data.cluster_id;
            // meta.boid_type = BoidType::Ghost;

            if self.ghosts_on
                && time_diff_abs
                    > match self.replay_direction {
                        ReplayDirection::Forward => 0,
                        ReplayDirection::Backward => 1,
                    }
            {
                meta.boid_type = BoidType::Ghost
            }
        }

        // println!("@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
        // println!("finished: {:?}", finished);
        // println!("starting: {:?}", starting);
        // println!("timeset: {:?}", timeset.iter().sorted());
        // println!("toskip: {:?}", to_skip);
        // println!("playhead before: {:?}", self.buff_playhead);

        // move buffer position
        if !finished {
            self.buff_playhead = match self.replay_direction {
                ReplayDirection::Forward => self.buff_playhead + init_boids,
                ReplayDirection::Backward => self.buff_playhead - init_boids,
            };
        }
        // println!("playhead after: {:?}", self.buff_playhead);
        // println!("bugg_boid len: {:?}", self.buff_boid.len());
        // println!("buff_meta len: {:?}", self.buff_meta.len());
        Ok(finished)
    }
}

impl Tracker for ReplayTracker {
    fn new(_: &[boid::Boid], _: &RunOptions) -> Self
    where
        Self: Sized,
    {
        todo!()
    }

    fn update(&mut self, run_options: &RunOptions) {
        if !(run_options.stop_movement || self.paused) {
            _ = self.load_buffer();
        }
    }

    fn insert_single(&mut self, _: boid::Boid, _: &RunOptions) {}

    fn restart(&mut self, _: &[boid::Boid]) {
        self.buff_playhead =
            ReplayTracker::get_start_playhead(self.buff_bumper, self.run_options.init_boids);
    }

    fn delete_last(&mut self) -> Option<boid::Boid> {
        todo!()
    }

    fn delete_multiple(&mut self, _: &[usize], _: &RunOptions) {
        todo!()
    }

    fn insert_multiple(&mut self, _: &[boid::Boid], _: &RunOptions) {
        todo!()
    }

    fn view<'a>(
        &'a self,
    ) -> Box<dyn Iterator<Item = (&'a boid::Boid, &'a boid::BoidMetadata)> + 'a> {
        // Box::new(self.buff_boid.iter().map(|e| (e, &self.buff_meta[e.id])))
        Box::new(self.buff_boid.iter().zip(self.buff_meta.iter()))
    }

    fn get_neighbours<'a>(&'a self, _: &boid::Boid, _: &RunOptions) -> Vec<&'a boid::Boid> {
        todo!()
    }

    fn signal(&mut self, signal: TrackerSignal) {
        match signal {
            TrackerSignal::ReplayTickForward => {
                if self.replay_direction == ReplayDirection::Backward {
                    self.buff_playhead += 1;
                }
                if !self.playhead_at_end() {
                    let direction = self.replay_direction;
                    self.replay_direction = ReplayDirection::Forward;
                    _ = self.load_buffer();
                    self.replay_direction = direction;
                }
                if self.replay_direction == ReplayDirection::Backward {
                    self.buff_playhead -= 1;
                }
            }
            TrackerSignal::ReplayTickBackward => {
                if !self.playhead_at_start() {
                    if self.replay_direction == ReplayDirection::Forward {
                        self.buff_playhead -= 1;
                    }
                    let direction = self.replay_direction;
                    self.replay_direction = ReplayDirection::Backward;
                    _ = self.load_buffer();
                    self.replay_direction = direction;
                    if self.replay_direction == ReplayDirection::Forward {
                        self.buff_playhead += 1;
                    }
                }
            }
            TrackerSignal::ReplayUpdatePlayhead(playhead) => {
                if playhead % self.run_options.init_boids != 0 {
                    // error
                }

                // >= 0 implicit by type
                if playhead < self.locations.len() {
                    self.buff_playhead = playhead;
                }
            }
            TrackerSignal::ReplaySwitchDirection => {
                self.replay_direction = match self.replay_direction {
                    ReplayDirection::Forward => ReplayDirection::Backward,
                    ReplayDirection::Backward => ReplayDirection::Forward,
                }
            }
            TrackerSignal::ReplayPlayPause => {
                self.paused = !self.paused;
            }
        }
    }

    fn get_no_entities(&self) -> usize {
        todo!()
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

#[cfg(test)]
mod test {
    use glam::{vec2, Mat2};

    use crate::{flock::tracker::Tracker, options::RunOptions};

    use super::ReplayTracker;

    #[test]
    fn matrix_scaling() {
        let rescale_mat = 3. * Mat2::IDENTITY;
        let heading = vec2(3., 4.);

        let expected = vec2(9., 12.);
        let result = rescale_mat * heading;

        assert_eq!(expected.x, result.x);
        assert_eq!(expected.y, result.y);
    }

    #[test]
    fn initialise() {
        let rp = "../boids-data.csv".to_owned();
        let ro = RunOptions::default();
        let _ = ReplayTracker::from_path(&rp, &ro);
        // assert!(rt.locations.len() != 0);
    }

    #[test]
    fn progress_forward() {
        let rp = "../boids-data.csv".to_owned();
        let ro = RunOptions::default();
        let rr = ReplayTracker::from_path(&rp, &ro);

        match rr {
            Err(err) => panic!("test, replay tracker returned an error: {:?}", err),
            Ok(mut rt) => {
                {
                    let first = rt.view().next();

                    match first {
                        None => panic!("failed to retrieve first pair"),
                        Some(tuple) => assert!(tuple.0.id == 0),
                    }
                }

                let res = rt.load_buffer();

                if let Err(rte) = res {
                    panic!("rte error: {}", rte)
                }

                let second = rt.view().next();

                match second {
                    None => panic!("failed to retrieve first pair"),
                    Some(tuple) => assert!(tuple.0.id == 0),
                }
            }
        }
    }
}
