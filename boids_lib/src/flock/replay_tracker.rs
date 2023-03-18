use std::path::Path;

use super::tracker::{Tracker, self, TrackerSignal};
use crate::{
    birdwatcher::BoidData,
    boid::{Boid, BoidMetadata},
    options::{RunOptions, get_window_size}, math_helpers::tor_vec_p,
};
use csv::Reader;
use glam::vec2;
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

pub enum ReplayDirection {
    Forward,
    Backward
}
pub struct ReplayTracker {
    pub locations: Vec<BoidData>,
    // pub buffer: Vec<&BoidData>,
    pub buff_boid: Vec<Boid>,
    pub buff_meta: Vec<BoidMetadata>,
    run_options: RunOptions,
    buff_playhead: usize,
    replay_direction: ReplayDirection

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
            .sorted_unstable_by_key(|e| (e.time, e.id))
            .collect::<Vec<BoidData>>();

        let mut ro = run_options.to_owned();

        let (max_x, max_y) = locations.iter().fold::<(f32, f32), _>((0., 0.), |acc: (f32, f32), cur| (cur.x.max(acc.0), cur.y.max(acc.1)));

        ro.window = get_window_size(max_x.round() as u32, max_y.round() as u32);
        // scan the items until one repetition is found, assume the number of 
        // items skipped is the number of boids per iteration
        let no_boids = locations.iter()
        .enumerate()
        .find_or_first(|(index, d)| index != &0 && d.id == locations[0].id)
        .expect("")
        .0;

        // overwrite defaults
        ro.init_boids = no_boids;
        
        let mut ret = ReplayTracker {
            locations: locations,
            // buffer: Vec::new(),
            buff_boid: vec![Boid::default(); ro.init_boids],
            buff_meta: vec![BoidMetadata::default(); ro.init_boids],
            run_options: ro,
            buff_playhead: 0,
            replay_direction: ReplayDirection::Forward
        };

        let res = ret.load_buffer();

        if let Err(err) = res {
            return Err(err)
        }

        if ret.buff_boid.len() % ret.run_options.init_boids != 0 {
            return Err(ReplayTrackerError::DataBlockError(
                "failed to initialise".to_owned(),
            ));
        }

        return Ok(ret)
    }

    fn playhead_at_end(&self) -> bool {
        match self.replay_direction {
            ReplayDirection::Forward => self.buff_playhead == self.locations.len() - self.run_options.init_boids,
            ReplayDirection::Backward => self.buff_playhead == 0, 
        }
    }

    fn playhead_at_start(&self) -> bool {
        match self.replay_direction {
            ReplayDirection::Forward => self.buff_playhead == 0,
            ReplayDirection::Backward => self.buff_playhead == self.locations.len() - self.run_options.init_boids,
        }
    }
    fn load_buffer(&mut self) -> Result<bool, ReplayTrackerError> {
        let RunOptions { 
            init_boids,
            window,
             .. } = self.run_options;

        let finished = self.playhead_at_end(); 
        let starting = self.playhead_at_start();

        // does the next step have enough data points?
        if self.locations.len() - self.buff_playhead < init_boids {
            return Err(ReplayTrackerError::DataBlockError(
                "next step does not have enough data points".to_owned(),
            ));
        }

        let time = self.locations[self.buff_playhead].time;

        for (index, data) in self
            .locations
            .iter()
            .skip(self.buff_playhead)
            .take(init_boids)
            .enumerate()
        {
            if index != data.id {
                return Err(ReplayTrackerError::DataBlockError(
                    "id missalignment".to_owned(),
                ));
            }

            if data.time != time {
                return Err(ReplayTrackerError::DataBlockError(
                    "time missalignment".to_owned(),
                ));
            }

            let boid = self.buff_boid.get_mut(data.id).unwrap();
            let meta = self.buff_meta.get_mut(data.id).unwrap();

            boid.id = data.id;
            boid.position.x = data.x;
            boid.position.y = data.y;
            boid.velocity = vec2(2., 2.);

            // depending on replay direction and on whether or not we are in a starting position,
            // set the velocity by taking a lag(-1) of the difference of positions
            boid.velocity = match (&self.replay_direction, starting) {
                (_, true) => {
                    vec2(0., 0.)
                },
                (ReplayDirection::Forward, false) => {
                    let last_index = (time as usize - 1)  * init_boids - init_boids + data.id;
                    let last_location = self.locations[last_index];
                    let vec = tor_vec_p(last_location.x, data.x, last_location.y, data.y, &window);
                    vec2(vec.0, vec.1)
                },
                (ReplayDirection::Backward, false) => {
                    let last_index = (time as usize) * init_boids + data.id;
                    let last_location = self.locations[last_index];
                    let vec = tor_vec_p(last_location.x, data.x, last_location.y, data.y, &window);
                    vec2(vec.0, vec.1)
                },
            };
            meta.id = data.id;
            meta.n_neighbours = data.n_neighbours;
            meta.cluster_id = data.cluster_id;
        }

        // move buffer position
        if !finished {
            self.buff_playhead = match self.replay_direction {
                ReplayDirection::Forward => self.buff_playhead + init_boids,
                ReplayDirection::Backward => self.buff_playhead - init_boids,
            };
        }

        Ok(finished)
    }
}

impl Tracker for ReplayTracker {
    fn new(entities: &[crate::boid::Boid], run_options: &RunOptions) -> Self
    where
        Self: Sized,
    {
        todo!()
    }

    fn update(&mut self, run_options: &RunOptions) {
        let res = self.load_buffer();
    }

    fn insert_single(&mut self, _: crate::boid::Boid, _: &RunOptions) {}

    fn restart(&mut self, _: &[crate::boid::Boid]) {
        self.buff_playhead = 0; // todo: this will change depending on direction
    }

    fn delete_last(&mut self) -> Option<crate::boid::Boid> {
        todo!()
    }

    fn delete_multiple(&mut self, _: &[usize], _: &RunOptions) {
        todo!()
    }

    fn insert_multiple(&mut self, _: &[crate::boid::Boid], _: &RunOptions) {
        todo!()
    }

    fn view(&self) -> (&Vec<crate::boid::Boid>, &Vec<crate::boid::BoidMetadata>) {
        (&self.buff_boid, &self.buff_meta)
    }

    fn view2<'a>(
        &'a self,
    ) -> Box<dyn Iterator<Item = (&'a crate::boid::Boid, &'a crate::boid::BoidMetadata)> + 'a> {
        Box::new(self.buff_boid.iter().map(|e| (e, &self.buff_meta[e.id])))
    }

    fn get_neighbours<'a>(
        &'a self,
        _: &crate::boid::Boid,
        _: &RunOptions,
    ) -> Vec<&'a crate::boid::Boid> {
        todo!()
    }

    fn signal(&mut self, signal: TrackerSignal) {
        match signal {
            TrackerSignal::ReplayTickForward => {
                if !self.playhead_at_end() {
                    self.buff_playhead += self.run_options.init_boids;
                }
            },
            TrackerSignal::ReplayTickBackward => {
                if !self.playhead_at_start() {
                    self.buff_playhead -= self.run_options.init_boids;
                }
            },
            TrackerSignal::ReplayUpdatePlayhead(playhead) => {
                if playhead % self.run_options.init_boids != 0 {
                    // error
                }

                // >= 0 implicit by type
                if playhead < self.locations.len() {
                    self.buff_playhead = playhead;
                }
            },
        }
    }
}

#[cfg(test)]
mod test {
    use crate::{flock::tracker::Tracker, options::RunOptions};

    use super::ReplayTracker;

    #[test]
    fn initialise() {
        let rp = "../boids-data.csv".to_owned();
        let ro = RunOptions::default();
        let rr = ReplayTracker::from_path(&rp, &ro);
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
                    let first = rt.view2().next();

                    match first {
                        None => panic!("failed to retrieve first pair"),
                        Some(tuple) => assert!(tuple.0.id == 0),
                    }
                }

                let res = rt.load_buffer();

                if let Err(rte) = res {
                    panic!("rte error: {}", rte)
                }

                let second = rt.view2().next();

                match second {
                    None => panic!("failed to retrieve first pair"),
                    Some(tuple) => assert!(tuple.0.id == 0),
                }
            }
        }
    }
}
