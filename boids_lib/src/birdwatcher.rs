use std::{fs::OpenOptions, mem};
use chrono::{DateTime, Utc};
use serde::{Serialize, Deserialize};
use crate::{flock::Flock, options::SaveOptions};

// so right now, this is more of a bird data acummulator than a birdwatcher
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct BoidData {
    pub id: usize,
    pub x: f32,
    pub y: f32,
    pub cluster_id: usize,
    pub n_neighbours: usize,
    pub time: u64,
}

pub struct Birdwatcher {
    locations: Vec<BoidData>,
    render_ticker: u64,
    sample_rate: u64,
    // ghost_queue: Option<CircularQueue<(Boid, BoidMetadata)>>
    // flock: &'a Flock,
}

const PREFIX: &'static str = "boids-data";

impl Birdwatcher {
    pub fn new(sample_rate: u16) -> Self {
        Birdwatcher {
            locations: Vec::new(),
            render_ticker: 0,
            sample_rate: sample_rate as u64,
            // ghost_queue: None,
            // flock
        }
    }

    /// Triggers data collection
    pub fn watch(&mut self, flock: &Flock) -> () {
        if !self.should_sample() {
            return;
        }

        // if let None = self.ghost_queue {
        //     self.ghost_queue =
        // }

        self.locations.extend(flock.view2().map(|(e, m)| {
            if e.id != m.id {
                panic!("mismatch of boid and metadata id");
            }

            BoidData {
                id: e.id,
                x: e.position.x,
                y: e.position.y,
                cluster_id: m.cluster_id,
                n_neighbours: m.n_neighbours,
                time: self.render_ticker / self.sample_rate,
            }
        }));
        // let mut current_locations = flock.view2().map(|(e, m)| {
        //     if e.id != m.id {
        //         panic!("mismatch of boid and metadata id");
        //     }

        //     BoidData {
        //         id: e.id,
        //         x: e.position.x,
        //         y: e.position.y,
        //         cluster_id: m.cluster_id,
        //         n_neighbours: m.n_neighbours,
        //         time: self.render_ticker / self.sample_rate,
        //     }
        // }).collect_vec();

        // self.locations.append(&mut current_locations);
    }

    // pub fn get_ghosts<'a>(&'a  mut self, no_ghosts: usize, run_options: &RunOptions) 
    // -> Option<Box<dyn Iterator<Item = (&'a Boid, &'a BoidMetadata)> + 'a> > {
    //     let RunOptions { 
    //         init_boids, 
    //         ..
    //     } = run_options;

    //     let time = self.locations.len() / *init_boids;

    //     if time < no_ghosts { return None }

    //     if let Some(gq) = &mut self.ghost_queue {

    //     } else {
    //         let mut gq  = CircularQueue::with_capacity(no_ghosts * init_boids);

    //         let skip = self.locations.len() - (no_ghosts + 1) * init_boids;
 
            
    //         self.ghost_queue = Some(gq);
    //     }

    //     // assuming number of agents is frozen since the beginning of the run (after start or reset)

        


    //     // let to_skip = (init_boids * (time - 1)) - (no_ghosts * init_boids);

    //     todo!();
    // }

    pub fn restart(&mut self) {
        self.locations.clear();
        self.render_ticker = 0;
    }

    pub fn pop_data(&mut self) -> Vec<BoidData> {
        mem::take(&mut self.locations)
    }

    /// Saves the latest data in CSV format, then returns it while emptying the birdwatcher's memory
    ///
    /// Depending on save options, either attempts to overwrite the current file or write's a new timestamped file
    pub fn pop_data_save(&mut self, save_options: &SaveOptions) -> Vec<BoidData> {
        let data = self.pop_data();
        // let data = self.pop_data();
        // // sort by id for convenience
        // data.sort_by(|a, b| a.id.cmp(&b.id));
        // data.sort_by(|a, b| a.t.cmp(&b.t));

        if !save_options.save_locations {
            return data;
        }

        if let Some(path) = &save_options.save_locations_path {
            let file_path = format!(
                "{path}{file_name}",
                file_name = Birdwatcher::get_dataset_name(save_options, Utc::now())
            );

            // open file
            let file = OpenOptions::new()
                .write(true)
                .create(true)
                .append(false)
                .open(file_path)
                .expect("Can't open file");
            let mut wtr = csv::Writer::from_writer(file);
            // let mut wtr = csv::Writer::from_writer(io::stdout());

            // write data points
            data.iter().for_each(|b| {
                wtr.serialize(b).expect("Can't serialize data point");
            });
            wtr.flush().expect("Can't write data file");
        }

        data
    }

    fn get_dataset_name(save_options: &SaveOptions, now: DateTime<Utc>) -> String {
        match save_options.save_locations_timestamp {
            true => {
                let datetime_part = now.timestamp_millis();
                format!(
                    "{prefix}_{datetime}.csv",
                    prefix = PREFIX,
                    datetime = datetime_part
                )
            }
            false => format!("{prefix}.csv", prefix = PREFIX),
        }
    }

    fn should_sample(&mut self) -> bool {
        self.render_ticker += 1;

        if self.render_ticker % self.sample_rate == 0 {
            true
        } else {
            false
        }
    }
}
#[cfg(test)]
#[allow(deprecated)]

mod tests {
    use crate::{birdwatcher::Birdwatcher, options::SaveOptions};
    use chrono::prelude::*;
    use chrono::Utc;

    #[test]
    fn test_name_timestamped() {
        let expected = "boids-data_1668038059490.csv";
        let save_options = SaveOptions {
            save_locations: true,
            save_locations_path: Some("".to_owned()),
            save_locations_timestamp: true,
        };
        let dt = Utc.ymd(2022, 11, 09).and_hms_milli_opt(23, 54, 19, 490);
        let actual = Birdwatcher::get_dataset_name(&save_options, dt.unwrap());

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_name_overwrite() {
        let expected = "boids-data.csv";
        let save_options = SaveOptions {
            save_locations: true,
            save_locations_path: Some("".to_owned()),
            save_locations_timestamp: false,
        };
        let dt = Utc.ymd(2022, 11, 09).and_hms_milli_opt(23, 54, 19, 490);
        let actual = Birdwatcher::get_dataset_name(&save_options, dt.unwrap());

        assert_eq!(actual, expected)
    }
}
