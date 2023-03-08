use std::{fs::OpenOptions, mem};

use chrono::{DateTime, Utc};
use serde::Serialize;

use crate::{flock::Flock, options::SaveOptions};

// todo: it would be cool to have bird watcher store a reference to the flock, both would have to have the same lifetime anotations
// the problem is that the nannou Model would also have to be annotated

// so right now, this is more of a bird data acummulator than a birdwatcher
#[derive(Serialize, Debug, Clone, Copy)]
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
    // flock: &'a Flock,
}

const PREFIX: &'static str = "boids-data";

impl Birdwatcher {
    pub fn new(sample_rate: u64) -> Self {
        Birdwatcher {
            locations: Vec::new(),
            render_ticker: 0,
            sample_rate: sample_rate,
            // flock
        }
    }

    /// Triggers data collection
    pub fn watch(&mut self, flock: &Flock) -> () {
        if !self.should_sample() {
            return;
        }

        let mut current_locations: Vec<BoidData> = flock
            .view()
            .0
            .iter()
            .zip(flock.view().1.iter())
            .map(|(e, m)| BoidData {
                id: e.id,
                x: e.position.x,
                y: e.position.y,
                cluster_id: m.cluster_id,
                n_neighbours: m.n_neighbours,
                time: self.render_ticker / self.sample_rate,
            })
            .collect();

        self.locations.append(&mut current_locations);
    }

    pub fn restart(&mut self) {
        self.locations.clear();
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
