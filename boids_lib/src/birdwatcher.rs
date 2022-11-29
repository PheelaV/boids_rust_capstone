use std::{mem, fs::OpenOptions};

use chrono::{Utc, DateTime};
use serde::Serialize;

use crate::{flock::Flock, options::SaveOptions};

// todo: it would be cool to have bird watcher store a reference to the flock, both would have to have the same lifetime anotations
// the problem is that the nannou Model would also have to be annotated 

// so right now, this is more of a bird data acummulator than a birdwatcher
#[derive(Serialize, Debug, Clone, Copy)]
pub struct BoidData {
    pub id: u32,
    pub x: f32,
    pub y: f32,
    pub cluster_id: usize,
}

pub struct Birdwatcher{
    locations: Vec<BoidData>,
    render_ticker: u128,
    sample_rate: u32,
    // flock: &'a Flock,
}

const PREFIX: &'static str = "boids-data"; 

impl Birdwatcher {
    pub fn new(sample_rate: u32) -> Self {

        Birdwatcher
        { 
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

        let mut current_locations : Vec<BoidData> = 
            flock.view().iter()
            .map(|b| {
                BoidData{
                    id: b.id, 
                    x: b.position.x, 
                    y:b.position.y,
                    cluster_id: b.cluster_id
                }
            })
            .collect();
     
        self.locations.append(&mut current_locations);
    }

    pub fn pop_data(&mut self) -> Vec<BoidData> {
        mem::take(&mut self.locations)
    }

    /// Saves the latest data in CSV format, then returns it while emptying the birdwatcher's memory
    /// 
    /// Depending on save options, either attempts to overwrite the current file or write's a new timestamped file
    pub fn pop_data_save(&mut self, save_options: &SaveOptions) -> Vec<BoidData> {
        let mut data = self.pop_data();
        // sort by id for convenience
        data.sort_by(|a, b| a.id.cmp(&b.id));

        if !save_options.save_locations {
            return data
        }

        let path = match &save_options.save_locations_path {
            Some(path) => path,
            None => "./",
        };

        let file_path = format!("{path}{file_name}", file_name = Birdwatcher::get_dataset_name(save_options, Utc::now())); 

        // open file
        let file = OpenOptions::new()
            .write(true)
            .create(true)
            .open(file_path)
            .expect("Can't open file");
        let mut wtr = csv::Writer::from_writer(file);
        // let mut wtr = csv::Writer::from_writer(io::stdout());
  
        // write data points
        data
        .iter()
        .for_each(|b|{
                wtr.serialize(b).expect("Can't serialize data point");
        });
        wtr.flush().expect("Can't write data file");
    
        data
    }

    fn get_dataset_name(save_options: &SaveOptions, now: DateTime<Utc>) -> String {
        match save_options.save_locations_timestamp {
            true => {
                let datetime_part = now.timestamp_millis();
                format!("{prefix}_{datetime}.csv", prefix= PREFIX, datetime = datetime_part)
            },
            false => format!("{prefix}.csv", prefix= PREFIX)
        }
    }

    fn should_sample(&mut self) -> bool {
        self.render_ticker += 1;

        if self.render_ticker % self.sample_rate as u128 == 0 {
            true
        } else {
            false
        }
    }
}
#[cfg(test)]
mod tests {
    use chrono::Utc;
    use chrono::prelude::*;
    use crate::{birdwatcher::Birdwatcher, options::SaveOptions};

    #[test]
    fn test_name_timestamped() {
        let expected = "boids-data_1668038059490.csv";
        let save_options = SaveOptions {save_locations: true, save_locations_path: Some("".to_owned()), save_locations_timestamp: true};
        let dt =  Utc.ymd(2022, 11, 09).and_hms_milli(23, 54, 19, 490);
        let actual = Birdwatcher::get_dataset_name(&save_options, dt);

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_name_overwrite() {
        let expected = "boids-data.csv";
        let save_options = SaveOptions {save_locations: true, save_locations_path: Some("".to_owned()), save_locations_timestamp: false};
        let dt =  Utc.ymd(2022, 11, 09).and_hms_milli(23, 54, 19, 490);
        let actual = Birdwatcher::get_dataset_name(&save_options, dt);
        
        assert_eq!(actual, expected)
    }
}
