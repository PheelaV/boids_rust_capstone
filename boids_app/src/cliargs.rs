// TODO: provide config file argument
// https://docs.rs/clap-serde/latest/clap_serde/#
// https://stackoverflow.com/questions/55133351/is-there-a-way-to-get-clap-to-use-default-values-from-a-file
// use clap::{Parser, arg};
use clap_serde_derive::{
    clap::{self, Parser},
    serde::Serialize,
    ClapSerde,
};

// #[derive(Parser, Debug)]
#[derive(Parser)]
#[derive(ClapSerde)]
#[command(author="PheelaV", version, about, long_about = None)]
/// Implementation and visualisation of the Boids (Reynolds '86) algorithms.
pub struct Args{
    /// Config file
    #[arg(short, long = "config", default_value = "config.toml")]
    pub config_path: std::path::PathBuf,

    /// Rest of arguments
    #[command(flatten)]
    pub config: <Config as ClapSerde>::Opt,
}

#[derive(ClapSerde, Serialize)]
/// Programatic configuration
/// 
/// Uses defaults, which can be overwritten by specifying a filepath for the `-c` or `--config` arg option
pub struct Config{
     
    #[default(128)]
    #[arg(short = 'n', long)]
    /// number of boids
    pub no_boids: usize,

    #[default(4)]
    #[arg(short = 'r', long)]
    /// ratio of renders/sample_rate, e,g, 4 = sample every 4th render
    pub sample_rate: u64,

    #[default(true)]
    #[arg(short = 's', long)]
    pub save: bool,
    
    #[default(false)]
    #[arg(short = 't', long)]
    pub save_timestamp: bool,
   
    #[default(1350)]
    #[arg(short = 'x', long)]
    pub init_width: u32,
    
    #[default(900)]
    #[arg(short = 'y', long)]
    pub init_height: u32,

    #[default(60.)]
    #[arg(long = "sens_dist")]
    pub sensory_distance: f32,

    #[default(0.02)]
    #[arg(long = "ali_coef")]
    pub allignment_coefficient: f32,
    #[default(0.002)]
    #[arg(long = "coh_coef")]
    pub cohesion_coefficient: f32,
    #[default(4.1)]
    #[arg(long = "sep_coef")]
    pub separation_coefficient: f32,

    #[default(1.15)]
    #[arg(long = "ali_trs_coef")]
    pub allignment_treshold_coefficient: f32,
    #[default(0.95)]
    #[arg(long = "coh_trs_coef")]
    pub cohesion_treshold_coefficient: f32,
    #[default(0.35)]
    #[arg(long = "sep_trs_coef")]
    pub separation_treshold_coefficient: f32,

    #[default(0.65)]
    #[arg(long = "min_speed")]
    pub min_speed: f32,
    #[default(2.35)]
    #[arg(long = "max_speed")]
    pub max_speed: f32,
    #[default(0.65)]
    #[arg(long = "max_steering")]
    pub max_steering: f32,
    #[default(180.)]
    #[arg(long = "fov")]
    pub field_of_vision: f32,
    #[default(false)]
    #[arg(long = "dbscan")]
    pub dbscan_flock_clustering_on: bool,
    #[default(false)]
    #[arg(long = "wander")]
    pub wander_on: bool, 
    #[default(0.03)]
    #[arg(long = "wrate")]
    pub wander_rate: f32, 
    #[default(5.2)]
    #[arg(long = "wradius")]
    pub wander_radius: f32, 
    #[default(21.5)]
    #[arg(long = "wdistance")]
    pub wander_distance: f32, 
    #[default(0.4)]
    #[arg(long = "wcoeff")]
    pub wander_coefficient: f32,
    #[default(8.)]
    #[arg(long = "size")]
    /// number of boids
    pub size: f32, 
}