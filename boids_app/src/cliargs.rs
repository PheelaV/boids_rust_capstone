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
#[derive(Parser, ClapSerde)]
#[command(author="PheelaV", version, about, long_about = None)]
/// Implementation and visualisation of the Boids (Reynolds '86) algorithms.
pub struct Args {
    /// Config file
    #[arg(short, long = "config", default_value = "config.toml")]
    pub config_path: std::path::PathBuf,

    /// Rest of arguments
    #[command(flatten)]
    pub config: <Config as ClapSerde>::Opt,
}

/// Programatic configuration
///
/// Uses defaults, which can be overwritten by specifying a filepath for the `-c` or `--config` arg option
#[derive(ClapSerde, Serialize)]
pub struct Config {
    #[default(1024)]
    #[arg(short = 'n', long)]
    /// number of boids
    pub no_boids: usize,

    #[default(64)]
    #[arg(short = 'r', long)]
    /// ratio of renders/sample_rate, e,g, 4 = sample every 4th render
    pub sample_rate: u16,

    #[default(false)]
    #[arg(short = 's', long)]
    pub save: bool,

    #[default(false)]
    #[arg(short = 't', long)]
    pub save_timestamp: bool,

    // #[default(1350)]
    #[default(800)]
    #[arg(short = 'x', long)]
    pub init_width: u32,

    #[default(800)]
    #[arg(short = 'y', long)]
    pub init_height: u32,

    #[default(35.)]
    #[arg(long = "sens_dist")]
    pub sensory_distance: f32,

    #[default(0.2)]
    #[arg(long = "ali_coef")]
    pub alignment_coefficient: f32,
    #[default(0.1)]
    #[arg(long = "coh_coef")]
    pub cohesion_coefficient: f32,
    #[default(2.5)]
    #[arg(long = "sep_coef")]
    pub separation_coefficient: f32,

    #[default(1.)]
    #[arg(long = "ali_trs_coef")]
    pub alignment_treshold_coefficient: f32,
    #[default(0.8)]
    #[arg(long = "coh_trs_coef")]
    pub cohesion_treshold_coefficient: f32,
    #[default(0.3)]
    #[arg(long = "sep_trs_coef")]
    pub separation_treshold_coefficient: f32,

    #[default(120.0)]
    #[arg(long = "max_speed")]
    /// Maximum speed in pixels/second
    pub max_speed: f32,
    #[default(40.0)]
    #[arg(long = "max_steering")]
    /// Maximum steering acceleration in pixels/second²
    pub max_steering: f32,
    #[default(false)]
    #[arg(long = "steer")]
    pub agent_steering: bool,
    #[default(280.)]
    #[arg(long = "fov")]
    pub field_of_vision: f32,
    #[default(true)]
    #[arg(long = "dbscan")]
    pub dbscan_flock_clustering_on: bool,
    #[default(false)]
    #[arg(long = "wander")]
    pub wander_on: bool,
    #[default(false)]
    #[arg(long = "wander_random")]
    pub wander_random: bool,
    #[default(0.4)]
    #[arg(long = "wrate")]
    pub wander_rate: f32,
    #[default(11.5)]
    #[arg(long = "wradius")]
    pub wander_radius: f32,
    #[default(13.)]
    #[arg(long = "wdistance")]
    pub wander_distance: f32,
    #[default(1.0)]
    #[arg(long = "wcoeff")]
    pub wander_coefficient: f32,
    #[default(5.)]
    #[arg(long = "size")]
    /// number of boids
    pub size: f32,
    #[default(false)]
    #[arg(long = "rules_impl")]
    pub rules_impl: bool,

    // Physics: Natural motion
    #[default(0.02)]
    #[arg(long = "vel_retention")]
    /// Velocity retained per second (0.02 = 2% remains, 0.9 = 90% remains)
    pub velocity_retention: f32,

    #[default(50.0)]
    #[arg(long = "forward_drive")]
    /// Constant forward acceleration in heading direction
    pub forward_drive: f32,

    #[default(60.0)]
    #[arg(long = "target_fps")]
    /// Target FPS for fixed timestep mode
    pub target_fps: f32,

    #[default(true)]
    #[arg(long = "fixed_timestep")]
    /// Use fixed timestep vs variable dt
    pub fixed_timestep: bool,

    // SI Units
    #[default(0.01)]
    #[arg(long = "meters_per_pixel")]
    /// Spatial scale: meters per pixel (0.01 = 1px = 1cm)
    pub meters_per_pixel: f32,
}
