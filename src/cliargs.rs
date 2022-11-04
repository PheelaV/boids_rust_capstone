use clap::{Parser, arg};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args{
    #[arg(long, default_value_t = 256)]
    pub no_boids: u32
}