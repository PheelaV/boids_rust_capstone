`mod fv_capstone;`
# Capstone project: Boids algorithms


Implementatioin of boids algorithms, with integrated modified version of the original Vicsek Model, utilisinng spatial hashing, with integration into R.

![cover](/docs/mellow_wandering.png)

Boid algorithms: Crag Reynolds 1987, [see](https://dl.acm.org/doi/10.1145/37402.37406) [see](http://www.red3d.com/cwr/0)
Vicsek model: Tamás Vicsek 1995, [see](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.75.1226)

Implementation inspired by:
  - Daniel Shiffman, Nature of Code, [see](https://natureofcode.com/book/)
  - Nannou, creative coding library for Rust inspired by Processing, [see](https://github.com/nannou-org/nannou)
  - Nature of Code, [see](https://github.com/nannou-org/nannou/tree/master/nature_of_code)
  - OpenSteer [see](https://sourceforge.net/projects/opensteer/)
  - Eisen Daniel, [see](https://github.com/eisendaniel/boids/blob/master/src/main.rs)
  - Craig Bester, [see](https://github.com/cycraig/rust-boids-wasm)
  - Hastings, Erin J. et al. “Optimization of Large-Scale , Real-Time Simulations by Spatial Hashing.” (2007)
  - Rust/R integration enabled by project extendr [see](https://github.com/extendr/extendr)


## How to run/compile:
### Boids simulation
1. if you do not have Rust installed, go to [rustup](https://rustup.rs/) and follow instructions for your platform
2. run `cargo run --release` in the root directory or `cargo run --releae -p boids_app`

### Boidranalysis
1. Ensure you have R installed, [ref](https://www.r-project.org/)
2. `boidr` is a package that wraps ```boids_lib``` and provides an interface into R for simulation execution and data retrieval. As of now it is unreleased and needs to be loaded for every new R session, an example is in the [experiment2](./boidranalysis/R/experiment2.R) at the very top bellow library imports.

## Keyboard controls
C - controls toggle
R - restart simulation
I - double the number of agents
D - delete half the agents
X - deselect agent
G - toggle tails (can't update number of agents whilist active)
F - clustering toggle
X - deselect clicked agent
F12 - density based colouring
[mouse click] - select an agent
Spacebar - pause simulation
Escape - close
Num1 - alignment toggle
Num2 - cohesion toggle
Num3 - separation toggle
Num4 - wander toggle
A - steering toggle
M - freeze movement, not updates
Return - switch replay direction
Right key - step forward
Left key - step backward
F5 - separation bias toggle (if anyone is within separation distance, will ignore cohesion)
F8 - diagnostics
F9 - agent labels
F10 - debug distance
F11 - debug grid

## Replay tracker (limited set)

C - controls toggle
Spacebar - pause replay
M - pause (separate from space)
X - deselect clicked agent
Escape - close
Return - switch replay direction
Right key - step forward
Left key - step backward
F - clustering toggle
F8 - diagnostics
F9 - agent labels
F10 - debug distance
F11 - debug grid
F12 - density based colours

![cover](/docs/arch.jpg)