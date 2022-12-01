use std::fmt;

use boids_lib::flock_base;
use boids_lib::options;
use boids_lib::options::RunOptions;
use boids_lib::options::SaveOptions;
use criterion::criterion_group;
use criterion::criterion_main;
use criterion::BenchmarkId;
use criterion::Criterion;
use criterion::SamplingMode;
use criterion::Throughput;

#[derive(Debug)]
struct FlockBench {
    no_boids: usize,
    no_iter: u64,
}

impl fmt::Display for FlockBench {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "(no_boids: {}, no_iter: {})",
            self.no_boids, self.no_iter
        )
    }
}

fn get_bench_params() -> Vec<FlockBench> {
    static NO_BOIDS: u64 = 64;

    let mut result: Vec<FlockBench> = Vec::new();

    // [128, 1024, 4096].iter().for_each(|no_boids| {
        [128, 512].iter().for_each(|no_boids|{
        [NO_BOIDS, 2 * NO_BOIDS, 4 * NO_BOIDS]
            .iter()
            .for_each(|no_iter| {
                result.push(FlockBench {
                    no_iter: *no_iter,
                    no_boids: *no_boids,
                })
            })
    });
    result
}

// benchmarks different no_iter and no_boids configurations
// fn from_elem(c: &mut Criterion) {
//     let mut group = c.benchmark_group("test_matrix");

//     for nn in get_bench_params().iter() {
//         // group.measurement_time(Duration::from_secs(40));
//         group.sampling_mode(SamplingMode::Flat);
//         // group.plot_config(Plot)
//         group.throughput(Throughput::Elements(
//             (nn.no_boids as u64 * nn.no_iter) as u64,
//         ));
//         group.bench_with_input(BenchmarkId::from_parameter(nn), nn, |b, n| {
//             b.iter(|| {
//                 flock_base(n.no_iter, {
//                     let mut ro: RunOptions = Default::default();
//                     ro.init_boids = n.no_boids;
//                     ro.window = options::get_window_size(5000., 5000.);
//                     ro.save_options = SaveOptions {
//                         save_locations: false,
//                         save_locations_path: None,
//                         save_locations_timestamp: false,
//                     };
//                     ro.sensory_distance = 20.;
//                     ro
//                 })
//             });
//         });
//     }
//     group.finish();
// }

// benchmarks the effect of FOV on or off
// fn from_elem2(c: &mut Criterion) {
//     let mut group = c.benchmark_group("test_FOV");

//     let no_boids: usize = 8096;
//     let no_iter: u64 = 1024;

//     for fov in [true, false] {
//         // group.measurement_time(Duration::from_secs(40));
//         group.sampling_mode(SamplingMode::Flat);
//         // group.plot_config(Plot)
//         group.throughput(Throughput::Elements(
//             no_boids as u64 * no_iter
//         ));
//         group.bench_with_input(BenchmarkId::from_parameter(fov), &fov, |b, fov| {
//             b.iter(|| {
//                 flock_base(no_iter, {
//                     let mut ro: RunOptions = Default::default();
//                     ro.init_boids = no_boids;
//                     ro.window = options::get_window_size(5000., 5000.);
//                     ro.save_options = SaveOptions {
//                         save_locations: false,
//                         save_locations_path: None,
//                         save_locations_timestamp: false,
//                     };
//                     ro.sensory_distance = 20.;
//                     ro.field_of_vision_on = *fov;
//                     ro.neighbours_cosidered = 30;
//                     ro
//                 })
//             });
//         });
//     }
//     group.finish();
// }
// benchmarks the effect no neighbours to take into consideration
fn from_elem3(c: &mut Criterion) {
    let mut group = c.benchmark_group("test_FOV");

    let no_boids: usize = 8192;
    let no_iter: u64 = 1024;

    for no_neighbours_considered in [0, 4, 8, 16, 32] {
        // group.measurement_time(Duration::from_secs(40));
        group.sampling_mode(SamplingMode::Flat);
        // group.plot_config(Plot)
        group.throughput(Throughput::Elements(
            no_boids as u64 * no_iter
        ));
        group.bench_with_input(BenchmarkId::from_parameter(no_neighbours_considered), &no_neighbours_considered, |b, nnc| {
            b.iter(|| {
                flock_base(no_iter, {
                    let mut ro: RunOptions = Default::default();
                    ro.init_boids = no_boids;
                    ro.window = options::get_window_size(5000., 5000.);
                    ro.save_options = SaveOptions {
                        save_locations: false,
                        save_locations_path: None,
                        save_locations_timestamp: false,
                    };
                    ro.sensory_distance = 20.;
                    ro.field_of_vision_on = true;
                    ro.neighbours_cosidered = *nnc;
                    ro
                })
            });
        });
    }
    group.finish();
}

criterion_group!(benches2, from_elem3);
// criterion_group!(benches, from_elem);
// criterion_main!(benches, benches2);
criterion_main!(benches2);
