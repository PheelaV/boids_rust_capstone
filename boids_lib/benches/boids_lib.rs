use std::fmt;

use boids_lib::flock_base;
use boids_lib::options::RunOptions;
use criterion::BenchmarkId;
use criterion::Criterion;
use criterion::SamplingMode;
use criterion::Throughput;
use criterion::criterion_group;
use criterion::criterion_main;

#[derive(Debug)]
struct flock_bench {
    no_boids: u32,
    no_iter: u32
}

impl fmt::Display for flock_bench {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "(no_boids: {}, no_iter: {})", self.no_boids, self.no_iter)
    }
}


fn get_bench_params() -> Vec<flock_bench> {
    static NO_BOIDS: u32 = 64;

    let mut result : Vec<flock_bench> = Vec::new();

    // [128, 1024, 4096].iter().for_each(|no_boids|{
    [128, 512].iter().for_each(|no_boids|{
        [NO_BOIDS, 2 * NO_BOIDS, 4 * NO_BOIDS].iter().for_each(|no_iter|{
            result.push(flock_bench{
                no_iter: *no_iter ,
                no_boids: *no_boids
            })
        }) 
    });
    result
}
fn from_elem(c: &mut Criterion) {
    let mut group = c.benchmark_group("test_matrix");

    for nn in get_bench_params().iter() {
        // group.measurement_time(Duration::from_secs(40));
        group.sampling_mode(SamplingMode::Flat);
        // group.plot_config(Plot)
        group.throughput(Throughput::Elements((nn.no_boids * nn.no_iter) as u64));
        group.bench_with_input(BenchmarkId::from_parameter(nn), nn, |b, n| {
            b.iter(|| flock_base(n.no_iter, {
                let mut ro: RunOptions = Default::default();
                ro.init_boids = n.no_boids;
                ro
            }));
        });
    }
    group.finish();
}

criterion_group!(benches, from_elem);
criterion_main!(benches);
