use criterion::{black_box, criterion_group, criterion_main, Criterion};

// fn fibonacci(n: u64) -> u64 {
//     match n {
//         0 => 1,
//         1 => 1,
//         n => fibonacci(n-1) + fibonacci(n-2),
//     }
// }

fn fibonacci(n: u64) -> u64 {
    let mut a = 0;
    let mut b = 1;

    match n {
        0 => b,
        _ => {
            for _ in 0..n {
                let c = a + b;
                a = b;
                b = c;
            }
            b
        }
    }
}


fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("fib1 20", |b| b.iter(|| fibonacci(black_box(20))));
    // c.bench_function("fib2 20", |b| b.iter(|| fibonacci2(black_box(20))));

}


criterion_group!(benches, criterion_benchmark);
// criterion_group!(benches, criterion_benchmark2);
criterion_main!(benches);