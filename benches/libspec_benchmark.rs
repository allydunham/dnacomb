//! Benchmark LibSpec/Library lookup performance
use criterion::{criterion_group, criterion_main, Criterion};
use dnacomb::lib_spec;

fn bench_lookup(c: &mut Criterion) {
    c.bench_function("lookup_hamming_partial", |b| {
        b.iter(|| {
            let _ = 1;
        });
    });
}

criterion_group!(benches, bench_lookup);
criterion_main!(benches);

