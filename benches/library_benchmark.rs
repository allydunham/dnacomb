//! Benchmark LibSpec/Library lookup performance
use criterion::{Criterion, criterion_group, criterion_main};
use dnacomb::{
    interning::{SeqHandle, region_id_from_str, seq_from_bytes},
    library,
};
use std::collections::HashMap;

fn bench_lookup(c: &mut Criterion) {
    let lib = library::SubLibrary::from_file("config/pegrna.tsv", HashMap::new(), 3, None)
        .expect("Expect library to load correctly");

    let exact_match: SeqHandle = seq_from_bytes(&vec![
        b'C', b'T', b'T', b'A', b'A', b'T', b'G', b'T', b'T', b'G', b'A', b'C', b'T', b'T', b'C',
        b'T', b'T', b'A', b'C', b'G', b'A', b'C', b'G', b'A', b'A',
    ]);

    let partial_match: SeqHandle = seq_from_bytes(&vec![
        b'C', b'T', b'T', b'A', b'G', b'T', b'G', b'G', b'T', b'G', b'A', b'C', b'T', b'T', b'C',
        b'T', b'T', b'A', b'A', b'G', b'C', b'C', b'G', b'A', b'A',
    ]);

    let region = region_id_from_str("extension");

    c.bench_function("lookup_hamming_exact", |b| {
        b.iter(|| {
            let _ = lib.lookup(
                &region,
                exact_match,
                library::DistanceMetric::Hamming,
                library::PartialMatching::Full,
            );
        });
    });

    c.bench_function("lookup_levenshtein_exact", |b| {
        b.iter(|| {
            let _ = lib.lookup(
                &region,
                exact_match,
                library::DistanceMetric::Levenshtein,
                library::PartialMatching::Full,
            );
        });
    });

    c.bench_function("lookup_bounded_levenshtein_exact", |b| {
        b.iter(|| {
            let _ = lib.lookup(
                &region,
                exact_match,
                library::DistanceMetric::BoundedLevenshtein,
                library::PartialMatching::Full,
            );
        });
    });

    c.bench_function("lookup_hamming_mismatch", |b| {
        b.iter(|| {
            let _ = lib.lookup(
                &region,
                partial_match,
                library::DistanceMetric::Hamming,
                library::PartialMatching::Full,
            );
        });
    });

    c.bench_function("lookup_levenshtein_mismatch", |b| {
        b.iter(|| {
            let _ = lib.lookup(
                &region,
                partial_match,
                library::DistanceMetric::Levenshtein,
                library::PartialMatching::Full,
            );
        });
    });

    c.bench_function("lookup_bounded_levenshtein_mismatch", |b| {
        b.iter(|| {
            let _ = lib.lookup(
                &region,
                partial_match,
                library::DistanceMetric::BoundedLevenshtein,
                library::PartialMatching::Full,
            );
        });
    });
}

criterion_group!(benches, bench_lookup);
criterion_main!(benches);
