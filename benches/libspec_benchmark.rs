//! Benchmark LibSpec/Library lookup performance
use bio::bio_types::sequence::Sequence;
use criterion::{Criterion, criterion_group, criterion_main};
use dnacomb::lib_spec;
use std::collections::HashMap;

fn bench_lookup(c: &mut Criterion) {
    let lib = lib_spec::Library::from_file("config/pegrna.tsv", HashMap::new(), 3)
        .expect("Expect library to load correctly");

    let exact_match: Sequence = vec![
        b'C', b'T', b'T', b'A', b'A', b'T', b'G', b'T', b'T', b'G', b'A', b'C', b'T', b'T', b'C',
        b'T', b'T', b'A', b'C', b'G', b'A', b'C', b'G', b'A', b'A',
    ];

    let partial_match: Sequence = vec![
        b'C', b'T', b'T', b'A', b'G', b'T', b'G', b'G', b'T', b'G', b'A', b'C', b'T', b'T', b'C',
        b'T', b'T', b'A', b'A', b'G', b'C', b'C', b'G', b'A', b'A',
    ];

    c.bench_function("lookup_hamming_exact", |b| {
        b.iter(|| {
            let _ = lib.lookup(
                "extension",
                &exact_match,
                lib_spec::DistanceMetric::Hamming,
                lib_spec::PartialMatching::Full,
            );
        });
    });

    c.bench_function("lookup_levenshtein_exact", |b| {
        b.iter(|| {
            let _ = lib.lookup(
                "extension",
                &exact_match,
                lib_spec::DistanceMetric::Levenshtein,
                lib_spec::PartialMatching::Full,
            );
        });
    });

    c.bench_function("lookup_bounded_levenshtein_exact", |b| {
        b.iter(|| {
            let _ = lib.lookup(
                "extension",
                &exact_match,
                lib_spec::DistanceMetric::BoundedLevenshtein,
                lib_spec::PartialMatching::Full,
            );
        });
    });

    c.bench_function("lookup_hamming_mismatch", |b| {
        b.iter(|| {
            let _ = lib.lookup(
                "extension",
                &partial_match,
                lib_spec::DistanceMetric::Hamming,
                lib_spec::PartialMatching::Full,
            );
        });
    });

    c.bench_function("lookup_levenshtein_mismatch", |b| {
        b.iter(|| {
            let _ = lib.lookup(
                "extension",
                &partial_match,
                lib_spec::DistanceMetric::Levenshtein,
                lib_spec::PartialMatching::Full,
            );
        });
    });

    c.bench_function("lookup_bounded_levenshtein_mismatch", |b| {
        b.iter(|| {
            let _ = lib.lookup(
                "extension",
                &partial_match,
                lib_spec::DistanceMetric::BoundedLevenshtein,
                lib_spec::PartialMatching::Full,
            );
        });
    });
}

criterion_group!(benches, bench_lookup);
criterion_main!(benches);
