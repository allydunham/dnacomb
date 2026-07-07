//! Benchmark library lookup
use criterion::{BenchmarkId, Criterion, Throughput, black_box, criterion_group, criterion_main};
use dnacomb::interning::{region_id_from_str, seq_from_bytes};
use dnacomb::library::{DistanceMetric, PartialMatching, SubLibrary};
use std::collections::HashMap;

mod bench_support;
use bench_support::*;

#[derive(Clone)]
struct LookupFixture {
    lib: SubLibrary,
    query: Vec<u8>,
    region: dnacomb::interning::RegionID,
    name: String,
}

impl LookupFixture {
    fn new(lib: &SubLibrary, query: &[u8], name: &str) -> Self {
        Self {
            lib: lib.clone(),
            query: query.to_vec(),
            region: region_id_from_str("r1"),
            name: name.to_string(),
        }
    }
}

fn build_sublib_from_seqs(seqs: Vec<Vec<u8>>, max_dist: u64) -> SubLibrary {
    let region = region_id_from_str("r1");
    let lib = HashMap::from([(region, seqs)]);
    SubLibrary::new(lib, None, HashMap::new(), max_dist, None).unwrap()
}

fn unique_library(size: usize, len: usize, max_dist: u64) -> (SubLibrary, Vec<u8>) {
    let seqs = (0..size).map(|i| make_base_seq(i, len)).collect::<Vec<_>>();
    let base = seqs[seqs.len() / 2].clone();
    (build_sublib_from_seqs(seqs, max_dist), base)
}

fn run_fixture(
    group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
    fixture: &LookupFixture,
    metric: DistanceMetric,
    partial: PartialMatching,
    function: &str,
    name: &str,
) {
    let query = seq_from_bytes(&fixture.query);

    group.bench_with_input(BenchmarkId::new(function, name), fixture, |b, f| {
        b.iter(|| {
            let out = f.lib.lookup(
                black_box(&f.region),
                black_box(&query),
                black_box(metric),
                black_box(partial),
            );
            let _ = black_box(out);
        })
    });
}

fn bench_lookup_metric(c: &mut Criterion) {
    let mut group = c.benchmark_group("library_metric");

    let (lib, base) = unique_library(64, 32, 2);

    let exact_hit = LookupFixture::new(&lib, &base, "exact");
    let sub_query = LookupFixture::new(&lib, &apply_n_subs(&base, 2), "subs");
    let indel_query = LookupFixture::new(&lib, &apply_mixed_edit(&base), "indels");

    let matrix = [
        (&exact_hit, DistanceMetric::Exact, "exact"),
        (&exact_hit, DistanceMetric::Hamming, "hamming"),
        (&sub_query, DistanceMetric::Hamming, "hamming"),
        (
            &exact_hit,
            DistanceMetric::BoundedLevenshtein,
            "bounded_levenshtein",
        ),
        (
            &sub_query,
            DistanceMetric::BoundedLevenshtein,
            "bounded_levenshtein",
        ),
        (
            &indel_query,
            DistanceMetric::BoundedLevenshtein,
            "bounded_levenshtein",
        ),
        (&exact_hit, DistanceMetric::Levenshtein, "levenshtein"),
        (&sub_query, DistanceMetric::Levenshtein, "levenshtein"),
        (&indel_query, DistanceMetric::Levenshtein, "levenshtein"),
    ];

    for (fix, metric, name) in matrix {
        run_fixture(
            &mut group,
            fix,
            metric,
            PartialMatching::Full,
            name,
            &fix.name,
        );
    }

    group.finish();
}

fn bench_lookup_lib_size(c: &mut Criterion) {
    let mut group = c.benchmark_group("library_size");

    for size in &[32usize, 128, 512, 1024] {
        group.throughput(Throughput::Elements(*size as u64));

        let (lib, base) = unique_library(*size, 32, 2);
        let fix = LookupFixture::new(&lib, &apply_n_subs(&base, 2), "subs");

        run_fixture(
            &mut group,
            &fix,
            DistanceMetric::Hamming,
            PartialMatching::Full,
            "lookup",
            &format!("{}", size),
        );
    }

    group.finish();
}

fn bench_lookup_query_size(c: &mut Criterion) {
    let mut group = c.benchmark_group("query_size");

    for size in &[16usize, 32, 128] {
        group.throughput(Throughput::Elements(*size as u64));

        let (lib, base) = unique_library(128, *size, 2);
        let fix = LookupFixture::new(&lib, &apply_n_subs(&base, 2), "subs");

        run_fixture(
            &mut group,
            &fix,
            DistanceMetric::Hamming,
            PartialMatching::Full,
            "hamming",
            &format!("{}", size),
        );

        run_fixture(
            &mut group,
            &fix,
            DistanceMetric::BoundedLevenshtein,
            PartialMatching::Full,
            "bounded_levenshtein",
            &format!("{}", size),
        );
    }

    group.finish();
}

fn bench_lookup_partial(c: &mut Criterion) {
    let mut group = c.benchmark_group("partial");

    let (lib, base) = unique_library(64, 32, 2);

    let mut query = apply_n_subs(&base, 2);
    query.truncate(24);

    let sub_query = LookupFixture::new(&lib, &query, "subs");

    let matrix = [
        (
            DistanceMetric::Hamming,
            PartialMatching::Full,
            "hamming/full",
        ),
        (
            DistanceMetric::Hamming,
            PartialMatching::FivePrimeOnly,
            "hamming/partial3",
        ),
        (
            DistanceMetric::Hamming,
            PartialMatching::ThreePrimeOnly,
            "hamming/partial5",
        ),
        (
            DistanceMetric::BoundedLevenshtein,
            PartialMatching::Full,
            "bounded_levenshtein/full",
        ),
        (
            DistanceMetric::BoundedLevenshtein,
            PartialMatching::FivePrimeOnly,
            "bounded_levenshtein/partial3",
        ),
        (
            DistanceMetric::BoundedLevenshtein,
            PartialMatching::ThreePrimeOnly,
            "bounded_levenshtein/partial5",
        ),
        (
            DistanceMetric::Levenshtein,
            PartialMatching::Full,
            "levenshtein/full",
        ),
        (
            DistanceMetric::Levenshtein,
            PartialMatching::FivePrimeOnly,
            "levenshtein/partial3",
        ),
        (
            DistanceMetric::Levenshtein,
            PartialMatching::ThreePrimeOnly,
            "levenshtein/partial5",
        ),
    ];

    for (metric, partial, name) in matrix {
        run_fixture(&mut group, &sub_query, metric, partial, "lookup", name);
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_lookup_metric,
    bench_lookup_lib_size,
    bench_lookup_query_size,
    bench_lookup_partial
);
criterion_main!(benches);
