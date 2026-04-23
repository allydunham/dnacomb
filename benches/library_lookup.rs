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

fn neighbor_library_subs(len: usize, n_neighbors: usize, max_dist: u64) -> (SubLibrary, Vec<u8>) {
    let mut seqs = (0..(n_neighbors + 9))
        .map(|i| make_base_seq(i, len))
        .collect::<Vec<_>>();

    let mid = seqs.len() / 2;
    let base = seqs[mid].clone();

    for i in 0..n_neighbors {
        let pos = (i * 5 + len / 4) % len;
        let neighbor = mutate_sub(base.clone(), pos);

        let insert_at = (mid + i + 1).min(seqs.len());
        seqs.insert(insert_at, neighbor);
    }

    (build_sublib_from_seqs(seqs, max_dist), base)
}

fn neighbor_library_subs_and_indels(
    len: usize,
    n_sub_neighbors: usize,
    n_indel_neighbors: usize,
    max_dist: u64,
) -> (SubLibrary, Vec<u8>) {
    let mut seqs = (0..(n_sub_neighbors + n_indel_neighbors + 9))
        .map(|i| make_base_seq(i, len))
        .collect::<Vec<_>>();

    let mid = seqs.len() / 2;
    let base = seqs[mid].clone();

    for i in 0..n_sub_neighbors {
        let pos = (i * 5 + len / 4) % len;
        let neighbor = mutate_sub(base.clone(), pos);
        let insert_at = (mid + i + 1).min(seqs.len());
        seqs.insert(insert_at, neighbor);
    }

    for i in 0..n_indel_neighbors {
        let pos = usize::min((i * 3 + len / 3) % len, len.saturating_sub(1));
        let neighbor = if i % 2 == 0 {
            insert_base(&base, pos, b'T')
        } else if len > 1 {
            delete_base(&base, pos)
        } else {
            insert_base(&base, 0, b'T')
        };

        let insert_at = (mid + n_sub_neighbors + i + 1).min(seqs.len());
        seqs.insert(insert_at, neighbor);
    }

    (build_sublib_from_seqs(seqs, max_dist), base)
}

fn prefix_tie_library(
    size: usize,
    len: usize,
    shared_prefix_len: usize,
    max_dist: u64,
) -> (SubLibrary, Vec<u8>) {
    let prefix = vec![b'A'; shared_prefix_len];
    let mut seqs = Vec::with_capacity(size);
    for i in 0..size {
        let mut suffix = make_base_seq(i + 100, len - shared_prefix_len);
        let mut seq = prefix.clone();
        seq.append(&mut suffix);
        seqs.push(seq);
    }
    let query = vec![b'A'; shared_prefix_len];
    (build_sublib_from_seqs(seqs, max_dist), query)
}

fn suffix_tie_library(
    size: usize,
    len: usize,
    shared_suffix_len: usize,
    max_dist: u64,
) -> (SubLibrary, Vec<u8>) {
    let suffix = vec![b'T'; shared_suffix_len];
    let mut seqs = Vec::with_capacity(size);
    for i in 0..size {
        let mut prefix = make_base_seq(i + 200, len - shared_suffix_len);
        prefix.extend_from_slice(&suffix);
        seqs.push(prefix);
    }
    let query = vec![b'T'; shared_suffix_len];
    (build_sublib_from_seqs(seqs, max_dist), query)
}

fn fixture_unique_exact_hit(size: usize, len: usize, max_dist: u64) -> LookupFixture {
    let (lib, base) = unique_library(size, len, max_dist);
    LookupFixture {
        lib,
        query: base,
        region: region_id_from_str("r1"),
        name: format!("unique_exact_hit_lib{size}_len{len}_k{max_dist}"),
    }
}

fn fixture_unique_sub_query(
    size: usize,
    len: usize,
    n_subs: usize,
    max_dist: u64,
) -> LookupFixture {
    let (lib, base) = unique_library(size, len, max_dist);
    LookupFixture {
        lib,
        query: apply_n_subs(&base, n_subs),
        region: region_id_from_str("r1"),
        name: format!("unique_sub{n_subs}_query_lib{size}_len{len}_k{max_dist}"),
    }
}

fn fixture_unique_indel_query(size: usize, len: usize, max_dist: u64) -> LookupFixture {
    let (lib, base) = unique_library(size, len, max_dist);
    LookupFixture {
        lib,
        query: insert_base(&base, len / 2, b'G'),
        region: region_id_from_str("r1"),
        name: format!("unique_indel_query_lib{size}_len{len}_k{max_dist}"),
    }
}

fn fixture_sub_neighbor_tie(len: usize, n_neighbors: usize, max_dist: u64) -> LookupFixture {
    let (lib, base) = neighbor_library_subs(len, n_neighbors, max_dist);
    LookupFixture {
        lib,
        query: base,
        region: region_id_from_str("r1"),
        name: format!("sub_neighbor_tie_n{n_neighbors}_len{len}_k{max_dist}"),
    }
}

fn fixture_subs_and_indels_query(
    len: usize,
    n_sub_neighbors: usize,
    n_indel_neighbors: usize,
    max_dist: u64,
) -> LookupFixture {
    let (lib, base) =
        neighbor_library_subs_and_indels(len, n_sub_neighbors, n_indel_neighbors, max_dist);
    LookupFixture {
        lib,
        query: apply_n_subs(&base, 1),
        region: region_id_from_str("r1"),
        name: format!(
            "subs_and_indels_query_subn{n_sub_neighbors}_indeln{n_indel_neighbors}_len{len}_k{max_dist}"
        ),
    }
}

fn fixture_prefix_tie(size: usize, len: usize, prefix_len: usize, max_dist: u64) -> LookupFixture {
    let (lib, query) = prefix_tie_library(size, len, prefix_len, max_dist);
    LookupFixture {
        lib,
        query,
        region: region_id_from_str("r1"),
        name: format!("prefix_tie_lib{size}_len{len}_prefix{prefix_len}_k{max_dist}"),
    }
}

fn fixture_suffix_tie(size: usize, len: usize, suffix_len: usize, max_dist: u64) -> LookupFixture {
    let (lib, query) = suffix_tie_library(size, len, suffix_len, max_dist);
    LookupFixture {
        lib,
        query,
        region: region_id_from_str("r1"),
        name: format!("suffix_tie_lib{size}_len{len}_suffix{suffix_len}_k{max_dist}"),
    }
}

fn run_fixture(
    group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
    fixture: &LookupFixture,
    metric: DistanceMetric,
    partial: PartialMatching,
    label: &str,
) {
    let query = seq_from_bytes(&fixture.query);

    group.bench_with_input(BenchmarkId::new(label, &fixture.name), fixture, |b, f| {
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

fn bench_lookup_matrix(c: &mut Criterion) {
    let mut group = c.benchmark_group("library_lookup_matrix");

    group.sample_size(30);

    for &size in &[16usize, 64usize, 256usize] {
        group.throughput(Throughput::Elements(size as u64));

        let exact_hit = fixture_unique_exact_hit(size, 24, 2);
        let sub_query = fixture_unique_sub_query(size, 24, 1, 2);
        let sub_query_far = fixture_unique_sub_query(size, 24, 3, 1);
        let indel_query = fixture_unique_indel_query(size, 24, 2);

        run_fixture(
            &mut group,
            &exact_hit,
            DistanceMetric::Exact,
            PartialMatching::Full,
            "exact_full",
        );
        run_fixture(
            &mut group,
            &sub_query,
            DistanceMetric::Exact,
            PartialMatching::Full,
            "exact_full_sub_query",
        );
        run_fixture(
            &mut group,
            &indel_query,
            DistanceMetric::Exact,
            PartialMatching::Full,
            "exact_full_indel_query",
        );

        run_fixture(
            &mut group,
            &exact_hit,
            DistanceMetric::Hamming,
            PartialMatching::Full,
            "hamming_full_exact_hit",
        );
        run_fixture(
            &mut group,
            &sub_query,
            DistanceMetric::Hamming,
            PartialMatching::Full,
            "hamming_full_sub_query",
        );
        run_fixture(
            &mut group,
            &sub_query_far,
            DistanceMetric::Hamming,
            PartialMatching::Full,
            "hamming_full_sub_query_miss",
        );
        run_fixture(
            &mut group,
            &indel_query,
            DistanceMetric::Hamming,
            PartialMatching::Full,
            "hamming_full_indel_query",
        );

        run_fixture(
            &mut group,
            &exact_hit,
            DistanceMetric::Levenshtein,
            PartialMatching::Full,
            "lev_full_exact_hit",
        );
        run_fixture(
            &mut group,
            &sub_query,
            DistanceMetric::Levenshtein,
            PartialMatching::Full,
            "lev_full_sub_query",
        );
        run_fixture(
            &mut group,
            &indel_query,
            DistanceMetric::Levenshtein,
            PartialMatching::Full,
            "lev_full_indel_query",
        );

        run_fixture(
            &mut group,
            &exact_hit,
            DistanceMetric::BoundedLevenshtein,
            PartialMatching::Full,
            "bounded_lev_full_exact_hit",
        );
        run_fixture(
            &mut group,
            &sub_query,
            DistanceMetric::BoundedLevenshtein,
            PartialMatching::Full,
            "bounded_lev_full_sub_query",
        );
        run_fixture(
            &mut group,
            &indel_query,
            DistanceMetric::BoundedLevenshtein,
            PartialMatching::Full,
            "bounded_lev_full_indel_query",
        );
    }

    let sub_tie = fixture_sub_neighbor_tie(24, 4, 2);
    run_fixture(
        &mut group,
        &sub_tie,
        DistanceMetric::Hamming,
        PartialMatching::Full,
        "hamming_full_sub_tie",
    );
    run_fixture(
        &mut group,
        &sub_tie,
        DistanceMetric::Levenshtein,
        PartialMatching::Full,
        "lev_full_sub_tie",
    );
    run_fixture(
        &mut group,
        &sub_tie,
        DistanceMetric::BoundedLevenshtein,
        PartialMatching::Full,
        "bounded_lev_full_sub_tie",
    );

    let mixed_neighbors = fixture_subs_and_indels_query(24, 4, 4, 2);
    run_fixture(
        &mut group,
        &mixed_neighbors,
        DistanceMetric::Levenshtein,
        PartialMatching::Full,
        "lev_full_mixed_neighbors",
    );
    run_fixture(
        &mut group,
        &mixed_neighbors,
        DistanceMetric::BoundedLevenshtein,
        PartialMatching::Full,
        "bounded_lev_full_mixed_neighbors",
    );

    for &size in &[16usize, 64usize, 256usize] {
        let prefix_tie = fixture_prefix_tie(size, 24, 12, 2);
        run_fixture(
            &mut group,
            &prefix_tie,
            DistanceMetric::Exact,
            PartialMatching::FivePrimeOnly,
            "exact_5p_tie",
        );
        run_fixture(
            &mut group,
            &prefix_tie,
            DistanceMetric::Hamming,
            PartialMatching::FivePrimeOnly,
            "hamming_5p_tie",
        );
        run_fixture(
            &mut group,
            &prefix_tie,
            DistanceMetric::Levenshtein,
            PartialMatching::FivePrimeOnly,
            "lev_5p_tie",
        );
        run_fixture(
            &mut group,
            &prefix_tie,
            DistanceMetric::BoundedLevenshtein,
            PartialMatching::FivePrimeOnly,
            "bounded_lev_5p_tie",
        );

        let suffix_tie = fixture_suffix_tie(size, 24, 12, 2);
        run_fixture(
            &mut group,
            &suffix_tie,
            DistanceMetric::Exact,
            PartialMatching::ThreePrimeOnly,
            "exact_3p_tie",
        );
        run_fixture(
            &mut group,
            &suffix_tie,
            DistanceMetric::Hamming,
            PartialMatching::ThreePrimeOnly,
            "hamming_3p_tie",
        );
        run_fixture(
            &mut group,
            &suffix_tie,
            DistanceMetric::Levenshtein,
            PartialMatching::ThreePrimeOnly,
            "lev_3p_tie",
        );
        run_fixture(
            &mut group,
            &suffix_tie,
            DistanceMetric::BoundedLevenshtein,
            PartialMatching::ThreePrimeOnly,
            "bounded_lev_3p_tie",
        );
    }

    group.finish();
}

criterion_group!(benches, bench_lookup_matrix);
criterion_main!(benches);
