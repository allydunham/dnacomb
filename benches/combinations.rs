//! Benchmark combination performance
use bio::io::fastq;
use criterion::{
    BatchSize, BenchmarkId, Criterion, Throughput, black_box, criterion_group, criterion_main,
};
use dnacomb::combination::CombinationKey;
use dnacomb::combinations::{CacheHit, ObservedCombinations};
use dnacomb::filters::FilterReason;
use dnacomb::groups::ReadGroup;
use dnacomb::interning::{region_id_from_str, seq_from_bytes};
use dnacomb::library::{DistanceMetric, Library, SubLibrary};
use dnacomb::region::{RegionCompleteness, RegionKey};
use dnacomb::seqs::ReadPair;
use std::collections::HashMap;

mod bench_support;
use bench_support::{apply_n_subs, make_base_seq};

use crate::bench_support::no_filter;

fn make_counts(region_ids: &[&str]) -> ObservedCombinations {
    ObservedCombinations::new(
        region_ids.iter().map(|x| region_id_from_str(x)).collect(),
        no_filter(),
    )
}

fn reg(id: &str, seq: &[u8], completeness: RegionCompleteness) -> RegionKey {
    RegionKey::new(region_id_from_str(id), seq_from_bytes(seq), completeness)
}

fn comb_key(regs: Vec<RegionKey>) -> CombinationKey {
    CombinationKey::new(None, regs)
}

fn readpair(fwd: &[u8], rev: Option<&[u8]>, group: ReadGroup) -> ReadPair {
    let f_qual = vec![b'I'; fwd.len()];
    let forward = fastq::Record::with_attrs("f", None, fwd, &f_qual);
    let reverse = rev.map(|r| {
        let r_qual = vec![b'I'; r.len()];
        fastq::Record::with_attrs("r", None, r, &r_qual)
    });

    ReadPair {
        forward,
        reverse,
        group,
    }
}

fn make_repeated_key() -> CombinationKey {
    comb_key(vec![
        reg("r1", b"AAAA", RegionCompleteness::Complete),
        reg("r2", b"CCCC", RegionCompleteness::Complete),
    ])
}

fn make_unique_keys(n: usize) -> Vec<CombinationKey> {
    (0..n)
        .map(|i| {
            let s1 = make_base_seq(i, 12);
            let s2 = make_base_seq(i + 100_000, 12);
            comb_key(vec![
                reg("r1", &s1, RegionCompleteness::Complete),
                reg("r2", &s2, RegionCompleteness::Complete),
            ])
        })
        .collect()
}

fn make_shared_region_keys(n: usize) -> Vec<CombinationKey> {
    let shared = make_base_seq(42, 12);
    (0..n)
        .map(|i| {
            let s2 = make_base_seq(i + 500_000, 12);
            comb_key(vec![
                reg("r1", &shared, RegionCompleteness::Complete),
                reg("r2", &s2, RegionCompleteness::Complete),
            ])
        })
        .collect()
}

fn make_grouped_repeated_keys(n: usize, n_groups: usize) -> Vec<(CombinationKey, ReadGroup)> {
    let key = make_repeated_key();
    (0..n)
        .map(|i| {
            (
                key.clone(),
                ReadGroup::grouped(&format!("g{}", i % n_groups)),
            )
        })
        .collect()
}

fn make_shifted_unique_keys(n: usize) -> Vec<CombinationKey> {
    (0..n)
        .map(|i| {
            let s1 = make_base_seq(i + 9_000_000, 12);
            let s2 = make_base_seq(i + 10_000_000, 12);
            comb_key(vec![
                reg("r1", &s1, RegionCompleteness::Complete),
                reg("r2", &s2, RegionCompleteness::Complete),
            ])
        })
        .collect()
}

fn make_library() -> Library {
    let mut map: HashMap<_, Vec<Vec<u8>>> = HashMap::new();
    map.insert(
        region_id_from_str("r1"),
        vec![
            b"AAAA".to_vec(),
            b"AAAT".to_vec(),
            b"GGGG".to_vec(),
            b"TTTT".to_vec(),
        ],
    );
    map.insert(
        region_id_from_str("r2"),
        vec![
            b"CCCC".to_vec(),
            b"CCCC".to_vec(),
            b"TTTT".to_vec(),
            b"GGGG".to_vec(),
        ],
    );

    let ids = Some(vec![
        "seq1".to_string(),
        "seq2".to_string(),
        "seq3".to_string(),
        "seq4".to_string(),
    ]);

    let sub = SubLibrary::new(map, ids, HashMap::new(), 2, None).unwrap();
    Library::new(vec![sub]).unwrap()
}

fn populate_counts_from_keys(keys: &[CombinationKey]) -> ObservedCombinations {
    let mut counts = make_counts(&["r1", "r2"]);
    for key in keys {
        counts
            .add_or_increment_combination(key, ReadGroup::ungrouped())
            .unwrap();
    }
    counts
}

fn populate_compare_counts(n: usize) -> ObservedCombinations {
    let mut counts = make_counts(&["r1", "r2"]);

    for i in 0..n {
        let s1 = if i % 4 == 0 {
            b"AAAA".to_vec()
        } else {
            apply_n_subs(b"AAAA", 1)
        };

        let s2 = if i % 3 == 0 {
            b"CCCC".to_vec()
        } else {
            b"TTTT".to_vec()
        };

        let key = comb_key(vec![
            reg("r1", &s1, RegionCompleteness::Complete),
            reg("r2", &s2, RegionCompleteness::Complete),
        ]);

        counts
            .add_or_increment_combination(&key, ReadGroup::ungrouped())
            .unwrap();
    }

    counts
}

fn bench_add_or_increment(c: &mut Criterion) {
    let mut group = c.benchmark_group("combinations");

    group.throughput(Throughput::Elements(1000 as u64));

    let repeated = make_repeated_key();
    group.bench_with_input(
        BenchmarkId::new("add_or_increment", "repeated_same_key"),
        &repeated,
        |b, key| {
            b.iter_batched(
                || make_counts(&["r1", "r2"]),
                |mut counts| {
                    for _ in 0..1000 {
                        counts
                            .add_or_increment_combination(
                                black_box(key),
                                ReadGroup::ungrouped(),
                            )
                            .unwrap();
                    }
                    black_box(counts.len());
                },
                BatchSize::SmallInput,
            )
        }
    );

    let unique = make_unique_keys(1000);
    group.bench_with_input(
        BenchmarkId::new("add_or_increment", "all_unique_keys"),
        &unique,
        |b, keys| {
            b.iter_batched(
                || make_counts(&["r1", "r2"]),
                |mut counts| {
                    for key in keys {
                        counts
                            .add_or_increment_combination(
                                black_box(key),
                                ReadGroup::ungrouped(),
                            )
                            .unwrap();
                    }
                    black_box(counts.len());
                },
                BatchSize::SmallInput,
            )
        },
    );

    let shared = make_shared_region_keys(1000);
    group.bench_with_input(
        BenchmarkId::new("add_or_increment", "shared_region_keys"),
        &shared,
        |b, keys| {
            b.iter_batched(
                || make_counts(&["r1", "r2"]),
                |mut counts| {
                    for key in keys {
                        counts
                            .add_or_increment_combination(
                                black_box(key),
                                ReadGroup::ungrouped(),
                            )
                            .unwrap();
                    }
                    black_box(counts.len());
                },
                BatchSize::SmallInput,
            )
        },
    );

    let grouped = make_grouped_repeated_keys(1000, 32);
    group.bench_with_input(
        BenchmarkId::new("add_or_increment", "repeated_key_grouped_counts"),
        &grouped,
        |b, items| {
            b.iter_batched(
                || make_counts(&["r1", "r2"]),
                |mut counts| {
                    for (key, grp) in items {
                        counts
                            .add_or_increment_combination(
                                black_box(key),
                                black_box(grp.clone()),
                            )
                            .unwrap();
                    }
                    black_box(counts.len());
                },
                BatchSize::SmallInput,
            )
        },
    );

    group.finish();
}

fn bench_cache_replay(c: &mut Criterion) {
    let mut group = c.benchmark_group("combinations_cache");

    let key = make_repeated_key();
    let read = readpair(b"AAAA", None, ReadGroup::grouped("g1"));
    let filtered_read = readpair(b"TTTT", None, ReadGroup::grouped("g1"));

    group.bench_function("cache_combination_replay", |b| {
        b.iter_batched(
            || {
                let mut counts = make_counts(&["r1", "r2"]);
                counts.cache(read.key(), CacheHit::Comb(key.clone()));
                counts
            },
            |mut counts| {
                let out = counts.check_cache(black_box(&read), true).unwrap();
                black_box(out);
                black_box(counts.len());
            },
            BatchSize::SmallInput,
        )
    });

    group.bench_function("cache_filter_replay", |b| {
        b.iter_batched(
            || {
                let mut counts = make_counts(&["r1", "r2"]);
                counts.cache(
                    filtered_read.key(),
                    CacheHit::Filter(FilterReason::ShortRead),
                );
                counts
            },
            |mut counts| {
                let out = counts.check_cache(black_box(&filtered_read), true).unwrap();
                black_box(out);
                black_box(counts.total_filtered());
            },
            BatchSize::SmallInput,
        )
    });

    group.finish();
}

fn bench_merge(c: &mut Criterion) {
    let mut group = c.benchmark_group("combinations");

    group.throughput(Throughput::Elements(1000 as u64));

    let left_overlap = populate_counts_from_keys(&make_shared_region_keys(1000));
    let right_overlap = populate_counts_from_keys(&make_shared_region_keys(1000));

    group.bench_function(BenchmarkId::new("merge", "mostly_overlapping"), |b, | {
        b.iter_batched(
            || (left_overlap.clone(), right_overlap.clone()),
            |(mut a, b_counts)| {
                a.merge(black_box(b_counts)).unwrap();
                black_box(a.len());
            },
            BatchSize::SmallInput,
        )
    });

    let left_disjoint = populate_counts_from_keys(&make_unique_keys(1000));
    let right_disjoint = populate_counts_from_keys(&make_shifted_unique_keys(1000));

    group.bench_function(BenchmarkId::new("merge", "mostly_disjoint"), |b| {
        b.iter_batched(
            || (left_disjoint.clone(), right_disjoint.clone()),
            |(mut a, b_counts)| {
                a.merge(black_box(b_counts)).unwrap();
                black_box(a.len());
            },
            BatchSize::SmallInput,
        )
    });

    group.finish();
}

fn bench_compare_and_summary(c: &mut Criterion) {
    let mut group = c.benchmark_group("combinations");

    group.throughput(Throughput::Elements(1000 as u64));

    let lib = make_library();
    let compare_counts = populate_compare_counts(1000);

    group.bench_function(
        BenchmarkId::new("compare", "Hamming"),
        |b| {
            b.iter_batched(
                || compare_counts.clone(),
                |mut counts| {
                    counts
                        .compare_to_library(lib.clone(), None, DistanceMetric::Hamming, 10, 1)
                        .unwrap();
                    black_box(counts.len());
                },
                BatchSize::SmallInput,
            )
        },
    );

    let mut compared_counts = populate_compare_counts(1000);
    compared_counts
        .compare_to_library(lib.clone(), None, DistanceMetric::Hamming, 10, 1)
        .unwrap();

    group.bench_function(
        BenchmarkId::new("summarise", ""),
        |b| {
            b.iter(|| {
                let s = compared_counts.summarise();
                black_box(s.total());
            })
        },
    );

    let unique_counts = populate_counts_from_keys(&make_unique_keys(1000));
    group.bench_function(BenchmarkId::new("to_vector", "unsorted"), |b| {
        b.iter(|| {
            let v = unique_counts.to_vector(false);
            black_box(v.len());
        })
    });

    group.bench_function(BenchmarkId::new("to_vector", "sorted"), |b| {
        b.iter(|| {
            let v = unique_counts.to_vector(true);
            black_box(v.len());
        })
    });

    let shared_counts = populate_counts_from_keys(&make_shared_region_keys(1000));
    let mut shared_compared = shared_counts.clone();
    shared_compared
        .compare_to_library(lib.clone(), None, DistanceMetric::Hamming, 10, 1)
        .unwrap();

    group.bench_function(
        BenchmarkId::new("to_library_vector", "unsorted"),
        |b| {
            b.iter(|| {
                let v = shared_compared.to_library_vector(false).unwrap();
                black_box(v.len());
            })
        },
    );

    group.bench_function(
        BenchmarkId::new("to_library_vector", "sorted"),
        |b| {
            b.iter(|| {
                let v = shared_compared.to_library_vector(true).unwrap();
                black_box(v.len());
            })
        },
    );


    group.finish();
}

criterion_group!(
    benches,
    bench_add_or_increment,
    bench_cache_replay,
    bench_merge,
    bench_compare_and_summary
);
criterion_main!(benches);
