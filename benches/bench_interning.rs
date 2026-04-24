//! Benchmark the interner

#[cfg(feature = "interning")]
use criterion::BatchSize;
use criterion::{BenchmarkId, Criterion, Throughput, black_box, criterion_group, criterion_main};
use dnacomb::groups::ReadGroup;
use dnacomb::interning::{
    GroupID, LibraryID, RegionID, group_id_from_str, group_id_to_str, library_id_from_str,
    library_id_to_str, region_id_from_str, region_id_to_str,
};

#[cfg(feature = "interning")]
use dnacomb::interning::SeqInterner;

fn make_unique_seqs(n: usize, len: usize, seed: usize) -> Vec<Vec<u8>> {
    let alphabet = [b'A', b'C', b'G', b'T'];
    (0..n)
        .map(|i| {
            let mut v = Vec::with_capacity(len);
            let mut x = i ^ seed;
            for _ in 0..len {
                v.push(alphabet[x & 0b11]);
                x = x.rotate_left(3) ^ 0x9E37_79B9usize;
            }
            v
        })
        .collect()
}

fn make_duplicate_seqs(n: usize, len: usize, distinct: usize, seed: usize) -> Vec<Vec<u8>> {
    let seed_set = make_unique_seqs(distinct.max(1), len, seed);
    (0..n)
        .map(|i| seed_set[i % seed_set.len()].clone())
        .collect()
}

fn make_unique_names(prefix: &str, n: usize) -> Vec<String> {
    (0..n).map(|i| format!("{prefix}_{i}")).collect()
}

fn make_duplicate_names(prefix: &str, n: usize, distinct: usize) -> Vec<String> {
    let names = make_unique_names(prefix, distinct.max(1));
    (0..n).map(|i| names[i % names.len()].clone()).collect()
}

#[cfg(feature = "interning")]
fn bench_seq_interning_isolated(c: &mut Criterion) {
    let mut group = c.benchmark_group("interning");

    for &len in &[20, 100] {
        let unique = make_unique_seqs(10000, len, 0x1111);
        let duplicate_heavy = make_duplicate_seqs(10000, len, 16, 0x2222);

        group.throughput(Throughput::Elements(10000 as u64));

        group.bench_with_input(
            BenchmarkId::new("seq_intern", format!("unique_{len}bp")),
            &unique,
            |b, seqs| {
                b.iter_batched(
                    || SeqInterner::default(),
                    |interner| {
                        for s in seqs {
                            let h = interner.intern(black_box(s));
                            black_box(h);
                        }
                        black_box(interner.num_interned_forward());
                        black_box(interner.num_interned_reverse());
                    },
                    BatchSize::SmallInput,
                )
            },
        );

        group.bench_with_input(
            BenchmarkId::new("seq_intern", format!("duplicate_{len}bp")),
            &duplicate_heavy,
            |b, seqs| {
                b.iter(|| {
                    let interner = SeqInterner::default();
                    for s in seqs {
                        let h = interner.intern(black_box(s));
                        black_box(h);
                    }
                    black_box(interner.num_interned_forward());
                    black_box(interner.num_interned_reverse());
                })
            },
        );
    }

    group.finish();
}

#[cfg(feature = "interning")]
fn bench_seq_resolve_warm_local(c: &mut Criterion) {
    let mut group = c.benchmark_group("interning");

    let seqs = make_duplicate_seqs(10_000, 40, 64, 0x3333);
    let interner = SeqInterner::default();
    let handles: Vec<_> = seqs.iter().map(|s| interner.intern(s)).collect();

    group.throughput(Throughput::Elements(handles.len() as u64));

    group.bench_function(
        BenchmarkId::new("seq_resolve", format!("duplicate_40bp")),
        |b| {
            b.iter(|| {
                for h in &handles {
                    let bytes = interner.resolve(black_box(h));
                    black_box(bytes);
                }
            })
        },
    );

    group.finish();
}

fn bench_group_string_interner_warm(c: &mut Criterion) {
    let mut group = c.benchmark_group("interning");

    let unique = make_unique_names("group_unique_only", 10000);
    let duplicate_heavy = make_duplicate_names("group_dup_only", 10000, 32);

    group.throughput(Throughput::Elements(10000 as u64));

    group.bench_with_input(
        BenchmarkId::new("group_intern", "unique"),
        &unique,
        |b, names| {
            b.iter(|| {
                for name in names {
                    let id = group_id_from_str(black_box(name));
                    black_box(id);
                }
            })
        },
    );

    group.bench_with_input(
        BenchmarkId::new("group_intern", "duplicates"),
        &duplicate_heavy,
        |b, names| {
            b.iter(|| {
                for name in names {
                    let id = group_id_from_str(black_box(name));
                    black_box(id);
                }
            })
        },
    );

    let unique_ids: Vec<GroupID> = unique.iter().map(|x| group_id_from_str(x)).collect();
    let duplicate_ids: Vec<GroupID> = duplicate_heavy
        .iter()
        .map(|x| group_id_from_str(x))
        .collect();

    group.bench_with_input(
        BenchmarkId::new("group_resolve", "unique"),
        &unique_ids,
        |b, ids| {
            b.iter(|| {
                for id in ids {
                    let id = group_id_to_str(black_box(id));
                    black_box(id);
                }
            })
        },
    );

    group.bench_with_input(
        BenchmarkId::new("group_resolve", "duplicates"),
        &duplicate_ids,
        |b, ids| {
            b.iter(|| {
                for id in ids {
                    let id = group_id_to_str(black_box(id));
                    black_box(id);
                }
            })
        },
    );

    group.finish();
}

fn bench_region_string_interner_warm(c: &mut Criterion) {
    let mut group = c.benchmark_group("interning");

    let unique = make_unique_names("region_unique_only", 100);
    let duplicate_heavy = make_duplicate_names("region_dup_only", 100, 32);

    group.throughput(Throughput::Elements(100 as u64));

    group.bench_with_input(
        BenchmarkId::new("region_intern", "unique"),
        &unique,
        |b, names| {
            b.iter(|| {
                for name in names {
                    let id = region_id_from_str(black_box(name));
                    black_box(id);
                }
            })
        },
    );

    group.bench_with_input(
        BenchmarkId::new("region_intern", "duplicates"),
        &duplicate_heavy,
        |b, names| {
            b.iter(|| {
                for name in names {
                    let id = region_id_from_str(black_box(name));
                    black_box(id);
                }
            })
        },
    );

    let unique_ids: Vec<RegionID> = unique.iter().map(|x| region_id_from_str(x)).collect();
    let duplicate_ids: Vec<RegionID> = duplicate_heavy
        .iter()
        .map(|x| region_id_from_str(x))
        .collect();

    group.bench_with_input(
        BenchmarkId::new("region_resolve", "unique"),
        &unique_ids,
        |b, ids| {
            b.iter(|| {
                for id in ids {
                    let id = region_id_to_str(black_box(id));
                    black_box(id);
                }
            })
        },
    );

    group.bench_with_input(
        BenchmarkId::new("region_resolve", "duplicates"),
        &duplicate_ids,
        |b, ids| {
            b.iter(|| {
                for id in ids {
                    let id = region_id_to_str(black_box(id));
                    black_box(id);
                }
            })
        },
    );

    group.finish();
}

fn bench_library_string_interner_warm(c: &mut Criterion) {
    let mut group = c.benchmark_group("interning");

    let unique = make_unique_names("library_unique_only", 1000);
    let duplicate_heavy = make_duplicate_names("library_dup_only", 1000, 32);

    group.throughput(Throughput::Elements(1000 as u64));

    group.bench_with_input(
        BenchmarkId::new("library_intern", "unique"),
        &unique,
        |b, names| {
            b.iter(|| {
                for name in names {
                    let id = library_id_from_str(black_box(name));
                    black_box(id);
                }
            })
        },
    );

    group.bench_with_input(
        BenchmarkId::new("library_intern", "duplicates"),
        &duplicate_heavy,
        |b, names| {
            b.iter(|| {
                for name in names {
                    let id = library_id_from_str(black_box(name));
                    black_box(id);
                }
            })
        },
    );

    let unique_ids: Vec<LibraryID> = unique.iter().map(|x| library_id_from_str(x)).collect();
    let duplicate_ids: Vec<LibraryID> = duplicate_heavy
        .iter()
        .map(|x| library_id_from_str(x))
        .collect();

    group.bench_with_input(
        BenchmarkId::new("library_resolve", "unique"),
        &unique_ids,
        |b, ids| {
            b.iter(|| {
                for id in ids {
                    let id = library_id_to_str(black_box(id));
                    black_box(id);
                }
            })
        },
    );

    group.bench_with_input(
        BenchmarkId::new("library_resolve", "duplicates"),
        &duplicate_ids,
        |b, ids| {
            b.iter(|| {
                for id in ids {
                    let id = library_id_to_str(black_box(id));
                    black_box(id);
                }
            })
        },
    );

    group.finish();
}

fn bench_readgroup_surface_warm(c: &mut Criterion) {
    let mut group = c.benchmark_group("interning");

    let unique_names = make_unique_names("readgroup_unique", 10000);
    let duplicate_names = make_duplicate_names("readgroup_duplicates", 10000, 128);

    group.throughput(Throughput::Elements(10_000));

    group.bench_function(BenchmarkId::new("readgroup", "unique"), |b| {
        b.iter(|| {
            for name in &unique_names {
                let g = ReadGroup::grouped(black_box(name));
                black_box(g);
            }
        })
    });

    group.bench_function(BenchmarkId::new("readgroup", "duplcates"), |b| {
        b.iter(|| {
            for name in &duplicate_names {
                let g = ReadGroup::grouped(black_box(name));
                black_box(g);
            }
        })
    });

    group.finish();
}

#[cfg(feature = "interning")]
criterion_group!(
    benches,
    bench_seq_interning_isolated,
    bench_seq_resolve_warm_local,
    bench_group_string_interner_warm,
    bench_region_string_interner_warm,
    bench_library_string_interner_warm,
    bench_readgroup_surface_warm
);

#[cfg(not(feature = "interning"))]
criterion_group!(
    benches,
    bench_group_string_interner_warm,
    bench_region_string_interner_warm,
    bench_library_string_interner_warm,
    bench_readgroup_surface_warm
);

criterion_main!(benches);
