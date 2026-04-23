//! Benchmark the interner

use criterion::{BenchmarkId, Criterion, Throughput, black_box, criterion_group, criterion_main};
use dnacomb::groups::ReadGroup;
use dnacomb::interning::{group_id_from_str, library_id_from_str, region_id_from_str};

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
    let mut group = c.benchmark_group("seq_interning_isolated");

    for &(n, len) in &[(1_000usize, 20usize), (10_000, 20), (10_000, 100)] {
        let unique = make_unique_seqs(n, len, 0x1111);
        let duplicate_heavy = make_duplicate_seqs(n, len, 16, 0x2222);

        group.throughput(Throughput::Elements(n as u64));

        group.bench_with_input(
            BenchmarkId::new("unique_insert_cold", format!("n{n}_len{len}")),
            &unique,
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

        group.bench_with_input(
            BenchmarkId::new("duplicate_insert_cold", format!("n{n}_len{len}")),
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

        group.bench_with_input(
            BenchmarkId::new("unique_insert_then_resolve_cold", format!("n{n}_len{len}")),
            &unique,
            |b, seqs| {
                b.iter(|| {
                    let interner = SeqInterner::default();
                    let handles: Vec<_> =
                        seqs.iter().map(|s| interner.intern(black_box(s))).collect();

                    for h in &handles {
                        let bytes = interner.resolve(black_box(h));
                        black_box(bytes);
                    }
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new(
                "duplicate_insert_then_resolve_cold",
                format!("n{n}_len{len}"),
            ),
            &duplicate_heavy,
            |b, seqs| {
                b.iter(|| {
                    let interner = SeqInterner::default();
                    let handles: Vec<_> =
                        seqs.iter().map(|s| interner.intern(black_box(s))).collect();

                    for h in &handles {
                        let bytes = interner.resolve(black_box(h));
                        black_box(bytes);
                    }
                })
            },
        );
    }

    group.finish();
}

#[cfg(feature = "interning")]
fn bench_seq_resolve_warm_local(c: &mut Criterion) {
    let mut group = c.benchmark_group("seq_resolve_local_warm");

    let seqs = make_duplicate_seqs(10_000, 40, 64, 0x3333);
    let interner = SeqInterner::default();
    let handles: Vec<_> = seqs.iter().map(|s| interner.intern(s)).collect();

    group.throughput(Throughput::Elements(handles.len() as u64));

    group.bench_function("resolve_duplicate_heavy_local_warm", |b| {
        b.iter(|| {
            for h in &handles {
                let bytes = interner.resolve(black_box(h));
                black_box(bytes);
            }
        })
    });

    group.finish();
}

fn bench_group_string_interner_warm(c: &mut Criterion) {
    let mut group = c.benchmark_group("group_string_interner_warm");

    for &n in &[1_000usize, 10_000usize] {
        // Distinct prefixes so these names do not collide with other bench families
        let unique = make_unique_names("bench_group_unique_only", n);
        let duplicate_heavy = make_duplicate_names("bench_group_dup_only", n, 32);

        group.throughput(Throughput::Elements(n as u64));

        group.bench_with_input(BenchmarkId::new("unique", n), &unique, |b, names| {
            b.iter(|| {
                for name in names {
                    let id = group_id_from_str(black_box(name));
                    black_box(id);
                }
            })
        });

        group.bench_with_input(
            BenchmarkId::new("duplicate_heavy", n),
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
    }

    group.finish();
}

fn bench_region_string_interner_warm(c: &mut Criterion) {
    let mut group = c.benchmark_group("region_string_interner_warm");

    for &n in &[1_000usize, 10_000usize] {
        let unique = make_unique_names("bench_region_unique_only", n);
        let duplicate_heavy = make_duplicate_names("bench_region_dup_only", n, 32);

        group.throughput(Throughput::Elements(n as u64));

        group.bench_with_input(BenchmarkId::new("unique", n), &unique, |b, names| {
            b.iter(|| {
                for name in names {
                    let id = region_id_from_str(black_box(name));
                    black_box(id);
                }
            })
        });

        group.bench_with_input(
            BenchmarkId::new("duplicate_heavy", n),
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
    }

    group.finish();
}

fn bench_library_string_interner_warm(c: &mut Criterion) {
    let mut group = c.benchmark_group("library_string_interner_warm");

    for &n in &[1_000usize, 10_000usize] {
        let unique = make_unique_names("bench_library_unique_only", n);
        let duplicate_heavy = make_duplicate_names("bench_library_dup_only", n, 32);

        group.throughput(Throughput::Elements(n as u64));

        group.bench_with_input(BenchmarkId::new("unique", n), &unique, |b, names| {
            b.iter(|| {
                for name in names {
                    let id = library_id_from_str(black_box(name));
                    black_box(id);
                }
            })
        });

        group.bench_with_input(
            BenchmarkId::new("duplicate_heavy", n),
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
    }

    group.finish();
}

fn bench_readgroup_surface_warm(c: &mut Criterion) {
    let mut group = c.benchmark_group("readgroup_surface_warm");

    let unique_names = make_unique_names("bench_readgroup_unique_only", 10_000);
    let duplicate_names = make_duplicate_names("bench_readgroup_surface_only", 10_000, 128);

    group.throughput(Throughput::Elements(10_000));

    group.bench_function("grouped_unique", |b| {
        b.iter(|| {
            for name in &unique_names {
                let g = ReadGroup::grouped(black_box(name));
                black_box(g);
            }
        })
    });

    group.bench_function("grouped_duplicate_heavy", |b| {
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
