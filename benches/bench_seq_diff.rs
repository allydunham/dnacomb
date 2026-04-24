//! Benchmark SeqDiff
use criterion::{BenchmarkId, Criterion, Throughput, black_box, criterion_group, criterion_main};
use dnacomb::seq_diff::SequenceDiff;

mod bench_support;
use bench_support::*;

#[derive(Clone, Copy, Debug)]
enum DiffCaseKind {
    Exact,
    SingleSub,
    MultiSubSparse,
    MultiSubClustered,
    SingleInsertion,
    SingleDeletion,
    MixedEdit,
}

#[derive(Clone, Copy, Debug)]
struct DiffCase {
    len: usize,
    kind: DiffCaseKind,
}

fn make_case(case: DiffCase) -> (Vec<u8>, Vec<u8>, String) {
    let expected = make_base_seq(12345, case.len);

    let observed = match case.kind {
        DiffCaseKind::Exact => expected.clone(),
        DiffCaseKind::SingleSub => apply_n_subs(&expected, 1),
        DiffCaseKind::MultiSubSparse => apply_n_subs(&expected, 4),
        DiffCaseKind::MultiSubClustered => {
            let start = expected.len().saturating_div(2).saturating_sub(2);
            apply_clustered_subs(&expected, start, 4)
        }
        DiffCaseKind::SingleInsertion => insert_base(&expected, expected.len() / 2, b'A'),
        DiffCaseKind::SingleDeletion => {
            if expected.is_empty() {
                expected.clone()
            } else {
                delete_base(&expected, expected.len() / 2)
            }
        }
        DiffCaseKind::MixedEdit => apply_mixed_edit(&expected),
    };

    let label = match case.kind {
        DiffCaseKind::Exact => "exact",
        DiffCaseKind::SingleSub => "single_sub",
        DiffCaseKind::MultiSubSparse => "multi_sub_sparse",
        DiffCaseKind::MultiSubClustered => "multi_sub_clustered",
        DiffCaseKind::SingleInsertion => "single_insertion",
        DiffCaseKind::SingleDeletion => "single_deletion",
        DiffCaseKind::MixedEdit => "mixed_edit",
    };

    (observed, expected, format!("{label}_len{}", case.len))
}

fn bench_seq_diff_compute(c: &mut Criterion) {
    let mut group = c.benchmark_group("seq_diff");

    for &len in &[16usize, 64usize, 128usize] {
        group.throughput(Throughput::Elements(len as u64));

        let cases = [
            DiffCase {
                len,
                kind: DiffCaseKind::Exact,
            },
            DiffCase {
                len,
                kind: DiffCaseKind::SingleSub,
            },
            DiffCase {
                len,
                kind: DiffCaseKind::MultiSubSparse,
            },
            DiffCase {
                len,
                kind: DiffCaseKind::MultiSubClustered,
            },
            DiffCase {
                len,
                kind: DiffCaseKind::SingleInsertion,
            },
            DiffCase {
                len,
                kind: DiffCaseKind::SingleDeletion,
            },
            DiffCase {
                len,
                kind: DiffCaseKind::MixedEdit,
            },
        ];

        for case in cases {
            let (observed, expected, label) = make_case(case);

            group.bench_with_input(BenchmarkId::new("compute", label), &case, |b, _| {
                b.iter(|| {
                    let diff = SequenceDiff::compute(black_box(&observed), black_box(&expected));
                    black_box(diff);
                })
            });
        }
    }

    group.finish();
}

criterion_group!(benches, bench_seq_diff_compute);
criterion_main!(benches);
