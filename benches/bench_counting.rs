//! Benchmark counting
use criterion::{
    BatchSize, BenchmarkId, Criterion, Throughput, black_box, criterion_group, criterion_main,
};
use dnacomb::counting::AlignmentScorer;
use dnacomb::filters::{AlignmentTolerance, FilterConfig};
use dnacomb::parsing::ReadPairProducer;
use dnacomb::{
    Compression, CountMode, LibrarySpec, ReadPairParser, SeqFormat, SeqPath, count_reads,
};
use std::path::{Path, PathBuf};

const PEGRNA: (&str, &str) = ("benches/pegrna.fq", "benches/pegrna.json");

const GRNA: (&str, &str) = ("benches/grna.fq", "benches/grna.json");

#[derive(Clone, Copy, Debug)]
struct FixtureCase {
    name: &'static str,
    reads: &'static str,
    libspec: &'static str,
    allow_inframe: bool,
}

#[derive(Clone, Copy, Debug)]
struct CountCase {
    fixture: FixtureCase,
    mode: CountMode,
}

fn repo_path(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

fn alignment_scorer() -> AlignmentScorer {
    AlignmentScorer::new(6, -2, -3, -10, -4)
}

fn load_libspec(path: &str) -> LibrarySpec {
    LibrarySpec::from_file(repo_path(path).to_str().unwrap(), None, None, None, None).unwrap()
}

fn make_parser(path: &str) -> ReadPairParser {
    let forward = SeqPath::new(
        repo_path(path).to_string_lossy().to_string(),
        SeqFormat::Auto,
        Compression::Auto,
    );

    ReadPairParser::from_paths(forward, None, None, 0, b'I').unwrap()
}

fn make_filter_config(
    spec: &LibrarySpec,
    mode: CountMode,
    scorer: AlignmentScorer,
    reader_has_reverse: bool,
) -> FilterConfig {
    let alignment_tolerance = if matches!(mode, CountMode::Align) {
        let rev = if reader_has_reverse {
            Some(spec.expected_reverse_read())
        } else {
            None
        };

        Some(
            AlignmentTolerance::from_expected_reads(
                &spec.expected_forward_read(),
                rev.as_ref(),
                &spec.template_sequence(),
                &scorer,
                0.8,
                false,
            )
            .unwrap(),
        )
    } else {
        None
    };

    FilterConfig::new(None, alignment_tolerance, None, None, true)
}

fn case_name(case: CountCase) -> String {
    let mode = match case.mode {
        CountMode::FullRead => "full_read",
        CountMode::Inframe => "inframe",
        CountMode::Pattern => "pattern",
        CountMode::Align => "align",
    };

    format!("{}/{}", case.fixture.name, mode,)
}

fn valid_case(case: CountCase) -> bool {
    if matches!(case.mode, CountMode::Inframe) && !case.fixture.allow_inframe {
        return false;
    }

    true
}

fn run_count_case(case: CountCase) {
    let reader = make_parser(case.fixture.reads);
    let spec = load_libspec(case.fixture.libspec);
    let scorer = alignment_scorer();
    let filter_config = make_filter_config(&spec, case.mode, scorer, reader.has_reverse());

    let counts = count_reads(
        reader,
        &Some(spec.clone()),
        case.mode,
        false,
        filter_config,
        Some(scorer),
        Some(12),
        Some(1),
        false,
        1,
        None,
    )
    .unwrap();

    black_box(counts.len());
    black_box(counts.total_filtered());
}

fn build_cases() -> Vec<CountCase> {
    let fixtures = [
        FixtureCase {
            name: "pegrna",
            reads: PEGRNA.0,
            libspec: PEGRNA.1,
            allow_inframe: false,
        },
        FixtureCase {
            name: "grna",
            reads: GRNA.0,
            libspec: GRNA.1,
            allow_inframe: true,
        },
    ];

    let modes = [
        CountMode::FullRead,
        CountMode::Inframe,
        CountMode::Pattern,
        CountMode::Align,
    ];

    let mut out = Vec::new();

    for fixture in fixtures {
        for mode in modes {
            let case = CountCase { fixture, mode };

            if valid_case(case) {
                out.push(case);
            }
        }
    }

    out
}

fn bench_counting_modes(c: &mut Criterion) {
    let mut group = c.benchmark_group("counting_modes");

    group.sample_size(100);

    for case in build_cases() {
        group.throughput(Throughput::Elements(1000));

        let name = case_name(case);
        group.bench_with_input(BenchmarkId::new("count", name), &case, |b, &case| {
            b.iter_batched(|| (), |_| run_count_case(case), BatchSize::PerIteration)
        });
    }

    group.finish();
}

criterion_group!(benches, bench_counting_modes);
criterion_main!(benches);
