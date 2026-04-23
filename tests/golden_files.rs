//! Integration test comparing observed and expected output
//! for a few small exact cases
use std::fs;
use std::fs::File;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use dnacomb::filters::AlignmentTolerance;
use dnacomb::library::Library;
use dnacomb::parsing::ReadPairProducer;
use dnacomb::{
    Compression, CountMode, DistanceMetric, FilterConfig, LibrarySpec, ProgressStyle,
    ReadPairParser, SeqFormat, SeqPath, count_reads, write_counts, write_filter_summary,
    write_library_counts, write_summary,
};
use regex::Regex;

/// One golden integration-test case.
///
/// Fixture directory convention:
/// - reads.fastq or reads.fasta
/// - reads_R1.fastq / reads_R2.fastq for paired-end
/// - libspec.json when structured counting is used
/// - library.tsv or multiple library TSVs when library comparison is used
/// - expected.counts.tsv
/// - expected.summary.tsv
/// - expected.filtered.tsv
/// - expected.library_counts.tsv (only when compare_to_library = true)
#[derive(Debug, Clone)]
struct GoldenCase {
    /// Relative to CARGO_MANIFEST_DIR, e.g. "tests/full_read"
    dir: &'static str,

    mode: CountMode,
    paired: bool,
    full_seq: bool,

    /// Optional grouping regex applied to forward read names
    group_regex: Option<&'static str>,

    /// Optional LibSpec loading
    use_libspec: bool,

    /// Read filtering
    mean_quality_threshold: Option<f32>,
    minimum_read_length: Option<usize>,
    maximum_read_length: Option<usize>,
    filter_empty: bool,
    alignment_tolerance: Option<f32>,

    /// Pattern-mode settings
    pattern_length: Option<usize>,
    pattern_tolerance: Option<u64>,

    /// Runtime
    cache: bool,
    threads: usize,

    /// Optional library comparison
    compare_to_library: bool,
    library_files: &'static [&'static str],
    distance_metric: DistanceMetric,
    max_distance: u64,
    max_matches: usize,
}

impl GoldenCase {
    fn dir(&self) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join(self.dir)
    }

    fn forward_input_path(&self) -> PathBuf {
        let dir = self.dir();
        for candidate in ["forward.fq", "forward.fa", "forward.fq.gz", "forward.fa.gz"] {
            let p = dir.join(candidate);
            if p.exists() {
                return p;
            }
        }
        panic!("No forward reads file found in {}", dir.display());
    }

    fn reverse_input_path(&self) -> Option<PathBuf> {
        if !self.paired {
            return None;
        }

        let dir = self.dir();
        for candidate in ["reverse.fq", "reverse.fa", "reverse.fq.gz", "reverse.fa.gz"] {
            let p = dir.join(candidate);
            if p.exists() {
                return Some(p);
            }
        }
        panic!(
            "paired=true but no reverse reads file found in {}",
            dir.display()
        );
    }

    fn libspec_path(&self) -> Option<PathBuf> {
        if self.use_libspec {
            let p = self.dir().join("libspec.json");
            assert!(p.exists(), "Missing libspec.json in {}", self.dir);
            Some(p)
        } else {
            None
        }
    }

    fn library_paths(&self) -> Vec<PathBuf> {
        self.library_files
            .iter()
            .map(|name| self.dir().join(name))
            .collect()
    }
}

fn read_tsv(path: PathBuf) -> Result<String> {
    fs::read_to_string(&path).with_context(|| format!("reading {}", path.display()))
}

/// Compare a produced file to `expected.<suffix>.tsv` inside the fixture directory.
/// If the expected file is missing, this fails with a helpful message so you can
/// inspect the produced tempdir and generate the golden file yourself.
fn assert_matches_expected(expected: &Path, observed: &Path) -> Result<()> {
    if !observed.exists() {
        anyhow::bail!("Missing observed file: {}", observed.display());
    }

    if !expected.exists() {
        anyhow::bail!(
            "Missing expected golden file: {}. Create it from {} and rerun.",
            expected.display(),
            observed.display(),
        );
    }

    let expected_str = read_tsv(expected.into())?;
    let observed_str = read_tsv(observed.into())?;

    if expected_str == observed_str {
        let _ = fs::remove_file(observed);
        return Ok(());
    }

    anyhow::bail!(
        "Golden file mismatch:\n\
         Expected: {}\n\
         Observed: {}\n\
         \
         Inspect the observed file and update the golden if the new output is correct.",
        expected.display(),
        observed.display(),
    );
}

/// Run a golden case end-to-end and write outputs to a temp directory.
///
/// This:
/// - loads optional LibSpec from file
/// - opens reads via SeqPath + ReadPairParser
/// - runs count_reads
/// - optionally compiles expected library TSV(s) and compares observed reads
/// - writes counts/summary/filtered/(library_counts) TSVs
fn run_case(case: &GoldenCase) -> Result<()> {
    let group = case
        .group_regex
        .map(|pat| Regex::new(pat))
        .transpose()
        .context("compiling grouping regex")?;

    let forward = SeqPath::new(
        case.forward_input_path().to_string_lossy().to_string(),
        SeqFormat::Auto,
        Compression::Auto,
    );

    let reverse = case.reverse_input_path().map(|p| {
        SeqPath::new(
            p.to_string_lossy().to_string(),
            SeqFormat::Auto,
            Compression::Auto,
        )
    });

    let reader = ReadPairParser::from_paths(
        forward, reverse, group, 0, // max_reads
        b'I',
    )
    .context("creating read parser")?;

    let lib_spec = match case.libspec_path() {
        Some(p) => Some(LibrarySpec::from_file(
            &p.to_string_lossy(),
            None,
            None,
            None,
            None,
        )?),
        None => None,
    };

    let alignment_scorer = dnacomb::AlignmentScorer::new(6, -2, -3, -10, -4);

    let alignment_tolerance = match (&lib_spec, case.alignment_tolerance) {
        (None, None) => None,
        (Some(_), None) => None,
        (None, Some(_)) => panic!("Alignment tolerance with no LibSpec"),
        (Some(l), Some(t)) => {
            let r = if reader.has_reverse() {
                Some(l.expected_reverse_read())
            } else {
                None
            };

            Some(
                AlignmentTolerance::from_expected_reads(
                    &l.expected_forward_read(),
                    r.as_ref(),
                    &l.template_sequence(),
                    &alignment_scorer,
                    t,
                    false,
                )
                .expect("AlignmentTolerance build failed"),
            )
        }
    };

    let filter_config = FilterConfig::new(
        case.mean_quality_threshold,
        alignment_tolerance,
        case.minimum_read_length,
        case.maximum_read_length,
        case.filter_empty,
    );

    let progress_style = ProgressStyle::new(None, false);

    let mut counts = count_reads(
        reader,
        &lib_spec,
        case.mode,
        case.full_seq,
        filter_config,
        Some(alignment_scorer),
        case.pattern_length,
        case.pattern_tolerance,
        case.cache,
        case.threads,
        Some(&progress_style),
    )
    .context("counting reads")?;

    if case.compare_to_library {
        let spec = lib_spec
            .as_ref()
            .context("compare_to_library=true requires use_libspec=true")?;

        let library_paths = case.library_paths();
        if library_paths.is_empty() {
            anyhow::bail!("compare_to_library=true but no library_files configured");
        }

        let library_path_strings: Vec<String> = library_paths
            .iter()
            .map(|p| p.to_string_lossy().to_string())
            .collect();

        let library = Library::from_files(&library_path_strings, spec, case.max_distance)?;

        counts.compare_to_library(
            library,
            Some(&progress_style),
            case.distance_metric,
            case.max_matches,
            case.threads,
        )?;
    }

    // Write output
    let dir = case.dir();
    let counts_tsv = dir.join("observed.counts.tsv");
    let summary_tsv = dir.join("observed.summary.tsv");
    let filtered_tsv = dir.join("observed.filtered.tsv");
    let library_counts_tsv = if case.compare_to_library {
        Some(dir.join("observed.library_counts.tsv"))
    } else {
        None
    };

    write_counts(&counts, File::create(&counts_tsv)?, true)?;

    if case.compare_to_library {
        write_library_counts(
            &counts,
            File::create(library_counts_tsv.as_ref().unwrap())?,
            true,
        )?;
    }

    let summary = counts.summarise();
    write_summary(&summary, File::create(&summary_tsv)?).context("writing summary.tsv")?;

    write_filter_summary(&counts.filtered_reads(), File::create(&filtered_tsv)?, true)?;

    // Check outputs
    assert_matches_expected(&dir.join("expected.counts.tsv"), &counts_tsv)?;
    assert_matches_expected(&dir.join("expected.summary.tsv"), &summary_tsv)?;
    assert_matches_expected(&dir.join("expected.filtered.tsv"), &filtered_tsv)?;

    match (case.compare_to_library, library_counts_tsv) {
        (true, None) => anyhow::bail!("Expected library counts TSV missing"),
        (false, Some(_)) => anyhow::bail!("Unexpected library counts TSV present"),
        (false, None) => {}
        (true, Some(path)) => {
            assert_matches_expected(&dir.join("expected.library_counts.tsv"), &path)?;
        }
    }

    Ok(())
}

// Tests
#[test]
fn full_read_single() -> Result<()> {
    let case = GoldenCase {
        dir: "tests/full_read_single",
        mode: CountMode::FullRead,
        paired: false,
        full_seq: true,
        group_regex: None,
        use_libspec: false,
        mean_quality_threshold: None,
        minimum_read_length: None,
        maximum_read_length: None,
        alignment_tolerance: None,
        filter_empty: true,
        pattern_length: None,
        pattern_tolerance: None,
        cache: true,
        threads: 1,
        compare_to_library: false,
        library_files: &[],
        distance_metric: DistanceMetric::Hamming,
        max_distance: 3,
        max_matches: 10,
    };

    run_case(&case)
}

#[test]
fn full_read_paired() -> Result<()> {
    let case = GoldenCase {
        dir: "tests/full_read_paired",
        mode: CountMode::FullRead,
        paired: true,
        full_seq: true,
        group_regex: None,
        use_libspec: false,
        mean_quality_threshold: None,
        minimum_read_length: None,
        maximum_read_length: None,
        alignment_tolerance: None,
        filter_empty: true,
        pattern_length: None,
        pattern_tolerance: None,
        cache: true,
        threads: 1,
        compare_to_library: false,
        library_files: &[],
        distance_metric: DistanceMetric::Hamming,
        max_distance: 3,
        max_matches: 10,
    };

    run_case(&case)
}

#[test]
fn inframe() -> Result<()> {
    let case = GoldenCase {
        dir: "tests/inframe",
        mode: CountMode::Inframe,
        paired: true,
        full_seq: true,
        group_regex: None,
        use_libspec: true,
        mean_quality_threshold: None,
        minimum_read_length: None,
        maximum_read_length: None,
        alignment_tolerance: None,
        filter_empty: true,
        pattern_length: None,
        pattern_tolerance: None,
        cache: true,
        threads: 1,
        compare_to_library: true,
        library_files: &["library.tsv"],
        distance_metric: DistanceMetric::Hamming,
        max_distance: 3,
        max_matches: 10,
    };

    run_case(&case)
}

#[test]
fn filtering() -> Result<()> {
    let case = GoldenCase {
        dir: "tests/filtering",
        mode: CountMode::Align,
        paired: false,
        full_seq: true,
        group_regex: None,
        use_libspec: true,
        mean_quality_threshold: Some(30.0),
        minimum_read_length: Some(10),
        maximum_read_length: Some(31),
        alignment_tolerance: Some(0.8),
        filter_empty: true,
        pattern_length: None,
        pattern_tolerance: None,
        cache: true,
        threads: 1,
        compare_to_library: true,
        library_files: &["library.tsv"],
        distance_metric: DistanceMetric::Hamming,
        max_distance: 3,
        max_matches: 10,
    };

    run_case(&case)
}
