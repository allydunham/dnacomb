//! CLI entry point for the DNAComb structured-read counting tool.
//!
//! This binary wires together argument parsing, input validation, read parsing,
//! counting, optional library comparison, and TSV output generation.
use anyhow::Error;
use clap::{ArgAction, Parser};
use log::{self, LevelFilter, debug, error, info, warn};
use regex::Regex;
use std::fs;
use std::process::exit;
use std::str;

use dnacomb::counting::{AlignmentScorer, CountMode, count_reads};
use dnacomb::filters::{AlignmentTolerance, FilterConfig};
use dnacomb::lib_spec::LibrarySpec;
use dnacomb::library::{DistanceMetric, Library};
use dnacomb::logging::ProgressStyle;
use dnacomb::parsing::{Compression, ReadPairParser, ReadPairProducer, SeqFormat, SeqPath};
use dnacomb::{
    ObservedCombinations, write_counts, write_filter_summary, write_library_counts, write_summary,
};

/// Fast general-purpose read counter for structured sequencing reads.
///
/// DNAComb extracts variable regions from single-end or paired-end reads using
/// one of several region-extraction modes, optionally compares those regions to
/// one or more expected libraries, and writes TSV outputs summarising observed
/// combinations, inferred library assignments, read-category summaries, and
/// filtered reads.
///
/// Supported counting modes:
/// - `align`: semi-global alignment to the LibSpec template; most robust, slowest
/// - `pattern`: flank-pattern matching; faster, less robust to flank mutations
/// - `inframe`: positional extraction from expected read coordinates; fastest structured mode
/// - `full-read`: count complete read sequences without structured extraction
///
/// Outputs written:
/// - `{prefix}.counts.tsv`
/// - `{prefix}.library_counts.tsv` (when library comparison is performed)
/// - `{prefix}.summary.tsv`
/// - `{prefix}.filtered.tsv`
#[derive(Parser, Debug)]
#[command(author, version)]
struct Cli {
    /// Forward reads in Fasta or Fastq format
    forward: String,

    /// Reverse reads in Fasta or Fastq format
    reverse: Option<String>,

    /// Path to JSON file defining the library specification
    #[arg(short = 'l', long, help_heading = "Input")]
    library_spec: Option<String>,

    /// Input file(s) format
    #[arg(short = 'f', long, value_enum, default_value_t = SeqFormat::Auto, help_heading = "Input")]
    format: SeqFormat,

    /// Input file(s) compression. Currently only supports Gzip.
    #[arg(short = 'z', long, value_enum, default_value_t = Compression::Auto, help_heading = "Input")]
    compression: Compression,

    /// Override fordard start region
    #[arg(long, help_heading = "Input")]
    forward_start: Option<String>,

    /// Override fordard read length
    #[arg(long, help_heading = "Input")]
    forward_length: Option<u32>,

    /// Override reverse start region
    #[arg(long, help_heading = "Input")]
    reverse_start: Option<String>,

    /// Override reverse read length
    #[arg(long, help_heading = "Input")]
    reverse_length: Option<u32>,

    /// Prefix for output TSV files
    #[arg(
        short = 'o',
        long,
        default_value = "read_counts",
        help_heading = "Output"
    )]
    output: String,

    /// Sort output by descending count. Considers the total count across read groups
    /// when using grouping.
    #[arg(short = 's', long, action, help_heading = "Output")]
    sort: bool,

    /// Overwrite output if it exists
    #[arg(short = 'w', long, action, help_heading = "Output")]
    overwrite: bool,

    /// Verbose logging. More fine grained control available with RUST_LOG environment
    /// variable. -v is equivalent to RUST_LOG=info. Progress bars with timings and similar progress
    /// indicators are output at INFO level. Otherwise only warnings and errors will be displayed.
    #[arg(short = 'v', long, action, help_heading = "Output")]
    verbose: bool,

    /// Method for extracting regions of interest
    #[arg(short = 'm', long, value_enum, default_value_t = CountMode::Align, help_heading = "Counting")]
    mode: CountMode,

    /// Store full read sequence(s) alongside each observed combination.
    /// Useful for debugging and inspection, but increases output size because
    /// many observed combinations correspond to multiple underlying full reads.
    #[arg(short = 'F', long, action, help_heading = "Counting")]
    full_seq: bool,

    /// Group counts by applying this capture-group regex to forward read names
    /// and using the first capture as the group label; reads without a match are
    /// assigned to `_unmatched_`
    #[arg(short = 'g', long, help_heading = "Counting")]
    group: Option<String>,

    /// Compare observed regions to one or more expected library TSVs and write
    /// an additional library-summary output table
    #[arg(short = 'c', long, num_args = 1.., value_delimiter = ' ', help_heading = "Library Comparison")]
    library: Option<Vec<String>>,

    /// Distance metric to use for library comparison.
    ///
    /// Hamming counts substitutions only, while Levenshtein allows substitutions,
    /// insertions, and deletions. Bounded Levenshtein applies the configured
    /// maximum distance threshold during lookup and is usually faster while giving
    /// the same accepted matches.
    #[arg(
        short = 'd',
        long,
        value_enum,
        default_value_t = DistanceMetric::Hamming,
        help_heading = "Library Comparison"
    )]
    distance_metric: DistanceMetric,

    /// Max distance to consider for library comparisons when not specified for that
    /// region in LibSpec
    #[arg(
        short = 'x',
        long,
        default_value_t = 3,
        help_heading = "Library Comparison"
    )]
    max_distance: u64,

    /// Max number of matches to consider for library comparisons. More than this many equivalent
    /// matches will be considered indeterminant.
    #[arg(
        short = 't',
        long,
        default_value_t = 10,
        help_heading = "Library Comparison"
    )]
    max_matches: usize,

    /// Filter reads with mean Phred score below this threshold
    #[arg(short = 'q', long, help_heading = "Filtering")]
    mean_quality_threshold: Option<f32>,

    /// Minimum proportion of expected alignment score to keep. Only used in alignment mode.
    #[arg(short = 'r', long, help_heading = "Filtering")]
    alignment_tolerance: Option<f32>,

    /// Filter reads shorter than this length. Both F and R must meet the threshold.
    /// For technical reasons reads were either F or R is empty are always filtered.
    #[arg(short = 'L', long, help_heading = "Filtering")]
    minimum_read_length: Option<usize>,

    /// Filter reads longer than this length. Both F and R must meet the threshold.
    #[arg(short = 'M', long, help_heading = "Filtering")]
    maximum_read_length: Option<usize>,

    /// Length of flanking pattern to use (where possible)
    #[arg(long, default_value_t = 10, help_heading = "Pattern Matching")]
    pattern_length: usize,

    /// Number of mismatches to accept while matching flanking patterns. Only using in pattern mode.
    #[arg(long, default_value_t = 1, help_heading = "Pattern Matching")]
    pattern_tolerance: u64,

    /// Match score for alignment
    #[arg(
        long,
        default_value_t = 6,
        allow_hyphen_values = true,
        help_heading = "Alignment"
    )]
    match_score: i32,

    /// Match score against Ns in alignment
    #[arg(long, default_value_t = -2, allow_hyphen_values = true, help_heading = "Alignment")]
    n_match_score: i32,

    /// Misatch score for alignment
    #[arg(long, default_value_t = -3, allow_hyphen_values = true, help_heading = "Alignment")]
    mismatch_score: i32,

    /// Gap open score for alignment
    #[arg(long, default_value_t = -10, allow_hyphen_values = true, help_heading = "Alignment")]
    gap_open_score: i32,

    /// Gap extension score for alignment
    #[arg(long, default_value_t = -4, allow_hyphen_values = true, help_heading = "Alignment")]
    gap_extend_score: i32,

    /// Disable read-level caching in align mode, trading lower memory usage for
    /// slower repeated processing of duplicate read sequences
    #[arg(long, action = ArgAction::SetTrue, help_heading = "Technical")]
    no_cache: bool,

    /// Only process the first n reads for testing. Set to 0 (default) for all reads
    #[arg(long, default_value_t = 0, help_heading = "Technical")]
    max_reads: u64,

    /// Phred byte to assign when reading FASTA input, which lacks quality scores. Only meaningful
    /// if comparing Fasta and Fastq.
    #[arg(long, default_value_t = b'I', help_heading = "Technical")]
    default_phred: u8,

    /// Number of threads to use
    #[arg(short = 'T', long, default_value_t = 1, help_heading = "Technical")]
    threads: usize,
}

/// Parse CLI arguments, initialise logging, and run the main pipeline.
fn main() -> Result<(), Error> {
    // Process arguments
    let args: Cli = Cli::parse();

    if args.verbose {
        unsafe { std::env::set_var("RUST_LOG", "info") };
    }

    env_logger::Builder::from_default_env()
        .format_target(false)
        .format_indent(None)
        .init();

    match run(args) {
        Ok(_) => {}
        Err(e) => {
            error!("{}", e);
            exit(1)
        }
    }
    Ok(())
}

/// Execute the full CLI workflow.
///
/// This function:
/// - validates output paths and options,
/// - initialises read parsing and grouping,
/// - loads and validates the LibSpec,
/// - optionally compiles expected library TSVs,
/// - constructs alignment/filtering configuration,
/// - runs read counting,
/// - optionally performs library comparison,
/// - and writes all requested TSV outputs.
fn run(args: Cli) -> Result<(), Error> {
    info!("Using options: {:#?}", args);

    let progress_style: ProgressStyle = match log::max_level() {
        LevelFilter::Off | LevelFilter::Error | LevelFilter::Warn => {
            ProgressStyle::new(None, false)
        }
        LevelFilter::Info | LevelFilter::Debug | LevelFilter::Trace => ProgressStyle::default(),
    };

    // Check if output files exist
    let count_path = format!("{}.counts.tsv", args.output);
    if !args.overwrite && fs::exists(&count_path)? {
        error!("File exists \"{count_path}\" with --overwrite disabled. Exiting");
        exit(1);
    }

    let library_summary_path = format!("{}.library_counts.tsv", args.output);
    if !args.overwrite && fs::exists(&library_summary_path)? {
        error!("File exists \"{library_summary_path}\" with --overwrite disabled. Exiting");
        exit(1);
    }

    let read_summary_path = format!("{}.summary.tsv", args.output);
    if !args.overwrite && fs::exists(&read_summary_path)? {
        error!("File exists \"{read_summary_path}\" with --overwrite disabled. Exiting");
        exit(1);
    }

    let filtered_path: String = format!("{}.filtered.tsv", args.output);
    if !args.overwrite && fs::exists(&filtered_path)? {
        error!("File exists \"{filtered_path}\" with --overwrite disabled. Exiting");
        exit(1);
    }

    // Check SIMD status
    check_simd_features();

    // Process grouping
    let group = match args.group {
        None => None,
        Some(x) => Some(Regex::new(&x)?),
    };
    match group {
        None => {}
        Some(ref x) => info!("Grouping reads using regex: {}", x),
    };

    // Initialise Parser
    let forward = SeqPath::new(args.forward, args.format, args.compression);
    let reverse = match args.reverse {
        Some(x) => Some(SeqPath::new(x, args.format, args.compression)),
        None => None,
    };

    match &reverse {
        None => info!("Reading seqs from F: {}, R: None", forward),
        Some(r) => info!("Reading seqs from F: {}, R: {}", forward, r),
    }

    let reader: ReadPairParser =
        ReadPairParser::from_paths(forward, reverse, group, args.max_reads, args.default_phred)?;

    // Load library specification
    let lib_spec: Option<LibrarySpec> = match &args.library_spec {
        Some(lib_spec) => Some(LibrarySpec::from_file(
            lib_spec,
            args.forward_start,
            args.forward_length,
            args.reverse_start,
            args.reverse_length,
        )?),
        None => None,
    };

    if lib_spec.is_some() {
        info!("Read LibSpec: {}", &args.library_spec.unwrap());
        debug!("Parsed LibSpec: {:?}", lib_spec);
    }

    // Validate alignment mode is suitable
    if let Some(l) = &lib_spec {
        if matches!(args.mode, CountMode::Inframe) && l.variable_length_regions() > 0 {
            error!("Can't use inframe matching with variable length regions. Exiting");
            exit(1)
        }
    }

    // Compare observed combinations to library
    let library: Option<Library>;
    match (args.library, &lib_spec) {
        (None, _) => {
            library = None;
            info!("No library counts requested");
        }
        (Some(_), None) => {
            error!(
                "Library TSV(s) supplied but no LibSpec, which is required for library counts. Exiting"
            );
            exit(1)
        }
        (Some(libs), Some(spec)) => {
            library = Some(Library::from_files(&libs, spec, args.max_distance)?);
            info!("Compiled library from files: {:?}", libs.join(", "));
        }
    };

    // Setup aligner - established with default
    let alignment_scorer = AlignmentScorer::new(
        args.match_score,
        args.n_match_score,
        args.mismatch_score,
        args.gap_open_score,
        args.gap_extend_score,
    );

    // Setup read filtering
    if args.alignment_tolerance.is_some() && !matches!(args.mode, CountMode::Align) {
        warn!("--alignment-tolerance passed without --mode align. Tolerance has no effect.")
    }

    // First check expected alignment score for filter threshold
    let alignment_tolerance = match (&lib_spec, args.alignment_tolerance) {
        (None, None) => None,
        (Some(_), None) => None,
        (None, Some(_)) => {
            warn!(
                "--alignment-tolerance passed without a LibSpec, meaning no alignment will be done."
            );
            None
        }
        (Some(l), Some(t)) => {
            let r = if reader.has_reverse() {
                Some(l.expected_reverse_read())
            } else {
                None
            };

            Some(AlignmentTolerance::from_expected_reads(
            &l.expected_forward_read(),
            r.as_ref(),
            &l.template_sequence(),
            &alignment_scorer,
            t,
            true
        )?)},
    };

    let filter_config = FilterConfig::new(
        args.mean_quality_threshold,
        alignment_tolerance,
        args.minimum_read_length,
        args.maximum_read_length,
        true,
    );

    info!("Filtering reads with config: {:?}", filter_config);

    // Enumerate region combinations
    let mut counts: ObservedCombinations = count_reads(
        reader,
        &lib_spec,
        args.mode,
        args.full_seq,
        filter_config,
        Some(alignment_scorer),
        Some(args.pattern_length),
        Some(args.pattern_tolerance),
        !args.no_cache,
        args.threads,
        Some(&progress_style),
    )?;

    if counts.is_empty() && counts.total_filtered() == 0 {
        error!("No observed combinations or filtered reads (empty fasta?). Exiting");
        debug!("ObservedCombinations:\n{:?}", counts);
        exit(1);
    }

    if counts.is_empty() && counts.total_filtered() > 0 {
        warn!("All reads filtered. Check input files and filter settings.");
    }

    match library {
        None => (),
        Some(x) => {
            info!("Comparing observed regions to library");
            counts.compare_to_library(
                x,
                Some(&progress_style),
                args.distance_metric,
                args.max_matches,
                args.threads,
            )?;
        }
    }

    // Write Overall count table
    {
        info!("Writing full count TSV to: {}", count_path);
        let count_file = match args.overwrite {
            true => fs::File::create(count_path)?,
            false => fs::File::create_new(count_path)?,
        };
        write_counts(&counts, count_file, args.sort)?;
    }

    // Write library count table if applicable
    if counts.is_compared_to_library() {
        info!(
            "Writing library count summary TSV to: {}",
            library_summary_path
        );
        let library_summary_file = match args.overwrite {
            true => fs::File::create(library_summary_path)?,
            false => fs::File::create_new(library_summary_path)?,
        };
        write_library_counts(&counts, library_summary_file, args.sort)?;
    }

    // Calculate and output summary statistics
    let read_summary = counts.summarise();
    {
        info!("Writing read count summary TSV to: {}", read_summary_path);
        let summary_file = match args.overwrite {
            true => fs::File::create(read_summary_path)?,
            false => fs::File::create_new(read_summary_path)?,
        };
        write_summary(&read_summary, summary_file)?;
    }

    // Write filter count table
    {
        info!("Writing filtered reads counts TSV to: {}", filtered_path);
        let filtered_file = match args.overwrite {
            true => fs::File::create(filtered_path)?,
            false => fs::File::create_new(filtered_path)?,
        };
        write_filter_summary(counts.filtered_reads(), filtered_file, args.sort)?;
    }

    Ok(())
}

/// Log whether SIMD-accelerated distance calculations are available and enabled
/// for this build/runtime combination.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
fn check_simd_features() {
    match (
        cfg!(all(target_feature = "avx2", target_feature = "sse4.1")),
        std::is_x86_feature_detected!("avx2") && std::is_x86_feature_detected!("sse4.1"),
    ) {
        (true, true) => info!("Using SIMD instructions to speed up distance calculations."),
        (true, false) => info!(
            "Compiled with SIMD instructions but AVX2/SSE4.1 not available so \
                                falling back to regular distance metrics."
        ),
        (false, true) => info!(
            "Compiled without SIMD instructions but AVX2 & SSE4.1 are \
                                available, consider re-compiling to benefit from SIMD \
                                speed increases."
        ),
        (false, false) => info!("Compiled without SIMD instructions and AVX2/SSE4.1 unavailable."),
    }
}

/// Log that SIMD-accelerated distance calculations are unavailable on this architecture.
#[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
fn check_simd_features() {
    log::info!("SIMD features are not available on this architecture.");
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    // ---------- CLI PARSING ----------
    #[test]
    fn cli_parse_minimal_args() {
        let args = Cli::try_parse_from(vec!["prog", "forward.fq"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert_eq!(cli.forward, "forward.fq");
        assert!(cli.reverse.is_none());
        assert_eq!(cli.output, "read_counts");
    }

    #[test]
    fn cli_parse_with_reverse() {
        let args = Cli::try_parse_from(vec!["prog", "forward.fq", "reverse.fq"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert_eq!(cli.forward, "forward.fq");
        assert_eq!(cli.reverse.as_ref().unwrap(), "reverse.fq");
    }

    #[test]
    fn cli_parse_with_library_spec() {
        let args = Cli::try_parse_from(vec!["prog", "-l", "spec.json", "forward.fq"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert_eq!(cli.library_spec.unwrap(), "spec.json");
    }

    #[test]
    fn cli_parse_with_output_prefix() {
        let args = Cli::try_parse_from(vec!["prog", "-o", "my_output", "forward.fq"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert_eq!(cli.output, "my_output");
    }

    #[test]
    fn cli_parse_with_format() {
        let args = Cli::try_parse_from(vec!["prog", "-f", "Fasta", "forward.fa"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert_eq!(cli.format, SeqFormat::Fasta);
    }

    #[test]
    fn cli_parse_with_compression() {
        let args = Cli::try_parse_from(vec!["prog", "-z", "Gzip", "forward.fq.gz"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert_eq!(cli.gzip, Compression::Gzip);
    }

    #[test]
    fn cli_parse_with_mode() {
        let args = Cli::try_parse_from(vec!["prog", "-m", "pattern", "forward.fq"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert_eq!(cli.mode, CountMode::Pattern);
    }

    #[test]
    fn cli_parse_with_distance_metric() {
        let args = Cli::try_parse_from(vec!["prog", "-d", "Levenshtein", "forward.fq"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert_eq!(cli.distance_metric, DistanceMetric::Levenshtein);
    }

    #[test]
    fn cli_parse_sorting_flag() {
        let args = Cli::try_parse_from(vec!["prog", "-s", "forward.fq"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert!(cli.sort);
    }

    #[test]
    fn cli_parse_verbose_flag() {
        let args = Cli::try_parse_from(vec!["prog", "-v", "forward.fq"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert!(cli.verbose);
    }

    #[test]
    fn cli_parse_overwrite_flag() {
        let args = Cli::try_parse_from(vec!["prog", "-w", "forward.fq"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert!(cli.overwrite);
    }

    #[test]
    fn cli_parse_full_seq_flag() {
        let args = Cli::try_parse_from(vec!["prog", "-F", "forward.fq"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert!(cli.full_seq);
    }

    #[test]
    fn cli_parse_grouping_regex() {
        let args = Cli::try_parse_from(vec!["prog", "-g", r"([A-Z]+)_", "forward.fq"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert_eq!(cli.group.as_ref().unwrap(), r"([A-Z]+)_");
    }

    #[test]
    fn cli_parse_quality_threshold() {
        let args = Cli::try_parse_from(vec!["prog", "-q", "20.5", "forward.fq"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert_eq!(cli.mean_quality_threshold.unwrap(), 20.5);
    }

    #[test]
    fn cli_parse_alignment_tolerance() {
        let args = Cli::try_parse_from(vec!["prog", "-r", "0.9", "forward.fq"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert_eq!(cli.alignment_tolerance.unwrap(), 0.9);
    }

    #[test]
    fn cli_parse_read_length_filters() {
        let args = Cli::try_parse_from(vec!["prog", "-L", "50", "-M", "200", "forward.fq"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert_eq!(cli.minimum_read_length.unwrap(), 50);
        assert_eq!(cli.maximum_read_length.unwrap(), 200);
    }

    #[test]
    fn cli_parse_pattern_length_and_tolerance() {
        let args = Cli::try_parse_from(vec![
            "prog",
            "--pattern-length",
            "15",
            "--pattern-tolerance",
            "2",
            "forward.fq",
        ]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert_eq!(cli.pattern_length, 15);
        assert_eq!(cli.pattern_tolerance, 2);
    }

    #[test]
    fn cli_parse_alignment_scores() {
        let args = Cli::try_parse_from(vec![
            "prog",
            "--match-score",
            "5",
            "--mismatch-score",
            "-4",
            "--gap-open-score",
            "-8",
            "forward.fq",
        ]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert_eq!(cli.match_score, 5);
        assert_eq!(cli.mismatch_score, -4);
        assert_eq!(cli.gap_open_score, -8);
    }

    #[test]
    fn cli_parse_threads() {
        let args = Cli::try_parse_from(vec!["prog", "-T", "8", "forward.fq"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert_eq!(cli.threads, 8);
    }

    #[test]
    fn cli_parse_max_reads() {
        let args = Cli::try_parse_from(vec!["prog", "--max-reads", "1000", "forward.fq"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert_eq!(cli.max_reads, 1000);
    }

    #[test]
    fn cli_parse_default_phred() {
        let args = Cli::try_parse_from(vec!["prog", "--default-phred", "35", "forward.fq"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert_eq!(cli.default_phred, 35);
    }

    #[test]
    fn cli_parse_no_cache_flag() {
        let args = Cli::try_parse_from(vec!["prog", "--no-cache", "forward.fq"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        assert!(cli.no_cache);
    }

    #[test]
    fn cli_parse_library_files() {
        let args = Cli::try_parse_from(vec!["prog", "-c", "lib1.tsv", "lib2.tsv", "forward.fq"]);
        assert!(args.is_ok());
        let cli = args.unwrap();
        let libs = cli.library.unwrap();
        assert_eq!(libs.len(), 2);
        assert_eq!(libs[0], "lib1.tsv");
        assert_eq!(libs[1], "lib2.tsv");
    }

    // Simd detection
    #[test]
    fn check_simd_features_runs() {
        // Just verify the function runs without panicking
        check_simd_features();
    }
}
