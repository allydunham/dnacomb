//! CLI interface to the DNAComb read counting and library comparison tool.
use anyhow::Error;
use bio::alignment::pairwise::Aligner;
use clap::{ArgAction, Parser};
use log::{self, LevelFilter, debug, error, info, warn};
use regex::Regex;
use std::fs;
use std::process::exit;
use std::str;

use dnacomb::ObservedCombinations;
use dnacomb::counting::{AlignmentScorer, CountMode, count_reads};
use dnacomb::filters::{AlignmentTolerance, FilterConfig};
use dnacomb::lib_spec::{DistanceMetric, Library, LibrarySpec};
use dnacomb::logging::ProgressStyle;
use dnacomb::parsing::{Compression, ReadPairParser, ReadPairProducer, SeqFormat, SeqPath};

/// Fast general purpose read counter supporting complex structured reads
///
/// Count occurances of structured sequence reads, extracting regions of interest from
/// a library specification and comparing observations to a library of expected combinations.
/// Supports several alignment modes appropriate for different read structures and data sizes:
///
/// Alignment - Align reads to the template via semi-global alignment. Most thorough but
/// slowest, allows variable region lengths [Default]
///
/// Pattern matching - Use flanking regions to identify regions. Faster than alignment
/// while allowing variable region length but less robust against variation.
///
/// Inframe - Assume regions occur at the correct position in reads (for instance after
/// using cutadapt). Fastest structured read counting but can't handle variation.
///
/// Raw - Count full length sequences, fastest but unstructured
///
/// Outputs 1-3 files and logs to stderr:
/// * {prefix}.counts.tsv - full count table
/// * {prefix}.library_counts.tsv - summarised counts of library matches only
/// * {prefix}.summary.tsv - summary read counts for e.g. matches, recombinations, mismatches
/// * {prefix}.filtered.tsv - filtered read counts for e.g. bad alignments, low quality, short reads
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

    /// Additionally store the full read sequence associated with each combination. This is useful for debugging
    /// computation and experiments but also generally each combination is associated with many sequences so it
    /// doesn't compress results as much.
    #[arg(short = 'F', long, action, help_heading = "Counting")]
    full_seq: bool,

    /// Group counts by applying this capture group regex to forward read names and
    /// extracting the first capture group match
    #[arg(short = 'g', long, help_heading = "Counting")]
    group: Option<String>,

    /// Calculate similarity to oligo library and output an additional table of library counts
    #[arg(short = 'c', long, action, help_heading = "Library Comparison")]
    library_counts: bool,

    /// Distance metric to use for library comparison. Hamming counts the number of mismatches
    /// and levenshtein the number of subs/insertion/deletions required to go from A to B.
    /// Bounded levenshtein only considers up to some distance and should be faster than unbounded.
    /// Since we only consider matches with a max distance bounded should be prefered unless a
    /// special case applies.
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

    /// Minimum proportion of expected alignment score to keep
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

    /// Number of mismatches to accept while matching flanking patterns
    #[arg(long, default_value_t = 1, help_heading = "Pattern Matching")]
    pattern_tolerance: u64,

    /// Match score for alignment
    #[arg(long, default_value_t = 6, allow_hyphen_values = true, help_heading = "Alignment")]
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

    /// Don't cache reads when in align mode, trading lower memory usage for lower speed
    #[arg(long, action = ArgAction::SetTrue, help_heading = "Technical")]
    no_cache: bool,

    /// Only process the first n reads for testing. Set to 0 (default) for all reads
    #[arg(long, default_value_t = 0, help_heading = "Technical")]
    max_reads: u64,

    /// Phred value to assume for Fasta files. Only matters when comparing to Fastq.
    #[arg(long, default_value_t = b'I', help_heading = "Technical")]
    default_phred: u8,

    /// Number of threads to use
    #[arg(short = 'T', long, default_value_t = 1, help_heading = "Technical")]
    threads: usize,
}

/// Main function
///
/// Runs the pipeline and captures errors
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

/// Run the CLI pipeline
///
/// Manages control flow for CLI tool, dispatching to the appropriate functions
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
        Some(lib_spec) => Some(LibrarySpec::from_file(lib_spec)?),
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
    match (args.library_counts, &lib_spec) {
        (false, _) => {
            library = None;
            info!("No library counts requested, library ignored");
        }
        (true, None) => {
            library = None;
            warn!(
                "Library counts requested but no library path in LibSpec, skipping library counts"
            );
        }
        (true, Some(spec)) => {
            library = Library::from_lib_spec(spec, args.max_distance)?;
            info!("Compiled library {:?}", spec.library.as_ref().unwrap());
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
        (Some(l), Some(t)) => calculate_alignment_tolerance(l, &alignment_scorer, &reader, t)?,
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

    if args.library_counts && library.is_none() {
        warn!(
            "Can't calculate library counts without a library specification that includes a \
             library file. Skipping library comparison."
        )
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
        counts.write_tsv(count_file, args.sort)?;
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
        counts.write_summary_tsv(library_summary_file, args.sort)?;
    }

    // Calculate and output summary statistics
    let read_summary = counts.summarise();
    {
        info!("Writing read count summary TSV to: {}", read_summary_path);
        let summary_file = match args.overwrite {
            true => fs::File::create(read_summary_path)?,
            false => fs::File::create_new(read_summary_path)?,
        };
        read_summary.write_tsv(summary_file)?;
    }

    // Write filter count table
    {
        info!("Writing filtered reads counts TSV to: {}", filtered_path);
        let filtered_file = match args.overwrite {
            true => fs::File::create(filtered_path)?,
            false => fs::File::create_new(filtered_path)?,
        };
        counts.write_filtered_tsv(filtered_file, args.sort)?;
    }

    Ok(())
}

/// Calculate alignment tolerance based on the expected alignment scores
fn calculate_alignment_tolerance(
    lib_spec: &LibrarySpec,
    alignment_scorer: &AlignmentScorer,
    reader: &ReadPairParser,
    tolerance: f32,
) -> Result<Option<AlignmentTolerance>, anyhow::Error> {
    let exp_f_read = lib_spec.expected_forward_read();

    // Initialise aligner
    let scoring = alignment_scorer.get_scoring();
    let mut aligner = Aligner::with_capacity_and_scoring(400, 150, scoring);
    let template = lib_spec.template_sequence();

    let f_alignment = aligner.semiglobal(&exp_f_read, &template);

    info!(
        "Expected Fwd Alignment:\nScore: {}, Cigar: {}\n{}",
        f_alignment.score,
        f_alignment.cigar(false),
        f_alignment.pretty(&exp_f_read, &template, 100),
    );

    let mut r_score = 0;
    if reader.has_reverse() {
        let exp_r_read = lib_spec.expected_reverse_read();
        let r_alignment = aligner.semiglobal(&exp_r_read, &template);
        r_score = r_alignment.score;

        info!(
            "Expected Rev Alignment:\nScore: {}, Cigar: {}\n{}",
            r_alignment.score,
            r_alignment.cigar(false),
            r_alignment.pretty(&exp_r_read, &template, 100),
        );
    }

    Ok(Some(AlignmentTolerance::new(
        tolerance,
        f_alignment.score,
        r_score,
    )?))
}

/// Log SIMD feature presence and whether the tool is compiled
/// to use them.
///
/// This version is for x86/x86_64 where SIMD may be available
/// via AVX2/SSE4.1
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

/// Log SIMD feature presence and whether the tool is compiled
/// to use them.
///
/// This version is for other architectures where SIMD will not
/// be available.
#[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
fn check_simd_features() {
    log::info!("SIMD features are not available on this architecture.");
}
