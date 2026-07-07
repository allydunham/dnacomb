//! Library supporting the [DNAComb] structured sequence read processing tool.
//! The [GitHub] and [Crates.io] page give documentation for the CLI tool with
//! these pages documenting the library,
//!
//! DNAComb provides functionality for:
//! - defining expected construct structure with `LibrarySpec`,
//! - reading and filtering single-end or paired-end sequencing reads,
//! - extracting variable regions from structured reads,
//! - counting distinct observed region combinations,
//! - and optionally comparing observed regions to expected library sequences.
//!
//! The crate is designed primarily to support the `dnacomb` CLI, but the core
//! counting and library-comparison logic can also be used programmatically from
//! Rust.
//!
//! The main entry points are:
//! - [`LibrarySpec`] for describing construct structure,
//! - [`count_reads`] for extracting and counting observed read forms,
//! - [`ObservedCombinations`] for storing and post-processing counts,
//! - [`SubLibrary`] for compiling expected library sequence sets,
//! - and the `write_*` functions in [`output`] for generating TSV outputs.
//!
//! For CLI usage, see the project README and `dnacomb --help`.
//!
//! # Example
//!
//! ```no_run
//! use std::fs::File;
//!
//! use dnacomb::{
//!     count_reads, write_counts, Compression, CountMode, FilterConfig, LibrarySpec,
//!     ReadPairParser, SeqFormat, SeqPath,
//! };
//!
//! # fn run() -> Result<(), Box<dyn std::error::Error>> {
//! // Load a LibSpec
//! let lib_spec = LibrarySpec::from_file(
//!     "config/example_libspec.json",
//!     None,
//!     None,
//!     None,
//!     None,
//! )?;
//!
//! // Create input paths and parser
//! let forward = SeqPath::new(
//!     "reads_R1.fastq.gz".to_string(),
//!     SeqFormat::Auto,
//!     Compression::Auto,
//! );
//! let reverse = Some(SeqPath::new(
//!     "reads_R2.fastq.gz".to_string(),
//!     SeqFormat::Auto,
//!     Compression::Auto,
//! ));
//!
//! let reader = ReadPairParser::from_paths(
//!     forward,
//!     reverse,
//!     None,
//!     0,
//!     b'I',
//! )?;
//!
//! // Configure filtering
//! let filter_config = FilterConfig::new(
//!     Some(20.0),
//!     None,
//!     Some(50),
//!     Some(300),
//!     true,
//! );
//!
//! // Count reads
//! let counts = count_reads(
//!     reader,
//!     &Some(lib_spec),
//!     CountMode::Pattern,
//!     false,
//!     filter_config,
//!     None,
//!     Some(12),
//!     Some(2),
//!     true,
//!     1,
//!     None,
//! )?;
//!
//! // Write full count table
//! let out = File::create("example.counts.tsv")?;
//! write_counts(&counts, out, true, false)?;
//!
//! # Ok(())
//! # }
//! ```
//!
//! [DNAComb]: https://github.com/allydunham/dnacomb
//! [GitHub]: https://github.com/allydunham/dnacomb
//! [Crates.io]: https://crates.io/crates/dnacomb
pub mod combination;
pub mod combinations;
pub mod counting;
pub mod errors;
pub mod filters;
pub mod groups;
pub mod interning;
pub mod lib_spec;
pub mod library;
pub mod library_combination;
pub mod logging;
pub mod output;
pub mod parsing;
pub mod region;
pub mod seq_diff;
pub mod seqs;
pub mod utils;

/// Parse single-end or paired-end sequence files into `ReadPair`s.
pub use parsing::ReadPairParser;

/// Path plus format/compression metadata for sequence-file input.
pub use parsing::SeqPath;

/// Sequence file format selection.
pub use parsing::SeqFormat;

/// Compression mode for sequence-file input.
pub use parsing::Compression;

/// Specification of construct structure and expected read layout.
pub use lib_spec::LibrarySpec;

/// Compiled expected sequence library for one sublibrary.
pub use library::SubLibrary;

/// Distance metric used for library comparison.
pub use library::DistanceMetric;

/// Filter reads with various configuration options.
pub use filters::FilterConfig;

/// Count observed read forms.
pub use counting::count_reads;

/// Supported modes for extracting regions from reads
pub use counting::CountMode;

/// Alignment scoring scheme for template-based region extraction.
pub use counting::AlignmentScorer;

/// Container for observed region combinations and their counts.
pub use combinations::ObservedCombinations;

/// Write the full observed-combination count table as TSV.
pub use output::write_counts;

/// Write the filtered-read summary table as TSV.
pub use output::write_filter_summary;

/// Write the library-combination summary table as TSV.
pub use output::write_library_counts;

/// Write the high-level read summary table as TSV.
pub use output::write_summary;

/// Configure loggins
pub use logging::ProgressStyle;

#[cfg(not(any(target_pointer_width = "64", target_pointer_width = "32")))]
compile_error!(
    "Only 32bit and 64 bit targets are supported as usize is assumed to be >= u32 in various places"
);
