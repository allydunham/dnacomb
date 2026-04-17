//! Error and diagnostic types used throughout DNAComb.
//!
//! This module defines the main error enums used for counting, LibSpec parsing,
//! library import, and sequence-file parsing, as well as helper utilities for
//! converting sequence bytes into displayable strings for logs and output.
use bio::bio_types::alignment::Alignment;
use bio::bio_types::sequence::Sequence;
use bio::io::fastq;
use log::warn;
use std::fmt;
use std::io;

use crate::interning::RegionID;
use crate::interning::region_id_to_str;
use crate::region::RegionCompleteness;

/// Convert a Sequence `Vec<u8>` to a UTF-8 string for display/output.
///
/// If conversion fails, a warning is logged and an empty string is returned
/// rather than aborting processing. This is intended for diagnostics and TSV
/// writing, where best-effort output is preferable to panicking on unexpected
/// non-UTF-8 sequence content. UTF-8 errors should be very rare for normal input
/// and will generally be caught earlier.
pub fn seq_to_string_or_log(seq: &Sequence) -> String {
    match std::str::from_utf8(seq) {
        Ok(i) => i.into(),
        Err(_) => {
            warn!(
                "Error converting Vec<u8> Sequence {:?} to String via UTF-8",
                seq
            );
            "".to_string()
        }
    }
}

/// Error type for read counting and region extraction.
///
/// This is the main operational error type used during counting. It covers
/// invalid region structure, filter-configuration problems, unexpected alignment
/// failures, and generic counting errors.
#[derive(Debug)]
pub enum ReadCountError {
    UnexpectedRegion { region: RegionID },
    FilterConfigError { desc: String },
    BadAlignment { alignment: Box<AlignmentInfo> },
    Error { desc: String },
}

/// Detailed debugging information for an alignment/extraction failure.
///
/// This is attached to `ReadCountError::BadAlignment` to help diagnose cases
/// where alignment succeeded but region extraction from the alignment path
/// produced inconsistent or invalid coordinates.
#[derive(Debug)]
pub struct AlignmentInfo {
    /// Read name
    pub read_id: String,

    /// Read number in input file
    pub read_number: usize,

    /// Alignment string
    pub pretty_alignment: String,

    /// Alignment object
    pub alignment: Alignment,

    /// Vector of region names being matched
    pub region_ids: Vec<RegionID>,

    /// Positions of regions in the template sequence
    pub region_positions: Vec<(usize, usize)>,

    /// Vector of identified region positions positions in the input read
    pub mapped_positions: Vec<Option<(usize, usize, RegionCompleteness)>>,
}

impl fmt::Display for ReadCountError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ReadCountError::UnexpectedRegion { region } => {
                write!(
                    f,
                    "Added combination contains an unexpected region: {}",
                    region_id_to_str(*region)
                )
            }
            ReadCountError::BadAlignment { alignment } => {
                write!(
                    f,
                    "Alignment or region extraction error\nRead {}, id: {}\nAlignment:\n{}\n\
                     Path:\n{:?}\n\nCigar: {}\nScore: {:?}\nRegions: {:?}\nRegion positions: {:?}\n\
                     Extracted positions: {:?}",
                    alignment.read_number,
                    alignment.read_id,
                    alignment.pretty_alignment,
                    alignment.alignment.path(),
                    alignment.alignment.cigar(false),
                    alignment.alignment.score,
                    alignment.region_ids,
                    alignment.region_positions,
                    alignment.mapped_positions
                )
            }
            ReadCountError::Error { desc } => {
                write!(f, "{}", desc)
            }
            ReadCountError::FilterConfigError { desc } => {
                write!(f, "{}", desc)
            }
        }
    }
}

impl std::error::Error for ReadCountError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None // No underlying error
    }
}

/// Error type for LibSpec parsing and validation.
///
/// Covers JSON parsing, file I/O, and logical validation errors in the sequence
/// specification, such as duplicate regions, invalid lengths, or unsupported
/// region layouts.
#[derive(Debug)]
pub enum LibSpecError {
    /// Generic LibSpec error
    LibSpec { desc: String },

    /// One or more errors invalidating a library, returned from .validate()
    InvalidLibSpec { errs: Vec<String> },

    /// A region has min length greater than max length
    MinGreaterThanMax {
        id: RegionID,
        min: usize,
        max: usize,
    },

    /// Duplicate regions in library
    DuplicateRegion { id: RegionID },

    /// Required region missing
    MissingRegion { id: RegionID },

    /// Two variable regions appear consecutively without a fixed anchor region.
    NeighbouringVariable { id: RegionID },

    /// IO errors
    IOError(io::Error),

    /// JSON Error
    ParsingError(serde_json::Error),
}

impl fmt::Display for LibSpecError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            LibSpecError::InvalidLibSpec { errs } => {
                writeln!(f, "Multiple LibSpec errors detected:")?;
                for err in errs {
                    writeln!(f, "{}", err)?;
                }
                Ok(())
            }
            LibSpecError::MinGreaterThanMax { id, min, max } => {
                write!(
                    f,
                    "Region {}: min_length ({}) cannot be greater than max_length ({})",
                    region_id_to_str(*id),
                    min,
                    max
                )
            }
            LibSpecError::DuplicateRegion { id } => {
                write!(
                    f,
                    "Duplciated region id {} in LibSpec",
                    region_id_to_str(*id)
                )
            }
            LibSpecError::MissingRegion { id } => {
                write!(
                    f,
                    "{} not found in LibSpec Region list",
                    region_id_to_str(*id)
                )
            }
            LibSpecError::LibSpec { desc } => write!(f, "{}", desc),
            LibSpecError::NeighbouringVariable { id } => {
                write!(
                    f,
                    "Variable region {} follows another variable region",
                    region_id_to_str(*id)
                )
            }
            LibSpecError::IOError(e) => write!(f, "Error reading LibSpec JSON file: {}", e),
            LibSpecError::ParsingError(e) => write!(f, "Error parsing LibSpec JSON: {}", e),
        }
    }
}

impl std::error::Error for LibSpecError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None // No underlying error
    }
}

impl From<io::Error> for LibSpecError {
    fn from(err: io::Error) -> LibSpecError {
        LibSpecError::IOError(err)
    }
}

impl From<serde_json::Error> for LibSpecError {
    fn from(err: serde_json::Error) -> LibSpecError {
        LibSpecError::ParsingError(err)
    }
}

/// Error type for expected-library import and lookup setup.
///
/// Covers malformed library TSV input, duplicate or missing region definitions,
/// and incompatibilities between a library TSV and the corresponding LibSpec.
#[derive(Debug)]
pub enum LibraryError {
    /// Generic Library error
    Library { desc: String },

    /// Duplicate regions in sub-library
    DuplicateSubLibraryRegion { id: RegionID },

    /// Duplicate regions in library
    DuplicateRegion { id: RegionID },

    /// Required region missing
    MissingRegion { id: RegionID },

    /// IO errors
    IOError(csv::Error),
}

impl fmt::Display for LibraryError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            LibraryError::DuplicateSubLibraryRegion { id } => {
                write!(
                    f,
                    "Region id {} is found in multiple Libraries",
                    region_id_to_str(*id)
                )
            }
            LibraryError::DuplicateRegion { id } => {
                write!(
                    f,
                    "Duplicated region id {} in Library",
                    region_id_to_str(*id)
                )
            }
            LibraryError::MissingRegion { id } => {
                write!(
                    f,
                    "{} not found in Library Region list",
                    region_id_to_str(*id)
                )
            }
            LibraryError::Library { desc } => write!(f, "{}", desc),
            LibraryError::IOError(e) => write!(f, "Error reading Library TSV file: {}", e),
        }
    }
}

impl std::error::Error for LibraryError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None // No underlying error
    }
}

impl From<csv::Error> for LibraryError {
    fn from(err: csv::Error) -> LibraryError {
        LibraryError::IOError(err)
    }
}

/// Error while reading an individual sequence record from FASTA or FASTQ input.
///
/// Think interface for Rust Bio errors.
#[derive(Debug)]
pub enum FastaError {
    Fasta(io::Error),
    Fastq(fastq::Error),
}

impl fmt::Display for FastaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FastaError::Fasta(e) => write!(f, "{}", e),
            FastaError::Fastq(e) => write!(f, "{}", e),
        }
    }
}

impl std::error::Error for FastaError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None // No underlying error
    }
}

impl From<io::Error> for FastaError {
    fn from(err: io::Error) -> FastaError {
        FastaError::Fasta(err)
    }
}

impl From<fastq::Error> for FastaError {
    fn from(err: fastq::Error) -> FastaError {
        FastaError::Fastq(err)
    }
}

/// Error while reading or constructing a forward/reverse read pair.
///
/// This includes per-record parsing failures, file-format problems, paired-file
/// synchronisation issues, and lower-level I/O errors.
#[derive(Debug)]
pub enum ReadPairError {
    /// Forward and/or reverse record parsing failed for the current pair.
    ReadPair {
        forward: Option<FastaError>,
        reverse: Option<FastaError>,
    },
    /// Input file format or auto-detection was invalid.
    Format { desc: String },
    /// One paired-end file ended before the other.
    EarlyExhastion { read: String },
    /// Generic file IO error.
    IO(io::Error),
}

impl fmt::Display for ReadPairError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ReadPairError::ReadPair { forward, reverse } => match (forward, reverse) {
                (Some(forward), Some(reverse)) => write!(
                    f,
                    "Error in both reads.\nForward: {}\nReverse: {}",
                    forward, reverse
                ),
                (Some(forward), None) => write!(f, "Error in forward read: {}", forward),
                (None, Some(reverse)) => write!(f, "Error in reverse read: {}", reverse),
                (None, None) => write!(f, "Unknown read parsing error"),
            },
            ReadPairError::Format { desc } => write!(f, "{}", desc),
            ReadPairError::EarlyExhastion { read } => {
                write!(f, "Paired reads out of sync: {} exhausted first", read)
            }
            ReadPairError::IO(e) => write!(f, "{}", e),
        }
    }
}

impl std::error::Error for ReadPairError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None // No underlying error
    }
}

impl From<io::Error> for ReadPairError {
    fn from(err: io::Error) -> ReadPairError {
        ReadPairError::IO(err)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seq_to_string() {
        let seq: Sequence = vec![b'A', b'C', b'G', b'T'];
        let string: String = "ACGT".to_string();
        assert_eq!(seq_to_string_or_log(&seq), string)
    }

    #[test]
    fn seq_to_string_empty() {
        let seq: Sequence = vec![];
        let string: String = "".to_string();
        assert_eq!(seq_to_string_or_log(&seq), string)
    }

    #[test]
    fn seq_to_string_warning() {
        let seq: Sequence = vec![b'A', b'C', b'G', 0xC0]; // Invalid UTF-8 byte
        let string: String = "".to_string();
        assert_eq!(seq_to_string_or_log(&seq), string)
    }
}
