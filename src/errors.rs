//! Custom error types for DNAComb
//!
//! Defines error types for read counting, LibSpec and read parsing
use bio::bio_types::alignment::Alignment;
use bio::bio_types::sequence::Sequence;
use bio::io::fastq;
use log::warn;
use std::fmt;
use std::io;

use crate::interning::RegionID;
use crate::interning::region_id_to_str;
use crate::region::RegionCompleteness;

/// Convert a `Vec<u8>` Sequence to a string, logging failure but not panicing
///
/// This is useful for writing output files so that bad UTF8 is flagged but doesn't
/// abort the whole write, meaning the user can more easily observe what has occured
/// in combination with the warnings. In theory this should rarely occur with good input
/// and bad input should be caught earlier.
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

/// Error type for read counting
///
/// Mostly the generic ReadCountError since the CLI doesn't need to differentiate much.
/// UnexpectedRegionError is included for ergonomics and clarity.
#[derive(Debug)]
pub enum ReadCountError {
    UnexpectedRegion { region: RegionID },
    FilterConfigError { desc: String },
    BadAlignment { alignment: Box<AlignmentInfo> },
    Error { desc: String },
}

#[derive(Debug)]
pub struct AlignmentInfo {
    pub read_id: String,
    pub read_number: usize,
    pub pretty_alignment: String,
    pub alignment: Alignment,
    pub region_ids: Vec<RegionID>,
    pub region_positions: Vec<(usize, usize)>,
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

/// Error type for LibSpec
///
/// Includes a range of possible errors and wraps downstream errors from other
/// modules.
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

    /// Required region missing
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

/// Error type for library
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

/// Error type for individual reads
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

/// Error type for read pairs
#[derive(Debug)]
pub enum ReadPairError {
    ReadPair {
        forward: Option<FastaError>,
        reverse: Option<FastaError>,
    },
    Format {
        desc: String,
    },
    EarlyExhastion {
        read: String,
    },
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
