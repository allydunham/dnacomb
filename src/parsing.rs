//! Parse Fasta and Fastq files with a consistent interface
//!
//! Provides a parser that can read single and paired-end fasta
//! and fastq files with a single interface, always returning a
//! standard ReadPair object.
use bio::io::{fasta, fastq};
use clap::ValueEnum;
use crossbeam::channel::Receiver;
use flate2::read::MultiGzDecoder;
use log::debug;
use regex::Regex;
use std::fmt;
use std::fs::File;
use std::io::{self, BufReader};
use std::str;

use crate::errors::{FastaError, ReadPairError};

pub type ReadKey = (Vec<u8>, Option<Vec<u8>>);

/// Pair of linked Fastq reads
#[derive(Debug)]
pub struct ReadPair {
    pub forward: fastq::Record,
    pub reverse: Option<fastq::Record>,
    pub group: ReadGroup,
}

impl ReadPair {
    /// Generate a key to identify unique read types
    pub fn key(&self) -> ReadKey {
        if self.reverse.is_some() {
            (
                self.forward.seq().to_vec(),
                Some(self.reverse.as_ref().unwrap().seq().to_vec()),
            )
        } else {
            (self.forward.seq().to_vec(), None)
        }
    }
}

/// Parser for Fastq files
///
/// Creates an iterator that yields Fastq records
struct FastqParser<R>
where
    R: Iterator<Item = Result<fastq::Record, fastq::Error>>,
{
    reads: R,
}

impl<R> FastqParser<R>
where
    R: Iterator<Item = Result<fastq::Record, fastq::Error>>,
{
    fn new(reads: R) -> Self {
        FastqParser { reads }
    }
}

impl<R> Iterator for FastqParser<R>
where
    R: Iterator<Item = Result<fastq::Record, fastq::Error>>,
{
    type Item = Result<fastq::Record, FastaError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reads.next() {
            Some(Ok(x)) => Some(Ok(x)),
            Some(Err(e)) => Some(Err(FastaError::Fastq(e))),
            None => None,
        }
    }
}

/// Parser for Fasta files
///
/// Creates an iterator yielding Fasta records
struct FastaParser<R>
where
    R: Iterator<Item = Result<fasta::Record, io::Error>>,
{
    reads: R,
    default_quality: u8,
}

impl<R> FastaParser<R>
where
    R: Iterator<Item = Result<fasta::Record, io::Error>>,
{
    fn new(reads: R, default_quality: u8) -> Self {
        FastaParser {
            reads,
            default_quality,
        }
    }
}

impl<R> Iterator for FastaParser<R>
where
    R: Iterator<Item = Result<fasta::Record, io::Error>>,
{
    type Item = Result<fastq::Record, FastaError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reads.next() {
            Some(Ok(x)) => Some(Ok(fasta_to_fastq(x, self.default_quality))),
            Some(Err(e)) => Some(Err(FastaError::Fasta(e))),
            None => None,
        }
    }
}

/// Create a mock Fastq record from a Fasta record
fn fasta_to_fastq(fasta_record: fasta::Record, default_quality: u8) -> fastq::Record {
    let id = fasta_record.id().to_string();
    let desc = fasta_record.desc().map(|d| d.to_string());
    let seq = fasta_record.seq().to_vec();
    let qual = vec![default_quality; seq.len()]; // Assign default quality for each base

    fastq::Record::with_attrs(&id, desc.as_deref(), &seq, &qual)
}

/// Group status of a read
#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub enum ReadGroup {
    Ungrouped,
    Unmatched,
    Match(String),
}

pub trait ReadPairProducer: Iterator<Item = Result<ReadPair, ReadPairError>> {
    fn has_reverse(&self) -> bool;
    fn group(&self) -> &Option<Regex>;
    fn max_reads(&self) -> u64;
    fn read_count(&self) -> u64;
}

/// Parser outputing ReadPair objects
///
/// Internally uses boxed forward and reverse parsers that may yield Fasta or Fastq
/// records then converts to Fastq format
pub struct ReadPairParser {
    /// Forward parser
    forward: Box<dyn Iterator<Item = Result<fastq::Record, FastaError>>>,

    /// Reverse parser
    reverse: Option<Box<dyn Iterator<Item = Result<fastq::Record, FastaError>>>>,

    /// Regex to process forward read names with to identify groups, for instance cells
    /// in single cell assays
    group: Option<Regex>,

    /// Maximum number of reads to process
    max_reads: u64,

    /// Number of reads processed
    read_count: u64,

    /// Implementation detail for repeated Regex search for read groups
    /// Instead of concating id/desc into a new string each time, reuse this
    /// buffer. id/desc are generally the same length so should quickly converge on
    /// a good capacity.
    group_haystack: String,
}

impl ReadPairParser {
    fn new(
        forward: Box<dyn Iterator<Item = Result<fastq::Record, FastaError>>>,
        reverse: Option<Box<dyn Iterator<Item = Result<fastq::Record, FastaError>>>>,
        group: Option<Regex>,
        max_reads: u64,
    ) -> Self {
        ReadPairParser {
            forward,
            reverse,
            group,
            max_reads,
            read_count: 0,
            group_haystack: String::with_capacity(200),
        }
    }

    /// Initialise a read parser from file paths
    pub fn from_paths(
        forward: SeqPath,
        reverse: Option<SeqPath>,
        group: Option<Regex>,
        max_reads: u64,
        default_quality: u8,
    ) -> Result<Self, ReadPairError> {
        let f_records = forward.get_records(default_quality)?;
        let r_records = match reverse {
            Some(x) => Some(x.get_records(default_quality)?),
            None => None,
        };

        Ok(ReadPairParser::new(f_records, r_records, group, max_reads))
    }

    /// Extract read group from a record
    fn read_group(&mut self, f_record: &fastq::Record) -> ReadGroup {
        let re = match &self.group {
            None => return ReadGroup::Ungrouped,
            Some(x) => x,
        };

        self.group_haystack.clear();
        self.group_haystack.push_str(f_record.id());
        match f_record.desc() {
            None => {}
            Some(x) => {
                self.group_haystack.push(' ');
                self.group_haystack.push_str(x);
            }
        }

        match re.captures(&self.group_haystack) {
            None => ReadGroup::Unmatched,
            Some(cap) => match cap.get(1) {
                None => ReadGroup::Unmatched,
                Some(x) => ReadGroup::Match(x.as_str().to_string()),
            },
        }
    }
}

impl ReadPairProducer for ReadPairParser {
    fn has_reverse(&self) -> bool {
        self.reverse.is_some()
    }

    fn group(&self) -> &Option<Regex> {
        &self.group
    }

    fn max_reads(&self) -> u64 {
        self.max_reads
    }

    fn read_count(&self) -> u64 {
        self.read_count
    }
}

impl Iterator for ReadPairParser {
    type Item = Result<ReadPair, ReadPairError>;

    fn next(&mut self) -> Option<Self::Item> {
        if (self.max_reads > 0) && (self.read_count == self.max_reads) {
            return None;
        } else {
            self.read_count += 1;
        }

        let f_record: Option<Result<fastq::Record, FastaError>> = self.forward.next();

        if self.reverse.is_none() {
            // Unpaired reads
            match f_record {
                Some(Ok(f)) => Some(Ok(ReadPair {
                    // Group first so can safely move f into ReadPair next
                    group: self.read_group(&f),
                    forward: f,
                    reverse: None,
                })),
                Some(Err(e)) => Some(Err(ReadPairError::ReadPair {
                    forward: Some(e),
                    reverse: None,
                })),
                None => None,
            }
        } else {
            // Paired reads
            let r_record: Option<Result<fastq::Record, FastaError>> =
                self.reverse.as_mut().unwrap().next();

            match (f_record, r_record) {
                // Expecgted read pair
                (Some(Ok(f)), Some(Ok(r))) => {
                    Some(Ok(ReadPair {
                        group: self.read_group(&f),
                        // Group first so can safely move f into ReadPair next
                        forward: f,
                        reverse: Some(r),
                    }))
                }

                // Files exhausted
                (None, None) => None,

                // Error combinations
                (Some(Ok(_)), Some(Err(r))) => Some(Err(ReadPairError::ReadPair {
                    forward: None,
                    reverse: Some(r),
                })),
                (Some(Err(f)), Some(Ok(_))) => Some(Err(ReadPairError::ReadPair {
                    forward: Some(f),
                    reverse: None,
                })),
                (Some(Err(f)), Some(Err(r))) => Some(Err(ReadPairError::ReadPair {
                    forward: Some(f),
                    reverse: Some(r),
                })),
                (Some(_), None) => Some(Err(ReadPairError::EarlyExhastion {
                    read: "Reverse".to_string(),
                })),
                (None, Some(_)) => Some(Err(ReadPairError::EarlyExhastion {
                    read: "Forward".to_string(),
                })),
            }
        }
    }
}

/// Multithreading wrapper for ReadPairParser
///
/// Recieves reads from a ReadPairParser on another thread to process them in parallel
pub struct ThreadedReadPairParser {
    // Channel receiving new reads
    rx: Receiver<Result<ReadPair, ReadPairError>>,

    // Whether the incoming read pairs have reverse reads
    rev_reads: bool,

    // Regex used to parse grouping structure
    group: Option<Regex>,

    /// Maximum number of reads potentially coming from the channel
    max_reads: u64,

    /// Number of reads received on this channel
    read_count: u64,
}

impl ThreadedReadPairParser {
    pub fn new(
        rx: Receiver<Result<ReadPair, ReadPairError>>,
        rev_reads: bool,
        group: Option<Regex>,
        max_reads: u64,
    ) -> Self {
        ThreadedReadPairParser {
            rx,
            rev_reads,
            group,
            max_reads,
            read_count: 0,
        }
    }
}

impl ReadPairProducer for ThreadedReadPairParser {
    fn has_reverse(&self) -> bool {
        self.rev_reads
    }

    fn group(&self) -> &Option<Regex> {
        &self.group
    }

    fn max_reads(&self) -> u64 {
        self.max_reads
    }

    fn read_count(&self) -> u64 {
        self.read_count
    }
}

impl Iterator for ThreadedReadPairParser {
    type Item = Result<ReadPair, ReadPairError>;

    fn next(&mut self) -> Option<Self::Item> {
        let next = self.rx.recv().ok();
        self.read_count += 1;
        next
    }
}

/// Path to a sequence file
///
/// This struct implements automatic detection of filetype
/// and compression status for more ergonomic parsing in ReadPairParser
#[derive(Debug)]
pub struct SeqPath {
    path: String,
    format: SeqFormat,
    gzip: Compression,
}

impl SeqPath {
    pub fn new(path: String, format: SeqFormat, gzip: Compression) -> Self {
        SeqPath { path, format, gzip }
    }

    /// Process path to produce the appropriate Fasta/Fastq parser
    fn get_records(
        &self,
        default_quality: u8,
    ) -> Result<Box<dyn Iterator<Item = Result<fastq::Record, FastaError>>>, ReadPairError> {
        let fmt = match self.format {
            SeqFormat::Auto => detect_seq_format(&self.path)?,
            SeqFormat::Fasta | SeqFormat::Fastq => self.format,
        };

        let gzip = match self.gzip {
            Compression::Auto => detect_gzip(&self.path),
            Compression::Gzip | Compression::None => self.gzip,
        };
        debug!(
            "Using format {} and compression {} for {}",
            fmt, gzip, self.path
        );

        let reader: BufReader<File> = BufReader::new(File::open(&self.path)?);

        match (gzip, fmt) {
            (Compression::Auto, _) | (_, SeqFormat::Auto) => Err(ReadPairError::Format {
                desc: "Gzip::Auto or SeqFormat::Auto remained after parsing".to_string(),
            }),
            (Compression::None, SeqFormat::Fasta) => {
                let fasta_reader = fasta::Reader::from_bufread(reader);
                Ok(Box::new(FastaParser::new(
                    fasta_reader.records(),
                    default_quality,
                )))
            }
            (Compression::None, SeqFormat::Fastq) => {
                let fastq_reader = fastq::Reader::from_bufread(reader);
                Ok(Box::new(FastqParser::new(fastq_reader.records())))
            }
            (Compression::Gzip, SeqFormat::Fasta) => {
                // We use MultiGzDecoder here as in some cases seq files have multiple blocks and in others
                // a single one, this protects against this, although could give weird results if you abuse it
                // with an strange multi-file gzipped seq.
                let gz_decoder = BufReader::new(MultiGzDecoder::new(reader));
                let fasta_reader = fasta::Reader::from_bufread(gz_decoder);
                Ok(Box::new(FastaParser::new(
                    fasta_reader.records(),
                    default_quality,
                )))
            }
            (Compression::Gzip, SeqFormat::Fastq) => {
                let gz_decoder = BufReader::new(MultiGzDecoder::new(reader));
                let fastq_reader = fastq::Reader::from_bufread(gz_decoder);
                Ok(Box::new(FastqParser::new(fastq_reader.records())))
            }
        }
    }
}

impl fmt::Display for SeqPath {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} ({}, {})", self.path, self.format, self.gzip)
    }
}

/// Sequence file format
#[derive(Clone, ValueEnum, Debug, Copy, PartialEq)]
pub enum SeqFormat {
    Auto,
    Fasta,
    Fastq,
}

impl fmt::Display for SeqFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Auto => write!(f, "Auto Format"),
            Self::Fasta => write!(f, "Fasta"),
            Self::Fastq => write!(f, "Fastq"),
        }
    }
}

/// Detect file format from a string path
fn detect_seq_format(path: &str) -> Result<SeqFormat, ReadPairError> {
    if str::ends_with(path, ".fa")
        || str::ends_with(path, ".fa.gz")
        || str::ends_with(path, ".fasta")
        || str::ends_with(path, ".fasta.gz")
    {
        Ok(SeqFormat::Fasta)
    } else if str::ends_with(path, ".fq")
        || str::ends_with(path, ".fq.gz")
        || str::ends_with(path, ".fastq")
        || str::ends_with(path, ".fastq.gz")
    {
        Ok(SeqFormat::Fastq)
    } else {
        Err(ReadPairError::Format {
            desc: format!(
                "Can't auto-detect format of {path} (assumes .fa/.fasta \
             or .fq/fastq ending with optional .gz)"
            ),
        })
    }
}

/// File compression status
#[derive(Clone, ValueEnum, Debug, Copy)]
pub enum Compression {
    Auto,
    Gzip,
    None,
}

impl fmt::Display for Compression {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Auto => write!(f, "Auto Compression"),
            Self::Gzip => write!(f, "Gzip"),
            Self::None => write!(f, "Fastq"),
        }
    }
}

/// Detect gzip status from a string path
fn detect_gzip(path: &str) -> Compression {
    if str::ends_with(path, ".gz") {
        Compression::Gzip
    } else {
        Compression::None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seq_detection() {
        assert_eq!(
            detect_seq_format("path/file.fa.gz").unwrap_or(SeqFormat::Auto),
            SeqFormat::Fasta
        );
        assert_eq!(
            detect_seq_format("path/file.fasta.gz").unwrap_or(SeqFormat::Auto),
            SeqFormat::Fasta
        );
        assert_eq!(
            detect_seq_format("path/file.fa").unwrap_or(SeqFormat::Auto),
            SeqFormat::Fasta
        );
        assert_eq!(
            detect_seq_format("path/file.fasta").unwrap_or(SeqFormat::Auto),
            SeqFormat::Fasta
        );
        assert_eq!(
            detect_seq_format("path/file.fq.gz").unwrap_or(SeqFormat::Auto),
            SeqFormat::Fastq
        );
        assert_eq!(
            detect_seq_format("path/file.fastq.gz").unwrap_or(SeqFormat::Auto),
            SeqFormat::Fastq
        );
        assert_eq!(
            detect_seq_format("path/file.fq").unwrap_or(SeqFormat::Auto),
            SeqFormat::Fastq
        );
        assert_eq!(
            detect_seq_format("path/file.fastq").unwrap_or(SeqFormat::Auto),
            SeqFormat::Fastq
        );
        assert!(detect_seq_format("path/file.not_fasta_ext").is_err());
    }
}
