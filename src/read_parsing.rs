//! Parse Fasta and Fastq files with a consistent interface
//!
//! Provides a parser that can read single and paired-end fasta
//! and fastq files with a single interface, always returning a
//! standard ReadPair object.
use bio::io::{fasta, fastq};
use clap::ValueEnum;
use flate2::read::MultiGzDecoder;
use log::debug;
use regex::Regex;
use std::fmt;
use std::fs::File;
use std::io::{self, BufReader};
use std::str;

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

/// Error type for sequences
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
pub enum Error {
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

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::ReadPair { forward, reverse } => match (forward, reverse) {
                (Some(forward), Some(reverse)) => write!(
                    f,
                    "Error in both reads.\nForward: {}\nReverse: {}",
                    forward, reverse
                ),
                (Some(forward), None) => write!(f, "Error in forward read: {}", forward),
                (None, Some(reverse)) => write!(f, "Error in reverse read: {}", reverse),
                (None, None) => write!(f, "Unknown read parsing error"),
            },
            Error::Format { desc } => write!(f, "{}", desc),
            Error::EarlyExhastion { read } => {
                write!(f, "Paired reads out of sync: {} exhausted first", read)
            }
            Error::IO(e) => write!(f, "{}", e),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None // No underlying error
    }
}

impl From<io::Error> for Error {
    fn from(err: io::Error) -> Error {
        Error::IO(err)
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
    pub group: Option<Regex>,

    /// Maximum number of reads to process
    pub max_reads: u64,

    /// Number of reads processed
    pub read_count: u64,

    /// Implementation detail for repeated Regex search for read groups
    /// Instead of concating id/desc into a new string each time, reuse this
    /// buffer. id/desc are generally the same length so should quickly converge on
    /// a good capacity.
    group_haystack: String,
}

/// Group status of a read
#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub enum ReadGroup {
    Ungrouped,
    Unmatched,
    Match(String),
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
    ) -> Result<Self, Error> {
        let f_records = forward.get_records(default_quality)?;
        let r_records = match reverse {
            Some(x) => Some(x.get_records(default_quality)?),
            None => None,
        };

        Ok(ReadPairParser::new(f_records, r_records, group, max_reads))
    }

    /// Whether the parser includes reverse reads
    pub fn has_reverse(&self) -> bool {
        self.reverse.is_some()
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

impl Iterator for ReadPairParser {
    type Item = Result<ReadPair, Error>;

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
                Some(Err(e)) => Some(Err(Error::ReadPair {
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
                (Some(Ok(_)), Some(Err(r))) => Some(Err(Error::ReadPair {
                    forward: None,
                    reverse: Some(r),
                })),
                (Some(Err(f)), Some(Ok(_))) => Some(Err(Error::ReadPair {
                    forward: Some(f),
                    reverse: None,
                })),
                (Some(Err(f)), Some(Err(r))) => Some(Err(Error::ReadPair {
                    forward: Some(f),
                    reverse: Some(r),
                })),
                (Some(_), None) => Some(Err(Error::EarlyExhastion {
                    read: "Reverse".to_string(),
                })),
                (None, Some(_)) => Some(Err(Error::EarlyExhastion {
                    read: "Forward".to_string(),
                })),
            }
        }
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
    ) -> Result<Box<dyn Iterator<Item = Result<fastq::Record, FastaError>>>, Error> {
        let fmt = match self.format {
            SeqFormat::Auto => detect_seq_format(&self.path)?,
            SeqFormat::Fasta | SeqFormat::Fastq => self.format,
        };

        let gzip = match self.gzip {
            Compression::Auto => detect_gzip(&self.path),
            Compression::Gzip | Compression::None => self.gzip,
        };
        debug!("Using format {} and compression {} for {}", fmt, gzip, self.path);

        let reader: BufReader<File> = BufReader::new(File::open(&self.path)?);

        match (gzip, fmt) {
            (Compression::Auto, _) | (_, SeqFormat::Auto) => Err(Error::Format {
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
#[derive(Clone, ValueEnum, Debug, Copy)]
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
fn detect_seq_format(path: &str) -> Result<SeqFormat, Error> {
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
        Err(Error::Format {
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
