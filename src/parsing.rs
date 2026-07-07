//! Unified parsing of single-end and paired-end FASTA/FASTQ inputs.
//!
//! This module normalises sequence-file input into a common `ReadPair` stream,
//! regardless of whether the source data are FASTA or FASTQ, compressed or
//! uncompressed, single-end or paired-end.
//!
//! FASTA input is converted into FASTQ-like records by assigning a default
//! quality score to each base so downstream code can operate on a single record type.
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
use crate::groups::ReadGroup;
use crate::seqs::ReadPair;

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

/// Convert a FASTA record into a synthetic FASTQ record.
///
/// Because downstream processing expects FASTQ-style records, FASTA input is
/// normalised by assigning the same default quality byte to every base.
fn fasta_to_fastq(fasta_record: fasta::Record, default_quality: u8) -> fastq::Record {
    let id = fasta_record.id().to_string();
    let desc = fasta_record.desc().map(|d| d.to_string());
    let seq = fasta_record.seq().to_vec();
    let qual = vec![default_quality; seq.len()]; // Assign default quality for each base

    fastq::Record::with_attrs(&id, desc.as_deref(), &seq, &qual)
}

/// Common interface for sources that yield `ReadPair` records.
///
/// This trait abstracts over direct file parsing and threaded read delivery so
/// counting code can consume either through the same interface.
pub trait ReadPairProducer: Iterator<Item = Result<ReadPair, ReadPairError>> {
    /// Return `true` if produced records may include reverse reads.
    fn has_reverse(&self) -> bool;

    /// Return the grouping regex used to assign read groups, if any.
    fn group(&self) -> &Option<Regex>;

    /// Maximum number of reads configured for this producer, or `0` for no limit.
    fn max_reads(&self) -> u64;

    /// Number of reads yielded so far.
    fn read_count(&self) -> u64;
}

/// Parser that yields normalised `ReadPair` records from one or two sequence files.
///
/// Internally, forward and reverse inputs may be FASTA or FASTQ and may be
/// compressed or uncompressed. All records are converted into FASTQ-style
/// records before being assembled into `ReadPair`s.
///
/// For paired-end input, forward and reverse files are consumed in lockstep and
/// must remain synchronised. Grouping is derived from the forward read header only.
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
    /// Construct a parser from already-initialised forward and optional reverse iterators.
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

    /// Construct a `ReadPairParser` from forward and optional reverse file paths.
    ///
    /// File format and compression may be explicitly specified or auto-detected
    /// from the path. FASTA input is converted into FASTQ-style records using
    /// `default_quality`.
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

    /// Assign a read group from the forward read header.
    ///
    /// If no grouping regex is configured, reads are assigned to the `ungrouped`
    /// sentinel group. If a regex is configured but no first capture group is
    /// found, reads are assigned to the `unmatched` sentinel group.
    fn read_group(&mut self, f_record: &fastq::Record) -> ReadGroup {
        let re = match &self.group {
            None => return ReadGroup::ungrouped(),
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
            None => ReadGroup::unmatched(),
            Some(cap) => match cap.get(1) {
                None => ReadGroup::unmatched(),
                Some(x) => ReadGroup::grouped(x.as_str()),
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

/// Yields one `ReadPair` at a time, enforcing paired-file synchronisation when
/// reverse reads are present.
impl Iterator for ReadPairParser {
    type Item = Result<ReadPair, ReadPairError>;

    fn next(&mut self) -> Option<Self::Item> {
        if (self.max_reads > 0) && (self.read_count == self.max_reads) {
            return None;
        }

        let f_record: Option<Result<fastq::Record, FastaError>> = self.forward.next();

        if f_record.is_some() {
            self.read_count += 1;
        }

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

/// `ReadPairProducer` implementation backed by a channel from another thread.
///
/// This is used for multithreaded counting, where one thread produces parsed
/// reads and worker threads consume them through the same `ReadPairProducer`
/// interface as direct file parsing.
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
    /// Construct a threaded read producer from a receiving channel and parser metadata.
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

        if next.is_some() {
            self.read_count += 1;
        }

        next
    }
}

/// Path plus parsing metadata for a sequence file.
///
/// A `SeqPath` stores the file path together with either explicit or auto-detected
/// format/compression settings, and can open the file as a normalised record iterator.
#[derive(Debug)]
pub struct SeqPath {
    path: String,
    format: SeqFormat,
    gzip: Compression,
}

impl SeqPath {
    /// Construct a new sequence-file descriptor.
    pub fn new(path: String, format: SeqFormat, gzip: Compression) -> Self {
        SeqPath { path, format, gzip }
    }

    /// Open the sequence file and return an iterator over normalised FASTQ records.
    ///
    /// Format and compression are either taken directly from the stored settings
    /// or auto-detected from the file path. FASTA input is converted into FASTQ
    /// records using `default_quality`.
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

/// Sequence file format.
#[derive(Clone, ValueEnum, Debug, Copy, PartialEq)]
pub enum SeqFormat {
    /// Detect format from the file extension.
    Auto,
    /// Parse as FASTA.
    Fasta,
    /// Parse as FASTQ.
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

/// Compression mode for sequence-file input.
#[derive(Clone, ValueEnum, Debug, Copy, PartialEq, Eq)]
pub enum Compression {
    /// Detect compression from the file extension.
    Auto,
    /// Treat the file as gzip-compressed.
    Gzip,
    /// Treat the file as uncompressed.
    None,
}

impl fmt::Display for Compression {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Auto => write!(f, "Auto detect Compression"),
            Self::Gzip => write!(f, "Gzip"),
            Self::None => write!(f, "None"),
        }
    }
}

/// Detect compression format from the file extension.
///
/// Only Gzip with extension ".gz" ia supported currently.
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
    use regex::Regex;

    // Format detection
    #[test]
    fn test_seq_detection_fasta_variants() {
        assert_eq!(detect_seq_format("file.fa").unwrap(), SeqFormat::Fasta);
        assert_eq!(detect_seq_format("file.fasta").unwrap(), SeqFormat::Fasta);
        assert_eq!(detect_seq_format("file.fa.gz").unwrap(), SeqFormat::Fasta);
        assert_eq!(
            detect_seq_format("file.fasta.gz").unwrap(),
            SeqFormat::Fasta
        );
    }

    #[test]
    fn test_seq_detection_fastq_variants() {
        assert_eq!(detect_seq_format("file.fq").unwrap(), SeqFormat::Fastq);
        assert_eq!(detect_seq_format("file.fastq").unwrap(), SeqFormat::Fastq);
        assert_eq!(detect_seq_format("file.fq.gz").unwrap(), SeqFormat::Fastq);
        assert_eq!(
            detect_seq_format("file.fastq.gz").unwrap(),
            SeqFormat::Fastq
        );
    }

    #[test]
    fn test_seq_detection_with_paths() {
        assert_eq!(
            detect_seq_format("path/to/file.fa").unwrap(),
            SeqFormat::Fasta
        );
        assert_eq!(
            detect_seq_format("/absolute/path/file.fastq").unwrap(),
            SeqFormat::Fastq
        );
        assert_eq!(
            detect_seq_format("../relative/file.fq.gz").unwrap(),
            SeqFormat::Fastq
        );
    }

    #[test]
    fn test_seq_detection_case_sensitive() {
        // Extensions are case-sensitive
        assert!(detect_seq_format("file.FA").is_err());
        assert!(detect_seq_format("file.Fasta").is_err());
        assert!(detect_seq_format("file.FQ").is_err());
    }

    #[test]
    fn test_seq_detection_invalid_extensions() {
        assert!(detect_seq_format("file.txt").is_err());
        assert!(detect_seq_format("file.seq").is_err());
        assert!(detect_seq_format("file.gz").is_err());
        assert!(detect_seq_format("file").is_err());
    }

    #[test]
    fn test_seq_detection_multiple_dots() {
        assert_eq!(
            detect_seq_format("file.backup.fa.gz").unwrap(),
            SeqFormat::Fasta
        );
    }

    // Compression detection
    #[test]
    fn test_gzip_detection_gzip() {
        assert_eq!(detect_gzip("file.fa.gz"), Compression::Gzip);
        assert_eq!(detect_gzip("file.txt.gz"), Compression::Gzip);
    }

    #[test]
    fn test_gzip_detection_none() {
        assert_eq!(detect_gzip("file.fa"), Compression::None);
        assert_eq!(detect_gzip("file.txt"), Compression::None);
    }

    #[test]
    fn test_gzip_detection_case_sensitive() {
        assert_eq!(detect_gzip("file.GZ"), Compression::None);
    }

    #[test]
    fn test_gzip_detection_multiple_gz() {
        assert_eq!(detect_gzip("file.fa.gz.gz"), Compression::Gzip);
    }

    // SeqPath
    #[test]
    fn seqpath_new_auto_format_auto_compression() {
        let sp = SeqPath::new("test.fa.gz".to_string(), SeqFormat::Auto, Compression::Auto);
        assert_eq!(sp.format, SeqFormat::Auto);
        assert_eq!(sp.gzip, Compression::Auto);
        assert_eq!(sp.path, "test.fa.gz");
    }

    #[test]
    fn seqpath_new_explicit_format_and_compression() {
        let sp = SeqPath::new("test.txt".to_string(), SeqFormat::Fastq, Compression::Gzip);
        assert_eq!(sp.format, SeqFormat::Fastq);
        assert_eq!(sp.gzip, Compression::Gzip);
    }

    #[test]
    fn seqpath_display() {
        let sp = SeqPath::new(
            "test.fa.gz".to_string(),
            SeqFormat::Fasta,
            Compression::Gzip,
        );
        let display = format!("{}", sp);
        assert!(display.contains("test.fa.gz"));
        assert!(display.contains("Fasta"));
        assert!(display.contains("Gzip"));
    }

    // Fasta to Fastq
    #[test]
    fn fasta_to_fastq_basic() {
        let fasta = fasta::Record::with_attrs("seq1", None, b"ACGT");
        let fastq = fasta_to_fastq(fasta, b'I');

        assert_eq!(fastq.id(), "seq1");
        assert_eq!(fastq.seq(), b"ACGT");
        assert_eq!(fastq.qual(), b"IIII");
    }

    #[test]
    fn fasta_to_fastq_with_description() {
        let fasta = fasta::Record::with_attrs("seq1", Some("description text"), b"ACGT");
        let fastq = fasta_to_fastq(fasta, b'I');

        assert_eq!(fastq.id(), "seq1");
        assert_eq!(fastq.seq(), b"ACGT");
        assert_eq!(fastq.qual(), b"IIII");
    }

    #[test]
    fn fasta_to_fastq_empty_sequence() {
        let fasta = fasta::Record::with_attrs("empty", None, b"");
        let fastq = fasta_to_fastq(fasta, b'I');

        assert_eq!(fastq.id(), "empty");
        assert_eq!(fastq.seq(), b"");
        assert_eq!(fastq.qual(), b"");
    }

    #[test]
    fn fasta_to_fastq_different_quality_scores() {
        let fasta = fasta::Record::with_attrs("seq", None, b"ACGTACGT");
        let fastq_low = fasta_to_fastq(fasta.clone(), b'!');
        let fastq_mid = fasta_to_fastq(fasta.clone(), b'I');
        let fastq_high = fasta_to_fastq(fasta, b'~');

        assert_eq!(fastq_low.qual(), b"!!!!!!!!");
        assert_eq!(fastq_mid.qual(), b"IIIIIIII");
        assert_eq!(fastq_high.qual(), b"~~~~~~~~");
    }

    #[test]
    fn fasta_to_fastq_long_sequence() {
        let seq = vec![b'A'; 10_000];
        let fasta = fasta::Record::with_attrs("long", None, &seq);
        let fastq = fasta_to_fastq(fasta, b'I');

        assert_eq!(fastq.seq().len(), 10_000);
        assert_eq!(fastq.qual().len(), 10_000);
        assert!(fastq.qual().iter().all(|&q| q == b'I'));
    }

    // ReadPair parser properties
    #[test]
    fn readpair_parser_properties() {
        let fwd_iter = vec![Ok(fastq::Record::with_attrs("r1", None, b"ACGT", b"IIII"))];
        let fwd: Box<dyn Iterator<Item = Result<fastq::Record, FastaError>>> =
            Box::new(fwd_iter.into_iter());

        let parser = ReadPairParser::new(fwd, None, None, 100);
        assert!(!parser.has_reverse());
        assert_eq!(parser.max_reads(), 100);
        assert_eq!(parser.read_count(), 0);
        assert!(parser.group().is_none());
    }

    #[test]
    fn readpair_parser_with_grouping() {
        let re = Regex::new(r"([A-Z0-9]+)").unwrap();
        let fwd_iter = vec![Ok(fastq::Record::with_attrs(
            "GROUPX_read1",
            None,
            b"ACGT",
            b"IIII",
        ))];
        let fwd: Box<dyn Iterator<Item = Result<fastq::Record, FastaError>>> =
            Box::new(fwd_iter.into_iter());

        let parser = ReadPairParser::new(fwd, None, Some(re.clone()), 100);
        assert!(parser.group().is_some());
        assert_eq!(parser.group().as_ref().unwrap().as_str(), re.as_str());
    }

    // ReadPair parser grouping
    #[test]
    fn readpair_parser_group_ungrouped() {
        let mut parser = ReadPairParser::new(Box::new(std::iter::empty()), None, None, 0);

        let record = fastq::Record::with_attrs("read1", None, b"ACGT", b"IIII");
        let group = parser.read_group(&record);
        assert!(group.is_ungrouped());
    }

    #[test]
    fn readpair_parser_group_with_regex_match() {
        let re = Regex::new(r"([A-Z]+)_").unwrap();
        let mut parser = ReadPairParser::new(Box::new(std::iter::empty()), None, Some(re), 0);

        let record = fastq::Record::with_attrs("CELLTYPE_read1", None, b"ACGT", b"IIII");
        let group = parser.read_group(&record);
        assert!(group.is_match());
    }

    #[test]
    fn readpair_parser_group_with_regex_no_match() {
        let re = Regex::new(r"([0-9]+)").unwrap();
        let mut parser = ReadPairParser::new(Box::new(std::iter::empty()), None, Some(re), 0);

        let record = fastq::Record::with_attrs("readABC", None, b"ACGT", b"IIII");
        let group = parser.read_group(&record);
        assert!(group.is_unmatched());
    }

    #[test]
    fn readpair_parser_group_with_description() {
        let re = Regex::new(r"group=([A-Za-z0-9]+)").unwrap();
        let mut parser = ReadPairParser::new(Box::new(std::iter::empty()), None, Some(re), 0);

        let record =
            fastq::Record::with_attrs("read1", Some("group=mygroup extra"), b"ACGT", b"IIII");
        let group = parser.read_group(&record);
        assert!(group.is_match());
    }

    #[test]
    fn readpair_parser_group_haystack_reuse() {
        let re = Regex::new(r"([A-Z]+)").unwrap();
        let mut parser = ReadPairParser::new(Box::new(std::iter::empty()), None, Some(re), 0);

        let r1 = fastq::Record::with_attrs("ABC_read", None, b"ACGT", b"IIII");
        let r2 = fastq::Record::with_attrs("XYZ_read", None, b"ACGT", b"IIII");

        parser.read_group(&r1);
        let after_first = parser.group_haystack.capacity();
        parser.read_group(&r2);
        let after_second = parser.group_haystack.capacity();

        // Haystack capacity should stabilize after first use
        assert_eq!(after_first, after_second);
    }

    // ReadPair Parser iteration
    #[test]
    fn readpair_parser_max_reads_limit() {
        let fwd_iter = vec![
            Ok(fastq::Record::with_attrs("r1", None, b"A", b"I")),
            Ok(fastq::Record::with_attrs("r2", None, b"C", b"I")),
            Ok(fastq::Record::with_attrs("r3", None, b"G", b"I")),
        ];
        let fwd: Box<dyn Iterator<Item = Result<fastq::Record, FastaError>>> =
            Box::new(fwd_iter.into_iter());

        let mut parser = ReadPairParser::new(fwd, None, None, 2);

        assert!(parser.next().is_some());
        assert!(parser.next().is_some());
        assert!(parser.next().is_none(), "should stop at max_reads");
    }

    #[test]
    fn readpair_parser_zero_max_reads() {
        let fwd_iter = vec![Ok(fastq::Record::with_attrs("r1", None, b"A", b"I"))];
        let fwd: Box<dyn Iterator<Item = Result<fastq::Record, FastaError>>> =
            Box::new(fwd_iter.into_iter());

        let mut parser = ReadPairParser::new(fwd, None, None, 0);
        assert!(parser.next().is_some(), "0 means no limit");
    }

    // Format and compression display
    #[test]
    fn seqformat_display() {
        assert_eq!(format!("{}", SeqFormat::Auto), "Auto Format");
        assert_eq!(format!("{}", SeqFormat::Fasta), "Fasta");
        assert_eq!(format!("{}", SeqFormat::Fastq), "Fastq");
    }

    #[test]
    fn compression_display() {
        assert_eq!(format!("{}", Compression::Auto), "Auto detect Compression");
        assert_eq!(format!("{}", Compression::Gzip), "Gzip");
        assert_eq!(format!("{}", Compression::None), "None");
    }

    // Threaded parser properties
    #[test]
    fn threaded_readpair_parser_properties() {
        use crossbeam::channel::bounded;

        let (_tx, rx) = bounded::<Result<ReadPair, ReadPairError>>(10);
        let parser = ThreadedReadPairParser::new(rx, true, None, 1000);

        assert!(parser.has_reverse());
        assert_eq!(parser.max_reads(), 1000);
        assert_eq!(parser.read_count(), 0);
        assert!(parser.group().is_none());
    }

    #[test]
    fn threaded_readpair_parser_with_group() {
        use crossbeam::channel::bounded;

        let re = Regex::new(r"(\w+)").unwrap();
        let (_tx, rx) = bounded::<Result<ReadPair, ReadPairError>>(10);
        let parser = ThreadedReadPairParser::new(rx, false, Some(re.clone()), 500);

        assert!(!parser.has_reverse());
        assert!(parser.group().is_some());
        assert_eq!(parser.max_reads(), 500);
    }

    // Fastq parser
    #[test]
    fn fastq_parser_iter() {
        let records = vec![
            Ok(fastq::Record::with_attrs("r1", None, b"ACGT", b"IIII")),
            Ok(fastq::Record::with_attrs("r2", None, b"GGGG", b"!!!!")),
        ];
        let mut parser = FastqParser::new(records.into_iter());

        let r1 = parser.next().unwrap().unwrap();
        assert_eq!(r1.id(), "r1");
        let r2 = parser.next().unwrap().unwrap();
        assert_eq!(r2.id(), "r2");
        assert!(parser.next().is_none());
    }

    // Fasta parser
    #[test]
    fn fasta_parser_iter_with_quality() {
        let records = vec![
            Ok(fasta::Record::with_attrs("s1", None, b"ACGT")),
            Ok(fasta::Record::with_attrs("s2", None, b"GGGG")),
        ];
        let mut parser = FastaParser::new(records.into_iter(), b'#');

        let r1 = parser.next().unwrap().unwrap();
        assert_eq!(r1.id(), "s1");
        assert_eq!(r1.qual(), b"####");

        let r2 = parser.next().unwrap().unwrap();
        assert_eq!(r2.id(), "s2");
        assert_eq!(r2.qual(), b"####");
    }
}
