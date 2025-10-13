//! Filters to identify anb count failing reads
//!
//! Provides a range of filtering critera that a function
//! processing reads can tap into via a simple config and
//! interface, easily filtering reads for a range of reasons.
//! Also makes it easy to keep track of filtered reads and
//! print a summary to file.
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

use crate::errors::{ReadCountError, seq_to_string_or_log};
use crate::parsing::{ReadKey, ReadPair};
use bio::bio_types::alignment::Alignment;

/// Calculate the mean of a fastq quality vector
pub fn mean_quality(qual: &[u8]) -> f32 {
    let total: u32 = qual.iter().fold(0, |a, e| a + *e as u32);
    total as f32 / qual.len() as f32 - 33.0 // Subtract 33 as Phred scores are shifted 33 in byte codepoints
}

#[derive(Debug, PartialEq)]
pub enum FilterReason {
    None,
    LowMeanQuality,
    BadAlignment,
    ShortRead,
    EmptyRead,
}

#[derive(Debug, Clone)]
struct FilteredCounts {
    low_mean_quality: u64,
    bad_alignment: u64,
    short_read: u64,
    empty_read: u64,
}

impl FilteredCounts {
    fn new() -> Self {
        Self {
            low_mean_quality: 0,
            bad_alignment: 0,
            short_read: 0,
            empty_read: 0,
        }
    }

    fn increment_count(&mut self, reason: &FilterReason) {
        match reason {
            FilterReason::None => {}
            FilterReason::LowMeanQuality => self.low_mean_quality += 1,
            FilterReason::BadAlignment => self.bad_alignment += 1,
            FilterReason::ShortRead => self.short_read += 1,
            FilterReason::EmptyRead => self.empty_read += 1,
        };
    }

    /// Total count of filtered reads
    pub fn total(&self) -> u64 {
        self.low_mean_quality + self.bad_alignment + self.empty_read + self.short_read
    }

    /// Merge counts from another FilteredCounts
    pub fn merge(&mut self, new_counts: FilteredCounts) {
        self.low_mean_quality += new_counts.low_mean_quality;
        self.bad_alignment += new_counts.bad_alignment;
        self.empty_read += new_counts.empty_read;
        self.short_read += new_counts.short_read;
    }

    /// Headers for the TSV line produced by to_tsv_line
    fn tsv_headers() -> String {
        "count\tproportion\tlow_mean_quality\tbad_alignment\tshort_read\tempty_read".to_string()
    }

    /// TSV line
    fn to_tsv_line(&self, total: f32) -> String {
        format!(
            "{}\t{:.4}\t{}\t{}\t{}\t{}",
            self.total(),
            self.total() as f32 / total,
            self.low_mean_quality,
            self.bad_alignment,
            self.short_read,
            self.empty_read
        )
    }
}

/// Container for filtered reads
///
/// Keeps track of counts for filtered reads during counting
#[derive(Debug, Clone)]
pub struct FilteredReads {
    config: FilterConfig,
    totals: FilteredCounts,
    counts: HashMap<ReadKey, FilteredCounts>,
}

impl FilteredReads {
    pub fn new(config: FilterConfig) -> Self {
        Self {
            config,
            totals: FilteredCounts::new(),
            counts: HashMap::new(),
        }
    }

    /// Update the count of a specifc read
    ///
    /// Counts are updated without checking the read meets the criteria. If using outside of
    /// filter_readpair/filter_alignment (for instance from a cache) be sure it is correct.
    pub fn increment_count(&mut self, read: &ReadKey, reason: &FilterReason) {
        self.totals.increment_count(reason);

        match self.counts.get_mut(read) {
            Some(counts) => counts.increment_count(reason),
            None => {
                let mut new_counts = FilteredCounts::new();
                new_counts.increment_count(reason);
                self.counts.insert(read.clone(), new_counts);
            }
        }
    }

    /// Determine if a readpair should be filtered based on the supplied config
    ///
    /// Checks whether the read should be filtered, adding it to the appropriate count if so, and
    /// returns a FilterReason determining why it was filtered.
    pub fn filter_readpair(&mut self, record: &ReadPair) -> FilterReason {
        let f_read = &record.forward;
        let r_read = match &record.reverse {
            Some(r) => Some(r),
            None => None,
        };

        // Check if reads are empty
        match (self.config.filter_empty, r_read) {
            (false, _) => {}
            (true, None) => {
                if f_read.seq().is_empty() {
                    self.increment_count(&record.key(), &FilterReason::EmptyRead);
                    return FilterReason::EmptyRead;
                }
            }
            (true, Some(r)) => {
                if f_read.seq().is_empty() || r.seq().is_empty() {
                    self.increment_count(&record.key(), &FilterReason::EmptyRead);
                    return FilterReason::EmptyRead;
                }
            }
        }

        // Check if reads reach the minimum length
        // Check if reads are empty
        match (self.config.minimum_length, r_read) {
            (None, _) => {}
            (Some(len), None) => {
                if f_read.seq().len() < len {
                    self.increment_count(&record.key(), &FilterReason::ShortRead);
                    return FilterReason::ShortRead;
                }
            }
            (Some(len), Some(r)) => {
                if f_read.seq().len() < len || r.seq().len() < len {
                    self.increment_count(&record.key(), &FilterReason::ShortRead);
                    return FilterReason::ShortRead;
                }
            }
        }

        // Check mean quality is high enough
        match (self.config.mean_quality_threshold, r_read) {
            (None, _) => {}
            (Some(q), None) => {
                if mean_quality(f_read.qual()) < q {
                    self.increment_count(&record.key(), &FilterReason::LowMeanQuality);
                    return FilterReason::LowMeanQuality;
                }
            }
            (Some(q), Some(r)) => {
                if mean_quality(f_read.qual()) < q || mean_quality(r.qual()) < q {
                    self.increment_count(&record.key(), &FilterReason::LowMeanQuality);
                    return FilterReason::LowMeanQuality;
                }
            }
        }

        FilterReason::None
    }

    /// Determine if an alignment should be filtered based on the supplied config
    ///
    /// Checks if an alignment should be filtered, adding it to the appropriate count if so, and
    /// returns a FilterReason determining why it was filtered
    pub fn filter_alignment(
        &mut self,
        record: &ReadPair,
        f_alignment: &Alignment,
        r_alignment: Option<&Alignment>,
    ) -> FilterReason {
        match &self.config.alignment_tolerance {
            None => {}
            Some(t) => {
                if f_alignment.score < t.minimum_f_score {
                    self.increment_count(&record.key(), &FilterReason::BadAlignment);
                    return FilterReason::BadAlignment;
                }

                if let Some(r_al) = r_alignment {
                    if r_al.score < t.minimum_r_score {
                        self.increment_count(&record.key(), &FilterReason::BadAlignment);
                        return FilterReason::BadAlignment;
                    }
                }
            }
        }

        FilterReason::None
    }

    /// Total count of filtered reads
    pub fn total(&self) -> u64 {
        self.totals.total()
    }

    /// Merge counts from another FilteredReads object
    pub fn merge(&mut self, new_reads: FilteredReads) -> Result<(), ReadCountError> {
        if !(self.config == new_reads.config) {
            return Err(ReadCountError::Error {
                desc: "Can't merge FilteredReads with different FilterConfigs".to_string(),
            });
        }

        self.totals.merge(new_reads.totals);

        for (key, new_counts) in new_reads.counts {
            match self.counts.get_mut(&key) {
                Some(counts) => counts.merge(new_counts),
                None => {
                    self.counts.insert(key, new_counts);
                }
            }
        }

        Ok(())
    }

    /// Write a TSV file listing the filteed reads
    ///
    /// Writes a TSV with columns for forward sequence, reverse sequence,
    /// total count, frequency, then one for each filter reason count.
    pub fn write_filter_tsv(&self, file: File, sort: bool) -> Result<(), anyhow::Error> {
        let total = self.total() as f32;

        let mut writer = BufWriter::new(file);
        let mut keys: Vec<(&ReadKey, u64)> =
            self.counts.iter().map(|x| (x.0, x.1.total())).collect();

        if sort {
            // Invert count to get desc order
            keys.sort_unstable_by_key(|x| 0 - i128::from(x.1));
        }

        // Write header
        writeln!(
            writer,
            "forward\treverse\t{}",
            FilteredCounts::tsv_headers()
        )?;

        for (key, _) in keys {
            let counts = self
                .counts
                .get(key)
                .expect("Count key from extracted key list missing from FilteredReads");

            writeln!(
                writer,
                "{}\t{}\t{}",
                seq_to_string_or_log(&key.0),
                match &key.1 {
                    Some(x) => seq_to_string_or_log(x),
                    None => "".to_string(),
                },
                counts.to_tsv_line(total),
            )?;
        }

        writer.flush()?;
        Ok(())
    }

    /// Generate a string of TSV lines representing the total filter counts
    ///
    /// Creates a TSV string with columns for filter reason, count and
    /// proportion of the supplied total read count. Mainly for use
    /// wihen outputing from ReadSummary
    pub fn to_summary_tsv_lines(&self, total: u64) -> String {
        let filtered_total = self.total();

        let mut out = String::with_capacity(300);

        out.push_str(&format!(
            "filtered\ttotal\t{}\t{:.4}\t1.000\n",
            filtered_total,
            filtered_total as f32 / total as f32,
        ));

        out.push_str(&format!(
            "filtered\tlow_mean_quality\t{}\t{:.4}\t{:.4}\n",
            self.totals.low_mean_quality,
            self.totals.low_mean_quality as f32 / total as f32,
            self.totals.low_mean_quality as f32 / filtered_total as f32,
        ));

        out.push_str(&format!(
            "filtered\tbad_alignment\t{}\t{:.4}\t{:.4}\n",
            self.totals.bad_alignment,
            self.totals.bad_alignment as f32 / total as f32,
            self.totals.bad_alignment as f32 / filtered_total as f32,
        ));

        out.push_str(&format!(
            "filtered\tempty_read\t{}\t{:.4}\t{:.4}\n",
            self.totals.empty_read,
            self.totals.empty_read as f32 / total as f32,
            self.totals.empty_read as f32 / filtered_total as f32,
        ));

        out.push_str(&format!(
            "filtered\tshort_read\t{}\t{:.4}\t{:.4}\n",
            self.totals.short_read,
            self.totals.short_read as f32 / total as f32,
            self.totals.short_read as f32 / filtered_total as f32,
        ));

        out
    }
}

/// Configuration for filtering
///
/// Instructions for how to filter reads
#[derive(Debug, Clone, PartialEq)]
pub struct FilterConfig {
    mean_quality_threshold: Option<f32>,
    alignment_tolerance: Option<AlignmentTolerance>,
    minimum_length: Option<usize>,
    filter_empty: bool,
}

impl FilterConfig {
    pub fn new(
        mean_quality_threshold: Option<f32>,
        alignment_tolerance: Option<AlignmentTolerance>,
        minimum_length: Option<usize>,
        allow_empty: bool,
    ) -> Self {
        Self {
            mean_quality_threshold,
            alignment_tolerance,
            minimum_length,
            filter_empty: allow_empty,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)]
pub struct AlignmentTolerance {
    tolerance: f32,
    expected_f_score: i32,
    expected_r_score: i32,
    minimum_f_score: i32,
    minimum_r_score: i32,
}

impl AlignmentTolerance {
    pub fn new(
        tolerance: f32,
        expected_f_score: i32,
        expected_r_score: i32,
    ) -> Result<Self, ReadCountError> {
        if !(0.0..=1.0).contains(&tolerance) {
            return Err(ReadCountError::FilterConfigError {
                desc: "Alignment tolerance must be between 0 and 1".to_string(),
            });
        };

        Ok(Self {
            tolerance,
            expected_f_score,
            expected_r_score,
            minimum_f_score: (expected_f_score as f32 * tolerance) as i32,
            minimum_r_score: (expected_r_score as f32 * tolerance) as i32,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parsing::{ReadGroup, ReadPair};
    use bio::bio_types::alignment::Alignment;

    #[test]
    fn test_mean_quality() {
        let qual = vec![b'F', b'F', b'H', b'H'];
        let mean: f32 = 38.0;
        assert_eq!(mean_quality(&qual), mean, "Mean quality incorrect")
    }

    // Test empty read filter
    #[test]
    fn test_empty_filter_single_end_denies() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("seq", None, b"", b""),
            reverse: None,
            group: ReadGroup::Ungrouped,
        };
        let key = readpair.key();

        let mut f = FilteredReads::new(FilterConfig::new(None, None, None, true));

        let out = f.filter_readpair(&readpair);

        assert_eq!(out, FilterReason::EmptyRead, "Empty read not filtered");
        assert_eq!(
            f.totals.empty_read, 1,
            "Inccorect filtered read count for empty read denied filter"
        );
        assert_eq!(
            f.total(),
            1,
            "Inccorect total filtered read count for empty read denied filter"
        );

        // counts map checks
        assert!(
            f.counts.contains_key(&key),
            "Counts hashmap missing read key for empty read denied filter"
        );
        let counts = f.counts.get(&key).unwrap();
        assert_eq!(
            counts.empty_read, 1,
            "Incorrect per-read empty_read count for empty read denied filter"
        );
        assert_eq!(
            counts.total(),
            1,
            "Incorrect per-read total count for empty read denied filter"
        );
    }

    #[test]
    fn test_empty_filter_single_end_allows() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("seq", None, b"", b""),
            reverse: None,
            group: ReadGroup::Ungrouped,
        };
        let key = readpair.key();

        let mut f = FilteredReads::new(FilterConfig::new(None, None, None, false));

        let out = f.filter_readpair(&readpair);

        assert_eq!(out, FilterReason::None, "Empty read not allowed");
        assert_eq!(
            f.totals.empty_read, 0,
            "Inccorect filtered read count for empty read allowed filter"
        );

        assert_eq!(
            f.total(),
            0,
            "Inccorect total filtered read count for empty read allowed filter"
        );

        // counts map should not get an entry for allowed read
        assert!(
            !f.counts.contains_key(&key),
            "Counts hashmap unexpectedly contains key for empty read allowed filter"
        );
    }

    #[test]
    fn test_empty_filter_paired_end_denies() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("seq", None, b"ACGT", b"FFFF"),
            reverse: Some(bio::io::fastq::Record::with_attrs("seq", None, b"", b"")),
            group: ReadGroup::Ungrouped,
        };
        let key = readpair.key();

        let mut f = FilteredReads::new(FilterConfig::new(None, None, None, true));

        let out = f.filter_readpair(&readpair);

        assert_eq!(
            out,
            FilterReason::EmptyRead,
            "Empty paired read not filtered"
        );
        assert_eq!(
            f.totals.empty_read, 1,
            "Inccorect filtered read count for empty paired read denied filter"
        );
        assert_eq!(
            f.total(),
            1,
            "Inccorect total filtered read count for empty paired read denied filter"
        );

        // counts map checks
        assert!(
            f.counts.contains_key(&key),
            "Counts hashmap missing read key for empty paired read denied filter"
        );
        let counts = f.counts.get(&key).unwrap();
        assert_eq!(
            counts.empty_read, 1,
            "Incorrect per-read empty_read count for empty paired read denied filter"
        );
        assert_eq!(
            counts.total(),
            1,
            "Incorrect per-read total count for empty paired read denied filter"
        );
    }

    #[test]
    fn test_empty_filter_paired_end_allows() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("seq", None, b"ACGT", b"FFFF"),
            reverse: Some(bio::io::fastq::Record::with_attrs("seq", None, b"", b"")),
            group: ReadGroup::Ungrouped,
        };
        let key = readpair.key();

        let mut f = FilteredReads::new(FilterConfig::new(None, None, None, false));

        let out = f.filter_readpair(&readpair);

        assert_eq!(out, FilterReason::None, "Empty paired read not allowed");
        assert_eq!(
            f.totals.empty_read, 0,
            "Inccorect filtered read count for empty paired read allowed filter"
        );
        assert_eq!(
            f.total(),
            0,
            "Inccorect total filtered read count for empty paired read allowed filter"
        );

        // counts map should not get an entry for allowed read
        assert!(
            !f.counts.contains_key(&key),
            "Counts hashmap unexpectedly contains key for empty paired read allowed filter"
        );
    }

    // Test short read filter
    #[test]
    fn test_short_filter_single_end_denies() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("seq", None, b"ACTG", b"FFFF"),
            reverse: None,
            group: ReadGroup::Ungrouped,
        };
        let key = readpair.key();

        let mut f = FilteredReads::new(FilterConfig::new(None, None, Some(10), false));

        let out = f.filter_readpair(&readpair);

        assert_eq!(out, FilterReason::ShortRead, "Short read not filtered");
        assert_eq!(
            f.totals.short_read, 1,
            "Inccorect filtered read count for short read denied filter"
        );
        assert_eq!(
            f.total(),
            1,
            "Inccorect total filtered read count for short read denied filter"
        );

        // counts map checks
        assert!(
            f.counts.contains_key(&key),
            "Counts hashmap missing read key for short read denied filter"
        );
        let counts = f.counts.get(&key).unwrap();
        assert_eq!(
            counts.short_read, 1,
            "Incorrect per-read short_read count for short read denied filter"
        );
        assert_eq!(
            counts.total(),
            1,
            "Incorrect per-read total count for short read denied filter"
        );
    }

    #[test]
    fn test_short_filter_single_end_allows() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("seq", None, b"ACTG", b"FFFF"),
            reverse: None,
            group: ReadGroup::Ungrouped,
        };
        let key = readpair.key();

        let mut f = FilteredReads::new(FilterConfig::new(None, None, None, false));

        let out = f.filter_readpair(&readpair);

        assert_eq!(out, FilterReason::None, "Short read not allowed");
        assert_eq!(
            f.totals.short_read, 0,
            "Inccorect filtered read count for short read allowed filter"
        );
        assert_eq!(
            f.total(),
            0,
            "Inccorect total filtered read count for short read allowed filter"
        );

        // counts map should not get an entry for allowed read
        assert!(
            !f.counts.contains_key(&key),
            "Counts hashmap unexpectedly contains key for short read allowed filter"
        );
    }

    #[test]
    fn test_short_filter_paired_end_denies() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs(
                "seq",
                None,
                b"AAACCCGGGTTT",
                b"FFFFFFFFFFFF",
            ),
            reverse: Some(bio::io::fastq::Record::with_attrs(
                "seq", None, b"ACTG", b"FFFF",
            )),
            group: ReadGroup::Ungrouped,
        };
        let key = readpair.key();

        let mut f = FilteredReads::new(FilterConfig::new(None, None, Some(10), false));

        let out = f.filter_readpair(&readpair);

        assert_eq!(
            out,
            FilterReason::ShortRead,
            "Short paired read not filtered"
        );
        assert_eq!(
            f.totals.short_read, 1,
            "Inccorect filtered read count for short paired read denied filter"
        );
        assert_eq!(
            f.total(),
            1,
            "Inccorect total filtered read count for short paired read denied filter"
        );

        // counts map checks
        assert!(
            f.counts.contains_key(&key),
            "Counts hashmap missing read key for short paired read denied filter"
        );
        let counts = f.counts.get(&key).unwrap();
        assert_eq!(
            counts.short_read, 1,
            "Incorrect per-read short_read count for short paired read denied filter"
        );
        assert_eq!(
            counts.total(),
            1,
            "Incorrect per-read total count for short paired read denied filter"
        );
    }

    #[test]
    fn test_short_filter_paired_end_allows() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs(
                "seq",
                None,
                b"AAACCCGGGTTT",
                b"FFFFFFFFFFFF",
            ),
            reverse: Some(bio::io::fastq::Record::with_attrs(
                "seq", None, b"ACTG", b"FFFF",
            )),
            group: ReadGroup::Ungrouped,
        };
        let key = readpair.key();

        let mut f = FilteredReads::new(FilterConfig::new(None, None, None, false));

        let out = f.filter_readpair(&readpair);

        assert_eq!(out, FilterReason::None, "Short paired read not allowed");
        assert_eq!(
            f.totals.short_read, 0,
            "Inccorect filtered read count for short paired read allowed filter"
        );
        assert_eq!(
            f.total(),
            0,
            "Inccorect total filtered read count for short paired read allowed filter"
        );

        // counts map should not get an entry for allowed read
        assert!(
            !f.counts.contains_key(&key),
            "Counts hashmap unexpectedly contains key for short paired read allowed filter"
        );
    }

    // Test low quality
    #[test]
    fn test_quality_filter_single_end_denies() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("seq", None, b"ACTG", b"AAAA"),
            reverse: None,
            group: ReadGroup::Ungrouped,
        };
        let key = readpair.key();

        let mut f = FilteredReads::new(FilterConfig::new(Some(40.0), None, None, false));

        let out = f.filter_readpair(&readpair);

        assert_eq!(
            out,
            FilterReason::LowMeanQuality,
            "Low quality read not filtered"
        );
        assert_eq!(
            f.totals.low_mean_quality, 1,
            "Inccorect filtered read count for low quality denied filter"
        );
        assert_eq!(
            f.total(),
            1,
            "Inccorect total filtered read count for low quality denied filter"
        );

        // counts map checks
        assert!(
            f.counts.contains_key(&key),
            "Counts hashmap missing read key for low quality denied filter"
        );
        let counts = f.counts.get(&key).unwrap();
        assert_eq!(
            counts.low_mean_quality, 1,
            "Incorrect per-read low_mean_quality count for low quality denied filter"
        );
        assert_eq!(
            counts.total(),
            1,
            "Incorrect per-read total count for low quality denied filter"
        );
    }

    #[test]
    fn test_quality_filter_single_end_allows() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("seq", None, b"ACTG", b"AAAA"),
            reverse: None,
            group: ReadGroup::Ungrouped,
        };
        let key = readpair.key();

        let mut f = FilteredReads::new(FilterConfig::new(None, None, None, false));

        let out = f.filter_readpair(&readpair);

        assert_eq!(out, FilterReason::None, "Low quality read not allowed");
        assert_eq!(
            f.totals.short_read, 0,
            "Inccorect filtered read count for low quality read allowed filter"
        );
        assert_eq!(
            f.total(),
            0,
            "Inccorect total filtered read count for low quality read allowed filter"
        );

        // counts map should not get an entry for allowed read
        assert!(
            !f.counts.contains_key(&key),
            "Counts hashmap unexpectedly contains key for low quality read allowed filter"
        );
    }

    #[test]
    fn test_quality_filter_paired_end_denies() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("seq", None, b"ACGT", b"KKKK"),
            reverse: Some(bio::io::fastq::Record::with_attrs(
                "seq", None, b"ACTG", b"AAAA",
            )),
            group: ReadGroup::Ungrouped,
        };
        let key = readpair.key();

        let mut f = FilteredReads::new(FilterConfig::new(Some(40.0), None, None, false));

        let out = f.filter_readpair(&readpair);

        assert_eq!(
            out,
            FilterReason::LowMeanQuality,
            "Low quality paired read not filtered"
        );
        assert_eq!(
            f.totals.low_mean_quality, 1,
            "Inccorect filtered read count for low quality paired read denied filter"
        );
        assert_eq!(
            f.total(),
            1,
            "Inccorect total filtered read count for low quality paired read denied filter"
        );

        // counts map checks
        assert!(
            f.counts.contains_key(&key),
            "Counts hashmap missing read key for low quality paired read denied filter"
        );
        let counts = f.counts.get(&key).unwrap();
        assert_eq!(
            counts.low_mean_quality, 1,
            "Incorrect per-read low_mean_quality count for low quality paired read denied filter"
        );
        assert_eq!(
            counts.total(),
            1,
            "Incorrect per-read total count for low quality paired read denied filter"
        );
    }

    #[test]
    fn test_quality_filter_paired_end_allows() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("seq", None, b"ACGT", b"KKKK"),
            reverse: Some(bio::io::fastq::Record::with_attrs(
                "seq", None, b"ACTG", b"AAAA",
            )),
            group: ReadGroup::Ungrouped,
        };
        let key = readpair.key();

        let mut f = FilteredReads::new(FilterConfig::new(None, None, None, false));

        let out = f.filter_readpair(&readpair);

        assert_eq!(
            out,
            FilterReason::None,
            "Low quality paired read not allowed"
        );
        assert_eq!(
            f.totals.short_read, 0,
            "Inccorect filtered read count for low quality paired read allowed filter"
        );
        assert_eq!(
            f.total(),
            0,
            "Inccorect total filtered read count for low quality paired read allowed filter"
        );

        // counts map should not get an entry for allowed read
        assert!(
            !f.counts.contains_key(&key),
            "Counts hashmap unexpectedly contains key for low quality paired read allowed filter"
        );
    }

    #[test]
    fn test_alignment_filter_single_end_denies() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("seq", None, b"ACGTACGT", b"FFFFFFFF"),
            reverse: None,
            group: ReadGroup::Ungrouped,
        };
        let key = readpair.key();

        // tolerance 0.8 of expected=100 -> minimum = 80
        let aln_tol =
            AlignmentTolerance::new(0.8, 100, 100).expect("AlignmentTolerance config failed");
        let mut f = FilteredReads::new(FilterConfig::new(None, Some(aln_tol), None, false));

        // forward below threshold
        let f_aln = Alignment {
            score: 79,
            ..Default::default()
        };
        let out = f.filter_alignment(&readpair, &f_aln, None);

        assert_eq!(
            out,
            FilterReason::BadAlignment,
            "Bad alignment single-end read not filtered"
        );
        assert_eq!(
            f.totals.bad_alignment, 1,
            "Inccorect filtered read count for bad alignment single-end denied filter"
        );

        // counts map checks
        assert!(
            f.counts.contains_key(&key),
            "Counts hashmap missing read key for bad alignment single-end denied filter"
        );
        let counts = f.counts.get(&key).unwrap();
        assert_eq!(
            counts.bad_alignment, 1,
            "Inccorect per-read bad_alignment count for bad alignment single-end denied filter"
        );
        assert_eq!(
            counts.total(),
            1,
            "Inccorect per-read total count for bad alignment single-end denied filter"
        );
    }

    #[test]
    fn test_alignment_filter_single_end_allows() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("seq", None, b"ACGTACGT", b"FFFFFFFF"),
            reverse: None,
            group: ReadGroup::Ungrouped,
        };
        let key = readpair.key();

        // threshold = 80 again
        let aln_tol =
            AlignmentTolerance::new(0.8, 100, 100).expect("AlignmentTolerance config failed");
        let mut f = FilteredReads::new(FilterConfig::new(None, Some(aln_tol), None, false));

        // forward meets threshold
        let f_aln = Alignment {
            score: 80,
            ..Default::default()
        };
        let out = f.filter_alignment(&readpair, &f_aln, None);

        assert_eq!(
            out,
            FilterReason::None,
            "Bad alignment single-end read not allowed"
        );
        assert_eq!(
            f.totals.bad_alignment, 0,
            "Inccorect filtered read count for bad alignment single-end allowed filter"
        );

        // counts map should not get an entry for allowed read
        assert!(
            !f.counts.contains_key(&key),
            "Counts hashmap unexpectedly contains key for bad alignment single-end allowed filter"
        );
    }

    #[test]
    fn test_alignment_filter_paired_end_denies() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs(
                "seq",
                None,
                b"AAACCCGGGTTT",
                b"FFFFFFFFFFFF",
            ),
            reverse: Some(bio::io::fastq::Record::with_attrs(
                "seq",
                None,
                b"ACTGACTG",
                b"FFFFFFFF",
            )),
            group: ReadGroup::Ungrouped,
        };
        let key = readpair.key();

        // threshold = 80
        let aln_tol =
            AlignmentTolerance::new(0.8, 100, 100).expect("AlignmentTolerance config failed");
        let mut f = FilteredReads::new(FilterConfig::new(None, Some(aln_tol), None, false));

        // forward fails, reverse passes
        let f_aln = Alignment {
            score: 79,
            ..Default::default()
        };
        let r_aln = Alignment {
            score: 95,
            ..Default::default()
        };
        let out = f.filter_alignment(&readpair, &f_aln, Some(&r_aln));

        assert_eq!(
            out,
            FilterReason::BadAlignment,
            "Bad alignment paired read not filtered"
        );
        assert_eq!(
            f.totals.bad_alignment, 1,
            "Inccorect filtered read count for bad alignment paired-end denied filter"
        );

        // counts map checks
        assert!(
            f.counts.contains_key(&key),
            "Counts hashmap missing read key for bad alignment paired-end denied filter"
        );
        let counts = f.counts.get(&key).unwrap();
        assert_eq!(
            counts.bad_alignment, 1,
            "Inccorect per-read bad_alignment count for bad alignment paired-end denied filter"
        );
        assert_eq!(
            counts.total(),
            1,
            "Inccorect per-read total count for bad alignment paired-end denied filter"
        );
    }

    #[test]
    fn test_alignment_filter_paired_end_allows() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs(
                "seq",
                None,
                b"AAACCCGGGTTT",
                b"FFFFFFFFFFFF",
            ),
            reverse: Some(bio::io::fastq::Record::with_attrs(
                "seq",
                None,
                b"ACTGACTG",
                b"FFFFFFFF",
            )),
            group: ReadGroup::Ungrouped,
        };
        let key = readpair.key();

        // threshold = 80
        let aln_tol =
            AlignmentTolerance::new(0.8, 100, 100).expect("AlignmentTolerance config failed");
        let mut f = FilteredReads::new(FilterConfig::new(None, Some(aln_tol), None, false));

        // both pass threshold
        let f_aln = Alignment {
            score: 100,
            ..Default::default()
        };
        let r_aln = Alignment {
            score: 80,
            ..Default::default()
        };
        let out = f.filter_alignment(&readpair, &f_aln, Some(&r_aln));

        assert_eq!(
            out,
            FilterReason::None,
            "Bad alignment paired read not allowed"
        );
        assert_eq!(
            f.totals.bad_alignment, 0,
            "Inccorect filtered read count for bad alignment paired-end allowed filter"
        );

        // counts map should not get an entry for allowed read
        assert!(
            !f.counts.contains_key(&key),
            "Counts hashmap unexpectedly contains key for bad alignment paired-end allowed filter"
        );
    }

    #[test]
    fn test_repeated_empty_filter_increments_counts() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("id", None, b"", b""),
            reverse: None,
            group: ReadGroup::Ungrouped,
        };
        let key = readpair.key();

        let mut f = FilteredReads::new(FilterConfig::new(None, None, None, true));

        // apply twice
        f.filter_readpair(&readpair);
        f.filter_readpair(&readpair);

        assert_eq!(
            f.totals.empty_read, 2,
            "Incorrect totals.empty_read after repeated empty filter"
        );
        assert_eq!(
            f.total(),
            2,
            "Incorrect FilteredReads.total() after repeated empty filter"
        );

        let counts = f.counts.get(&key).unwrap();
        assert_eq!(
            counts.empty_read, 2,
            "Incorrect per-read empty_read count after repeated empty filter"
        );
        assert_eq!(
            counts.total(),
            2,
            "Incorrect per-read total count after repeated empty filter"
        );
    }

    #[test]
    fn test_repeated_short_filter_increments_counts() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("id", None, b"ACTG", b"FFFF"),
            reverse: None,
            group: ReadGroup::Ungrouped,
        };
        let key = readpair.key();

        let mut f = FilteredReads::new(FilterConfig::new(None, None, Some(10), false));

        // apply twice
        f.filter_readpair(&readpair);
        f.filter_readpair(&readpair);

        assert_eq!(
            f.totals.short_read, 2,
            "Incorrect totals.short_read after repeated short filter"
        );
        assert_eq!(
            f.total(),
            2,
            "Incorrect FilteredReads.total() after repeated short filter"
        );

        let counts = f.counts.get(&key).unwrap();
        assert_eq!(
            counts.short_read, 2,
            "Incorrect per-read short_read count after repeated short filter"
        );
        assert_eq!(
            counts.total(),
            2,
            "Incorrect per-read total count after repeated short filter"
        );
    }

    #[test]
    fn test_repeated_quality_filter_increments_counts() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("id", None, b"ACTG", b"AAAA"),
            reverse: None,
            group: ReadGroup::Ungrouped,
        };
        let key = readpair.key();

        let mut f = FilteredReads::new(FilterConfig::new(Some(40.0), None, None, false));

        // apply twice
        f.filter_readpair(&readpair);
        f.filter_readpair(&readpair);

        assert_eq!(
            f.totals.low_mean_quality, 2,
            "Incorrect totals.low_mean_quality after repeated quality filter"
        );
        assert_eq!(
            f.total(),
            2,
            "Incorrect FilteredReads.total() after repeated quality filter"
        );

        let counts = f.counts.get(&key).unwrap();
        assert_eq!(
            counts.low_mean_quality, 2,
            "Incorrect per-read low_mean_quality count after repeated quality filter"
        );
        assert_eq!(
            counts.total(),
            2,
            "Incorrect per-read total count after repeated quality filter"
        );
    }

    #[test]
    fn test_repeated_alignment_filter_increments_counts() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("id", None, b"ACGTACGT", b"FFFFFFFF"),
            reverse: None,
            group: ReadGroup::Ungrouped,
        };
        let key = readpair.key();

        let aln_tol =
            AlignmentTolerance::new(0.8, 100, 100).expect("AlignmentTolerance config failed");
        let mut f = FilteredReads::new(FilterConfig::new(None, Some(aln_tol), None, false));

        let f_aln = Alignment {
            score: 79,
            ..Default::default()
        }; // below threshold

        // apply twice
        f.filter_alignment(&readpair, &f_aln, None);
        f.filter_alignment(&readpair, &f_aln, None);

        assert_eq!(
            f.totals.bad_alignment, 2,
            "Incorrect totals.bad_alignment after repeated alignment filter"
        );
        assert_eq!(
            f.total(),
            2,
            "Incorrect FilteredReads.total() after repeated alignment filter"
        );

        let counts = f.counts.get(&key).unwrap();
        assert_eq!(
            counts.bad_alignment, 2,
            "Incorrect per-read bad_alignment count after repeated alignment filter"
        );
        assert_eq!(
            counts.total(),
            2,
            "Incorrect per-read total count after repeated alignment filter"
        );
    }

    #[test]
    fn test_mixed_filters_same_sequence_two_readpairs() {
        let low_q_rp = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("id", None, b"ACTGACTG", b"AAAAAAAA"),
            reverse: None,
            group: ReadGroup::Ungrouped,
        };
        let high_q_rp = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("id", None, b"ACTGACTG", b"KKKKKKKK"),
            reverse: None,
            group: ReadGroup::Ungrouped,
        };

        let aln_tol =
            AlignmentTolerance::new(0.8, 100, 100).expect("AlignmentTolerance config failed");
        let mut f = FilteredReads::new(FilterConfig::new(Some(40.0), Some(aln_tol), None, false));
        let key = low_q_rp.key(); // same key for both readpairs (same forward seq; no reverse)

        // Filter_readpair on low-quality read => LowMeanQuality
        let out1 = f.filter_readpair(&low_q_rp);
        assert_eq!(
            out1,
            FilterReason::LowMeanQuality,
            "Low-quality read was not filtered as LowMeanQuality"
        );

        // Filter_readpair on high-quality read => None (doesn't trigger length or empty)
        let out2 = f.filter_readpair(&high_q_rp);
        assert_eq!(
            out2,
            FilterReason::None,
            "High-quality read unexpectedly filtered by filter_readpair"
        );

        // Mocked alignment below threshold on high-quality read => BadAlignment
        let bad_f = Alignment {
            score: 79,
            ..Default::default()
        }; // threshold is 80 (0.8 * 100)
        let out3 = f.filter_alignment(&high_q_rp, &bad_f, None);
        assert_eq!(
            out3,
            FilterReason::BadAlignment,
            "Bad-alignment case did not return BadAlignment"
        );

        // ---- totals checks ----
        assert_eq!(
            f.totals.low_mean_quality, 1,
            "Incorrect totals.low_mean_quality in mixed test"
        );
        assert_eq!(
            f.totals.bad_alignment, 1,
            "Incorrect totals.bad_alignment in mixed test"
        );
        assert_eq!(
            f.totals.short_read, 0,
            "Incorrect totals.short_read in mixed test"
        );
        assert_eq!(
            f.totals.empty_read, 0,
            "Incorrect totals.empty_read in mixed test (should be zero)"
        );
        assert_eq!(
            f.total(),
            2,
            "Incorrect FilteredReads.total() in mixed test"
        );

        // ---- per-read counts (same key for both readpairs) ----
        let counts = f
            .counts
            .get(&key)
            .expect("Counts hashmap missing key in mixed test");
        assert_eq!(
            counts.low_mean_quality, 1,
            "Incorrect per-read low_mean_quality in mixed test"
        );
        assert_eq!(
            counts.bad_alignment, 1,
            "Incorrect per-read bad_alignment in mixed test"
        );
        assert_eq!(
            counts.short_read, 0,
            "Incorrect per-read short_read in mixed test"
        );
        assert_eq!(
            counts.empty_read, 0,
            "Incorrect per-read empty_read in mixed test (should be zero)"
        );
        assert_eq!(counts.total(), 2, "Incorrect per-read total in mixed test");
    }
}
