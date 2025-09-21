//! Filters to identify anb count failing reads
//!
//! Provides a range of filtering critera that a function
//! processing reads can tap into via a simple config and
//! interface, easily filtering reads for a range of reasons.
//! Also makes it easy to keep track of filtered reads and
//! print a summary to file.
use std::{collections::HashMap, fs::File, io::BufWriter};

use crate::{errors::ReadCountError, parsing::{ReadKey, ReadPair}};
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
}

/// Container for filtered reads
///
/// Keeps track of counts for filtered reads during counting
#[derive(Debug, Clone)]
pub struct FilteredReads {
    config: FilterConfig,
    totals: FilteredCounts,
    counts: HashMap<ReadKey, FilteredCounts>
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
            },
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
                None => {self.counts.insert(key, new_counts);},
            }
        }

        Ok(())
    }

    /// Write a TSV file listing the filteed reads
    ///
    /// Writes a TSV with columns for forward sequence, reverse sequence,
    /// total count, frequency, then one for each filter reason count.
    pub fn write_filter_tsv(&self, file: File, sort: bool) -> Result<(), anyhow::Error> {
        let total = self.total();

        let mut writer = BufWriter::new(file);
        let mut keys: Vec<(&ReadKey, u64)> = self
            .counts
            .iter()
            .map(|x| (x.0, x.1.total()))
            .collect();

        if sort {
            // Invert count to get desc order
            keys.sort_unstable_by_key(|x| 0 - i128::from(x.1));
        }

        todo!();

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

        let mut f = FilteredReads::new(FilterConfig::new(None, None, None, true));

        let out = f.filter_readpair(&readpair);

        assert_eq!(out, FilterReason::EmptyRead, "Empty read not filtered");
        assert_eq!(
            f.totals.empty_read, 1,
            "Inccorect filtered read count for empty read denied filter"
        );
    }

    #[test]
    fn test_empty_filter_single_end_allows() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("seq", None, b"", b""),
            reverse: None,
            group: ReadGroup::Ungrouped,
        };

        let mut f = FilteredReads::new(FilterConfig::new(None, None, None, false));

        let out = f.filter_readpair(&readpair);

        assert_eq!(out, FilterReason::None, "Empty read not allowed");
        assert_eq!(
            f.totals.empty_read, 0,
            "Inccorect filtered read count for empty read allowed filter"
        );
    }

    #[test]
    fn test_empty_filter_paired_end_denies() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("seq", None, b"ACGT", b"FFFF"),
            reverse: Some(bio::io::fastq::Record::with_attrs("seq", None, b"", b"")),
            group: ReadGroup::Ungrouped,
        };

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
    }

    #[test]
    fn test_empty_filter_paired_end_allows() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("seq", None, b"ACGT", b"FFFF"),
            reverse: Some(bio::io::fastq::Record::with_attrs("seq", None, b"", b"")),
            group: ReadGroup::Ungrouped,
        };

        let mut f = FilteredReads::new(FilterConfig::new(None, None, None, false));

        let out = f.filter_readpair(&readpair);

        assert_eq!(out, FilterReason::None, "Empty paired read not allowed");
        assert_eq!(
            f.totals.empty_read, 0,
            "Inccorect filtered read count for empty paired read allowed filter"
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

        let mut f = FilteredReads::new(FilterConfig::new(None, None, Some(10), false));

        let out = f.filter_readpair(&readpair);

        assert_eq!(out, FilterReason::ShortRead, "Short read not filtered");
        assert_eq!(
            f.totals.short_read, 1,
            "Inccorect filtered read count for short read denied filter"
        );
    }

    #[test]
    fn test_short_filter_single_end_allows() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("seq", None, b"ACTG", b"FFFF"),
            reverse: None,
            group: ReadGroup::Ungrouped,
        };

        let mut f = FilteredReads::new(FilterConfig::new(None, None, None, false));

        let out = f.filter_readpair(&readpair);

        assert_eq!(out, FilterReason::None, "Short read not allowed");
        assert_eq!(
            f.totals.short_read, 0,
            "Inccorect filtered read count for short read allowed filter"
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

        let mut f = FilteredReads::new(FilterConfig::new(None, None, None, false));

        let out = f.filter_readpair(&readpair);

        assert_eq!(out, FilterReason::None, "Short paired read not allowed");
        assert_eq!(
            f.totals.short_read, 0,
            "Inccorect filtered read count for short paired read allowed filter"
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
    }

    #[test]
    fn test_quality_filter_single_end_allows() {
        let readpair = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("seq", None, b"ACTG", b"AAAA"),
            reverse: None,
            group: ReadGroup::Ungrouped,
        };

        let mut f = FilteredReads::new(FilterConfig::new(None, None, None, false));

        let out = f.filter_readpair(&readpair);

        assert_eq!(out, FilterReason::None, "Low quality read not allowed");
        assert_eq!(
            f.totals.short_read, 0,
            "Inccorect filtered read count for low quality read allowed filter"
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
    }
}
