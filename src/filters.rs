//! Filters to identify and count failing reads
//!
//! Provides a range of filtering criteria that a function
//! processing reads can tap into via a simple config and
//! interface, easily filtering reads for a range of reasons.
//! Also makes it easy to keep track of filtered reads and
//! print a summary to file.
use bio::bio_types::alignment::Alignment;
use bio::io::fastq::Record;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

use crate::errors::{ReadCountError, seq_to_string_or_log};
use crate::parsing::{ReadGroup, ReadKey, ReadPair};
use crate::utils::{div_or_zero, mean_quality};

/// Function filtering based on a read pair
type ReadFilter = fn(&Record, Option<&Record>, &FilterConfig) -> Option<FilterReason>;

/// Function filtering based on an alignment
type AlignmentFilter = fn(&Alignment, Option<&Alignment>, &FilterConfig) -> Option<FilterReason>;

/// Reasons for filtering a read
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FilterReason {
    LowMeanQuality,
    BadAlignment,
    LongRead,
    ShortRead,
    EmptyRead,
}

/// Metadata about filters
#[derive(Debug)]
pub struct FilterMeta {
    pub id: &'static str,
    pub label: &'static str,
}

impl FilterReason {
    const N_REASONS: usize = FilterReason::ALL_FILTERS.len();
    const ALL_FILTERS: &[FilterReason] = &[
        FilterReason::EmptyRead,
        FilterReason::ShortRead,
        FilterReason::LongRead,
        FilterReason::LowMeanQuality,
        FilterReason::BadAlignment,
    ];

    #[inline]
    pub const fn as_index(self) -> usize {
        self as usize
    }

    pub fn meta(self) -> FilterMeta {
        match self {
            FilterReason::EmptyRead => FilterMeta {
                id: "empty_read",
                label: "Read length = 0",
            },
            FilterReason::ShortRead => FilterMeta {
                id: "short_read",
                label: "Read length < minimum",
            },
            FilterReason::LongRead => FilterMeta {
                id: "long_read",
                label: "Read length > maximum",
            },
            FilterReason::LowMeanQuality => FilterMeta {
                id: "low_mean_quality",
                label: "Mean quality < minimum",
            },
            FilterReason::BadAlignment => FilterMeta {
                id: "bad_alignment",
                label: "Alignment quality < tolerance",
            },
        }
    }
}

/// Functions used to filter based on read pairs
const READPAIR_FILTERS: &[ReadFilter] = &[
    empty_read_filter,
    short_read_filter,
    long_read_filter,
    low_mean_quality_filter,
];

/// Functions used to filter based on alignments
const ALIGNMENT_FILTERS: &[AlignmentFilter] = &[bad_alignment_filter];

fn empty_read_filter(f: &Record, r: Option<&Record>, cfg: &FilterConfig) -> Option<FilterReason> {
    if !cfg.filter_empty {
        return None;
    }

    if f.seq().is_empty() {
        return Some(FilterReason::EmptyRead);
    }

    if let Some(r) = r {
        if r.seq().is_empty() {
            return Some(FilterReason::EmptyRead);
        }
    }

    None
}

fn short_read_filter(f: &Record, r: Option<&Record>, cfg: &FilterConfig) -> Option<FilterReason> {
    let min = cfg.minimum_length?;

    if f.seq().len() < min {
        return Some(FilterReason::ShortRead);
    }

    if let Some(r) = r {
        if r.seq().len() < min {
            return Some(FilterReason::ShortRead);
        }
    }

    None
}

fn long_read_filter(f: &Record, r: Option<&Record>, cfg: &FilterConfig) -> Option<FilterReason> {
    let max = cfg.maximum_length?;

    if f.seq().len() > max {
        return Some(FilterReason::LongRead);
    }

    if let Some(r) = r {
        if r.seq().len() > max {
            return Some(FilterReason::LongRead);
        }
    }

    None
}

fn low_mean_quality_filter(
    f: &Record,
    r: Option<&Record>,
    cfg: &FilterConfig,
) -> Option<FilterReason> {
    let q = cfg.mean_quality_threshold?;

    if mean_quality(f.qual()) < q {
        return Some(FilterReason::LowMeanQuality);
    }

    if let Some(r) = r {
        if mean_quality(r.qual()) < q {
            return Some(FilterReason::LowMeanQuality);
        }
    }

    None
}

// Alignment-based
fn bad_alignment_filter(
    f: &Alignment,
    r: Option<&Alignment>,
    cfg: &FilterConfig,
) -> Option<FilterReason> {
    let tol = &cfg.alignment_tolerance?;

    if f.score < tol.minimum_f_score {
        return Some(FilterReason::BadAlignment);
    }

    if let Some(r) = r {
        if r.score < tol.minimum_r_score {
            return Some(FilterReason::BadAlignment);
        }
    }

    None
}

/// Configuration for filtering
///
/// Instructions for how to filter reads
#[derive(Debug, Clone, PartialEq)]
pub struct FilterConfig {
    pub mean_quality_threshold: Option<f32>,
    pub alignment_tolerance: Option<AlignmentTolerance>,
    pub minimum_length: Option<usize>,
    pub maximum_length: Option<usize>,
    pub filter_empty: bool,
}

impl FilterConfig {
    pub fn new(
        mean_quality_threshold: Option<f32>,
        alignment_tolerance: Option<AlignmentTolerance>,
        minimum_length: Option<usize>,
        maximum_length: Option<usize>,
        filter_empty: bool,
    ) -> Self {
        Self {
            mean_quality_threshold,
            alignment_tolerance,
            minimum_length,
            maximum_length,
            filter_empty,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
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

/// Count of reads filtered for each reason, stored as an array indexed by FilterReason u8 values
#[derive(Debug, Clone)]
pub struct FilteredCounts([u64; FilterReason::N_REASONS]);

impl FilteredCounts {
    pub fn new() -> Self {
        Self([0; FilterReason::N_REASONS])
    }

    #[inline]
    pub fn get(&self, r: &FilterReason) -> u64 {
        self.0[r.as_index()]
    }

    fn increment_count(&mut self, r: FilterReason) {
        self.0[r.as_index()] += 1;
    }

    /// Total count of filtered reads
    pub fn total(&self) -> u64 {
        self.0.iter().sum()
    }

    /// Merge counts from another FilteredCounts
    pub fn merge(&mut self, new_counts: FilteredCounts) {
        for i in 0..FilterReason::N_REASONS {
            self.0[i] += new_counts.0[i]
        }
    }

    /// Headers for the TSV produced by to_tsv_line
    fn wide_tsv_headers() -> String {
        let ids: Vec<&str> = FilterReason::ALL_FILTERS
            .iter()
            .map(|r| r.meta().id)
            .collect();

        format!("count\tproportion\t{}", ids.join("\t"))
    }

    /// TSV line giving the total and counts of each reason for a given
    fn to_wide_tsv_line(&self, total: f32) -> String {
        let string_counts: Vec<String> = self.0.iter().map(|i| i.to_string()).collect();

        format!(
            "{}\t{:.4}\t{}",
            self.total(),
            self.total() as f32 / total,
            string_counts.join("\t"),
        )
    }

    /// Generate a string of TSV lines representing total counts
    ///
    /// Creates a TSV string with columns for filter reason, count and
    /// proportion of the supplied total read count. Mainly for use
    /// when outputting from ReadSummary
    pub fn to_long_tsv_lines(&self, total: u64) -> String {
        let filtered_total = self.total();

        let mut out = String::with_capacity(300);

        out.push_str(&format!(
            "filtered\ttotal\t{}\t{:.4}\t1.000\n",
            filtered_total,
            div_or_zero(filtered_total as f32, total as f32),
        ));

        for r in FilterReason::ALL_FILTERS {
            out.push_str(&format!(
                "filtered\t{}\t{}\t{:.4}\t{:.4}\n",
                r.meta().id,
                self.get(r),
                div_or_zero(self.get(r) as f32, total as f32),
                div_or_zero(self.get(r) as f32, filtered_total as f32),
            ));
        }

        out
    }
}

impl Default for FilteredCounts {
    fn default() -> Self {
        Self::new()
    }
}

/// Container for filtered reads
///
/// Keeps track of counts for filtered reads during counting
#[derive(Debug, Clone)]
pub struct FilteredReads {
    pub config: FilterConfig,
    pub totals: FilteredCounts,
    pub counts: HashMap<ReadKey, HashMap<ReadGroup, FilteredCounts>>,
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
    pub fn increment_count(&mut self, read: &ReadPair, reason: FilterReason) {
        self.totals.increment_count(reason);

        let key = read.key();
        let group = &read.group;

        match self.counts.get_mut(&key) {
            Some(groups) => match groups.get_mut(group) {
                Some(c) => c.increment_count(reason),
                None => {
                    let mut new_counts = FilteredCounts::new();
                    new_counts.increment_count(reason);

                    groups.insert(group.clone(), new_counts);
                }
            },
            None => {
                let mut new_counts = FilteredCounts::new();
                new_counts.increment_count(reason);

                let mut new_groups = HashMap::new();
                new_groups.insert(group.clone(), new_counts);

                self.counts.insert(key.clone(), new_groups);
            }
        }
    }

    /// Determine if a readpair should be filtered based on the supplied config
    ///
    /// Checks whether the read should be filtered, adding it to the appropriate count if so, and
    /// returns a FilterReason determining why it was filtered.
    pub fn filter_readpair(&mut self, record: &ReadPair) -> Option<FilterReason> {
        let f_read = &record.forward;
        let r_read = record.reverse.as_ref();

        for f in READPAIR_FILTERS {
            if let Some(r) = f(f_read, r_read, &self.config) {
                self.increment_count(record, r);
                return Some(r);
            }
        }

        None
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
    ) -> Option<FilterReason> {
        for f in ALIGNMENT_FILTERS {
            if let Some(r) = f(f_alignment, r_alignment, &self.config) {
                self.increment_count(record, r);
                return Some(r);
            }
        }

        None
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

        for (key, new_groups) in new_reads.counts {
            match self.counts.get_mut(&key) {
                Some(old_groups) => {
                    for (new_group, new_counts) in new_groups {
                        match old_groups.get_mut(&new_group) {
                            Some(old_counts) => old_counts.merge(new_counts),
                            None => {
                                old_groups.insert(new_group, new_counts);
                            }
                        }
                    }
                }
                None => {
                    self.counts.insert(key, new_groups);
                }
            }
        }

        Ok(())
    }

    /// Write a TSV file listing the filtered reads
    ///
    /// Writes a TSV with columns for forward sequence, reverse sequence,
    /// total count, frequency, then one for each filter reason count.
    pub fn write_filter_tsv(&self, file: File, sort: bool) -> Result<(), anyhow::Error> {
        let total = self.total() as f32;

        let mut writer = BufWriter::new(file);
        let mut keys: Vec<(&ReadKey, u64)> = self
            .counts
            .iter()
            .map(|x| (x.0, x.1.iter().map(|y| y.1.total()).sum()))
            .collect();

        if sort {
            // Invert count to get desc order
            keys.sort_unstable_by_key(|x| std::cmp::Reverse(x.1));
        }

        // Write header
        writeln!(
            writer,
            "group\tforward\treverse\t{}",
            FilteredCounts::wide_tsv_headers()
        )?;

        for (key, _) in keys {
            let groups = self
                .counts
                .get(key)
                .expect("Count key from extracted key list missing from FilteredReads");

            for (group, counts) in groups {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}",
                    group,
                    seq_to_string_or_log(&key.0),
                    match &key.1 {
                        Some(x) => seq_to_string_or_log(x),
                        None => "".to_string(),
                    },
                    counts.to_wide_tsv_line(total),
                )?;
            }
        }

        writer.flush()?;
        Ok(())
    }

    /// Generate a string of TSV lines representing the total filter counts
    ///
    /// Creates a TSV string with columns for filter reason, count and
    /// proportion of the supplied total read count. Mainly for use
    /// when outputting from ReadSummary
    pub fn to_summary_tsv_lines(&self, total: u64) -> String {
        let filtered_total = self.total();

        let mut out = String::with_capacity(300);

        out.push_str(&format!(
            "filtered\ttotal\t{}\t{:.4}\t1.000\n",
            filtered_total,
            filtered_total as f32 / total as f32,
        ));

        for r in FilterReason::ALL_FILTERS {
            out.push_str(&format!(
                "filtered\t{}\t{}\t{:.4}\t{:.4}\n",
                r.meta().id,
                self.totals.get(r),
                self.totals.get(r) as f32 / total as f32,
                self.totals.get(r) as f32 / filtered_total as f32,
            ));
        }

        out
    }
}

#[cfg(test)]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::parsing::{ReadGroup, ReadPair};
    use bio::bio_types::alignment::Alignment;

    /// Generate a read pair
    fn rp(f_seq: &str, f_qual: &str, r_seq: Option<&str>, r_qual: Option<&str>) -> ReadPair {
        let f = bio::io::fastq::Record::with_attrs("f", None, f_seq.as_bytes(), f_qual.as_bytes());
        let r = r_seq.map(|s| {
            bio::io::fastq::Record::with_attrs(
                "r",
                None,
                s.as_bytes(),
                r_qual.expect("r_qual required if r_seq is Some").as_bytes(),
            )
        });
        ReadPair {
            forward: f,
            reverse: r,
            group: ReadGroup::Ungrouped,
        }
    }

    /// Create alignment tolerance
    fn tol(frac: f32, exp_f: i32, exp_r: i32) -> AlignmentTolerance {
        AlignmentTolerance::new(frac, exp_f, exp_r).unwrap()
    }

    /// Generate mock alignment with score
    fn aln(score: i32) -> Alignment {
        Alignment {
            score,
            ..Default::default()
        }
    }

    // ---------- READPAIR FILTERS ----------
    #[test]
    fn readpair_filters() {
        #[derive(Debug)]
        struct Case {
            name: &'static str,
            cfg: FilterConfig,
            rp: ReadPair,
            expected: Option<FilterReason>,
        }

        let cases = vec![
            // Empty reads
            Case {
                name: "empty single-end, filter_empty=true",
                cfg: FilterConfig::new(None, None, None, None, true),
                rp: rp("", "", None, None),
                expected: Some(FilterReason::EmptyRead),
            },
            Case {
                name: "empty single-end, filter_empty=false",
                cfg: FilterConfig::new(None, None, None, None, false),
                rp: rp("", "", None, None),
                expected: None,
            },
            Case {
                name: "empty paired-end, filter_empty=true",
                cfg: FilterConfig::new(None, None, None, None, true),
                rp: rp("ACTG", "FFFF", Some(""), Some("")),
                expected: Some(FilterReason::EmptyRead),
            },
            Case {
                name: "empty paired-end, filter_empty=false",
                cfg: FilterConfig::new(None, None, None, None, false),
                rp: rp("ACTG", "FFFF", Some(""), Some("")),
                expected: None,
            },
            // Short reads
            Case {
                name: "short single-end, minimum_length=5",
                cfg: FilterConfig::new(None, None, Some(5), None, false),
                rp: rp("ACTG", "FFFF", None, None),
                expected: Some(FilterReason::ShortRead),
            },
            Case {
                name: "short single-end, minimum_length=None",
                cfg: FilterConfig::new(None, None, None, None, false),
                rp: rp("ACTG", "FFFF", None, None),
                expected: None,
            },
            Case {
                name: "short paired-end, minimum_length=5",
                cfg: FilterConfig::new(None, None, Some(5), None, false),
                rp: rp("ACTGACTG", "FFFFFFFF", Some("ACTG"), Some("FFFF")),
                expected: Some(FilterReason::ShortRead),
            },
            Case {
                name: "short paired-end, minimum_length=None",
                cfg: FilterConfig::new(None, None, None, None, false),
                rp: rp("ACTGACTG", "FFFFFFFF", Some("ACTG"), Some("FFFF")),
                expected: None,
            },
            // Long reads
            Case {
                name: "long single-end, maximum_length=5",
                cfg: FilterConfig::new(None, None, None, Some(5), false),
                rp: rp("ACTGACTG", "FFFFFFFF", None, None),
                expected: Some(FilterReason::LongRead),
            },
            Case {
                name: "long single-end, maximum_length=None",
                cfg: FilterConfig::new(None, None, None, None, false),
                rp: rp("ACTGACTG", "FFFFFFFF", None, None),
                expected: None,
            },
            Case {
                name: "long paired-end, maximum_length=5",
                cfg: FilterConfig::new(None, None, None, Some(5), false),
                rp: rp("ACTG", "FFFF", Some("ACTGACTG"), Some("FFFFFFFF")),
                expected: Some(FilterReason::LongRead),
            },
            Case {
                name: "long paired-end, maximum_length=None",
                cfg: FilterConfig::new(None, None, None, None, false),
                rp: rp("ACTG", "FFFF", Some("ACTGACTG"), Some("FFFFFFFF")),
                expected: None,
            },
            // Mean quality
            Case {
                name: "low mean quality single-end, threshold=40",
                cfg: FilterConfig::new(Some(40.0), None, None, None, false),
                rp: rp("ACTG", "AAAA", None, None), // 32 phred
                expected: Some(FilterReason::LowMeanQuality),
            },
            Case {
                name: "high mean quality single-end, threshold=40",
                cfg: FilterConfig::new(Some(40.0), None, None, None, false),
                rp: rp("ACTG", "KKKK", None, None), // 42 phred
                expected: None,
            },
            Case {
                name: "low mean quality paired-end, threshold=40",
                cfg: FilterConfig::new(Some(40.0), None, None, None, false),
                rp: rp("ACTG", "KKKK", Some("ACTG"), Some("AAAA")), // 32 phred
                expected: Some(FilterReason::LowMeanQuality),
            },
            Case {
                name: "high mean quality paired-end, threshold=40",
                cfg: FilterConfig::new(Some(40.0), None, None, None, false),
                rp: rp("ACTG", "KKKK", Some("ACTG"), Some("KKKK")), // 42 phred
                expected: None,
            },
        ];

        for c in cases {
            let mut fr = FilteredReads::new(c.cfg.clone());
            let key = c.rp.key();

            let got = fr.filter_readpair(&c.rp);
            assert_eq!(
                got, c.expected,
                "Unexpected filter output (case: {})",
                c.name
            );

            match c.expected {
                Some(reason) => {
                    // totals increased
                    assert_eq!(
                        fr.totals.get(&reason),
                        1,
                        "Total count not incremented (case: {})",
                        c.name
                    );
                    let read = fr.counts.get(&key).expect("Missing per-read counts");
                    let grp = read
                        .get(&ReadGroup::Ungrouped)
                        .expect("Missing group counts");
                    assert_eq!(grp.get(&reason), 1, "Read not tracked (case: {})", c.name);
                }
                None => {
                    assert!(
                        !fr.counts.contains_key(&key),
                        "Total count incorrectly incremented (case: {})",
                        c.name
                    );
                }
            }
        }
    }

    // ---------- ALIGNMENT FILTERS ----------
    #[test]
    fn alignment_filters() {
        #[derive(Debug)]
        struct Case {
            name: &'static str,
            cfg: FilterConfig,
            f: Alignment,
            r: Option<Alignment>,
            expected: Option<FilterReason>,
        }

        let rp = rp("ACGTACGT", "FFFFFFFF", Some("ACGTACGT"), Some("FFFFFFFF"));
        let t = tol(0.8, 100, 100); // min = 80

        let cases = vec![
            Case {
                name: "single-end below threshold",
                cfg: FilterConfig::new(None, Some(t.clone()), None, None, false),
                f: aln(79),
                r: None,
                expected: Some(FilterReason::BadAlignment),
            },
            Case {
                name: "single-end meets threshold",
                cfg: FilterConfig::new(None, Some(t.clone()), None, None, false),
                f: aln(80),
                r: None,
                expected: None,
            },
            Case {
                name: "paired-end forward fails, reverse fails",
                cfg: FilterConfig::new(None, Some(t.clone()), None, None, false),
                f: aln(79),
                r: Some(aln(79)),
                expected: Some(FilterReason::BadAlignment),
            },
            Case {
                name: "paired-end forward fails, reverse passes",
                cfg: FilterConfig::new(None, Some(t.clone()), None, None, false),
                f: aln(79),
                r: Some(aln(80)),
                expected: Some(FilterReason::BadAlignment),
            },
            Case {
                name: "paired-end forward passes, reverse fails",
                cfg: FilterConfig::new(None, Some(t.clone()), None, None, false),
                f: aln(80),
                r: Some(aln(79)),
                expected: Some(FilterReason::BadAlignment),
            },
            Case {
                name: "paired-end forward passes, reverse passes",
                cfg: FilterConfig::new(None, Some(t.clone()), None, None, false),
                f: aln(80),
                r: Some(aln(80)),
                expected: None,
            },
        ];

        for c in cases {
            let mut fr = FilteredReads::new(c.cfg.clone());
            let key = rp.key();

            let got = fr.filter_alignment(&rp, &c.f, c.r.as_ref());
            assert_eq!(
                got, c.expected,
                "Unexpected filter output (case: {})",
                c.name
            );

            match c.expected {
                Some(reason) => {
                    assert_eq!(
                        fr.totals.get(&reason),
                        1,
                        "Total count not incremented (case: {})",
                        c.name
                    );
                    let read = fr.counts.get(&key).expect("Missing per-read counts");
                    let grp = read
                        .get(&ReadGroup::Ungrouped)
                        .expect("Missing group counts");
                    assert_eq!(grp.get(&reason), 1, "Read not tracked (case: {})", c.name);
                }
                None => {
                    assert!(
                        fr.counts.get(&key).is_none(),
                        "Total count incorrectly incremented (case: {})",
                        c.name
                    );
                }
            }
        }
    }
}
