//! Read-filtering logic and filtered-read summaries.
//!
//! This module defines the configurable filters applied during counting, along
//! with the data structures used to track why reads were discarded and how often
//! each filtered read sequence was observed.
use bio::bio_types::alignment::Alignment;
use bio::io::fastq::Record;
use std::collections::HashMap;

use crate::errors::ReadCountError;
use crate::groups::ReadGroup;
use crate::seqs::{ReadPair, SeqPair};
use crate::utils::mean_quality;

/// Function filtering based on a read pair
type ReadFilter = fn(&Record, Option<&Record>, &FilterConfig) -> Option<FilterReason>;

/// Function filtering based on an alignment
type AlignmentFilter = fn(&Alignment, Option<&Alignment>, &FilterConfig) -> Option<FilterReason>;

/// Reason why a read or read pair was filtered out during counting.
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FilterReason {
    EmptyRead,
    ShortRead,
    LongRead,
    LowMeanQuality,
    BadAlignment,
}

/// Stable metadata describing a filter reason for display and TSV output.
#[derive(Debug)]
pub struct FilterMeta {
    pub id: &'static str,
    pub label: &'static str,
}

impl FilterReason {
    /// Number of filter reasons defined
    pub const N_REASONS: usize = FilterReason::ALL_FILTERS.len();

    /// List of all FilterReasons for iteration
    pub const ALL_FILTERS: &[FilterReason] = &[
        FilterReason::EmptyRead,
        FilterReason::ShortRead,
        FilterReason::LongRead,
        FilterReason::LowMeanQuality,
        FilterReason::BadAlignment,
    ];

    /// Get index of reason as enum integer and ALL_REASONS index
    #[inline]
    pub const fn as_index(self) -> usize {
        self as usize
    }

    /// Return the stable output identifier and human-readable label for this
    /// filter reason.
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

/// Configuration controlling which filters are applied during counting.
///
/// Filters are applied in a fixed order, and the first matching reason is
/// recorded for a read. Depending on counting mode, some filters act directly on
/// reads (`filter_readpair`) while alignment-quality filters act after alignment
/// (`filter_alignment`).
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

/// Precomputed alignment-score thresholds for filtering alignment-based modes.
///
/// The expected forward and reverse alignment scores are calculated from the
/// LibSpec template and alignment scoring scheme. The supplied tolerance is then
/// used to convert those expected scores into minimum acceptable scores for
/// filtering observed alignments.
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
    /// Create alignment-score thresholds from an expected-score baseline and a
    /// fractional tolerance in the range `[0.0, 1.0]`.
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

/// Counts of how many times a read was filtered for each `FilterReason`.
///
/// Internally this is stored as a fixed-size array indexed by `FilterReason`.
#[derive(Debug, Clone)]
pub struct FilteredCounts([u64; FilterReason::N_REASONS]);

impl FilteredCounts {
    pub fn new() -> Self {
        Self([0; FilterReason::N_REASONS])
    }

    /// Get the count for one specific filter reason.
    #[inline]
    pub fn get(&self, r: &FilterReason) -> u64 {
        self.0[r.as_index()]
    }

    /// Iterate over counts in the stable `FilterReason::ALL_FILTERS` order.
    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = &u64> {
        self.0.iter()
    }

    fn increment_count(&mut self, r: FilterReason) {
        self.0[r.as_index()] += 1;
    }

    /// Total number of filtered observations across all reasons.
    pub fn total(&self) -> u64 {
        self.0.iter().sum()
    }

    /// Merge counts from another FilteredCounts into this one
    pub fn merge(&mut self, new_counts: FilteredCounts) {
        for i in 0..FilterReason::N_REASONS {
            self.0[i] += new_counts.0[i]
        }
    }
}

impl Default for FilteredCounts {
    fn default() -> Self {
        Self::new()
    }
}

/// Summary of reads discarded by filtering.
///
/// This tracks:
/// - the filter configuration used,
/// - total counts by filter reason,
/// - and per-read/per-group filtered counts keyed by read sequence.
///
/// It is used both to decide whether a read should be discarded and to generate
/// summary/output tables of filtered reads.
#[derive(Debug, Clone)]
pub struct FilteredReads {
    pub config: FilterConfig,
    pub totals: FilteredCounts,
    pub counts: HashMap<SeqPair, HashMap<ReadGroup, FilteredCounts>>,
}

impl FilteredReads {
    /// Create an empty filtered-read container for a specific filter configuration.
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

    /// Apply read-level filters to a read pair.
    ///
    /// Filters are evaluated in the order defined by `READPAIR_FILTERS`. The first
    /// matching `FilterReason` is returned. If `increment` is true, that reason is
    /// also recorded in the filtered-read counts.
    pub fn filter_readpair(&mut self, record: &ReadPair, increment: bool) -> Option<FilterReason> {
        let f_read = &record.forward;
        let r_read = record.reverse.as_ref();

        for f in READPAIR_FILTERS {
            if let Some(r) = f(f_read, r_read, &self.config) {
                if increment {
                    self.increment_count(record, r);
                }
                return Some(r);
            }
        }

        None
    }

    /// Apply alignment-level filters to an aligned read pair.
    ///
    /// This is used after alignment-based extraction to discard reads whose
    /// alignment scores fall below the configured thresholds. If `increment` is
    /// true, the matching reason is also recorded in the filtered-read counts.
    pub fn filter_alignment(
        &mut self,
        record: &ReadPair,
        f_alignment: &Alignment,
        r_alignment: Option<&Alignment>,
        increment: bool,
    ) -> Option<FilterReason> {
        for f in ALIGNMENT_FILTERS {
            if let Some(r) = f(f_alignment, r_alignment, &self.config) {
                if increment {
                    self.increment_count(record, r);
                }
                return Some(r);
            }
        }

        None
    }

    /// Total count of filtered reads
    pub fn total(&self) -> u64 {
        self.totals.total()
    }

    /// Merge another `FilteredReads` into this one.
    ///
    /// This is primarily used when combining results from multiple counting
    /// threads. The filter configurations must be identical.
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

    /// Output filtered-read counts for iteration or output.
    ///
    /// Returns one entry per `(read sequence, read group)` combination together
    /// with its per-reason counts. If `sort` is true, rows are ordered by
    /// descending total filtered count per read sequence.
    ///
    /// This allocates a new vector of references.
    pub fn to_vector(&self, sort: bool) -> Vec<(&SeqPair, &ReadGroup, &FilteredCounts)> {
        let mut keys: Vec<(&SeqPair, u64)> = self
            .counts
            .iter()
            .map(|x| (x.0, x.1.iter().map(|y| y.1.total()).sum()))
            .collect();

        if sort {
            // Invert count to get desc order
            keys.sort_unstable_by_key(|x| std::cmp::Reverse(x.1));
        }

        keys.iter()
            .flat_map(|(k, _)| {
                let groups = self
                    .counts
                    .get(k)
                    .expect("Key missing despite coming from FilteredReads");

                groups.iter().map(move |(g, c)| (*k, g, c))
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seqs::ReadPair;
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
            group: ReadGroup::ungrouped(),
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

            let got = fr.filter_readpair(&c.rp, true);
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
                        .get(&ReadGroup::ungrouped())
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

            let got = fr.filter_alignment(&rp, &c.f, c.r.as_ref(), true);
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
                        .get(&ReadGroup::ungrouped())
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

    #[test]
    fn test_filter_reason_as_index() {
        assert_eq!(FilterReason::EmptyRead.as_index(), 0);
        assert_eq!(FilterReason::ShortRead.as_index(), 1);
        assert_eq!(FilterReason::LongRead.as_index(), 2);
        assert_eq!(FilterReason::LowMeanQuality.as_index(), 3);
        assert_eq!(FilterReason::BadAlignment.as_index(), 4);
    }

    // Smoke test that meta-data is consistent - this is in outputs so want to
    // be delibrate about changes
    #[test]
    fn test_filter_reason_meta_empty_read() {
        let meta = FilterReason::EmptyRead.meta();
        assert_eq!(meta.id, "empty_read");
        assert_eq!(meta.label, "Read length = 0");
    }

    #[test]
    fn test_filter_reason_meta_short_read() {
        let meta = FilterReason::ShortRead.meta();
        assert_eq!(meta.id, "short_read");
        assert_eq!(meta.label, "Read length < minimum");
    }

    #[test]
    fn test_filter_reason_meta_long_read() {
        let meta = FilterReason::LongRead.meta();
        assert_eq!(meta.id, "long_read");
        assert_eq!(meta.label, "Read length > maximum");
    }

    #[test]
    fn test_filter_reason_meta_low_mean_quality() {
        let meta = FilterReason::LowMeanQuality.meta();
        assert_eq!(meta.id, "low_mean_quality");
        assert_eq!(meta.label, "Mean quality < minimum");
    }

    #[test]
    fn test_filter_reason_meta_bad_alignment() {
        let meta = FilterReason::BadAlignment.meta();
        assert_eq!(meta.id, "bad_alignment");
        assert_eq!(meta.label, "Alignment quality < tolerance");
    }

    #[test]
    fn test_filter_reason_n_reasons() {
        assert_eq!(FilterReason::N_REASONS, 5);
    }

    #[test]
    fn test_filter_reason_all_filters_len() {
        assert_eq!(FilterReason::ALL_FILTERS.len(), 5);
    }

    // ---------- ALIGNMENT TOLERANCE ----------
    #[test]
    fn test_alignment_tolerance_valid() {
        let tol = AlignmentTolerance::new(0.8, 100, 120).unwrap();
        assert_eq!(tol.tolerance, 0.8);
        assert_eq!(tol.expected_f_score, 100);
        assert_eq!(tol.expected_r_score, 120);
        assert_eq!(tol.minimum_f_score, 80);
        assert_eq!(tol.minimum_r_score, 96);
    }

    #[test]
    fn test_alignment_tolerance_zero() {
        let tol = AlignmentTolerance::new(0.0, 100, 100).unwrap();
        assert_eq!(tol.minimum_f_score, 0);
        assert_eq!(tol.minimum_r_score, 0);
    }

    #[test]
    fn test_alignment_tolerance_one() {
        let tol = AlignmentTolerance::new(1.0, 100, 100).unwrap();
        assert_eq!(tol.minimum_f_score, 100);
        assert_eq!(tol.minimum_r_score, 100);
    }

    #[test]
    fn test_alignment_tolerance_below_range() {
        let result = AlignmentTolerance::new(-0.1, 100, 100);
        assert!(result.is_err());
        match result {
            Err(ReadCountError::FilterConfigError { desc }) => {
                assert!(desc.contains("between 0 and 1"));
            }
            _ => panic!("Expected FilterConfigError"),
        }
    }

    #[test]
    fn test_alignment_tolerance_above_range() {
        let result = AlignmentTolerance::new(1.1, 100, 100);
        assert!(result.is_err());
        match result {
            Err(ReadCountError::FilterConfigError { desc }) => {
                assert!(desc.contains("between 0 and 1"));
            }
            _ => panic!("Expected FilterConfigError"),
        }
    }

    #[test]
    fn test_alignment_tolerance_fractional() {
        let tol = AlignmentTolerance::new(0.5, 100, 100).unwrap();
        assert_eq!(tol.minimum_f_score, 50);
        assert_eq!(tol.minimum_r_score, 50);
    }

    // ---------- FILTERED COUNTS ----------
    #[test]
    fn test_filtered_counts_new() {
        let counts = FilteredCounts::new();
        assert_eq!(counts.total(), 0);
        for reason in FilterReason::ALL_FILTERS {
            assert_eq!(counts.get(reason), 0);
        }
    }

    #[test]
    fn test_filtered_counts_get() {
        let mut counts = FilteredCounts::new();
        counts.increment_count(FilterReason::EmptyRead);
        counts.increment_count(FilterReason::EmptyRead);
        counts.increment_count(FilterReason::ShortRead);

        assert_eq!(counts.get(&FilterReason::EmptyRead), 2);
        assert_eq!(counts.get(&FilterReason::ShortRead), 1);
        assert_eq!(counts.get(&FilterReason::LongRead), 0);
    }

    #[test]
    fn test_filtered_counts_total() {
        let mut counts = FilteredCounts::new();
        assert_eq!(counts.total(), 0);

        counts.increment_count(FilterReason::EmptyRead);
        assert_eq!(counts.total(), 1);

        counts.increment_count(FilterReason::ShortRead);
        counts.increment_count(FilterReason::LongRead);
        assert_eq!(counts.total(), 3);
    }

    #[test]
    fn test_filtered_counts_iter() {
        let mut counts = FilteredCounts::new();
        counts.increment_count(FilterReason::EmptyRead);
        counts.increment_count(FilterReason::ShortRead);

        let vec: Vec<u64> = counts.iter().copied().collect();
        assert_eq!(vec[0], 1); // empty
        assert_eq!(vec[1], 1); // short
        assert_eq!(vec[2], 0); // long
        assert_eq!(vec[3], 0); // quality
        assert_eq!(vec[4], 0); // alignment
    }

    #[test]
    fn test_filtered_counts_merge() {
        let mut counts1 = FilteredCounts::new();
        counts1.increment_count(FilterReason::EmptyRead);
        counts1.increment_count(FilterReason::ShortRead);

        let mut counts2 = FilteredCounts::new();
        counts2.increment_count(FilterReason::EmptyRead);
        counts2.increment_count(FilterReason::LongRead);

        counts1.merge(counts2);
        assert_eq!(counts1.get(&FilterReason::EmptyRead), 2);
        assert_eq!(counts1.get(&FilterReason::ShortRead), 1);
        assert_eq!(counts1.get(&FilterReason::LongRead), 1);
    }

    #[test]
    fn test_filtered_counts_default() {
        let counts = FilteredCounts::default();
        assert_eq!(counts.total(), 0);
    }

    // ---------- FILTERED READS ----------
    #[test]
    fn test_filtered_reads_new() {
        let cfg = FilterConfig::new(None, None, None, None, false);
        let fr = FilteredReads::new(cfg.clone());
        assert_eq!(fr.config, cfg);
        assert_eq!(fr.total(), 0);
        assert!(fr.counts.is_empty());
    }

    #[test]
    fn test_filtered_reads_increment_new_read() {
        let mut fr = FilteredReads::new(FilterConfig::new(None, None, None, None, false));
        let rp = rp("ACTG", "FFFF", None, None);

        fr.increment_count(&rp, FilterReason::EmptyRead);

        assert_eq!(fr.total(), 1);
        let key = rp.key();
        assert!(fr.counts.contains_key(&key));
    }

    #[test]
    fn test_filtered_reads_increment_existing_read() {
        let mut fr = FilteredReads::new(FilterConfig::new(None, None, None, None, false));
        let rp = rp("ACTG", "FFFF", None, None);

        fr.increment_count(&rp, FilterReason::EmptyRead);
        fr.increment_count(&rp, FilterReason::EmptyRead);
        fr.increment_count(&rp, FilterReason::ShortRead);

        assert_eq!(fr.total(), 3);
        let key = rp.key();
        let groups = fr.counts.get(&key).unwrap();
        let counts = groups.get(&ReadGroup::ungrouped()).unwrap();
        assert_eq!(counts.get(&FilterReason::EmptyRead), 2);
        assert_eq!(counts.get(&FilterReason::ShortRead), 1);
    }

    #[test]
    fn test_filtered_reads_increment_different_groups() {
        let mut fr = FilteredReads::new(FilterConfig::new(None, None, None, None, false));
        let rp1 = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("f", None, b"ACTG", b"FFFF"),
            reverse: None,
            group: ReadGroup::grouped("0"),
        };
        let rp2 = ReadPair {
            forward: bio::io::fastq::Record::with_attrs("f", None, b"ACTG", b"FFFF"),
            reverse: None,
            group: ReadGroup::grouped("1"),
        };

        fr.increment_count(&rp1, FilterReason::EmptyRead);
        fr.increment_count(&rp2, FilterReason::EmptyRead);

        assert_eq!(fr.total(), 2);
        let key = rp1.key();
        let groups = fr.counts.get(&key).unwrap();
        assert_eq!(groups.len(), 2);
    }

    #[test]
    fn test_filtered_reads_filter_readpair_no_increment() {
        let mut fr = FilteredReads::new(FilterConfig::new(None, None, None, None, true));
        let rp = rp("", "", None, None);

        let result = fr.filter_readpair(&rp, false);
        assert_eq!(result, Some(FilterReason::EmptyRead));
        assert_eq!(fr.total(), 0); // Not incremented
    }

    #[test]
    fn test_filtered_reads_filter_alignment_no_increment() {
        let cfg = FilterConfig::new(
            None,
            Some(AlignmentTolerance::new(0.8, 100, 100).unwrap()),
            None,
            None,
            false,
        );
        let mut fr = FilteredReads::new(cfg);
        let rp = rp("ACTG", "FFFF", None, None);
        let aln = Alignment {
            score: 50,
            ..Default::default()
        };

        let result = fr.filter_alignment(&rp, &aln, None, false);
        assert_eq!(result, Some(FilterReason::BadAlignment));
        assert_eq!(fr.total(), 0); // Not incremented
    }

    #[test]
    fn test_filtered_reads_merge_same_config() {
        let cfg = FilterConfig::new(None, None, None, None, false);
        let mut fr1 = FilteredReads::new(cfg.clone());
        let mut fr2 = FilteredReads::new(cfg.clone());

        let rp = rp("ACTG", "FFFF", None, None);
        fr1.increment_count(&rp, FilterReason::EmptyRead);
        fr2.increment_count(&rp, FilterReason::ShortRead);

        let result = fr1.merge(fr2);
        assert!(result.is_ok());
        assert_eq!(fr1.total(), 2);
    }

    #[test]
    fn test_filtered_reads_merge_different_config() {
        let cfg1 = FilterConfig::new(None, None, None, None, false);
        let cfg2 = FilterConfig::new(Some(20.0), None, None, None, false);

        let mut fr1 = FilteredReads::new(cfg1);
        let fr2 = FilteredReads::new(cfg2);

        let result = fr1.merge(fr2);
        assert!(result.is_err());
    }

    #[test]
    fn test_filtered_reads_to_vector_empty() {
        let fr = FilteredReads::new(FilterConfig::new(None, None, None, None, false));
        let vec = fr.to_vector(false);
        assert!(vec.is_empty());
    }

    #[test]
    fn test_filtered_reads_to_vector_unsorted() {
        let mut fr = FilteredReads::new(FilterConfig::new(None, None, None, None, false));
        let rp1 = rp("ACTG", "FFFF", None, None);
        let rp2 = rp("GGGG", "FFFF", None, None);

        fr.increment_count(&rp1, FilterReason::EmptyRead);
        fr.increment_count(&rp1, FilterReason::EmptyRead);
        fr.increment_count(&rp2, FilterReason::ShortRead);

        let vec = fr.to_vector(false);
        assert_eq!(vec.len(), 2);
    }

    #[test]
    fn test_filtered_reads_to_vector_sorted() {
        let mut fr = FilteredReads::new(FilterConfig::new(None, None, None, None, false));
        let rp1 = rp("ACTG", "FFFF", None, None);
        let rp2 = rp("GGGG", "FFFF", None, None);

        fr.increment_count(&rp1, FilterReason::EmptyRead);
        fr.increment_count(&rp1, FilterReason::ShortRead);
        fr.increment_count(&rp2, FilterReason::LongRead);

        let vec = fr.to_vector(true);
        assert_eq!(vec.len(), 2);
        // First entry should be rp1 with total 2, second rp2 with total 1
        assert_eq!(vec[0].2.total(), 2);
        assert_eq!(vec[1].2.total(), 1);
    }

    #[test]
    fn test_filter_config_equality() {
        let cfg1 = FilterConfig::new(Some(20.0), None, Some(50), Some(300), true);
        let cfg2 = FilterConfig::new(Some(20.0), None, Some(50), Some(300), true);
        assert_eq!(cfg1, cfg2);
    }

    #[test]
    fn test_filter_config_inequality() {
        let cfg1 = FilterConfig::new(Some(20.0), None, Some(50), Some(300), true);
        let cfg2 = FilterConfig::new(Some(30.0), None, Some(50), Some(300), true);
        assert_ne!(cfg1, cfg2);
    }
}
