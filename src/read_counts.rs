//! Counting the occurance of different reads in sequence files
//!
//! Contains structs and methods for counting combinations of expected
//! regions in DNA sequence input as well as comparing them to expected
//! library members as defiend in a LibSpec object. Supports multiple
//! approaches for extracting regions of interest from the input sequence:
//! alignment, pattern matching, inframe position matching and full read
//! counting.
use anyhow::{self};
use bio::alignment::AlignmentOperation;
use bio::alignment::pairwise::{Aligner, MatchFunc, Scoring};
use bio::alphabets::dna::revcomp;
use bio::bio_types::{alignment::Alignment, sequence::Sequence};
use clap::ValueEnum;
use itertools::{Itertools, izip};
use log::{info, debug, warn};
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::iter::zip;
use std::rc::Rc;
use std::sync::Once;

use crate::lib_spec::{
    DistanceMetric, LibSpecError, Library, LibraryError, LibraryRegion, LibrarySpec,
    PartialMatching, merge_matches,
};
use crate::log_progress::{Progress, ProgressStyle};
use crate::read_parsing::{ReadGroup, ReadKey, ReadPair, ReadPairParser};

/// Position in an alignment where a region is found
///
/// Has the format (start, stop, RegionCompleteness status)
pub type AlignmentPosition = (usize, usize, RegionCompleteness);

static WARN_MERGE: Once = Once::new();

// Module consts
/// Message for logging region extraction progress
const PROG_MSG: &str = "Extracting regions:";
/// Message emitted at the end of region extraction
const FINAL_MSG: &str = "Extracted regions:";

/// Single end align mode message
const ALIGN_START_SINGLE_MSG: &str = "Extracting regions from single end reads by alignment";
/// Paired end align mode message
const ALIGN_START_PAIRED_MSG: &str = "Extracting regions from paired end reads by alignment";
/// Logging interval when aligning
#[cfg(debug_assertions)]
const ALIGN_LOG_INTERVAL: u64 = 1000;
#[cfg(not(debug_assertions))]
const ALIGN_LOG_INTERVAL: u64 = 250000;

/// Single end pattern match mode message
const PATTERN_START_SINGLE_MSG: &str =
    "Extracting regions from single end reads by pattern matching";
/// Paired end pattern match mode message
const PATTERN_START_PAIRED_MSG: &str =
    "Extracting regions from paired end reads by pattern matching";
/// Logging interval when pattern matching
#[cfg(debug_assertions)]
const PATTERN_LOG_INTERVAL: u64 = 1000;
#[cfg(not(debug_assertions))]
const PATTERN_LOG_INTERVAL: u64 = 2000000;

/// Single end inframe mode message
const INFRAME_START_SINGLE_MSG: &str = "Extracting regions from single end reads in-frame";
/// Paired end inframe mode message
const INFRAME_START_PAIRED_MSG: &str = "Extracting regions from paired end reads in-frame";
/// Logging interval in inframe mode
#[cfg(debug_assertions)]
const INFRAME_LOG_INTERVAL: u64 = 1000;
#[cfg(not(debug_assertions))]
const INFRAME_LOG_INTERVAL: u64 = 5000000;

/// Single end raw mode message
const RAW_START_SINGLE_MSG: &str = "Counting raw single end reads";
/// Paired end raw mode message
const RAW_START_PAIRED_MSG: &str = "Counting raw paired end reads";
/// Logging interval in raw mode
#[cfg(debug_assertions)]
const RAW_LOG_INTERVAL: u64 = 1000;
#[cfg(not(debug_assertions))]
const RAW_LOG_INTERVAL: u64 = 10000000;

/// Customisable alignment scoring scheme allowing Ns
#[derive(Debug, Copy, Clone)]
pub struct AlignmentScorer {
    std_match: i32,
    n_match: i32,
    mismatch: i32,
    gap_open: i32,
    gap_extend: i32,
}

impl AlignmentScorer {
    pub fn new(
        std_match: i32,
        n_match: i32,
        mismatch: i32,
        gap_open: i32,
        gap_extend: i32,
    ) -> Self {
        Self {
            std_match,
            n_match,
            mismatch,
            gap_open,
            gap_extend,
        }
    }

    /// Generate a matching Scoring object to use with Rust Bio alignment
    pub fn get_scoring(self) -> Scoring<AlignmentScorer> {
        Scoring::new(self.gap_open, self.gap_extend, self)
    }
}

impl MatchFunc for AlignmentScorer {
    /// Alignment match scores allowing Ns
    ///
    /// Return a (mis)match score that allows alignment of anything against N
    /// with a moderate penalty. Penalty is greater than a mismatch but less
    /// than a gap to account for possible sequencing errors before variable
    /// regions, which otherwise get shunted into the N region
    fn score(&self, a: u8, b: u8) -> i32 {
        if a == b'N' || b == b'N' {
            self.n_match
        } else if a == b {
            self.std_match
        } else {
            self.mismatch
        }
    }
}

/// Extract the regions of a query sequence that match template sections by walking
/// an alignment path and region position vector together
///
/// Returns a vector of AlignmentPostions the same length as region_positions where each entry gives
/// the corresponding position in the query sequence plus a RegionCompleteness status
fn regions_from_alignment_path(
    region_positions: &[(usize, usize)],
    alignment_path: &[(usize, usize, AlignmentOperation)],
) -> Result<Vec<Option<AlignmentPosition>>, ReadCountError> {
    // Walk the alignment path and region list, adding each region when it completes
    let mut aln_idx: usize = 0;
    let mut seq_start: usize = 1; // AlignmentPath is 1 based
    let mut out_regions: Vec<Option<AlignmentPosition>> = vec![None; region_positions.len()];

    let mut reg_idx: usize = match region_positions
        .iter()
        .position(|x| x.1 >= alignment_path[0].1)
    {
        Some(x) => x,
        None => {
            // abort read if no region overlaps start of alignment, incrementing empty comb
            return Ok(out_regions);
        }
    };

    // If we start strictly inside a region the first is incomplete, otherwise complete
    let mut current_completeness = if alignment_path[0].1 > region_positions[reg_idx].0 {
        RegionCompleteness::Partial5Prime
    } else {
        RegionCompleteness::Complete
    };

    loop {
        // Operate based on current position in alignment/regions
        if alignment_path[aln_idx].1 < region_positions[reg_idx].0 {
            // Before current region
            aln_idx += 1;
            match alignment_path[aln_idx].2 {
                // Increment seq_start for match/subst
                AlignmentOperation::Match | AlignmentOperation::Subst => {
                    seq_start = alignment_path[aln_idx].0
                }
                // But not for del (it doesn't change) and ins (we want to include the start of
                // the insert in any putitive region)
                AlignmentOperation::Del | AlignmentOperation::Ins => {}
                AlignmentOperation::Xclip(_) | AlignmentOperation::Yclip(_) => {
                    return Err(ReadCountError::Error {
                        desc: "Unexpected end clipping in alignment".to_string(),
                    });
                }
            };
        } else if (alignment_path[aln_idx].1 >= region_positions[reg_idx].0)
            && (alignment_path[aln_idx].1 < region_positions[reg_idx].1)
        {
            // A del at the start of a region has a matching template position but preceeding
            // observed position => increment it
            if (alignment_path[aln_idx].1 == region_positions[reg_idx].0)
                && matches!(alignment_path[aln_idx].2, AlignmentOperation::Del)
            {
                seq_start += 1
            }
            // Inside current region
            aln_idx += 1;
        } else if alignment_path[aln_idx].1 >= region_positions[reg_idx].1 {
            // After current region
            if seq_start < alignment_path[aln_idx].0 {
                // If this isn't the case it indicates the whole region is
                // deleted in the alignment and so should be left as None
                out_regions[reg_idx] =
                    Some((seq_start, alignment_path[aln_idx].0, current_completeness));
            }
            current_completeness = RegionCompleteness::Complete;
            aln_idx += 1;
            reg_idx += 1;
        } else {
            return Err(ReadCountError::Error {
                desc: "Region/alignment out of sync while walking (this should \
                       be impossible...)"
                    .to_string(),
            });
        }

        // TODO - need to work out partial matches properly (both alignment end and region expected length?) and need to deal with different alignment opperations

        // Check break conditions
        if reg_idx == region_positions.len() {
            // Exhausted regions => break
            break;
        } else if aln_idx == alignment_path.len() - 1 {
            // Exhausted alignment with possibly open region =>
            // conclude possibly partial region and break
            if alignment_path[aln_idx].1 >= region_positions[reg_idx].0 {
                // Started current region
                if (alignment_path[aln_idx].1 - seq_start)
                    == (region_positions[reg_idx].1 - region_positions[reg_idx].0)
                {
                    // Length is correct, so complete
                    out_regions[reg_idx] = Some((
                        seq_start,
                        alignment_path[aln_idx].0,
                        RegionCompleteness::Complete,
                    ));
                } else {
                    // Shorter than expected to incomplete
                    out_regions[reg_idx] = Some((
                        seq_start,
                        alignment_path[aln_idx].0,
                        RegionCompleteness::Partial3Prime,
                    ));
                }
            }
            break;
        }
    }

    Ok(out_regions)
}

/// Convert a `Vec<u8>` Sequence to a string, logging failure but not panicing
///
/// This is useful for writing output files so that bad UTF8 is flagged but doesn't
/// abort the whole write, meaning the user can more easily observe what has occured
/// in combination with the warnings. In theory this should rarely occur with good input
/// and bad input should be caught earlier.
fn seq_to_string_or_log(seq: &Sequence) -> String {
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

/// Calculate the mean of a fastq quality vector
fn mean_quality(qual: &[u8]) -> f32 {
    let total: u32 = qual.iter().fold(0, |a, e| a + *e as u32);
    total as f32 / qual.len() as f32
}

/// Merge a forward and reverse sequence in a region
///
/// Looks at each position in turn and takes the highest quality option, if any.
/// Args are tuples with sequence, quality vector and a completeness tag, plus the
/// expected region length. This only works for fixed length, will need something
/// more complete if there is variable overlap (plus should really warn to use a read
/// merger at that point).
fn merge_seqs(
    fwd: Option<(Sequence, Vec<u8>, RegionCompleteness)>,
    rev: Option<(Sequence, Vec<u8>, RegionCompleteness)>,
    len: usize,
) -> Result<Option<(Sequence, RegionCompleteness)>, anyhow::Error> {
    // Deal with simple cases first
    let (f_reg, r_reg) = match (fwd, rev) {
        // Neither present
        (None, None) => return Ok(None),

        // Only rev
        (None, Some(r)) => return Ok(Some((r.0, r.2))),

        // Only fwd
        (Some(f), None) => return Ok(Some((f.0, f.2))),

        // Both, must decide which or possibly combine
        (Some(f), Some(r)) => (f, r),
    };

    // Choose path forward based on completeness flags
    // Overlap/missing centre shouldn't occur here as they imply it has already been merged
    match (f_reg.2, r_reg.2) {
        // Error options
        (RegionCompleteness::MissingCenter { .. }, _)
        | (_, RegionCompleteness::MissingCenter { .. }) => Err(ReadCountError::Error {
            desc: "merge_seqs passed MissingCenter regions, implying already merged".to_string(),
        }
        .into()),

        (RegionCompleteness::Overlapping { .. }, _)
        | (_, RegionCompleteness::Overlapping { .. }) => Err(ReadCountError::Error {
            desc: "merge_seqs passed Overlapping regions, implying already merged".to_string(),
        }
        .into()),

        // Both complete - use highest quality
        (RegionCompleteness::Complete, RegionCompleteness::Complete) => {
            WARN_MERGE.call_once(|| {
                log::warn!(
                    "F/R reads overlap. Consider merging reads first (e.g. using Pear) as this will compare and combine them much more robustly than the simple approach used here."
                );
            });

            if mean_quality(&f_reg.1) >= mean_quality(&r_reg.1) {
                Ok(Some((f_reg.0, RegionCompleteness::Complete)))
            } else {
                Ok(Some((r_reg.0, RegionCompleteness::Complete)))
            }
        }

        // One complete - use the complete one
        (RegionCompleteness::Complete, _) => Ok(Some((f_reg.0, f_reg.2))),
        (_, RegionCompleteness::Complete) => Ok(Some((r_reg.0, r_reg.2))),

        // Both partial in same way - use highest quality
        // Weird situation but some read trimming could lead her potentially
        (RegionCompleteness::Partial5Prime, RegionCompleteness::Partial5Prime) => {
            WARN_MERGE.call_once(|| {
                log::warn!(
                    "F/R reads overlap. Consider merging reads first (e.g. using Pear) as this will compare and combine them much more robustly than the simple approach used here."
                );
            });

            if mean_quality(&f_reg.1) >= mean_quality(&r_reg.1) {
                Ok(Some((f_reg.0, RegionCompleteness::Partial5Prime)))
            } else {
                Ok(Some((r_reg.0, RegionCompleteness::Partial5Prime)))
            }
        }
        (RegionCompleteness::Partial3Prime, RegionCompleteness::Partial3Prime) => {
            WARN_MERGE.call_once(|| {
                log::warn!(
                    "F/R reads overlap. Consider merging reads first (e.g. using Pear) as this will compare and combine them much more robustly than the simple approach used here."
                );
            });

            if mean_quality(&f_reg.1) >= mean_quality(&r_reg.1) {
                Ok(Some((f_reg.0, RegionCompleteness::Partial3Prime)))
            } else {
                Ok(Some((r_reg.0, RegionCompleteness::Partial3Prime)))
            }
        }

        // Both partial with gap - missing centre or overlapping based simply on max region length
        (RegionCompleteness::Partial5Prime, RegionCompleteness::Partial3Prime) => {
            let f_end = f_reg.0.len(); // `f` covers [0..f_end)
            let r_start = len - r_reg.0.len(); // `r` covers [r_start..seq_len)

            let mut out_seq = Vec::with_capacity(f_reg.0.len() + 1 + r_reg.0.len());
            out_seq.extend_from_slice(&f_reg.0);
            out_seq.push(b'/');
            out_seq.extend_from_slice(&r_reg.0);

            if r_start >= f_end {
                // No overlap
                Ok(Some((
                    out_seq,
                    RegionCompleteness::MissingCenter {
                        split_ind: f_reg.0.len(),
                    },
                )))
            } else {
                // Overlap
                WARN_MERGE.call_once(|| {
                    log::warn!(
                        "F/R reads overlap. Consider merging reads first (e.g. using Pear) as this will compare and combine them much more robustly than the simple approach used here."
                    );
                });

                Ok(Some((
                    out_seq,
                    RegionCompleteness::Overlapping {
                        split_ind: f_reg.0.len(),
                    },
                )))
            }
        }
        (RegionCompleteness::Partial3Prime, RegionCompleteness::Partial5Prime) => {
            // Reversed, should rarely see this case but some odd trimming could create in theory
            // Basically the reverse of above
            let r_end = r_reg.0.len(); // `r` covers [0..r_end)
            let f_start = len - f_reg.0.len(); // `f` covers [f_start..seq_len)

            let mut out_seq = Vec::with_capacity(r_reg.0.len() + 1 + f_reg.0.len());
            out_seq.extend_from_slice(&r_reg.0);
            out_seq.push(b'/');
            out_seq.extend_from_slice(&f_reg.0);

            if f_start >= r_end {
                // No overlap
                Ok(Some((
                    out_seq,
                    RegionCompleteness::MissingCenter {
                        split_ind: r_reg.0.len(),
                    },
                )))
            } else {
                // Overlap
                WARN_MERGE.call_once(|| {
                    log::warn!(
                        "F/R reads overlap. Consider merging reads first (e.g. using Pear) as this will compare and combine them much more robustly than the simple approach used here."
                    );
                });

                Ok(Some((
                    out_seq,
                    RegionCompleteness::Overlapping {
                        split_ind: r_reg.0.len(),
                    },
                )))
            }
        }
    }
}

/// Error type for read counting
///
/// Mostly the generic ReadCountError since the CLI doesn't need to differentiate much.
/// UnexpectedRegionError is included for ergonomics and clarity.
#[derive(Debug)]
pub enum ReadCountError {
    UnexpectedRegion { region: String },
    FilterConfigError { desc: String },
    BadAlignment { alignment: Box<AlignmentInfo> },
    Error { desc: String },
}

#[derive(Debug)]
pub struct AlignmentInfo {
    read_id: String,
    read_number: usize,
    pretty_alignment: String,
    alignment: Alignment,
    region_ids: Vec<String>,
    region_positions: Vec<(usize, usize)>,
    mapped_positions: Vec<Option<(usize, usize, RegionCompleteness)>>,
}

impl fmt::Display for ReadCountError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ReadCountError::UnexpectedRegion { region } => {
                write!(
                    f,
                    "Added combination contains an unexpected region: {}",
                    region
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

// Key types
/// Region keys identify via the region name, the observed sequence and completeness status
type RegionKey = (String, Sequence, RegionCompleteness);

/// Combination keys identify via a vector of RegionKeys
type CombinationKey = Vec<RegionKey>;

/// Library keys identify a particular library match type
type LibraryKey = Vec<(String, RegionMatch)>;

/// HashMap cache of observed reads and which combination they map to
pub type ObservedReads = HashMap<ReadKey, CacheHit>;

/// Options to cache for each identified read
///
/// The cache operates at a sequence level only, so you shouldn't cache filtering related
/// to quality
pub enum CacheHit {
    Comb(CombinationKey),
    Filter(FilterReason),
}

/// Container for ObservedCombinations
///
/// The core is a HashMap of the ObservedRegions seen so far and a HashMap
/// of ObservedCombination objects that link to the contained regions.
/// An additional HashMap plus the library is added when library comparison
/// is run. It also carries the region ids to be considered in order.
#[derive(Debug)]
pub struct ObservedCombinations {
    region_ids: Vec<String>,
    regions: HashMap<RegionKey, Rc<RefCell<ObservedRegion>>>,
    combinations: HashMap<CombinationKey, ObservedCombination>,
    library: Option<Library>,
    library_combinations: Option<HashMap<LibraryKey, LibraryCombination>>,
    filtered_reads: FilteredReads,
}

impl ObservedCombinations {
    fn new(region_ids: Vec<String>, filter_config: FilterConfig) -> Self {
        Self {
            region_ids,
            regions: HashMap::new(),
            combinations: HashMap::new(),
            library: None,
            library_combinations: None,
            filtered_reads: FilteredReads::new(filter_config),
        }
    }

    /// Number of distinct observed combination types
    pub fn len(&self) -> usize {
        self.combinations.len()
    }

    /// Increment a combination count or add a new combination if it hasn't been seen yet
    fn add_or_increment_combination(
        &mut self,
        key: &CombinationKey,
        group: ReadGroup,
    ) -> Result<(), anyhow::Error> {
        match self.combinations.get_mut(key) {
            Some(comb) => comb.increment_count(group),
            None => {
                let mut reg_map = HashMap::new();

                for reg_key in key {
                    if !self.region_ids.contains(&reg_key.0) {
                        return Err(ReadCountError::UnexpectedRegion {
                            region: reg_key.0.clone(),
                        }
                        .into());
                    }

                    match self.regions.get(reg_key) {
                        None => {
                            let new_reg = Rc::new(RefCell::new(ObservedRegion::new(
                                reg_key.0.clone(),
                                &reg_key.1,
                                reg_key.2,
                            )));
                            self.regions.insert(reg_key.clone(), new_reg.clone());
                            reg_map.insert(reg_key.0.clone(), new_reg.clone());
                        }
                        Some(r) => {
                            reg_map.insert(r.borrow().id.clone(), r.clone());
                        }
                    }
                }

                let mut comb = ObservedCombination::new(reg_map);
                comb.increment_count(group);
                self.combinations.insert(key.clone(), comb);
            }
        }

        Ok(())
    }

    /// Update filter counts without checking the read
    ///
    /// Passes through to self.filtered_reads.update_count, useful when using
    /// cached FilterReasons to prevent needing to re-align.
    fn update_filter_count(&mut self, reason: &FilterReason) {
        self.filtered_reads.update_count(reason)
    }

    /// Determine if a read should be filtered
    ///
    /// Pass through to self.filtered_reads.filter_read, which Checks whether the read should
    /// be filtered, adding it to the appropriate count if is, and returns a bool determining
    /// if it was filtered.
    fn filter_readpair(&mut self, record: &ReadPair) -> FilterReason {
        self.filtered_reads.filter_readpair(record)
    }

    /// Determine if an alignment should be filtered
    ///
    /// Passes through to self.filtered_reads.filter_alignment, which checks if an alignment should
    /// be filtered, adding it to the appropriate count if so, and returns a bool determining if it
    /// was filtered
    fn filter_alignment(
        &mut self,
        f_alignment: &Alignment,
        r_alignment: Option<&Alignment>,
    ) -> FilterReason {
        self.filtered_reads
            .filter_alignment(f_alignment, r_alignment)
    }

    /// Compare observed combinations to those expected in a Library
    pub fn compare_to_library(
        &mut self,
        library: Library,
        progress_style: Option<&ProgressStyle>,
        distance_metric: DistanceMetric,
        max_matches: usize,
    ) -> Result<(), LibraryError> {
        let n_regs = self.regions.len() as u64;
        let n_combs = self.combinations.len() as u64;

        // Compare each region to the library
        let mut reg_progress: Progress = Progress::from_style(
            progress_style.unwrap_or(&ProgressStyle::new(None)),
            "Matching regions:",
            "Matched regions:",
            Some(n_regs),
            match distance_metric {
                DistanceMetric::Hamming | DistanceMetric::Exact => 250000,
                DistanceMetric::BoundedLevenshtein => std::cmp::max(n_regs / 10, 50000),
                DistanceMetric::Levenshtein => std::cmp::max(n_regs / 20, 10000),
            },
        );

        for value in self.regions.values() {
            let val = value
                .borrow()
                .compare_to_library(&library, distance_metric, max_matches);
            value.borrow_mut().nearest_matches = val;
            reg_progress.inc(1);
        }
        reg_progress.finish();

        // Compare the combinations to the libary
        let mut comb_progress: Progress = Progress::from_style(
            progress_style.unwrap_or(&ProgressStyle::new(None)),
            "Comparing combinations:",
            "Compared combinations:",
            Some(n_combs),
            2000000,
        );

        for value in self.combinations.values_mut() {
            value.library_matches =
                value.compare_to_library(&self.region_ids, &library, distance_metric, max_matches);
            comb_progress.inc(1);
        }
        comb_progress.finish();

        // Make summary counts of observed library combinations
        let mut lib_summary_progress: Progress = Progress::from_style(
            progress_style.unwrap_or(&ProgressStyle::new(None)),
            "Summarising library matches:",
            "Summarised library matches:",
            Some(n_combs),
            std::cmp::max(n_combs / 4, 250000),
        );

        let mut lib_combs: HashMap<Vec<(String, RegionMatch)>, LibraryCombination> = HashMap::new();

        for comb in self.combinations.values_mut() {
            let mut key: LibraryKey = Vec::with_capacity(5);

            for reg in &self.region_ids {
                // Add the appropriate library match to the key. Simplify over distance,
                // match count, etc. to summarise the library and store NoLibrary
                // region seqs to capture e.g. barcodes.
                match comb.regions.get(reg) {
                    None => key.push((reg.to_string(), RegionMatch::Unmatched)),
                    Some(x) => key.push((
                        reg.to_string(),
                        match &x.borrow().nearest_matches {
                            RegionMatch::Unmatched => RegionMatch::Unmatched,
                            RegionMatch::Overmatched { .. } => RegionMatch::Overmatched {
                                distance: 0,
                                matches: 0,
                            },
                            RegionMatch::NoLibrary { .. } => RegionMatch::NoLibrary {
                                seq: Some(x.borrow().seq.clone()),
                            },
                            RegionMatch::Uncompared => RegionMatch::Uncompared,
                            RegionMatch::Match { seq_match, .. } => RegionMatch::Match {
                                seq_match: seq_match.clone(),
                                distance: 0,
                            },
                            RegionMatch::MultiMatch { seq_matches, .. } => {
                                RegionMatch::MultiMatch {
                                    seq_matches: seq_matches.to_vec(),
                                    distance: 0,
                                }
                            }
                        },
                    )),
                }
            }

            match lib_combs.get_mut(&key) {
                Some(x) => {
                    for (group, count) in &comb.counts {
                        x.increment_count(group, *count);
                    }
                }
                None => {
                    let mut x: LibraryCombination = LibraryCombination::new(
                        HashMap::from_iter(key.clone()),
                        match &comb.library_matches {
                            CombinationMatch::Uncompared => CombinationMatch::Uncompared,
                            CombinationMatch::Match { ind, .. } => CombinationMatch::Match {
                                ind: *ind,
                                distance: 0,
                            },
                            CombinationMatch::MultiMatch { inds, .. } => {
                                CombinationMatch::MultiMatch {
                                    inds: inds.clone(),
                                    distance: 0,
                                }
                            }
                            CombinationMatch::Recombination { .. } => {
                                CombinationMatch::Recombination { distance: 0 }
                            }
                            CombinationMatch::Mismatch => CombinationMatch::Mismatch,
                            CombinationMatch::Nonmatch => CombinationMatch::Nonmatch,
                        },
                    );
                    for (group, count) in &comb.counts {
                        x.increment_count(group, *count);
                    }
                    lib_combs.insert(key, x);
                }
            }

            lib_summary_progress.inc(1);
        }
        lib_summary_progress.finish();

        self.library = Some(library);
        self.library_combinations = Some(lib_combs);

        Ok(())
    }

    /// Check if library comparison has occured
    pub fn is_compared_to_library(&self) -> bool {
        self.library.is_some()
    }

    /// Summarise the observed  read count categories
    ///
    /// Counts the occurance of each CombinationMatch, so makes little sense if library comparison
    /// hasn't occured first as it will just sum the total reads
    pub fn summarise(&self) -> ReadSummary {
        let mut read_summary = ReadSummary::empty();

        read_summary.filtered_reads = Some(self.filtered_reads.clone());

        for comb in self.combinations.values() {
            let count: u64 = comb.total_count() as u64;

            match comb.library_matches {
                CombinationMatch::Uncompared => read_summary.uncompared += count,
                CombinationMatch::Match { distance, .. } => {
                    if distance == 0 {
                        read_summary.exact_match += count
                    } else {
                        read_summary.nearest_match += count
                    }
                }
                CombinationMatch::MultiMatch { .. } => read_summary.multimatch += count,
                CombinationMatch::Recombination { distance } => {
                    if distance == 0 {
                        read_summary.exact_recombination += count
                    } else {
                        read_summary.nearest_recombination += count
                    }
                }
                CombinationMatch::Mismatch => read_summary.mismatch += count,
                CombinationMatch::Nonmatch => read_summary.nonmatch += count,
            }
        }

        read_summary
    }

    /// Write all counts to file
    pub fn write_tsv(&self, file: File, sort: bool) -> Result<(), anyhow::Error> {
        let mut count_writer = BufWriter::new(file);
        let mut keys: Vec<(&CombinationKey, u32)> = self
            .combinations
            .iter()
            .map(|x| (x.0, x.1.total_count()))
            .collect();

        if sort {
            // Invert count to get desc order
            keys.sort_unstable_by_key(|x| 0 - i64::from(x.1));
        }

        // Write header
        write!(count_writer, "group\t")?;
        for r in &self.region_ids {
            write!(
                count_writer,
                "{}\t{}_nearest\t{}_distance\t{}_n_matches\t",
                r, r, r, r
            )?;
        }
        writeln!(
            count_writer,
            "combination_status\tcombination_distance\tcombinations_in_library\tcombination_indexes\tcount"
        )?;

        for (key, _) in keys {
            let combination = self.combinations.get(key).expect(
                "Combination key from extracted key list missing from ObservedCombinations",
            );
            write!(count_writer, "{}", combination.to_tsv(&self.region_ids))?;
        }

        count_writer.flush()?;
        Ok(())
    }

    /// Write matched library combination counts to file
    pub fn write_summary_tsv(&self, file: File, sort: bool) -> Result<(), anyhow::Error> {
        let combs = match &self.library_combinations {
            None => {
                return Err(ReadCountError::Error {
                    desc: "Combinations uncompared, compare before summarising".to_string(),
                }
                .into());
            }
            Some(x) => x,
        };

        let mut writer = BufWriter::new(file);
        let mut keys: Vec<(&LibraryKey, u32)> =
            combs.iter().map(|x| (x.0, x.1.total_count())).collect();

        // Sort
        if sort {
            // Invert count to get desc order
            keys.sort_unstable_by_key(|x| 0 - i64::from(x.1));
        }

        // Write
        write!(writer, "group\t")?;
        for r in &self.region_ids {
            write!(writer, "{}\t", r)?;
        }
        writeln!(
            writer,
            "combination_status\tcombinations_in_library\tcombination_indexes\tcount"
        )?;

        for (key, _) in keys {
            let combination = combs.get(key).expect(
                "Combination key from extracted key list missing from ObservedCombinations",
            );
            write!(writer, "{}", combination.to_tsv(&self.region_ids))?;
        }

        writer.flush()?;
        Ok(())
    }
}

/// Container for filtered reads
///
/// Keeps track of counts for filtered reads during counting
#[derive(Debug, Clone)]
pub struct FilteredReads {
    config: FilterConfig,
    low_mean_quality: u64,
    bad_alignment: u64,
}

pub enum FilterReason {
    None,
    LowMeanQuality,
    BadAlignment,
}

impl FilteredReads {
    fn new(config: FilterConfig) -> Self {
        Self {
            config,
            low_mean_quality: 0,
            bad_alignment: 0,
        }
    }

    /// Update the count without checking a read
    ///
    /// This lets you update counts if manually checking or when
    /// using a cached result
    pub fn update_count(&mut self, reason: &FilterReason) {
        match reason {
            FilterReason::None => {}
            FilterReason::LowMeanQuality => self.low_mean_quality += 1,
            FilterReason::BadAlignment => self.bad_alignment += 1,
        };
    }

    /// Determine if a readpair should be filtered based on the supplied config
    ///
    /// Checks whether the read should be filtered, adding it to the appropriate count if so, and
    /// returns a bool determining if it was filtered.
    fn filter_readpair(&mut self, record: &ReadPair) -> FilterReason {
        let f_read = &record.forward;
        let r_read = match &record.reverse {
            Some(r) => Some(r),
            None => None,
        };

        // Check mean quality is high enough
        match (self.config.mean_quality_threshold, r_read) {
            (None, _) => {}
            (Some(q), None) => {
                if mean_quality(f_read.qual()) < q {
                    self.low_mean_quality += 1;
                    return FilterReason::LowMeanQuality;
                }
            }
            (Some(q), Some(r)) => {
                if mean_quality(f_read.qual()) < q || mean_quality(r.qual()) < q {
                    self.low_mean_quality += 1;
                    return FilterReason::LowMeanQuality;
                }
            }
        }

        FilterReason::None
    }

    /// Determine if an alignment should be filtered based on the supplied config
    ///
    /// Checks if an alignment should be filtered, adding it to the appropriate count if so, and
    /// returns a bool determining if it was filtered
    fn filter_alignment(
        &mut self,
        f_alignment: &Alignment,
        r_alignment: Option<&Alignment>,
    ) -> FilterReason {
        match &self.config.alignment_tolerance {
            None => {}
            Some(t) => {
                if f_alignment.score < t.minimum_f_score {
                    self.bad_alignment += 1;
                    return FilterReason::BadAlignment;
                }

                if let Some(r_al) = r_alignment {
                    if r_al.score < t.minimum_r_score {
                        self.bad_alignment += 1;
                        return FilterReason::BadAlignment;
                    }
                }
            }
        }

        FilterReason::None
    }

    /// Total count of filtered reads
    fn total(&self) -> u64 {
        self.low_mean_quality + self.bad_alignment
    }

    /// Generate a string of TSV lines representing the filter counts
    ///
    /// Creates a TSV string with columns for filter reason, count and
    /// proportion of the supplied total read count. Mainly for use
    /// wihen outputing from ReadSummary
    fn to_tsv(&self, total: u64) -> String {
        let filtered_total = self.total();

        let mut out = String::with_capacity(300);

        out.push_str(&format!(
            "filtered\ttotal\t{}\t{:.4}\t1.000\n",
            filtered_total,
            filtered_total as f32 / total as f32,
        ));

        out.push_str(&format!(
            "filtered\tlow_mean_quality\t{}\t{:.4}\t{:.4}\n",
            self.low_mean_quality,
            self.low_mean_quality as f32 / total as f32,
            self.low_mean_quality as f32 / filtered_total as f32,
        ));

        out.push_str(&format!(
            "filtered\tbad_alignment\t{}\t{:.4}\t{:.4}\n",
            self.bad_alignment,
            self.bad_alignment as f32 / total as f32,
            self.bad_alignment as f32 / filtered_total as f32,
        ));

        out
    }
}

/// Configuration for filtering
///
/// Instructions for how to filter reads
#[derive(Debug, Clone)]
pub struct FilterConfig {
    mean_quality_threshold: Option<f32>,
    alignment_tolerance: Option<AlignmentTolerance>,
}

impl FilterConfig {
    pub fn new(
        mean_quality_threshold: Option<f32>,
        alignment_tolerance: Option<AlignmentTolerance>,
    ) -> Self {
        Self {
            mean_quality_threshold,
            alignment_tolerance,
        }
    }
}

#[derive(Debug, Clone)]
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

/// Combination of ObservedRegions seen in sequence reads
///
/// A set of observed regions determining the "type" of read, as defined in the LibSpec.
/// Also includes a count, the read grouping (for instance for different cells in single
/// cell studies) and whether it matches an expected library member.
#[derive(Debug)]
pub struct ObservedCombination {
    /// Count of observations for each read group. Ungrouped reads are stored in None
    counts: HashMap<ReadGroup, u32>,

    /// ObservedRegions defining the sequence form. References to ObservedRegion which
    /// should be stored in the parent ObservedCombinations object.
    regions: HashMap<String, Rc<RefCell<ObservedRegion>>>,

    /// Status and result of comparison with the expected library of sequences
    library_matches: CombinationMatch,
}

/// Status of the match between ObservedCombination and a Library
///
/// Includes the match status and the indeces of matches in the libray
/// plus the distance to the library.
#[derive(Debug, Clone)]
pub enum CombinationMatch {
    /// Comparison hasn't occured
    Uncompared,

    /// Full match with a specific library member and the distance
    Match { ind: usize, distance: u64 },

    /// Fully matches multiple library members and the distance
    MultiMatch { inds: HashSet<usize>, distance: u64 },

    /// Partially matches multiple library members and the total distance
    Recombination { distance: u64 },

    /// Regions exist but at least one cannot be assigned to the library
    Mismatch,

    /// Not all regions exist
    Nonmatch,
}

impl CombinationMatch {
    /// Output a TSV chunk for the combination match status
    ///
    /// Has the \t separated format:
    /// status combination_distance combinations_in_library combination_indexes
    fn to_tsv_chunk(&self) -> String {
        match self {
            CombinationMatch::Uncompared => "uncompared\t\t\t\t".to_string(),
            CombinationMatch::Match { ind, distance } => format!("match\t{distance}\t1\t{ind}\t",),
            CombinationMatch::MultiMatch { inds, distance } => format!(
                "match\t{}\t{}\t{}\t",
                distance,
                inds.len(),
                inds.iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            ),
            CombinationMatch::Recombination { distance } => format!("recombination\t{distance}\t0\t\t",),
            CombinationMatch::Mismatch => "mismatch\t\t0\t\t".to_string(),
            CombinationMatch::Nonmatch => "nonmatch\t\t0\t\t".to_string(),
        }
    }

    /// Output a summary TSV chunk for the combination match status
    ///
    /// Has the \t separated format:
    /// status combinations_in_library combination_indexes
    fn to_summary_tsv_chunk(&self) -> String {
        match self {
            CombinationMatch::Uncompared => "uncompared\t\t\t".to_string(),
            CombinationMatch::Match { ind, .. } => format!("match\t1\t{ind}\t",),
            CombinationMatch::MultiMatch { inds, .. } => format!(
                "match\t{}\t{}\t",
                inds.len(),
                inds.iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            ),
            CombinationMatch::Recombination { .. } => "recombination\t0\t\t".to_string(),
            CombinationMatch::Mismatch => "mismatch\t0\t\t".to_string(),
            CombinationMatch::Nonmatch => "nonmatch\t0\t\t".to_string(),
        }
    }
}

impl ObservedCombination {
    fn new(regions: HashMap<String, Rc<RefCell<ObservedRegion>>>) -> Self {
        Self {
            counts: HashMap::new(),
            regions,
            library_matches: CombinationMatch::Uncompared,
        }
    }

    /// Generate a key uniquely identifying a combination
    fn key(regions: Vec<RegionKey>) -> CombinationKey {
        regions
    }

    /// Total count across all read groups
    fn total_count(&self) -> u32 {
        self.counts.values().sum()
    }

    /// Incremment count for a read group
    fn increment_count(&mut self, group: ReadGroup) {
        match self.counts.get_mut(&group) {
            Some(x) => *x += 1,
            None => {
                self.counts.insert(group, 1);
            }
        }
    }

    /// Compare the combination to expected combinations in the library
    ///
    /// Returns a CombinationMatch object which can also be added to the ObservedCombination
    /// library matches field. Looks at each region in turn and identifies which library
    /// combinations are possible overall matches.
    fn compare_to_library(
        &self,
        region_ids: &Vec<String>,
        library: &Library,
        distance_metric: DistanceMetric,
        max_matches: usize,
    ) -> CombinationMatch {
        let mut comb_dist: u64 = 0;
        let mut candidate_matches: Option<HashSet<usize>> = None;

        for reg_id in region_ids {
            let reg = match self.regions.get(reg_id) {
                None => return CombinationMatch::Nonmatch,
                Some(x) => {
                    // Should never need this with the implementation in ObservedCombinations, but here as a back-up as otherwise could panic later. Do all regions first as slightly more efficient and easier to follow in log
                    if !x.borrow().is_compared_to_library() {
                        let val =
                            x.borrow()
                                .compare_to_library(library, distance_metric, max_matches);
                        x.borrow_mut().nearest_matches = val;
                    }
                    x
                }
            };

            match &reg.borrow().nearest_matches {
                RegionMatch::Uncompared => panic!("Region uncompared despite just comparing"),
                RegionMatch::Unmatched | RegionMatch::Overmatched { .. } => {
                    // An unmatched region or an indeterminate one means the
                    // combination cannot be assigned
                    return CombinationMatch::Mismatch;
                }
                RegionMatch::NoLibrary { .. } => {}
                RegionMatch::Match {
                    seq_match,
                    distance,
                } => {
                    comb_dist += distance;
                    match candidate_matches {
                        None => candidate_matches = Some(seq_match.inds.clone()),
                        Some(ref mut x) => {
                            x.retain(|i| seq_match.inds.contains(i));
                        }
                    };
                }
                RegionMatch::MultiMatch {
                    seq_matches,
                    distance,
                } => {
                    comb_dist += distance;
                    for mat in seq_matches {
                        match candidate_matches {
                            None => candidate_matches = Some(mat.inds.clone()),
                            Some(ref mut x) => {
                                x.retain(|i| mat.inds.contains(i));
                            }
                        };
                    }
                }
            }
        }

        match candidate_matches {
            None => CombinationMatch::Nonmatch, // Only occurs if no regions (e.g. empty reads)
            Some(x) => {
                if x.len() == 1 {
                    CombinationMatch::Match {
                        // Can unwrap because we know x.len() == 1
                        ind: *x.iter().next().unwrap(),
                        distance: comb_dist,
                    }
                } else if x.is_empty() {
                    CombinationMatch::Recombination {
                        distance: comb_dist,
                    }
                } else {
                    CombinationMatch::MultiMatch {
                        inds: x,
                        distance: comb_dist,
                    }
                }
            }
        }
    }

    /// Generate tsv line(s) corresponding to this combination. Each read group
    /// the combination is observed is given a separate line
    fn to_tsv(&self, region_ids: &Vec<String>) -> String {
        // Line has \t separated format:
        // group [{region} {region}_nearest {region}_distance {region}_n_matches for each region] status combination_distance combinations_in_library combination_indexes count

        let mut output = String::with_capacity(100 * self.counts.len());

        for (group, count) in self.counts.iter() {
            // Read group
            match group {
                ReadGroup::Ungrouped => output.push('\t'),
                ReadGroup::Unmatched => output.push_str("_unmatched_\t"),
                ReadGroup::Match(x) => {
                    output.push_str(x);
                    output.push('\t');
                }
            };

            // Region seq/nearest match(s)/distance per region
            for reg_id in region_ids {
                let region = self.regions.get(reg_id);

                match region {
                    None => output.push_str("\t\t\t\t"), // Missing regions 4 blanks
                    Some(r) => {
                        output.push_str(&r.borrow().to_tsv_chunk());
                        output.push('\t');
                    }
                }
            }

            output.push_str(&self.library_matches.to_tsv_chunk());
            output.push_str(&count.to_string());
            output.push('\n');
        }

        output
    }
}

/// Observed sequence for a given library region
///
/// Each time this sequence is observed for this region it will be linked back to this object to
/// save library comparison overhead. The struct includes the observed sequence, a marker for
/// whether the whole region has necessarily been observed (i.e. if it is at the the start or end
///  of a read) and a RegionMatch object detailing the library match
#[derive(Debug)]
pub struct ObservedRegion {
    /// Name of the region from LibSpec
    id: String,

    /// Observed Sequence
    seq: Sequence,

    /// Completeness status of the region
    complete: RegionCompleteness,

    /// Library indeces of nearest matches (Vec of Vecs containing indeces from each matching seq)
    nearest_matches: RegionMatch,
}

/// Status of the match between an ObservedRegion and a Library
///
/// Includes the status plus reference(s) to the matched sequence in
/// the library and how distant it is.
#[derive(Debug, Eq, Hash, PartialEq, Clone)]
pub enum RegionMatch {
    /// No comparison has occured yet
    Uncompared,

    /// A single match (`Vec<u8>` sequence and library indeces) and associated distance
    Match {
        seq_match: Rc<LibraryRegion>,
        distance: u64,
    },

    /// Multiple equidistant matches and the distance
    MultiMatch {
        seq_matches: Vec<Rc<LibraryRegion>>,
        distance: u64,
    },

    /// Too many matches
    Overmatched { distance: u64, matches: usize },

    /// No match found
    Unmatched,

    /// Region not in library, with option to store the observed sequence
    NoLibrary { seq: Option<Sequence> },
}

impl RegionMatch {
    /// Get library Sequence, distance and number of matches as strings for output.
    ///
    /// NoLibrary matches are considered not to have a matching sequence
    fn to_str_fields(&self) -> (String, String, String) {
        match self {
            RegionMatch::Uncompared | RegionMatch::Unmatched | RegionMatch::NoLibrary { .. } => {
                ("".to_string(), "".to_string(), "0".to_string())
            }
            RegionMatch::Overmatched { distance, matches } => {
                ("".to_string(), distance.to_string(), matches.to_string())
            }
            RegionMatch::Match {
                seq_match,
                distance,
            } => (
                seq_to_string_or_log(&seq_match.sequence),
                distance.to_string(),
                "1".to_string(),
            ),
            RegionMatch::MultiMatch {
                seq_matches,
                distance,
            } => (
                seq_matches
                    .iter()
                    .map(|x| seq_to_string_or_log(&x.sequence))
                    .collect::<Vec<_>>()
                    .join(","),
                distance.to_string(),
                seq_matches.len().to_string(),
            ),
        }
    }

    /// Get matching Sequence as a string, including the passed through raw
    /// Sequence for NoLibrary matches.
    fn str_sequence(&self) -> String {
        match self {
            RegionMatch::Uncompared | RegionMatch::Unmatched | RegionMatch::Overmatched { .. } => {
                "".to_string()
            }
            RegionMatch::NoLibrary { seq } => match seq {
                Some(s) => seq_to_string_or_log(s),
                None => "".to_string(),
            },
            RegionMatch::Match { seq_match, .. } => seq_to_string_or_log(&seq_match.sequence),
            RegionMatch::MultiMatch { seq_matches, .. } => seq_matches
                .iter()
                .map(|x| seq_to_string_or_log(&x.sequence))
                .collect::<Vec<_>>()
                .join(","),
        }
    }
}

/// Whether an observed region is complete
///
/// Complete regions are known to be full, either being flanked by neighbours on
/// both side in the read or being the full expected length. Partial regions might
/// have more sequence before/after because this is not the case.
#[derive(Debug, Eq, Hash, PartialEq, Clone, Copy)]
pub enum RegionCompleteness {
    /// A full region
    Complete,

    /// 5 prime end is potentially incomplete
    Partial5Prime,

    /// 3 prime end is potentially incomplete
    Partial3Prime,

    /// Reads come from both ends without meeting. Stored as
    /// both seqs divided by a / at split_ind
    MissingCenter { split_ind: usize },

    /// Reads come from both ends and overlap. Stored as
    /// both seqs divided by a / at split_ind
    Overlapping { split_ind: usize },
}

impl ObservedRegion {
    /// Create a new ObservedRegion
    fn new(id: String, seq: &[u8], complete: RegionCompleteness) -> Self {
        Self {
            id,
            seq: seq.to_vec(),
            complete,
            nearest_matches: RegionMatch::Uncompared,
        }
    }

    /// Generate a key for a region for mapping it in an ObservedRegions HashMap
    fn key(name: String, seq: &[u8], complete: RegionCompleteness) -> RegionKey {
        (name, seq.to_vec(), complete)
    }

    /// Check if library comparison has been performed
    fn is_compared_to_library(&self) -> bool {
        !matches!(self.nearest_matches, RegionMatch::Uncompared)
    }

    /// Compare the region sequence to a library of oligo options, dispatching to the appropriate implementation
    ///
    /// The distance to each is calculated and all nearest match indeces below the threshold
    /// supplied to library are returned.
    /// If no applicable matches are found nearest_matches will be an empty.
    /// If over max_matches hits are returned the region will be considered overmatched and
    /// indeterminant.
    fn compare_to_library(
        &self,
        library: &Library,
        distance_metric: DistanceMetric,
        max_matches: usize,
    ) -> RegionMatch {
        let lib_match = match self.complete {
            RegionCompleteness::Complete => {
                match library.lookup(&self.id, &self.seq, distance_metric, PartialMatching::Full) {
                    // lookup can only return Err(LibraryError::MissingRegion) which implies the
                    // region isn't in the library. Needs changing if more errors are added to it.
                    Err(_) => return RegionMatch::NoLibrary { seq: None },
                    Ok(x) => x,
                }
            }
            RegionCompleteness::Partial5Prime => {
                match library.lookup(
                    &self.id,
                    &self.seq,
                    distance_metric,
                    PartialMatching::ThreePrimeOnly,
                ) {
                    // lookup can only return Err(LibraryError::MissingRegion) which implies the
                    // region isn't in the library. Needs changing if more errors are added to it.
                    Err(_) => return RegionMatch::NoLibrary { seq: None },
                    Ok(x) => x,
                }
            }
            RegionCompleteness::Partial3Prime => {
                match library.lookup(
                    &self.id,
                    &self.seq,
                    distance_metric,
                    PartialMatching::FivePrimeOnly,
                ) {
                    // lookup can only return Err(LibraryError::MissingRegion) which implies the
                    // region isn't in the library. Needs changing if more errors are added to it.
                    Err(_) => return RegionMatch::NoLibrary { seq: None },
                    Ok(x) => x,
                }
            }
            RegionCompleteness::MissingCenter { split_ind }
            | RegionCompleteness::Overlapping { split_ind } => {
                let left_match = library.lookup(
                    &self.id,
                    &self.seq[0..split_ind],
                    distance_metric,
                    PartialMatching::FivePrimeOnly,
                );

                let right_match = library.lookup(
                    &self.id,
                    &self.seq[split_ind + 1..self.seq.len()],
                    distance_metric,
                    PartialMatching::ThreePrimeOnly,
                );

                match (left_match, right_match) {
                    (Ok(left), Ok(right)) => merge_matches(left, right),

                    // lookup can only return Err(LibraryError::MissingRegion) which implies the
                    // region isn't in the library. Needs changing if more errors are added to it.
                    (_, Err(_)) | (Err(_), _) => return RegionMatch::NoLibrary { seq: None },
                }
            }
        };

        match lib_match {
            None => RegionMatch::Unmatched,
            Some(x) => {
                if x.matches.len() == 1 {
                    RegionMatch::Match {
                        seq_match: x.matches[0].clone(),
                        distance: x.distance,
                    }
                } else if x.matches.len() > max_matches {
                    RegionMatch::Overmatched {
                        distance: x.distance,
                        matches: x.matches.len(),
                    }
                } else {
                    RegionMatch::MultiMatch {
                        seq_matches: x.matches,
                        distance: x.distance,
                    }
                }
            }
        }
    }

    /// Generate output TSV chunk describing the region
    fn to_tsv_chunk(&self) -> String {
        // Generate "{region} {region}_nearest {region}_distance {region}_n_matches" String
        let mut seq = match String::from_utf8(self.seq.clone()) {
            Ok(x) => x,
            Err(_) => panic!("Error converting Vec<u8> sequence to string"),
        };

        // Add appropriate marker to incomplete regions
        match self.complete {
            RegionCompleteness::Complete
            | RegionCompleteness::MissingCenter { .. }
            | RegionCompleteness::Overlapping { .. } => {}
            RegionCompleteness::Partial5Prime => seq.insert(0, '^'),
            RegionCompleteness::Partial3Prime => seq.push('^'),
        };

        let (nearest, distance, matches) = self.nearest_matches.to_str_fields();

        format!("{}\t{}\t{}\t{}", seq, nearest, distance, matches)
    }
}

/// Summary version of ObservedCombination
///
/// This counts a particular form of match with the expected library instead of a particular
/// combination of observed reads and so uses a hash map of RegionMatches instead of ObservedRegions
#[derive(Debug)]
pub struct LibraryCombination {
    /// Count of observations for each read group. Ungrouped reads are stored in None
    counts: HashMap<ReadGroup, u32>,

    /// RegionMatches determine the connection to the library
    regions: HashMap<String, RegionMatch>,

    /// Status and result of comparison with the expected library of sequences
    library_matches: CombinationMatch,
}

impl LibraryCombination {
    fn new(regions: HashMap<String, RegionMatch>, library_matches: CombinationMatch) -> Self {
        Self {
            counts: HashMap::new(),
            regions,
            library_matches,
        }
    }

    /// Get the total count across all read groups
    fn total_count(&self) -> u32 {
        self.counts.values().sum()
    }

    /// Increment the count for the desired read group by n
    fn increment_count(&mut self, group: &ReadGroup, n: u32) {
        match self.counts.get_mut(group) {
            Some(x) => *x += n,
            None => {
                self.counts.insert(group.clone(), n);
            }
        }
    }

    /// Generate tsv line(s) corresponding to this combination. Each read group
    /// the combination is observed is given a separate line
    fn to_tsv(&self, region_ids: &Vec<String>) -> String {
        // Line has \t separated format:
        // group [{region} for each region] status combinations_in_library combination_indexes count

        let mut output = String::with_capacity(100 * self.counts.len());

        for (group, count) in self.counts.iter() {
            // Read group
            match group {
                ReadGroup::Ungrouped => output.push('\t'),
                ReadGroup::Unmatched => output.push_str("_unmatched_\t"),
                ReadGroup::Match(x) => {
                    output.push_str(x);
                    output.push('\t');
                }
            };

            // Region seq
            for reg_id in region_ids {
                let region = self.regions.get(reg_id);

                match region {
                    None => output.push('\t'), // Missing regions 1 blanks
                    Some(r) => {
                        output.push_str(&r.str_sequence());
                        output.push('\t');
                    }
                }
            }

            output.push_str(&self.library_matches.to_summary_tsv_chunk());
            output.push_str(&count.to_string());
            output.push('\n');
        }

        output
    }
}

/// Summary counts of read types
pub struct ReadSummary {
    /// Comparison hasn't occured
    pub uncompared: u64,

    /// Full match with a specific library member
    pub exact_match: u64,

    /// Nearest match with a specific library member
    pub nearest_match: u64,

    /// Fully matches multiple library members and the distance
    pub multimatch: u64,

    /// Exactly matches multiple library members but recombined
    pub exact_recombination: u64,

    /// Partially matches multiple library members but recombinaed
    pub nearest_recombination: u64,

    /// Regions exist but at least one cannot be assigned to the library
    pub mismatch: u64,

    /// Not all regions exist
    pub nonmatch: u64,

    /// Filtered Reads
    pub filtered_reads: Option<FilteredReads>,
}

impl ReadSummary {
    #[allow(dead_code)]
    fn new(
        uncompared: u64,
        exact_match: u64,
        nearest_match: u64,
        multimatch: u64,
        exact_recombination: u64,
        nearest_recombination: u64,
        mismatch: u64,
        nonmatch: u64,
        filtered_reads: Option<FilteredReads>,
    ) -> Self {
        Self {
            uncompared,
            exact_match,
            nearest_match,
            multimatch,
            exact_recombination,
            nearest_recombination,
            mismatch,
            nonmatch,
            filtered_reads,
        }
    }

    /// Inititalise an empty ReadSummary
    ///
    /// Useful shortcut for using it as a counter
    fn empty() -> Self {
        Self {
            uncompared: 0,
            exact_match: 0,
            nearest_match: 0,
            multimatch: 0,
            exact_recombination: 0,
            nearest_recombination: 0,
            mismatch: 0,
            nonmatch: 0,
            filtered_reads: None,
        }
    }

    /// Sum of unfiltered reads
    fn total_unfiltered(&self) -> u64 {
        self.uncompared
            + self.exact_match
            + self.nearest_match
            + self.multimatch
            + self.exact_recombination
            + self.nearest_recombination
            + self.mismatch
            + self.nonmatch
    }

    /// Total reads observed across categories
    pub fn total(&self) -> u64 {
        self.total_unfiltered()
            + match &self.filtered_reads {
                Some(f) => f.total(),
                None => 0,
            }
    }

    /// Write to file in TSV format
    pub fn write_tsv(self, file: File) -> Result<(), anyhow::Error> {
        let mut writer = BufWriter::new(file);
        let unfiltered_total = self.total_unfiltered();
        let total = self.total();

        // Write data
        writeln!(
            writer,
            "group\tmetric\tcount\toverall_proportion\tgroup_proportion"
        )?;
        writeln!(writer, "all\ttotal\t{}\t1.0000\t1.000", total)?;
        writeln!(
            writer,
            "unfiltered\ttotal\t{}\t{:.4}\t1.000",
            unfiltered_total,
            unfiltered_total as f32 / total as f32,
        )?;
        writeln!(
            writer,
            "unfiltered\tuncompared\t{}\t{:.4}\t{:.4}",
            self.uncompared,
            self.uncompared as f32 / total as f32,
            self.uncompared as f32 / unfiltered_total as f32,
        )?;
        writeln!(
            writer,
            "unfiltered\texact_match\t{}\t{:.4}\t{:.4}",
            self.exact_match,
            self.exact_match as f32 / total as f32,
            self.exact_match as f32 / unfiltered_total as f32,
        )?;
        writeln!(
            writer,
            "unfiltered\tnearest_match\t{}\t{:.4}\t{:.4}",
            self.nearest_match,
            self.nearest_match as f32 / total as f32,
            self.nearest_match as f32 / unfiltered_total as f32,
        )?;
        writeln!(
            writer,
            "unfiltered\tmultimatch\t{}\t{:.4}\t{:.4}",
            self.multimatch,
            self.multimatch as f32 / total as f32,
            self.multimatch as f32 / unfiltered_total as f32,
        )?;
        writeln!(
            writer,
            "unfiltered\texact_recombination\t{}\t{:.4}\t{:.4}",
            self.exact_recombination,
            self.exact_recombination as f32 / total as f32,
            self.exact_recombination as f32 / unfiltered_total as f32,
        )?;
        writeln!(
            writer,
            "unfiltered\tnearest_recombination\t{}\t{:.4}\t{:.4}",
            self.nearest_recombination,
            self.nearest_recombination as f32 / total as f32,
            self.nearest_recombination as f32 / unfiltered_total as f32,
        )?;
        writeln!(
            writer,
            "unfiltered\tmismatch\t{}\t{:.4}\t{:.4}",
            self.mismatch,
            self.mismatch as f32 / total as f32,
            self.mismatch as f32 / unfiltered_total as f32,
        )?;
        writeln!(
            writer,
            "unfiltered\tnonmatch\t{}\t{:.4}\t{:.4}",
            self.nonmatch,
            self.nonmatch as f32 / total as f32,
            self.nonmatch as f32 / unfiltered_total as f32,
        )?;

        match self.filtered_reads {
            Some(f) => {
                write!(writer, "{}", f.to_tsv(total))?;
            }
            None => {
                writeln!(writer, "filtered\ttotal\t0\t0.0000\t0.0000",)?;
            }
        }

        writer.flush()?;
        Ok(())
    }
}

/// Count algorithm to apply
#[derive(Clone, ValueEnum, Debug, Copy)]
pub enum CountMode {
    FullRead,
    Inframe,
    Pattern,
    Align,
}

/// Count the occurance of query regions in sequencing reads
///
/// Dispatches counting to the appropriate implementation based on CountMode
pub fn count_reads(
    reads: ReadPairParser,
    lib_spec: &Option<LibrarySpec>,
    mode: CountMode,
    filter_config: FilterConfig,
    alignment_scorer: Option<AlignmentScorer>,
    cache: bool,
    progress_style: Option<&ProgressStyle>,
) -> Result<ObservedCombinations, anyhow::Error> {
    let default_progress = ProgressStyle::new(None);
    let progress = progress_style.unwrap_or(&default_progress);

    // Determine type of matching desired and despatch as appropriate
    match (reads.has_reverse(), lib_spec, mode, alignment_scorer) {
        (_, _, CountMode::Align, None) => Err(ReadCountError::Error {
            desc: "Mode is 'align' but no AlignmentScorer passed".to_string(),
        }
        .into()),
        (false, Some(lib_spec), CountMode::Align, Some(a)) => {
            count_single_align(reads, lib_spec, filter_config, a, cache, progress)
        }
        (true, Some(lib_spec), CountMode::Align, Some(a)) => {
            count_paired_align(reads, lib_spec, filter_config, a, cache, progress)
        }
        (false, Some(lib_spec), CountMode::Pattern, _) => {
            count_single_pattern(reads, lib_spec, filter_config, progress)
        }
        (true, Some(lib_spec), CountMode::Pattern, _) => {
            count_paired_pattern(reads, lib_spec, filter_config, progress)
        }
        (false, Some(lib_spec), CountMode::Inframe, _) => {
            count_single_inframe(reads, lib_spec, filter_config, progress)
        }
        (true, Some(lib_spec), CountMode::Inframe, _) => {
            count_paired_inframe(reads, lib_spec, filter_config, progress)
        }
        (false, Some(_), CountMode::FullRead, _) => {
            count_single_raw(reads, filter_config, progress)
        }
        (true, Some(_), CountMode::FullRead, _) => count_paired_raw(reads, filter_config, progress),
        (false, None, _, _) => count_single_raw(reads, filter_config, progress),
        (true, None, _, _) => count_paired_raw(reads, filter_config, progress),
    }
}

/// Count single end reads by aligning to the library template
fn count_single_align(
    reads: ReadPairParser,
    lib_spec: &LibrarySpec,
    filter_config: FilterConfig,
    alignment_scorer: AlignmentScorer,
    cache: bool,
    progress_style: &ProgressStyle,
) -> Result<ObservedCombinations, anyhow::Error> {
    info!("{}", ALIGN_START_SINGLE_MSG);

    let mut progress: Progress = Progress::from_style(
        progress_style,
        PROG_MSG,
        FINAL_MSG,
        None,
        ALIGN_LOG_INTERVAL,
    );

    // Identify regions
    let regions = lib_spec.variable_regions();
    let region_positions: Vec<(usize, usize)> = regions
        .iter()
        .map(|x| lib_spec.template_position(x).expect("Region from LibSpec"))
        // Alignment path is 1 based, so offset originally 0 based positions
        .map(|x| (x.0 + 1, x.1 + 1))
        .collect();

    // Initialise counts
    let mut counts = ObservedCombinations::new(regions.clone(), filter_config);
    let mut observed_reads: ObservedReads = HashMap::new();

    // Initialise aligner
    let scoring = alignment_scorer.get_scoring();
    let mut aligner = Aligner::with_capacity_and_scoring(400, 150, scoring);
    let template = lib_spec.template_sequence();

    // Iterate over reads
    info!("Printing example alignments. Check these are reasonable to help diagnose problems");
    for (i, result) in reads.enumerate() {
        let record: ReadPair = result?;

        // Check if read should be filtered
        if !matches!(counts.filter_readpair(&record), FilterReason::None) {
            continue;
        }

        let record_key = record.key();

        // Check cached reads and only realign if needed
        if cache && observed_reads.contains_key(&record_key) {
            let cache_hit = observed_reads
                .get(&record_key)
                .expect("Just checked read key is present in observed reads");
            match cache_hit {
                CacheHit::Comb(k) => {
                    counts.add_or_increment_combination(k, record.group)?;
                }
                CacheHit::Filter(r) => {
                    counts.update_filter_count(r);
                }
            }
        } else {
            let read_alignment = aligner.semiglobal(record.forward.seq(), &template);
            let alignment_path = read_alignment.path();

            // Filter reads with too low scores or otherwise unmatched
            match counts.filter_alignment(&read_alignment, None) {
                FilterReason::None => {}
                other => {
                    if cache {
                        observed_reads.insert(record_key, CacheHit::Filter(other));
                    }
                    continue;
                }
            }

            // Extract positions
            let query_regions = regions_from_alignment_path(&region_positions, &alignment_path)?;

            if i < 10 {
                info!(
                    "Alignment {}:\nScore: {}, Cigar: {}\n{}\nExtracted regions: {:?}",
                    i + 1,
                    read_alignment.score,
                    read_alignment.cigar(false),
                    read_alignment.pretty(record.forward.seq(), &template, 100),
                    query_regions
                );
                debug!("{:?}", read_alignment.path())
            }

            let mut comb_key_vec: Vec<RegionKey> = Vec::with_capacity(query_regions.len());
            for (id, opt_pos) in zip(&regions, &query_regions) {
                if let Some(pos) = opt_pos {
                    match record.forward.seq().get((pos.0 - 1)..(pos.1 - 1)) {
                        Some(s) => {
                            comb_key_vec.push(ObservedRegion::key(
                                id.to_string(),
                                // Offset seq lookup - rust vec 0 based and AlignmentPath 1 based
                                s,
                                pos.2,
                            ));
                        }
                        None => {
                            return Err(ReadCountError::BadAlignment {
                                alignment: Box::new(AlignmentInfo {
                                    read_id: record.forward.id().to_string(),
                                    read_number: i,
                                    pretty_alignment: read_alignment.pretty(
                                        record.forward.seq(),
                                        &template,
                                        100,
                                    ),
                                    alignment: read_alignment,
                                    region_ids: regions,
                                    region_positions,
                                    mapped_positions: query_regions,
                                }),
                            }
                            .into());
                        }
                    }
                }
            }

            // Construct combination
            let comb_key: CombinationKey = ObservedCombination::key(comb_key_vec);
            //eprintln!("{:?}", comb_key);

            counts.add_or_increment_combination(&comb_key, record.group)?;

            if cache {
                observed_reads.insert(record_key, CacheHit::Comb(comb_key));
            }
        }
        progress.inc(1);
    }
    progress.finish();

    Ok(counts)
}

/// Count paired end reads by aligning to the library template
fn count_paired_align(
    reads: ReadPairParser,
    lib_spec: &LibrarySpec,
    filter_config: FilterConfig,
    alignment_scorer: AlignmentScorer,
    cache: bool,
    progress_style: &ProgressStyle,
) -> Result<ObservedCombinations, anyhow::Error> {
    info!("{}", ALIGN_START_PAIRED_MSG);

    let mut progress: Progress = Progress::from_style(
        progress_style,
        PROG_MSG,
        FINAL_MSG,
        None,
        ALIGN_LOG_INTERVAL,
    );

    // Identify regions
    let regions = lib_spec.variable_regions();

    let region_lengths: Vec<usize> = regions
        .iter()
        .map(|x| {
            lib_spec
                .get_region(x)
                .expect("Region should exist as comes from lib")
                .len()
        })
        .collect();

    let region_positions: Vec<(usize, usize)> = regions
        .iter()
        .map(|x| lib_spec.template_position(x).expect("Region from LibSpec"))
        // Alignment path is 1 based, so offset originally 0 based positions
        .map(|x| (x.0 + 1, x.1 + 1))
        .collect();

    // Initialise counts
    let mut counts = ObservedCombinations::new(regions.clone(), filter_config);
    let mut observed_reads: ObservedReads = HashMap::new();

    // Initialise aligner
    let scoring = alignment_scorer.get_scoring();
    let mut aligner = Aligner::with_capacity_and_scoring(400, 150, scoring);
    let template = lib_spec.template_sequence();

    // Iterate over reads
    info!("Printing example alignments. Check these are reasonable to help diagnose problems");
    for (i, result) in reads.enumerate() {
        let record: ReadPair = result?;

        // Check if read should be filtered
        if !matches!(counts.filter_readpair(&record), FilterReason::None) {
            continue;
        }

        let record_key = record.key();

        // Check cached reads and only realign if needed
        if cache && observed_reads.contains_key(&record_key) {
            let cache_hit = observed_reads
                .get(&record_key)
                .expect("Just checked read key is present in observed reads");
            match cache_hit {
                CacheHit::Comb(k) => {
                    counts.add_or_increment_combination(k, record.group)?;
                }
                CacheHit::Filter(r) => {
                    counts.update_filter_count(r);
                }
            }
        } else {
            let f_read = record.forward.seq();
            let r_read = match &record.reverse {
                Some(x) => revcomp(x.seq()),
                None => {
                    return Err(ReadCountError::Error {
                        desc: format!("No reverse read found at read {i}"),
                    }
                    .into());
                }
            };

            let f_qual = record.forward.qual();
            let r_qual: Vec<u8> = match &record.reverse {
                Some(x) => x.qual().iter().rev().cloned().collect(),
                None => {
                    return Err(ReadCountError::Error {
                        desc: format!("No reverse read found at read {i}"),
                    }
                    .into());
                }
            };

            // Fwd Read
            let f_alignment = aligner.semiglobal(f_read, &template);
            let f_path = f_alignment.path();

            // Rev Read
            let r_alignment = aligner.semiglobal(&r_read, &template);
            let r_path = r_alignment.path();

            // Filter reads with too low scores or otherwise unmatched
            match counts.filter_alignment(&f_alignment, Some(&r_alignment)) {
                FilterReason::None => {}
                other => {
                    if cache {
                        observed_reads.insert(record_key, CacheHit::Filter(other));
                    }
                    continue;
                }
            }

            // Extract regions
            let f_regions = regions_from_alignment_path(&region_positions, &f_path)?;
            let r_regions = regions_from_alignment_path(&region_positions, &r_path)?;

            // Print first 10 alignments
            if i < 10 {
                info!(
                    "Fwd Alignment {}:\nScore: {}, Cigar: {}\n{}\nExtracted regions: {:?}",
                    i + 1,
                    f_alignment.score,
                    f_alignment.cigar(false),
                    f_alignment.pretty(f_read, &template, 100),
                    f_regions
                );
                debug!("{:?}", f_alignment.path());

                info!(
                    "Rev Alignment {}:\nScore: {}, Cigar: {}\n{}\nExtracted regions: {:?}",
                    i + 1,
                    r_alignment.score,
                    r_alignment.cigar(false),
                    r_alignment.pretty(&r_read, &template, 100),
                    r_regions
                );
                debug!("{:?}", r_alignment.path());
            }

            let mut comb_key_vec: Vec<RegionKey> =
                Vec::with_capacity(std::cmp::max(f_regions.len(), r_regions.len()));

            for (id, len, f_pos, r_pos) in izip!(&regions, &region_lengths, &f_regions, &r_regions)
            {
                match (f_pos, r_pos) {
                    (None, None) => continue,
                    (None, Some(r)) => {
                        match r_read.get((r.0 - 1)..(r.1 - 1)) {
                            Some(s) => {
                                comb_key_vec.push(ObservedRegion::key(
                                    id.to_string(),
                                    // Offset seq lookup - rust vec 0 based and AlignmentPath 1 based
                                    s,
                                    r.2,
                                ));
                            }
                            None => {
                                return Err(ReadCountError::BadAlignment {
                                    alignment: Box::new(AlignmentInfo {
                                        read_id: record
                                            .reverse
                                            .expect("Read known to be present")
                                            .id()
                                            .to_string(),
                                        read_number: i,
                                        pretty_alignment: r_alignment
                                            .pretty(&r_read, &template, 100),
                                        alignment: r_alignment,
                                        region_ids: regions,
                                        region_positions,
                                        mapped_positions: r_regions,
                                    }),
                                }
                                .into());
                            }
                        }
                    }
                    (Some(f), None) => {
                        match f_read.get((f.0 - 1)..(f.1 - 1)) {
                            Some(s) => {
                                comb_key_vec.push(ObservedRegion::key(
                                    id.to_string(),
                                    // Offset seq lookup - rust vec 0 based and AlignmentPath 1 based
                                    s,
                                    f.2,
                                ));
                            }
                            None => {
                                return Err(ReadCountError::BadAlignment {
                                    alignment: Box::new(AlignmentInfo {
                                        read_id: record.forward.id().to_string(),
                                        read_number: i,
                                        pretty_alignment: f_alignment
                                            .pretty(f_read, &template, 100),
                                        alignment: f_alignment,
                                        region_ids: regions,
                                        region_positions,
                                        mapped_positions: f_regions,
                                    }),
                                }
                                .into());
                            }
                        }
                    }
                    (Some(f), Some(r)) => {
                        let f_reg_seq = match f_read.get((f.0 - 1)..(f.1 - 1)) {
                            Some(s) => s,
                            None => {
                                return Err(ReadCountError::BadAlignment {
                                    alignment: Box::new(AlignmentInfo {
                                        read_id: record.forward.id().to_string(),
                                        read_number: i,
                                        pretty_alignment: f_alignment
                                            .pretty(f_read, &template, 100),
                                        alignment: f_alignment,
                                        region_ids: regions,
                                        region_positions,
                                        mapped_positions: f_regions,
                                    }),
                                }
                                .into());
                            }
                        };

                        let f_reg_qual = f_qual
                            .get((f.0 - 1)..(f.1 - 1))
                            .expect("If seq extracted then qual should be too as the same length");

                        let r_reg_seq = match r_read.get((r.0 - 1)..(r.1 - 1)) {
                            Some(s) => s,
                            None => {
                                return Err(ReadCountError::BadAlignment {
                                    alignment: Box::new(AlignmentInfo {
                                        read_id: record
                                            .reverse
                                            .expect("Read known to be present")
                                            .id()
                                            .to_string(),
                                        read_number: i,
                                        pretty_alignment: r_alignment
                                            .pretty(&r_read, &template, 100),
                                        alignment: r_alignment,
                                        region_ids: regions,
                                        region_positions,
                                        mapped_positions: r_regions,
                                    }),
                                }
                                .into());
                            }
                        };

                        let r_reg_qual = r_qual
                            .get((r.0 - 1)..(r.1 - 1))
                            .expect("If seq extracted then qual should be too as the same length");

                        match merge_seqs(
                            Some((f_reg_seq.to_vec(), f_reg_qual.to_vec(), f.2)),
                            Some((r_reg_seq.to_vec(), r_reg_qual.to_vec(), r.2)),
                            *len,
                        )? {
                            Some((seq, comp)) => {
                                comb_key_vec.push(ObservedRegion::key(id.clone(), &seq, comp))
                            }
                            None => continue,
                        };
                    }
                }
            }

            // Construct combination
            let comb_key: CombinationKey = ObservedCombination::key(comb_key_vec);
            //eprintln!("{:?}", comb_key);

            counts.add_or_increment_combination(&comb_key, record.group)?;

            if cache {
                observed_reads.insert(record_key, CacheHit::Comb(comb_key));
            }
        }
        progress.inc(1);
    }
    progress.finish();

    Ok(counts)
}

/// Count single end reads using surrounding patterns from the library template
fn count_single_pattern(
    reads: ReadPairParser,
    lib_spec: &LibrarySpec,
    filter_config: FilterConfig,
    progress_style: &ProgressStyle,
) -> Result<ObservedCombinations, anyhow::Error> {
    info!("{}", PATTERN_START_SINGLE_MSG);

    let mut progress: Progress = Progress::from_style(
        progress_style,
        PROG_MSG,
        FINAL_MSG,
        None,
        PATTERN_LOG_INTERVAL,
    );

    let regions = lib_spec.variable_regions();

    let flank_regions: Vec<(Option<Sequence>, Option<Sequence>)> = match regions
        .iter()
        .map(|x| lib_spec.flanking_regions(x, 10))
        .collect::<Result<Vec<_>, _>>()
    {
        Ok(x) => x,
        Err(e) => return Err(e.into()),
    };

    info!(
        "Extracting regions using flank seqs: {:?}",
        flank_regions
            .iter()
            .map(|(x, y)| match (x, y) {
                (None, None) => ("None".to_string(), "None".to_string()),
                (None, Some(e)) => ("None".to_string(), seq_to_string_or_log(e)),
                (Some(s), None) => (seq_to_string_or_log(s), "None".to_string()),
                (Some(s), Some(e)) => (seq_to_string_or_log(s), seq_to_string_or_log(e)),
            })
            .collect::<Vec<(String, String)>>()
    );

    // Check no flanking regions are duplicates, in which case can't use pattern matching and
    // alignment is required
    if !flank_regions
        .iter()
        .filter_map(|x| x.0.clone())
        .all_unique()
    {
        return Err(ReadCountError::Error {
            desc: "Duplicate leading flanking sequences - alignment required".to_string(),
        }
        .into());
    }

    if !flank_regions
        .iter()
        .filter_map(|x| x.1.clone())
        .all_unique()
    {
        return Err(ReadCountError::Error {
            desc: "Duplicate trailing flanking sequences - alignment required".to_string(),
        }
        .into());
    }

    // Count reads
    let mut counts = ObservedCombinations::new(regions.clone(), filter_config);

    for result in reads {
        let record: ReadPair = result?;

        // Check if read should be filtered
        if !matches!(counts.filter_readpair(&record), FilterReason::None) {
            continue;
        }

        let mut comb_key_vec: Vec<RegionKey> = Vec::with_capacity(regions.len());
        let seq = record.forward.seq();
        let seq_len = seq.len();

        // Not most efficient - iterate once for each region, could improve
        for (id, flanks) in zip(&regions, &flank_regions) {
            let mut complete = RegionCompleteness::Complete;

            let start_pos = match &flanks.0 {
                None => {
                    complete = RegionCompleteness::Partial5Prime;
                    Some(0)
                }
                Some(x) => seq
                    .windows(x.len())
                    .position(|window| window == x)
                    .map(|i| i + x.len()),
            };

            let end_pos = match &flanks.1 {
                None => {
                    complete = RegionCompleteness::Partial3Prime;
                    Some(seq_len - 1)
                }
                Some(x) => seq.windows(x.len()).position(|window| window == x), // End pos doesn't need offsetting due to slice being exclusive
            };

            let reg_seq: Sequence = match (start_pos, end_pos) {
                // Region not found
                (None, _) | (_, None) => continue,

                // Region found
                (Some(start), Some(end)) => {
                    if start > end {
                        // Flank seqs overlap - no region contained within
                        continue;
                    }
                    seq.get(start..end)
                        .expect("Seq should contain region as extracted from within it")
                        .to_vec()
                }
            };

            comb_key_vec.push(ObservedRegion::key(id.clone(), &reg_seq, complete));
        }

        let comb_key: CombinationKey = ObservedCombination::key(comb_key_vec);

        counts.add_or_increment_combination(&comb_key, record.group)?;
        progress.inc(1);
    }
    progress.finish();

    Ok(counts)
}

/// Count paired end reads using surrounding patterns from the library template
fn count_paired_pattern(
    reads: ReadPairParser,
    lib_spec: &LibrarySpec,
    filter_config: FilterConfig,
    progress_style: &ProgressStyle,
) -> Result<ObservedCombinations, anyhow::Error> {
    info!("{}", PATTERN_START_PAIRED_MSG);

    let mut progress: Progress = Progress::from_style(
        progress_style,
        PROG_MSG,
        FINAL_MSG,
        None,
        PATTERN_LOG_INTERVAL,
    );

    let regions = lib_spec.variable_regions();

    let region_lengths: Vec<usize> = regions
        .iter()
        .map(|x| {
            lib_spec
                .get_region(x)
                .expect("Region should exist as comes from lib")
                .len()
        })
        .collect();

    let flank_regions: Vec<(Option<Sequence>, Option<Sequence>)> = match regions
        .iter()
        .map(|x| lib_spec.flanking_regions(x, 10))
        .collect::<Result<Vec<_>, _>>()
    {
        Ok(x) => x,
        Err(e) => return Err(e.into()),
    };

    info!(
        "Extracting regions using flank seqs: {:?}",
        flank_regions
            .iter()
            .map(|(x, y)| match (x, y) {
                (None, None) => ("None".to_string(), "None".to_string()),
                (None, Some(e)) => ("None".to_string(), seq_to_string_or_log(e)),
                (Some(s), None) => (seq_to_string_or_log(s), "None".to_string()),
                (Some(s), Some(e)) => (seq_to_string_or_log(s), seq_to_string_or_log(e)),
            })
            .collect::<Vec<(String, String)>>()
    );

    // Check no flanking regions are duplicates, in which case can't use pattern matching and
    // alignment is required
    if !flank_regions
        .iter()
        .filter_map(|x| x.0.clone())
        .all_unique()
    {
        return Err(ReadCountError::Error {
            desc: "Duplicate leading flanking sequences - alignment required".to_string(),
        }
        .into());
    }

    if !flank_regions
        .iter()
        .filter_map(|x| x.1.clone())
        .all_unique()
    {
        return Err(ReadCountError::Error {
            desc: "Duplicate trailing flanking sequences - alignment required".to_string(),
        }
        .into());
    }

    // Count reads
    let mut counts = ObservedCombinations::new(regions.clone(), filter_config);

    for (i, result) in reads.enumerate() {
        let record: ReadPair = result?;

        // Check if read should be filtered
        if !matches!(counts.filter_readpair(&record), FilterReason::None) {
            continue;
        }

        let mut comb_key_vec: Vec<RegionKey> = Vec::with_capacity(regions.len());

        let f_read = record.forward.seq();
        let r_read = match &record.reverse {
            Some(x) => revcomp(x.seq()),
            None => {
                return Err(ReadCountError::Error {
                    desc: format!("No reverse read found at read {i}"),
                }
                .into());
            }
        };

        let f_qual = record.forward.qual();
        let r_qual: Vec<u8> = match &record.reverse {
            Some(x) => x.qual().iter().rev().cloned().collect(),
            None => {
                return Err(ReadCountError::Error {
                    desc: format!("No reverse read found at read {i}"),
                }
                .into());
            }
        };

        let f_len = f_read.len();
        let r_len = r_read.len();

        // Not most efficient - iterate once for each region, could improve
        for (id, len, flanks) in izip!(&regions, &region_lengths, &flank_regions) {
            let mut complete = RegionCompleteness::Complete;

            // Match in fwd read
            let start_pos = match &flanks.0 {
                None => {
                    complete = RegionCompleteness::Partial5Prime;
                    Some(0)
                }
                Some(x) => f_read
                    .windows(x.len())
                    .position(|window| window == x)
                    .map(|i| i + x.len()),
            };

            let end_pos = match &flanks.1 {
                None => {
                    complete = RegionCompleteness::Partial3Prime;
                    Some(f_len - 1)
                }
                Some(x) => f_read.windows(x.len()).position(|window| window == x), // End pos doesn't need offsetting due to slice being exclusive
            };

            let fwd: Option<(Sequence, Vec<u8>, RegionCompleteness)> = match (start_pos, end_pos) {
                // Region not found
                (None, _) | (_, None) => None,

                // Region found
                (Some(start), Some(end)) => {
                    if start > end {
                        // Flank seqs overlap - no region contained within
                        None
                    } else {
                        Some((
                            f_read
                                .get(start..end)
                                .expect("Seq should contain region as extracted from within it")
                                .to_vec(),
                            f_qual
                                .get(start..end)
                                .expect("Qual should contain region as extracted from within it")
                                .to_vec(),
                            complete,
                        ))
                    }
                }
            };

            // Match in rev read
            let start_pos = match &flanks.0 {
                None => {
                    complete = RegionCompleteness::Partial5Prime;
                    Some(0)
                }
                Some(x) => r_read
                    .windows(x.len())
                    .position(|window| window == x)
                    .map(|i| i + x.len()),
            };

            let end_pos = match &flanks.1 {
                None => {
                    complete = RegionCompleteness::Partial3Prime;
                    Some(r_len - 1)
                }
                Some(x) => r_read.windows(x.len()).position(|window| window == x), // End pos doesn't need offsetting due to slice being exclusive
            };

            let rev: Option<(Sequence, Vec<u8>, RegionCompleteness)> = match (start_pos, end_pos) {
                // Region not found
                (None, _) | (_, None) => None,

                // Region found
                (Some(start), Some(end)) => {
                    if start > end {
                        // Flank seqs overlap - no region contained within
                        None
                    } else {
                        Some((
                            r_read
                                .get(start..end)
                                .expect("Seq should contain region as extracted from within it")
                                .to_vec(),
                            r_qual
                                .get(start..end)
                                .expect("Qual should contain region as extracted from within it")
                                .to_vec(),
                            complete,
                        ))
                    }
                }
            };

            // Determine which read to use
            match merge_seqs(fwd, rev, *len)? {
                Some((seq, comp)) => comb_key_vec.push(ObservedRegion::key(id.clone(), &seq, comp)),
                None => continue,
            }
        }

        let comb_key: CombinationKey = ObservedCombination::key(comb_key_vec);

        counts.add_or_increment_combination(&comb_key, record.group)?;
        progress.inc(1);
    }
    progress.finish();

    Ok(counts)
}

/// Count single end reads based on their in-frame position in the template
fn count_single_inframe(
    reads: ReadPairParser,
    lib_spec: &LibrarySpec,
    filter_config: FilterConfig,
    progress_style: &ProgressStyle,
) -> Result<ObservedCombinations, anyhow::Error> {
    info!("{}", INFRAME_START_SINGLE_MSG);

    let mut progress: Progress = Progress::from_style(
        progress_style,
        PROG_MSG,
        FINAL_MSG,
        None,
        INFRAME_LOG_INTERVAL,
    );

    let regions = lib_spec.variable_regions();

    let offset = match lib_spec.template_position(&lib_spec.forward_start_region) {
        Ok(x) => x.0,
        Err(_) => {
            return Err(LibSpecError::LibSpec {
                desc: "Forward start region missing from LibSpec".to_string(),
            }
            .into());
        }
    };
    debug!("F offset: {}", offset);

    let region_positions: Vec<(usize, usize)> = regions
        .iter()
        .map(|x| {
            lib_spec
                .template_position(x)
                .expect("The region should be in LibSpec")
        })
        .map(|x| (x.0 + offset, x.1 + offset))
        .collect();
    debug!("F region positions: {:?}", region_positions);

    let mut counts = ObservedCombinations::new(regions.clone(), filter_config);

    for result in reads {
        let record: ReadPair = result?;

        // Check if read should be filtered
        if !matches!(counts.filter_readpair(&record), FilterReason::None) {
            continue;
        }

        let mut comb_key_vec: Vec<RegionKey> = Vec::with_capacity(regions.len());
        let seq_len = record.forward.seq().len();

        for (id, pos) in zip(&regions, &region_positions) {
            let reg_seq: Sequence;
            let complete: RegionCompleteness;

            if (pos.0 < seq_len) && (pos.1 < seq_len) {
                // Region fully within read
                reg_seq = record
                    .forward
                    .seq()
                    .get(pos.0..pos.1)
                    .expect("Seq should contain region as just tested")
                    .to_vec();
                complete = RegionCompleteness::Complete;
            } else if (pos.0 < seq_len) && (pos.1 > seq_len - 1) {
                // Region partially within read
                reg_seq = record
                    .forward
                    .seq()
                    .get(pos.0..seq_len)
                    .expect("Seq should contain region as just tested")
                    .to_vec();
                complete = RegionCompleteness::Partial3Prime;
            } else {
                // Region outside read
                // Regions are fetched in order so all remaining will be outside to => break
                break;
            }

            comb_key_vec.push(ObservedRegion::key(id.clone(), &reg_seq, complete));
        }

        let comb_key: CombinationKey = ObservedCombination::key(comb_key_vec);

        counts.add_or_increment_combination(&comb_key, record.group)?;
        progress.inc(1);
    }
    progress.finish();

    Ok(counts)
}

/// Count paired end reads based on their in-frame position in the template
fn count_paired_inframe(
    reads: ReadPairParser,
    lib_spec: &LibrarySpec,
    filter_config: FilterConfig,
    progress_style: &ProgressStyle,
) -> Result<ObservedCombinations, anyhow::Error> {
    info!("{}", INFRAME_START_PAIRED_MSG);

    let mut progress: Progress = Progress::from_style(
        progress_style,
        PROG_MSG,
        FINAL_MSG,
        None,
        INFRAME_LOG_INTERVAL,
    );

    let regions = lib_spec.variable_regions();

    let region_lengths: Vec<usize> = regions
        .iter()
        .map(|x| {
            lib_spec
                .get_region(x)
                .expect("Region should exist as comes from lib")
                .len()
        })
        .collect();

    let f_offset = match lib_spec.template_position(&lib_spec.forward_start_region) {
        Ok(x) => x.0,
        Err(_) => {
            return Err(LibSpecError::LibSpec {
                desc: "Forward start region missing from LibSpec".to_string(),
            }
            .into());
        }
    };
    debug!("F offset: {}", f_offset);

    let f_region_positions: Vec<(usize, usize)> = regions
        .iter()
        .map(|x| {
            lib_spec
                .template_position(x)
                .expect("The region should be in LibSpec")
        })
        .map(|x| (x.0 + f_offset, x.1 + f_offset))
        .collect();
    debug!("F region positions: {:?}", f_region_positions);

    let r_offset = match lib_spec.template_position(&lib_spec.reverse_start_region) {
        Ok(x) => x.1,
        Err(_) => {
            return Err(LibSpecError::LibSpec {
                desc: "Reverse start region missing from LibSpec".to_string(),
            }
            .into());
        }
    };
    debug!("R offset: {}", r_offset);

    let r_region_positions: Vec<(usize, usize)> = regions
        .iter()
        .map(|x| {
            lib_spec
                .template_position(x)
                .expect("The region should be in LibSpec")
        })
        .map(|x| (x.0.abs_diff(r_offset), x.1.abs_diff(r_offset)))
        .collect();
    debug!("R region positions: {:?}", r_region_positions);

    let mut counts = ObservedCombinations::new(regions.clone(), filter_config);

    for (i, result) in reads.enumerate() {
        let record: ReadPair = result?;

        // Check if read needs to be filtered
        if !matches!(counts.filter_readpair(&record), FilterReason::None) {
            continue;
        }

        let mut comb_key_vec: Vec<RegionKey> = Vec::with_capacity(regions.len());

        let f_read = record.forward.seq();
        let r_read = match &record.reverse {
            Some(x) => revcomp(x.seq()),
            None => {
                return Err(ReadCountError::Error {
                    desc: format!("No reverse read found at read {i}"),
                }
                .into());
            }
        };

        let f_qual = record.forward.qual();
        let r_qual: Vec<u8> = match &record.reverse {
            Some(x) => x.qual().iter().rev().cloned().collect(),
            None => {
                return Err(ReadCountError::Error {
                    desc: format!("No reverse read found at read {i}"),
                }
                .into());
            }
        };

        let f_len = f_read.len();
        let r_len = r_read.len();

        for (id, len, f_pos, r_pos) in izip!(
            &regions,
            &region_lengths,
            &f_region_positions,
            &r_region_positions
        ) {
            let mut fwd: Option<(Sequence, Vec<u8>, RegionCompleteness)> = None;
            let mut rev: Option<(Sequence, Vec<u8>, RegionCompleteness)> = None;

            // Extract f_seq
            if (f_pos.0 < f_len) && (f_pos.1 < f_len) {
                // Region fully within read
                fwd = Some((
                    f_read
                        .get(f_pos.0..f_pos.1)
                        .expect("Seq should contain region as just tested")
                        .to_vec(),
                    f_qual
                        .get(f_pos.0..f_pos.1)
                        .expect("Seq should contain region as just tested")
                        .to_vec(),
                    RegionCompleteness::Complete,
                ));
            } else if (f_pos.0 < f_len) && (f_pos.1 > f_len - 1) {
                // Region partially within read
                fwd = Some((
                    f_read
                        .get(f_pos.0..f_len)
                        .expect("Seq should contain region as just tested")
                        .to_vec(),
                    f_qual
                        .get(f_pos.0..f_len)
                        .expect("Seq should contain region as just tested")
                        .to_vec(),
                    RegionCompleteness::Partial3Prime,
                ));
            }

            // Extract r_seq
            if (r_len > r_pos.0) && (r_len > r_pos.1) {
                // Region fully within read
                rev = Some((
                    r_read
                        .get((r_len - r_pos.0)..(r_len - r_pos.1))
                        .expect("Seq should contain region as just tested")
                        .to_vec(),
                    r_qual
                        .get((r_len - r_pos.0)..(r_len - r_pos.1))
                        .expect("Seq should contain region as just tested")
                        .to_vec(),
                    RegionCompleteness::Complete,
                ));
            } else if (r_len <= r_pos.0) && (r_len > r_pos.1) {
                // Region partially within read
                rev = Some((
                    r_read
                        .get(0..(r_len - r_pos.1))
                        .expect("Seq should contain region as just tested")
                        .to_vec(),
                    r_qual
                        .get(0..(r_len - r_pos.1))
                        .expect("Seq should contain region as just tested")
                        .to_vec(),
                    RegionCompleteness::Partial5Prime,
                ));
            }

            // Determine which read to use
            match merge_seqs(fwd, rev, *len)? {
                Some((seq, comp)) => comb_key_vec.push(ObservedRegion::key(id.clone(), &seq, comp)),
                None => continue,
            }
        }

        let comb_key: CombinationKey = ObservedCombination::key(comb_key_vec);

        counts.add_or_increment_combination(&comb_key, record.group)?;
        progress.inc(1);
    }
    progress.finish();

    Ok(counts)
}

/// Count entire single end reads
fn count_single_raw(
    reads: ReadPairParser,
    filter_config: FilterConfig,
    progress_style: &ProgressStyle,
) -> Result<ObservedCombinations, anyhow::Error> {
    info!("{}", RAW_START_SINGLE_MSG);

    let mut progress: Progress =
        Progress::from_style(progress_style, PROG_MSG, FINAL_MSG, None, RAW_LOG_INTERVAL);

    let mut counts = ObservedCombinations::new(vec!["seq".to_string()], filter_config);

    for result in reads {
        let record: ReadPair = result?;

        // Check if read should be filtered
        if !matches!(counts.filter_readpair(&record), FilterReason::None) {
            continue;
        }

        let reg_key: RegionKey = ObservedRegion::key(
            "seq".to_string(),
            record.forward.seq(),
            RegionCompleteness::Complete,
        );
        let comb_key: CombinationKey = ObservedCombination::key(vec![reg_key]);

        counts.add_or_increment_combination(&comb_key, record.group)?;
        progress.inc(1);
    }
    progress.finish();

    Ok(counts)
}

/// Count entire paired end reads
fn count_paired_raw(
    reads: ReadPairParser,
    filter_config: FilterConfig,
    progress_style: &ProgressStyle,
) -> Result<ObservedCombinations, anyhow::Error> {
    info!("{}", RAW_START_PAIRED_MSG);

    let mut progress: Progress =
        Progress::from_style(progress_style, PROG_MSG, FINAL_MSG, None, RAW_LOG_INTERVAL);

    let mut counts =
        ObservedCombinations::new(vec!["fwd".to_string(), "rev".to_string()], filter_config);

    for result in reads {
        let record: ReadPair = result?;

        // Check if read should be filtered
        if !matches!(counts.filter_readpair(&record), FilterReason::None) {
            continue;
        }

        let fwd_key: RegionKey = ObservedRegion::key(
            "fwd".to_string(),
            record.forward.seq(),
            RegionCompleteness::Complete,
        );

        // Add check for reverse read
        let rev_key: RegionKey = match &record.reverse {
            Some(rev) => ObservedRegion::key(
                "rev".to_string(),
                rev.seq(),
                RegionCompleteness::Complete,
            ),
            None => {
                return Err(ReadCountError::Error {
                    desc: format!("No reverse read found in paired raw mode")
                }.into());
            }
        };

        let comb_key: CombinationKey = ObservedCombination::key(vec![fwd_key, rev_key]);

        counts.add_or_increment_combination(&comb_key, record.group)?;
        progress.inc(1);
    }
    progress.finish();

    Ok(counts)
}
