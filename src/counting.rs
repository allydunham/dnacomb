//! Counting and region-extraction algorithms for structured sequencing reads.
//!
//! This module contains the core logic for converting sequencing reads into
//! `ObservedCombinations`. It supports several extraction strategies with
//! different robustness/performance tradeoffs:
//! - alignment to a full template,
//! - flanking-pattern matching,
//! - in-frame positional extraction,
//! - and raw full-read counting.
//!
//! It also contains helper logic for pairing forward/reverse evidence within
//! a region and for dispatching counting work across multiple threads. The
//! `count_reads` function is one of the key entry points into DNAComb.
use bio::alignment::AlignmentOperation;
use bio::alignment::distance::hamming;
use bio::alignment::pairwise::{Aligner, MatchFunc, Scoring};
use bio::alphabets::dna::revcomp;
use bio::bio_types::sequence::Sequence;
use clap::ValueEnum;
use itertools::izip;
use log::{debug, info};
use std::iter::{repeat_with, zip};
use std::sync::Once;

use crate::combination::CombinationKey;
use crate::combinations::{CacheHit, ObservedCombinations};
use crate::errors::{AlignmentInfo, LibSpecError, ReadCountError};
use crate::filters::FilterConfig;
use crate::interning::seq_from_bytes;
use crate::lib_spec::{FlankingSequences, LibrarySpec};
use crate::logging::{Progress, ProgressStyle};
use crate::parsing::{ReadPairProducer, ThreadedReadPairParser};
use crate::region::{RegionCompleteness, RegionKey};
use crate::seqs::{ReadPair, SeqPair};
use crate::utils::mean_quality;

/// Query interval and completeness assigned to one expected region:
/// `(start, end, completeness)`.
///
/// Coordinates are query-sequence positions using the 1-based convention derived
/// from Rust-Bio alignment paths.
pub type AlignmentPosition = (usize, usize, RegionCompleteness);

/// Observed region sequence, quality values, and completeness status extracted
/// from one read.e
type RegionMatch = (Sequence, Vec<u8>, RegionCompleteness);

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
/// Logging interval in raw mode
#[cfg(debug_assertions)]
const RAW_LOG_INTERVAL: u64 = 1000;
#[cfg(not(debug_assertions))]
const RAW_LOG_INTERVAL: u64 = 10000000;

/// Alignment scoring scheme for template-based region extraction.
///
/// This wraps the match/mismatch/gap parameters used for semi-global alignment
/// against the LibSpec template. Matches involving `N` use a separate score so
/// that variable regions represented by `N` in the template can align flexibly
/// without being treated as either full matches or full mismatches.
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

    /// Build a Rust-Bio `Scoring` object using this scorer for match evaluation.
    pub fn get_scoring(self) -> Scoring<AlignmentScorer> {
        Scoring::new(self.gap_open, self.gap_extend, self)
    }
}

impl MatchFunc for AlignmentScorer {
    /// Score a single aligned base pair.
    ///
    /// Return a (mis)match score for a pair of bytes. `N` is treated specially so
    /// that template positions representing variable regions can absorb sequence
    /// variation with an intermediate penalty.
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

/// Map template-region intervals onto query-sequence intervals by walking an
/// alignment path.
///
/// `region_positions` should contain template intervals in ascending order,
/// expressed using the same 1-based coordinate convention as Rust-Bio alignment
/// paths. The returned vector has the same length, with each element giving the
/// corresponding query interval and `RegionCompleteness` status if that region
/// could be located in the alignment.
///
/// Regions may be returned as:
/// - complete,
/// - truncated at the 5' or 3' end,
/// - or absent (`None`) if no sequence could be assigned.
///
/// This helper is used by alignment-based counting to translate template-relative
/// region definitions into observed query substrings.
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

/// This combines two region observations derived from opposite reads into a
/// single observed sequence plus a `RegionCompleteness` state.
///
/// Behaviour depends on completeness:
/// - if only one side is present, that side is used;
/// - if one side is complete, it takes precedence;
/// - if both sides are partial and non-overlapping, the result is marked
///   `MissingCenter`;
/// - if both sides are partial and overlapping, the result is marked
///   `Overlapping`;
/// - if both sides cover the same span, the higher-quality sequence is chosen.
///
/// This function assumes a fixed expected region length and uses that to infer
/// whether partial forward/reverse observations overlap or leave a gap. In general
/// using hte max length for variable length regions gives reasonable outcomes.
///
/// This function is intentionally quite a basic attempt to merge information and
/// it is recommended (including in a warning on use) that using a dedicated read merger
/// will be more robust when overlap is expected.
fn merge_seqs(
    fwd: Option<RegionMatch>,
    rev: Option<RegionMatch>,
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
        // Weird situation but some read trimming could lead here potentially
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
        (RegionCompleteness::Partial3Prime, RegionCompleteness::Partial5Prime) => {
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
        (RegionCompleteness::Partial5Prime, RegionCompleteness::Partial3Prime) => {
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

/// Extract variable regions by searching for flanking fixed-sequence patterns.
///
/// Each variable region is defined by one of:
/// - a start-of-read open flank,
/// - an end-of-read open flank,
/// - or fixed flanks on both sides.
///
/// Matching proceeds left-to-right through the read. Regions found form a
/// continuous subsequence of the expected region list: once a required flank is
/// missed, downstream regions are not recovered later in the read. This simple algorithm
/// is robust when patterns are unique per region but can produce unexpected results
/// when they are not sufficiently different as the wrong starting region can be identified.
/// However, it is important to also deal with messy/truncated reads rather than require
/// all regions in order.
///
/// `tolerance` allows a bounded number of mismatches in flank-pattern matching.
/// This improves robustness to sequencing errors in fixed regions, but increases
/// the risk of ambiguous or incorrect matches when flanking patterns are similar.
///
/// Returns a vector of hits, one per input region with None if the region is missing or
/// Some((Sequence, Phred Quality, RegionCompleteness)) tuple
fn match_flank_patterns(
    seq: &[u8],
    qual: &[u8],
    flanks: &[FlankingSequences],
    tolerance: u64,
) -> Result<Vec<Option<RegionMatch>>, ReadCountError> {
    let mut out: Vec<Option<RegionMatch>> = repeat_with(|| None).take(flanks.len()).collect();

    // If no flanking regions then find nothing
    if flanks.is_empty() {
        return Ok(out);
    }

    // Check flank sequence is valid
    match LibrarySpec::validate_flank_seqs(flanks) {
        Ok(_) => {}
        Err(e) => {
            return Err(ReadCountError::Error {
                desc: format!("Invalid flanking sequences: {}", e),
            });
        }
    };

    // Initialise seach space
    let mut pos: usize = 0; // Position in sequence to search for match
    let mut end: usize; // end of current flank seq to match
    let mut reg: usize = 0; // region being matched
    let mut reg_start: usize = 0; // Start point of region seq
    let mut flank_seq: &Sequence; // Sequence being searched for
    let mut open: bool = false; // whether the region start is found
    let mut dist: u64; // Distance to region

    // Find opening region to assign start point.
    //
    // Scan along the sequence looking for:
    // - OpenStart(end): end flank closes a 5' partial region.
    // - Internal(start, end): start flank opens an internal region, end flank closes a truncated partial 5' region
    // - OpenEnd(start): start flank opens a 3' partial terminal region.
    'outer: while pos < seq.len() {
        for (i, flank) in flanks.iter().enumerate() {
            match flank {
                FlankingSequences::Unflanked => unreachable!("Unflanked already checked"),

                FlankingSequences::OpenStart(close) => {
                    end = pos + close.len();
                    if end > seq.len() {
                        continue;
                    }

                    dist = hamming(close, &seq[pos..end]);
                    if dist <= tolerance {
                        out[i] = Some((
                            seq[0..pos].to_vec(),
                            qual[0..pos].to_vec(),
                            RegionCompleteness::Partial5Prime,
                        ));

                        reg = i + 1;
                        open = false;
                        break 'outer;
                    }
                }

                FlankingSequences::Internal(start, close) => {
                    // Prefer the normal start flank if both could match at this position.
                    end = pos + start.len();
                    if end <= seq.len() {
                        dist = hamming(start, &seq[pos..end]);
                        if dist <= tolerance {
                            reg = i;
                            open = true;
                            reg_start = end;
                            pos = end;
                            break 'outer;
                        }
                    }

                    // Fallback: first visible region started before the read,
                    // and we have found its closing flank.
                    end = pos + close.len();
                    if end <= seq.len() {
                        dist = hamming(close, &seq[pos..end]);
                        if dist <= tolerance {
                            out[i] = Some((
                                seq[0..pos].to_vec(),
                                qual[0..pos].to_vec(),
                                RegionCompleteness::Partial5Prime,
                            ));

                            reg = i + 1;
                            open = false;
                            break 'outer;
                        }
                    }
                }

                FlankingSequences::OpenEnd(start) => {
                    end = pos + start.len();
                    if end > seq.len() {
                        continue;
                    }

                    dist = hamming(start, &seq[pos..end]);
                    if dist <= tolerance {
                        out[i] = Some((
                            seq[end..seq.len()].to_vec(),
                            qual[end..seq.len()].to_vec(),
                            RegionCompleteness::Partial3Prime,
                        ));

                        return Ok(out);
                    }
                }
            }
        }

        pos += 1;
    }

    // Walk the remaining sequence and region list in parallel, adding each newly found region to output.
    // Before finding the first region must consider all options at each point
    'outer: while reg < flanks.len() && pos < seq.len() {
        // Get region sequence
        flank_seq = match (open, &flanks[reg]) {
            (_, FlankingSequences::Unflanked) => unreachable!("Unflanked already checked"),
            (_, FlankingSequences::OpenStart(..)) => {
                unreachable!("Open start can only be first and already processed")
            }
            (_, FlankingSequences::OpenEnd(start)) => start,
            (false, FlankingSequences::Internal(start, _)) => start,
            (true, FlankingSequences::Internal(_, end)) => end,
        };

        end = pos + flank_seq.len();

        // Loop forward to find then process
        'inner: while pos < seq.len() {
            // If seq runs off end without finding we've exhausted
            if end > seq.len() {
                break 'outer;
            }

            dist = hamming(flank_seq, &seq[pos..end]);

            // If not found, continue (do this way to save indent below)
            if dist > tolerance {
                pos += 1;
                end += 1;
                continue;
            }

            match (open, &flanks[reg]) {
                (_, FlankingSequences::Unflanked) => unreachable!("Unflanked already checked"),
                (_, FlankingSequences::OpenStart(..)) => {
                    unreachable!("Open start can only be first and already processed")
                }
                (true, FlankingSequences::OpenEnd(..)) => {
                    unreachable!("Open end only has a start and is then processed below")
                }
                (false, FlankingSequences::OpenEnd(..)) => {
                    out[reg] = Some((
                        seq[end..seq.len()].to_vec(),
                        qual[end..seq.len()].to_vec(),
                        RegionCompleteness::Partial3Prime,
                    ));

                    // Open end must finish the seq so break outer
                    break 'outer;
                }
                (true, FlankingSequences::Internal(..)) => {
                    // Found end - set out[reg] and move to next region
                    out[reg] = Some((
                        seq[reg_start..pos].to_vec(),
                        qual[reg_start..pos].to_vec(),
                        RegionCompleteness::Complete,
                    ));

                    open = false;
                    reg += 1;
                    break 'inner;
                }
                (false, FlankingSequences::Internal(..)) => {
                    // Found start - set open and search for end
                    pos = end;
                    open = true;
                    reg_start = end;
                    break 'inner;
                }
            }
        }
    }

    // If an Internal region is still open, close it as partial 3'
    if open && matches!(&flanks[reg], FlankingSequences::Internal(..)) {
        out[reg] = Some((
            seq[reg_start..seq.len()].to_vec(),
            qual[reg_start..seq.len()].to_vec(),
            RegionCompleteness::Partial3Prime,
        ));
    }

    Ok(out)
}

/// Extract ObservedCombinations from a worker thread JoinHandle
///
/// Returns the appropriate result or an error that can be raised via ?
fn join_observed_combinations(
    handle: Option<std::thread::JoinHandle<Result<ObservedCombinations, anyhow::Error>>>,
) -> Result<ObservedCombinations, anyhow::Error> {
    match handle {
        Some(j) => match j.join() {
            Ok(r) => r,
            Err(e) => {
                if let Some(msg) = e.downcast_ref::<&str>() {
                    Err(ReadCountError::Error {
                        desc: format!("Counting thread paniced: {msg}"),
                    }
                    .into())
                } else if let Some(msg) = e.downcast_ref::<String>() {
                    Err(ReadCountError::Error {
                        desc: format!("Counting thread paniced: {msg}"),
                    }
                    .into())
                } else {
                    Err(ReadCountError::Error {
                        desc: "Counting thread paniced with an unknown payload".to_string(),
                    }
                    .into())
                }
            }
        },
        None => Err(ReadCountError::Error {
            desc: "No counting threads returned objects".to_string(),
        }
        .into()),
    }
}

/// Region-extraction strategy to use during counting.
///
/// Different modes trade off robustness, assumptions about read structure, and
/// speed.
#[derive(Clone, ValueEnum, Debug, Copy, PartialEq, Eq)]
pub enum CountMode {
    /// Count complete read sequences without structured region extraction.
    FullRead,

    /// Extract regions from their expected in-read positions.
    ///
    /// Fast, but assumes reads are already in frame and region lengths are fixed.
    Inframe,

    /// Extract regions using fixed flanking sequences.
    ///
    /// Faster than alignment and supports variable-length regions, but relies on
    /// intact flanking sequence and only recovers a continuous block of regions.
    Pattern,

    /// Extract regions by semi-global alignment to the full template.
    ///
    /// Most robust and the default choice for complex or noisy data, but slowest.
    Align,
}

/// Count observed read forms from a sequencing dataset.
///
/// This is the main entry point for region extraction and counting. It dispatches
/// to one of the supported counting modes, applies configured read/alignment
/// filtering, optionally caches sequence-derived results, and can parallelise the
/// counting step across worker threads.
///
/// Behaviour depends on `mode`:
/// - `Align` requires a `LibrarySpec` and an `AlignmentScorer`,
/// - `Pattern` requires a `LibrarySpec`, `pattern_length`, and `pattern_tolerance`,
/// - `Inframe` requires a `LibrarySpec`,
/// - `FullRead` can operate without a `LibrarySpec`.
///
/// If `full_seq` is true, the full read sequence(s) are stored alongside each
/// observed combination; otherwise only extracted regions are tracked.
///
/// The returned `ObservedCombinations` contains unfiltered counts, filtered-read
/// summaries, and any read-level cache accumulated during counting. Library
/// comparison is not performed here.
pub fn count_reads<T: ReadPairProducer>(
    reads: T,
    lib_spec: &Option<LibrarySpec>,
    mode: CountMode,
    full_seq: bool,
    filter_config: FilterConfig,
    alignment_scorer: Option<AlignmentScorer>,
    pattern_length: Option<usize>,
    pattern_tolerance: Option<u64>,
    cache: bool,
    threads: usize,
    progress_style: Option<&ProgressStyle>,
) -> Result<ObservedCombinations, anyhow::Error> {
    let default_progress = ProgressStyle::new(None, false);
    let progress = progress_style.unwrap_or(&default_progress);

    match threads.cmp(&1) {
        std::cmp::Ordering::Less => Err(ReadCountError::Error {
            desc: "Threads must be >0".to_string(),
        }
        .into()),
        std::cmp::Ordering::Equal => {
            // Determine type of matching desired and despatch as appropriate
            match (
                reads.has_reverse(),
                lib_spec,
                mode,
                alignment_scorer,
                pattern_length,
                pattern_tolerance,
            ) {
                (_, _, CountMode::Align, None, _, _) => Err(ReadCountError::Error {
                    desc: "Mode is 'align' but no AlignmentScorer passed".to_string(),
                }
                .into()),
                (false, Some(lib_spec), CountMode::Align, Some(a), _, _) => {
                    count_single_align(reads, lib_spec, full_seq, filter_config, a, cache, progress)
                }
                (true, Some(lib_spec), CountMode::Align, Some(a), _, _) => {
                    count_paired_align(reads, lib_spec, full_seq, filter_config, a, cache, progress)
                }

                (_, _, CountMode::Pattern, _, None, _) | (_, _, CountMode::Pattern, _, _, None) => {
                    Err(ReadCountError::Error {
                        desc: "Mode is 'pattern' but pattern length and/or tolerance is missing"
                            .to_string(),
                    }
                    .into())
                }
                (false, Some(lib_spec), CountMode::Pattern, _, Some(len), Some(tol)) => {
                    count_single_pattern(
                        reads,
                        lib_spec,
                        full_seq,
                        filter_config,
                        len,
                        tol,
                        cache,
                        progress,
                    )
                }
                (true, Some(lib_spec), CountMode::Pattern, _, Some(len), Some(tol)) => {
                    count_paired_pattern(
                        reads,
                        lib_spec,
                        full_seq,
                        filter_config,
                        len,
                        tol,
                        cache,
                        progress,
                    )
                }

                (false, Some(lib_spec), CountMode::Inframe, _, _, _) => {
                    count_single_inframe(reads, lib_spec, full_seq, filter_config, cache, progress)
                }
                (true, Some(lib_spec), CountMode::Inframe, _, _, _) => {
                    count_paired_inframe(reads, lib_spec, full_seq, filter_config, cache, progress)
                }

                (_, Some(_), CountMode::FullRead, _, _, _) => {
                    count_raw(reads, filter_config, progress)
                }
                (_, None, _, _, _, _) => count_raw(reads, filter_config, progress),
            }
        }
        std::cmp::Ordering::Greater => {
            // Set up communication channel
            let (read_tx, read_rx) = crossbeam::channel::bounded(10 * threads);
            let mut handles = Vec::new();

            // Spin up worker threads ready to receive data
            for _ in 0..threads {
                // each thread gets its own Read Receiver and Count Sender clone
                let rx = read_rx.clone();
                let rev_reads = reads.has_reverse();
                let group = reads.group().clone();
                let max_reads = reads.max_reads();
                let lspec = lib_spec.as_ref().cloned();
                let fconf = filter_config.clone();
                let mut pstyle = progress_style.cloned();
                if let Some(ref mut p) = pstyle {
                    p.use_thread_id = true
                }

                handles.push(std::thread::spawn(move || {
                    let local_reads = ThreadedReadPairParser::new(rx, rev_reads, group, max_reads);
                    count_reads(
                        local_reads,
                        &lspec,
                        mode,
                        full_seq,
                        fconf,
                        alignment_scorer,
                        pattern_length,
                        pattern_tolerance,
                        cache,
                        1,
                        pstyle.as_ref(),
                    )
                }));
            }

            // Produce reads on main thread - as they are sent they will be processed on worker threads
            for read in reads {
                read_tx
                    .send(read)
                    .expect("Region extraction thread send failed");
            }
            drop(read_tx);
            drop(read_rx);

            // Unpack initial count object
            let mut final_counts: ObservedCombinations = join_observed_combinations(handles.pop())?;

            // Merge rest of the results
            for handle in handles {
                let new_counts = join_observed_combinations(Some(handle))?;
                final_counts.merge(new_counts)?
            }

            Ok(final_counts)
        }
    }
}

/// Count single-end reads by semi-global alignment to the LibSpec template.
///
/// Variable regions are located by mapping the alignment path back onto template
/// region intervals. This is the most robust structured counting mode, but also
/// the slowest.
fn count_single_align<T: ReadPairProducer>(
    reads: T,
    lib_spec: &LibrarySpec,
    full_seq: bool,
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

    // Initialise aligner
    let scoring = alignment_scorer.get_scoring();
    let mut aligner = Aligner::with_capacity_and_scoring(400, 150, scoring);
    let template = lib_spec.template_sequence();

    // Iterate over reads
    info!("Printing example alignments. Check these are reasonable to help diagnose problems");
    for (i, result) in reads.enumerate() {
        let record: ReadPair = result?;

        // Check if read should be filtered
        if counts.filter_readpair(&record, true).is_some() {
            progress.inc(1);
            continue;
        }

        // Check if the read has been cached
        if cache && counts.check_cache(&record, true)?.is_some() {
            progress.inc(1);
            continue;
        }

        let read_alignment = aligner.semiglobal(record.forward.seq(), &template);
        let alignment_path = read_alignment.path();

        // Filter reads with too low scores or otherwise unmatched
        match counts.filter_alignment(&record, &read_alignment, None, true) {
            None => {}
            Some(reason) => {
                if cache {
                    counts.cache(record.into_seqpair(), CacheHit::Filter(reason));
                }
                progress.inc(1);
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

        let mut comb_key: CombinationKey = CombinationKey::new(
            if full_seq {
                Some(SeqPair::from_readpair(&record))
            } else {
                None
            },
            Vec::with_capacity(query_regions.len()),
        );

        for (id, opt_pos) in zip(&regions, &query_regions) {
            if let Some(pos) = opt_pos {
                match record.forward.seq().get((pos.0 - 1)..(pos.1 - 1)) {
                    Some(s) => {
                        comb_key.regions.push(RegionKey::new(
                            id.clone(),
                            // Offset seq lookup - rust vec 0 based and AlignmentPath 1 based
                            seq_from_bytes(s),
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

        counts.add_or_increment_combination(&comb_key, record.group.clone())?;

        if cache {
            counts.cache(record.into_seqpair(), CacheHit::Comb(comb_key));
        }
        progress.inc(1);
    }
    progress.finish();

    Ok(counts)
}

/// Count paired-end reads by aligning forward and reverse reads independently
/// to the LibSpec template and then merging per-region evidence.
///
/// Reverse reads are reverse-complemented before alignment. Where both reads
/// contribute evidence for the same region, `merge_seqs` is used to reconcile
/// complete, partial, overlapping, or gapped observations.
fn count_paired_align<T: ReadPairProducer>(
    reads: T,
    lib_spec: &LibrarySpec,
    full_seq: bool,
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

    // Initialise aligner
    let scoring = alignment_scorer.get_scoring();
    let mut aligner = Aligner::with_capacity_and_scoring(400, 150, scoring);
    let template = lib_spec.template_sequence();

    // Iterate over reads
    info!("Printing example alignments. Check these are reasonable to help diagnose problems");
    for (i, result) in reads.enumerate() {
        let record: ReadPair = result?;

        // Check if read should be filtered
        if counts.filter_readpair(&record, true).is_some() {
            progress.inc(1);
            continue;
        }

        // Check if the read has been cached
        if cache && counts.check_cache(&record, true)?.is_some() {
            progress.inc(1);
            continue;
        }

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
        match counts.filter_alignment(&record, &f_alignment, Some(&r_alignment), true) {
            None => {}
            Some(reason) => {
                if cache {
                    counts.cache(record.into_seqpair(), CacheHit::Filter(reason));
                }
                progress.inc(1);
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

        let mut comb_key: CombinationKey = CombinationKey::new(
            if full_seq {
                Some(SeqPair::from_readpair(&record))
            } else {
                None
            },
            Vec::with_capacity(std::cmp::max(f_regions.len(), r_regions.len())),
        );

        for (id, len, f_pos, r_pos) in izip!(&regions, &region_lengths, &f_regions, &r_regions) {
            match (f_pos, r_pos) {
                (None, None) => continue,
                (None, Some(r)) => {
                    match r_read.get((r.0 - 1)..(r.1 - 1)) {
                        Some(s) => {
                            comb_key.regions.push(RegionKey::new(
                                id.clone(),
                                // Offset seq lookup - rust vec 0 based and AlignmentPath 1 based
                                seq_from_bytes(s),
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
                                    pretty_alignment: r_alignment.pretty(&r_read, &template, 100),
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
                            comb_key.regions.push(RegionKey::new(
                                id.clone(),
                                // Offset seq lookup - rust vec 0 based and AlignmentPath 1 based
                                seq_from_bytes(s),
                                f.2,
                            ));
                        }
                        None => {
                            return Err(ReadCountError::BadAlignment {
                                alignment: Box::new(AlignmentInfo {
                                    read_id: record.forward.id().to_string(),
                                    read_number: i,
                                    pretty_alignment: f_alignment.pretty(f_read, &template, 100),
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
                                    pretty_alignment: f_alignment.pretty(f_read, &template, 100),
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
                                    pretty_alignment: r_alignment.pretty(&r_read, &template, 100),
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
                        Some((seq, comp)) => comb_key.regions.push(RegionKey::new(
                            id.clone(),
                            seq_from_bytes(&seq),
                            comp,
                        )),
                        None => continue,
                    };
                }
            }
        }

        counts.add_or_increment_combination(&comb_key, record.group.clone())?;

        if cache {
            counts.cache(record.into_seqpair(), CacheHit::Comb(comb_key));
        }

        progress.inc(1);
    }
    progress.finish();

    Ok(counts)
}

/// Count single-end reads by identifying variable regions from flanking fixed-sequence patterns.
///
/// This mode is faster than alignment and still supports variable-length regions,
/// but depends on reliable flanking sequence and may fail to recover downstream
/// regions once an earlier flank is missed.
fn count_single_pattern<T: ReadPairProducer>(
    reads: T,
    lib_spec: &LibrarySpec,
    full_seq: bool,
    filter_config: FilterConfig,
    pattern_length: usize,
    pattern_tolerance: u64,
    cache: bool,
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

    let flank_regions = lib_spec.get_all_flanking_regions(pattern_length)?;

    info!(
        "Extracting regions using flank seqs: {:?}",
        flank_regions
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
    );

    // Count reads
    let mut counts = ObservedCombinations::new(regions.clone(), filter_config);

    for result in reads {
        let record: ReadPair = result?;

        // Check if read should be filtered
        if counts.filter_readpair(&record, true).is_some() {
            progress.inc(1);
            continue;
        }

        // Check if the read has been cached
        if cache && counts.check_cache(&record, true)?.is_some() {
            progress.inc(1);
            continue;
        }

        let region_matches = match_flank_patterns(
            record.forward.seq(),
            record.forward.qual(),
            &flank_regions,
            pattern_tolerance,
        )?;

        let comb_key: CombinationKey = CombinationKey::new(
            if full_seq {
                Some(SeqPair::from_readpair(&record))
            } else {
                None
            },
            zip(&regions, region_matches)
                .filter_map(|(id, reg)| match reg {
                    Some(r) => Some(RegionKey::new(id.clone(), seq_from_bytes(&r.0), r.2)),
                    None => None,
                })
                .collect(),
        );

        counts.add_or_increment_combination(&comb_key, record.group.clone())?;

        if cache {
            counts.cache(record.into_seqpair(), CacheHit::Comb(comb_key));
        }

        progress.inc(1);
    }
    progress.finish();

    Ok(counts)
}

/// Count paired-end reads by identifying variable regions from flanking fixed-sequence patterns.
///
/// This mode is faster than alignment and still supports variable-length regions,
/// but depends on reliable flanking sequence and may fail to recover downstream
/// regions once an earlier flank is missed.
fn count_paired_pattern<T: ReadPairProducer>(
    reads: T,
    lib_spec: &LibrarySpec,
    full_seq: bool,
    filter_config: FilterConfig,
    pattern_length: usize,
    pattern_tolerance: u64,
    cache: bool,
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
    let n_regions = regions.len();

    let region_lengths: Vec<usize> = regions
        .iter()
        .map(|x| {
            lib_spec
                .get_region(x)
                .expect("Region should exist as comes from lib")
                .len()
        })
        .collect();

    let flank_regions = lib_spec.get_all_flanking_regions(pattern_length)?;

    info!(
        "Extracting regions using flank seqs: {:?}",
        flank_regions
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
    );

    // Count reads
    let mut counts = ObservedCombinations::new(regions.clone(), filter_config);

    for (i, result) in reads.enumerate() {
        let record: ReadPair = result?;

        // Check if read should be filtered
        if counts.filter_readpair(&record, true).is_some() {
            progress.inc(1);
            continue;
        }

        // Check if the read has been cached
        if cache && counts.check_cache(&record, true)?.is_some() {
            progress.inc(1);
            continue;
        }

        let f_seq = record.forward.seq();
        let r_seq = match &record.reverse {
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

        let f_matches = match_flank_patterns(f_seq, f_qual, &flank_regions, pattern_tolerance)?;

        let r_matches = match_flank_patterns(&r_seq, &r_qual, &flank_regions, pattern_tolerance)?;

        let mut comb_key: CombinationKey = CombinationKey::new(
            if full_seq {
                Some(SeqPair::from_readpair(&record))
            } else {
                None
            },
            Vec::with_capacity(n_regions),
        );

        for (id, len, fwd, rev) in izip!(&regions, &region_lengths, f_matches, r_matches) {
            if let Some(merged) = merge_seqs(fwd, rev, *len)? {
                comb_key.regions.push(RegionKey::new(
                    id.clone(),
                    seq_from_bytes(&merged.0),
                    merged.1,
                ));
            }
        }

        counts.add_or_increment_combination(&comb_key, record.group.clone())?;

        if cache {
            counts.cache(record.into_seqpair(), CacheHit::Comb(comb_key));
        }

        progress.inc(1);
    }
    progress.finish();

    Ok(counts)
}

/// Count single-end reads by extracting regions from expected in-read positions.
///
/// This mode assumes reads begin at the configured LibSpec start region and that
/// region boundaries can be inferred directly from template coordinates. It is
/// therefore best suited to fixed-length, well-framed reads.
fn count_single_inframe<T: ReadPairProducer>(
    reads: T,
    lib_spec: &LibrarySpec,
    full_seq: bool,
    filter_config: FilterConfig,
    cache: bool,
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
        if counts.filter_readpair(&record, true).is_some() {
            progress.inc(1);
            continue;
        }

        // Check if the read has been cached
        if cache && counts.check_cache(&record, true)?.is_some() {
            progress.inc(1);
            continue;
        }

        let mut comb_key: CombinationKey = CombinationKey::new(
            if full_seq {
                Some(SeqPair::from_readpair(&record))
            } else {
                None
            },
            Vec::with_capacity(regions.len()),
        );

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

            comb_key.regions.push(RegionKey::new(
                id.clone(),
                seq_from_bytes(&reg_seq),
                complete,
            ));
        }

        counts.add_or_increment_combination(&comb_key, record.group.clone())?;

        if cache {
            counts.cache(record.into_seqpair(), CacheHit::Comb(comb_key));
        }

        progress.inc(1);
    }
    progress.finish();

    Ok(counts)
}

/// Count paired-end reads by extracting regions from expected in-read positions.
///
/// This mode assumes reads begin at the configured LibSpec start region and that
/// region boundaries can be inferred directly from template coordinates. It is
/// therefore best suited to fixed-length, well-framed reads.
fn count_paired_inframe<T: ReadPairProducer>(
    reads: T,
    lib_spec: &LibrarySpec,
    full_seq: bool,
    filter_config: FilterConfig,
    cache: bool,
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
        if counts.filter_readpair(&record, true).is_some() {
            progress.inc(1);
            continue;
        }

        // Check if the read has been cached
        if cache && counts.check_cache(&record, true)?.is_some() {
            progress.inc(1);
            continue;
        }

        let mut comb_key: CombinationKey = CombinationKey::new(
            if full_seq {
                Some(SeqPair::from_readpair(&record))
            } else {
                None
            },
            Vec::with_capacity(regions.len()),
        );

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
            let mut fwd: Option<RegionMatch> = None;
            let mut rev: Option<RegionMatch> = None;

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
                Some((seq, comp)) => {
                    comb_key
                        .regions
                        .push(RegionKey::new(id.clone(), seq_from_bytes(&seq), comp))
                }
                None => continue,
            }
        }

        counts.add_or_increment_combination(&comb_key, record.group.clone())?;

        if cache {
            counts.cache(record.into_seqpair(), CacheHit::Comb(comb_key));
        }

        progress.inc(1);
    }
    progress.finish();

    Ok(counts)
}

/// Count complete read sequences without structured region extraction.
///
/// This mode still applies read-level filtering, but does not use the LibSpec
/// region structure and stores each full read (or read pair) as a distinct
/// observed combination.
fn count_raw<T: ReadPairProducer>(
    reads: T,
    filter_config: FilterConfig,
    progress_style: &ProgressStyle,
) -> Result<ObservedCombinations, anyhow::Error> {
    info!("{}", RAW_START_SINGLE_MSG);

    let mut progress: Progress =
        Progress::from_style(progress_style, PROG_MSG, FINAL_MSG, None, RAW_LOG_INTERVAL);

    let mut counts = ObservedCombinations::new(vec![], filter_config);

    for result in reads {
        let record: ReadPair = result?;

        // Check if read should be filtered
        if counts.filter_readpair(&record, true).is_some() {
            progress.inc(1);
            continue;
        }

        let comb_key: CombinationKey =
            CombinationKey::new(Some(SeqPair::from_readpair(&record)), vec![]);

        counts.add_or_increment_combination(&comb_key, record.group)?;
        progress.inc(1);
    }
    progress.finish();

    Ok(counts)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::errors::ReadPairError;
    use crate::filters::FilterConfig;
    use crate::groups::ReadGroup;
    use crate::lib_spec::FlankingSequences;
    use crate::parsing::ReadPairProducer;
    use crate::region::RegionCompleteness;
    use bio::alignment::AlignmentOperation;
    use bio::io::fastq;
    use regex::Regex;

    struct MockProducer {
        items: std::vec::IntoIter<Result<ReadPair, ReadPairError>>,
        has_reverse: bool,
        group: Option<Regex>,
        max_reads: u64,
        read_count: u64,
    }

    impl MockProducer {
        fn new(items: Vec<Result<ReadPair, ReadPairError>>, has_reverse: bool) -> Self {
            Self {
                items: items.into_iter(),
                has_reverse,
                group: None,
                max_reads: 0,
                read_count: 0,
            }
        }
    }

    impl Iterator for MockProducer {
        type Item = Result<ReadPair, ReadPairError>;

        fn next(&mut self) -> Option<Self::Item> {
            let next = self.items.next();
            if next.is_some() {
                self.read_count += 1;
            }
            next
        }
    }

    impl ReadPairProducer for MockProducer {
        fn has_reverse(&self) -> bool {
            self.has_reverse
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

    fn make_record(id: &str, seq: &[u8], qual: &[u8]) -> fastq::Record {
        fastq::Record::with_attrs(id, None, seq, qual)
    }

    fn make_readpair(
        f_seq: &[u8],
        f_qual: &[u8],
        r_seq: Option<&[u8]>,
        r_qual: Option<&[u8]>,
    ) -> ReadPair {
        ReadPair {
            forward: make_record("f", f_seq, f_qual),
            reverse: r_seq.map(|seq| make_record("r", seq, r_qual.expect("reverse qual required"))),
            group: ReadGroup::ungrouped(),
        }
    }

    #[test]
    fn alignment_scorer_scores_standard_n_and_mismatch_cases() {
        let scorer = AlignmentScorer::new(6, -2, -3, -10, -4);

        assert_eq!(scorer.score(b'A', b'A'), 6);
        assert_eq!(scorer.score(b'N', b'A'), -2);
        assert_eq!(scorer.score(b'A', b'N'), -2);
        assert_eq!(scorer.score(b'A', b'T'), -3);
    }

    #[test]
    fn regions_from_alignment_path_maps_complete_region() {
        let region_positions = vec![(1, 5)];
        let path = vec![
            (1, 1, AlignmentOperation::Match),
            (2, 2, AlignmentOperation::Match),
            (3, 3, AlignmentOperation::Match),
            (4, 4, AlignmentOperation::Match),
            (5, 5, AlignmentOperation::Match),
        ];

        let out = regions_from_alignment_path(&region_positions, &path).unwrap();
        assert_eq!(out.len(), 1);
        assert_eq!(out[0], Some((1, 5, RegionCompleteness::Complete)));
    }

    #[test]
    fn regions_from_alignment_path_marks_partial_5prime_when_alignment_starts_inside_region() {
        let region_positions = vec![(2, 5)];
        let path = vec![
            (1, 3, AlignmentOperation::Match),
            (2, 4, AlignmentOperation::Match),
            (3, 5, AlignmentOperation::Match),
            (4, 6, AlignmentOperation::Match),
        ];

        let out = regions_from_alignment_path(&region_positions, &path).unwrap();
        assert_eq!(out[0], Some((1, 3, RegionCompleteness::Partial5Prime)));
    }

    #[test]
    fn regions_from_alignment_path_leaves_fully_deleted_region_unmapped() {
        let region_positions = vec![(1, 3)];
        let path = vec![
            (1, 1, AlignmentOperation::Del),
            (1, 2, AlignmentOperation::Del),
            (1, 3, AlignmentOperation::Match),
            (2, 4, AlignmentOperation::Match),
        ];

        let out = regions_from_alignment_path(&region_positions, &path).unwrap();
        assert_eq!(out, vec![None]);
    }

    #[test]
    fn merge_seqs_returns_present_side_when_other_missing() {
        let fwd = Some((
            b"ACGT".to_vec(),
            b"IIII".to_vec(),
            RegionCompleteness::Complete,
        ));
        let out = merge_seqs(fwd.clone(), None, 4).unwrap();
        assert_eq!(out, Some((b"ACGT".to_vec(), RegionCompleteness::Complete)));

        let rev = Some((
            b"TGCA".to_vec(),
            b"####".to_vec(),
            RegionCompleteness::Partial5Prime,
        ));
        let out = merge_seqs(None, rev, 4).unwrap();
        assert_eq!(
            out,
            Some((b"TGCA".to_vec(), RegionCompleteness::Partial5Prime))
        );
    }

    #[test]
    fn merge_seqs_prefers_higher_quality_for_duplicate_complete_observations() {
        let low = Some((
            b"AAAA".to_vec(),
            b"!!!!".to_vec(),
            RegionCompleteness::Complete,
        ));
        let high = Some((
            b"TTTT".to_vec(),
            b"IIII".to_vec(),
            RegionCompleteness::Complete,
        ));

        let out = merge_seqs(low, high, 4).unwrap();
        assert_eq!(out, Some((b"TTTT".to_vec(), RegionCompleteness::Complete)));
    }

    #[test]
    fn merge_seqs_combines_non_overlapping_partials_into_missing_center() {
        let fwd = Some((
            b"AAA".to_vec(),
            b"III".to_vec(),
            RegionCompleteness::Partial3Prime,
        ));
        let rev = Some((
            b"TT".to_vec(),
            b"II".to_vec(),
            RegionCompleteness::Partial5Prime,
        ));

        let out = merge_seqs(fwd, rev, 6).unwrap();
        assert_eq!(
            out,
            Some((
                b"AAA/TT".to_vec(),
                RegionCompleteness::MissingCenter { split_ind: 3 },
            ))
        );
    }

    #[test]
    fn merge_seqs_combines_overlapping_partials_into_overlapping() {
        let fwd = Some((
            b"AAAA".to_vec(),
            b"IIII".to_vec(),
            RegionCompleteness::Partial3Prime,
        ));
        let rev = Some((
            b"TTTT".to_vec(),
            b"IIII".to_vec(),
            RegionCompleteness::Partial5Prime,
        ));

        let out = merge_seqs(fwd, rev, 6).unwrap();
        assert_eq!(
            out,
            Some((
                b"AAAA/TTTT".to_vec(),
                RegionCompleteness::Overlapping { split_ind: 4 },
            ))
        );
    }

    #[test]
    fn merge_seqs_rejects_already_merged_inputs() {
        let bad1 = Some((
            b"AA/TT".to_vec(),
            b"IIIII".to_vec(),
            RegionCompleteness::MissingCenter { split_ind: 2 },
        ));

        let bad2 = Some((
            b"AA/TT".to_vec(),
            b"IIIII".to_vec(),
            RegionCompleteness::MissingCenter { split_ind: 2 },
        ));

        let err = merge_seqs(bad1, bad2, 4).unwrap_err().to_string();
        assert!(err.contains("already merged") || err.contains("MissingCenter"));
    }

    #[test]
    fn count_reads_rejects_zero_threads() {
        let reads = MockProducer::new(vec![], false);
        let err = count_reads(
            reads,
            &None,
            CountMode::FullRead,
            false,
            FilterConfig::new(None, None, None, None, true),
            None,
            None,
            None,
            false,
            0,
            None,
        )
        .unwrap_err()
        .to_string();

        assert!(err.contains("Threads must be >0"));
    }

    #[test]
    fn count_reads_align_mode_requires_alignment_scorer() {
        let reads = MockProducer::new(vec![], false);
        let err = count_reads(
            reads,
            &None,
            CountMode::Align,
            false,
            FilterConfig::new(None, None, None, None, true),
            None,
            None,
            None,
            false,
            1,
            None,
        )
        .unwrap_err()
        .to_string();

        assert!(err.contains("no AlignmentScorer passed"));
    }

    #[test]
    fn count_reads_pattern_mode_requires_length_and_tolerance() {
        let reads = MockProducer::new(vec![], false);
        let err = count_reads(
            reads,
            &None,
            CountMode::Pattern,
            false,
            FilterConfig::new(None, None, None, None, true),
            None,
            Some(10),
            None,
            false,
            1,
            None,
        )
        .unwrap_err()
        .to_string();

        assert!(err.contains("pattern length and/or tolerance is missing"));
    }

    #[test]
    fn count_reads_full_read_mode_counts_unfiltered_reads_without_libspec() {
        let read = make_readpair(b"ACGT", b"IIII", None, None);
        let reads = MockProducer::new(vec![Ok(read)], false);

        let counts = count_reads(
            reads,
            &None,
            CountMode::FullRead,
            false,
            FilterConfig::new(None, None, None, None, true),
            None,
            None,
            None,
            false,
            1,
            None,
        )
        .unwrap();

        assert_eq!(counts.len(), 1);
        assert_eq!(counts.total_filtered(), 0);
        assert_eq!(counts.to_vector(false)[0].total_count(), 1);
    }

    #[test]
    fn test_perfect_flank_matching() {
        let _ = env_logger::try_init();

        let seq = b"CCCCAATTGGGCCGGAAAAGGCCGGTATAGGGGATATGGGCGCGTTTT";
        let qual = b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
        let flanks = vec![
            FlankingSequences::OpenStart(b"AATT".to_vec()),
            FlankingSequences::Internal(b"CCGG".to_vec(), b"GGCC".to_vec()),
            FlankingSequences::Internal(b"TATA".to_vec(), b"ATAT".to_vec()),
            FlankingSequences::OpenEnd(b"CGCG".to_vec()),
        ];
        let tolerance: u64 = 0;

        let exp: Vec<Option<RegionMatch>> = vec![
            Some((
                b"CCCC".to_vec(),
                b"FFFF".to_vec(),
                RegionCompleteness::Partial5Prime,
            )),
            Some((
                b"AAAA".to_vec(),
                b"FFFF".to_vec(),
                RegionCompleteness::Complete,
            )),
            Some((
                b"GGGG".to_vec(),
                b"FFFF".to_vec(),
                RegionCompleteness::Complete,
            )),
            Some((
                b"TTTT".to_vec(),
                b"FFFF".to_vec(),
                RegionCompleteness::Partial3Prime,
            )),
        ];

        if let Ok(obs) = match_flank_patterns(seq, qual, &flanks, tolerance) {
            assert_eq!(obs, exp, "Observed regions don't match expected");
        } else {
            assert!(false, "match_flank_patterns returned Err(...)")
        }
    }
    #[test]
    fn test_flank_matching_with_mismatches() {
        let _ = env_logger::try_init();

        let seq = b"CCCCAATCTAGGCGGAAAAGGCCGGTATAGGGGATATAAACGTTTTTT";
        let qual = b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
        let flanks = vec![
            FlankingSequences::OpenStart(b"AATT".to_vec()), // 1 mismatch: AATC
            FlankingSequences::Internal(b"CCGG".to_vec(), b"GGCC".to_vec()), // 1 mismatch: GCGG
            FlankingSequences::Internal(b"TATA".to_vec(), b"ATAT".to_vec()),
            FlankingSequences::OpenEnd(b"CGCG".to_vec()), // 2 mismatches: CGTT
        ];
        let tolerance: u64 = 1;

        let exp: Vec<Option<RegionMatch>> = vec![
            Some((
                b"CCCC".to_vec(),
                b"FFFF".to_vec(),
                RegionCompleteness::Partial5Prime,
            )),
            Some((
                b"AAAA".to_vec(),
                b"FFFF".to_vec(),
                RegionCompleteness::Complete,
            )),
            Some((
                b"GGGG".to_vec(),
                b"FFFF".to_vec(),
                RegionCompleteness::Complete,
            )),
            None,
        ];

        let obs =
            match_flank_patterns(seq, qual, &flanks, tolerance).expect("Pattern match failed");
        assert_eq!(
            obs, exp,
            "Observed regions don't match expected with mismatch tolerance"
        );
    }

    #[test]
    fn test_flank_matching_partial_path() {
        let seq = b"GGGACCGGAAAAGGCCGGTATAGGGG"; // Starts at region 2
        let qual = b"FFFFFFFFFFFFFFFFFFFFFFFFFF";
        let flanks = vec![
            FlankingSequences::OpenStart(b"AATT".to_vec()), // missing
            FlankingSequences::Internal(b"CCGG".to_vec(), b"GGCC".to_vec()),
            FlankingSequences::Internal(b"TATA".to_vec(), b"ATAT".to_vec()),
            FlankingSequences::OpenEnd(b"CGCG".to_vec()), // missing
        ];
        let tolerance: u64 = 0;

        let exp: Vec<Option<RegionMatch>> = vec![
            None,
            Some((
                b"AAAA".to_vec(),
                b"FFFF".to_vec(),
                RegionCompleteness::Complete,
            )),
            Some((
                b"GGGG".to_vec(),
                b"FFFF".to_vec(),
                RegionCompleteness::Partial3Prime,
            )),
            None,
        ];

        let obs =
            match_flank_patterns(seq, qual, &flanks, tolerance).expect("Pattern match failed");
        assert_eq!(obs, exp);
    }

    #[test]
    fn test_flank_matching_with_gap_stops_scan() {
        let seq = b"GGGGAATTGGGGGGGGCGCG"; // Region 2 is missing
        let qual = b"FFFFFFFFFFFFFFFFFFFF";
        let flanks = vec![
            FlankingSequences::OpenStart(b"AATT".to_vec()),
            FlankingSequences::Internal(b"CCGG".to_vec(), b"GGCC".to_vec()), // missing
            FlankingSequences::Internal(b"TATA".to_vec(), b"ATAT".to_vec()), // would match if not skipped
            FlankingSequences::OpenEnd(b"CGCG".to_vec()),
        ];
        let tolerance: u64 = 0;

        let exp: Vec<Option<RegionMatch>> = vec![
            Some((
                b"GGGG".to_vec(),
                b"FFFF".to_vec(),
                RegionCompleteness::Partial5Prime,
            )),
            None,
            None,
            None,
        ];

        let obs =
            match_flank_patterns(seq, qual, &flanks, tolerance).expect("Pattern match failed");
        assert_eq!(obs, exp);
    }

    #[test]
    fn test_empty_flank_list() {
        let seq = b"ACGTACGT";
        let qual = b"FFFFFFFF";
        let flanks: Vec<FlankingSequences> = vec![];
        let tolerance = 0;

        let obs =
            match_flank_patterns(seq, qual, &flanks, tolerance).expect("Pattern match failed");
        assert!(obs.is_empty());
    }

    #[test]
    fn test_flank_matching_opening_scan_can_start_from_internal_close() {
        let seq = b"AAAAGGCCGGTATAGGGGATATCGCGTTTT";
        let qual = b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";

        let flanks = vec![
            FlankingSequences::OpenStart(b"AATT".to_vec()), // missing
            FlankingSequences::Internal(b"CCGG".to_vec(), b"GGCC".to_vec()), // start missing, close present
            FlankingSequences::Internal(b"TATA".to_vec(), b"ATAT".to_vec()),
            FlankingSequences::OpenEnd(b"CGCG".to_vec()),
        ];

        let exp: Vec<Option<RegionMatch>> = vec![
            None,
            Some((
                b"AAAA".to_vec(),
                b"FFFF".to_vec(),
                RegionCompleteness::Partial5Prime,
            )),
            Some((
                b"GGGG".to_vec(),
                b"FFFF".to_vec(),
                RegionCompleteness::Complete,
            )),
            Some((
                b"TTTT".to_vec(),
                b"FFFF".to_vec(),
                RegionCompleteness::Partial3Prime,
            )),
        ];

        let obs = match_flank_patterns(seq, qual, &flanks, 0).expect("Pattern match failed");
        assert_eq!(obs, exp);
    }
}
