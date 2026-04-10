//! Counting the occurance of different reads in sequence files
//!
//! Contains methods for counting combinations of expected
//! regions in DNA sequence input.
//! Supports multiple approaches for extracting regions of interest
//! from the input sequence: alignment, pattern matching, inframe
//! position matching and full read counting.
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
use crate::lib_spec::{FlankingSequences, LibrarySpec};
use crate::logging::{Progress, ProgressStyle};
use crate::parsing::{ReadPairProducer, ThreadedReadPairParser};
use crate::region::{RegionCompleteness, RegionKey};
use crate::seqs::{ReadPair, SeqPair};
use crate::utils::mean_quality;

/// Position in an alignment where a region is found
///
/// Has the format (start, stop, RegionCompleteness status)
pub type AlignmentPosition = (usize, usize, RegionCompleteness);

/// Region sequence identified while searching an input sequence, containing
/// the found sequence, quality and whether it is complete. Convenience type
/// for counting functions that is quickly processed into region keys and the
/// proper combination structs.
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

/// Merge a forward and reverse sequence in a region
///
/// Looks at each position in turn and takes the highest quality option, if any.
/// Args are tuples with sequence, quality vector and a completeness tag, plus the
/// expected region length. This only works for fixed length, will need something
/// more complete if there is variable overlap (plus should really warn to use a read
/// merger at that point).
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

/// Find region matches by pattern
///
/// Search an input sequence for regions flanked by the input patterns, in the order
/// they occur. Not all regions need to be identified, but those found will be a
/// continuous subsequence. For instance, it may return hits for regions 2-4 but
/// be missing regions 1 and 5. Where the first/last flank sequence is None it is
/// considered to start/end at the sequence start/end. A mismatch tolerance allows
/// close flank matches to be considered hits to correct for errors, but this requires
/// care where multiple flank sequences have similar sequences.
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

    // Find opening region to assign start point
    let start_regs: Vec<Sequence> = flanks
        .iter()
        .map(|r| match r {
            FlankingSequences::Unflanked => unreachable!("Unflanked already checked"),
            FlankingSequences::OpenStart(end) => Ok(end.clone()),
            FlankingSequences::Internal(start, ..) => Ok(start.clone()),
            FlankingSequences::OpenEnd(start) => Ok(start.clone()),
        })
        .collect::<Result<Vec<Sequence>, ReadCountError>>()?;

    'outer: while pos < seq.len() {
        for (i, r) in start_regs.iter().enumerate() {
            end = pos + r.len();
            if end > seq.len() {
                continue;
            }

            dist = hamming(r, &seq[pos..end]);

            if dist <= tolerance {
                match flanks[i] {
                    FlankingSequences::Unflanked => unreachable!("Unflanked already checked"),
                    FlankingSequences::OpenStart(..) => {
                        out[i] = Some((
                            seq[0..pos].to_vec(),
                            qual[0..pos].to_vec(),
                            RegionCompleteness::Partial5Prime,
                        ));

                        reg = i + 1;
                        open = false;
                        pos = end;
                    }
                    FlankingSequences::Internal(..) => {
                        reg = i;
                        open = true;
                        reg_start = end;
                        pos = end;
                    }
                    FlankingSequences::OpenEnd(..) => {
                        out[i] = Some((
                            seq[end..seq.len()].to_vec(),
                            qual[end..seq.len()].to_vec(),
                            RegionCompleteness::Partial3Prime,
                        ));

                        // An open end region must be at the end (checked in validation)
                        // so directly return
                        return Ok(out);
                    }
                }
                break 'outer;
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

                    pos = end;
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

/// Count single end reads by aligning to the library template
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
            continue;
        }

        // Check if the read has been cached
        if cache && counts.check_cache(&record, true)?.is_some() {
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
                            *id,
                            // Offset seq lookup - rust vec 0 based and AlignmentPath 1 based
                            s.to_vec(),
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

/// Count paired end reads by aligning to the library template
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
            continue;
        }

        // Check if the read has been cached
        if cache && counts.check_cache(&record, true)?.is_some() {
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
                                *id,
                                // Offset seq lookup - rust vec 0 based and AlignmentPath 1 based
                                s.to_vec(),
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
                                *id,
                                // Offset seq lookup - rust vec 0 based and AlignmentPath 1 based
                                s.to_vec(),
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
                        Some((seq, comp)) => comb_key.regions.push(RegionKey::new(*id, seq, comp)),
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

/// Count single end reads using surrounding patterns from the library template
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
            continue;
        }

        // Check if the read has been cached
        if cache && counts.check_cache(&record, true)?.is_some() {
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
                    Some(r) => Some(RegionKey::new(*id, r.0, r.2)),
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

/// Count paired end reads using surrounding patterns from the library template
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
            continue;
        }

        // Check if the read has been cached
        if cache && counts.check_cache(&record, true)?.is_some() {
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
                comb_key
                    .regions
                    .push(RegionKey::new(*id, merged.0, merged.1));
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

/// Count single end reads based on their in-frame position in the template
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
            continue;
        }

        // Check if the read has been cached
        if cache && counts.check_cache(&record, true)?.is_some() {
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

            comb_key
                .regions
                .push(RegionKey::new(*id, reg_seq, complete));
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

/// Count paired end reads based on their in-frame position in the template
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
            continue;
        }

        // Check if the read has been cached
        if cache && counts.check_cache(&record, true)?.is_some() {
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
                Some((seq, comp)) => comb_key.regions.push(RegionKey::new(*id, seq, comp)),
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

/// Count entire single end reads
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
    use crate::lib_spec::FlankingSequences;
    use crate::region::RegionCompleteness;

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

        let seq = b"CCCCAATCGGGGCGGAAAAGGCCGGTATAGGGGATATAAACGTTTTTT";
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
        let seq = b"GGGCCGGAAAAGGCCGGTATAGGGG"; // Starts at region 2
        let qual = b"FFFFFFFFFFFFFFFFFFFFFFFFF";
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
}
