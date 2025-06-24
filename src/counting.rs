//! Counting the occurance of different reads in sequence files
//!
//! Contains methods for counting combinations of expected
//! regions in DNA sequence input.
//! Supports multiple approaches for extracting regions of interest
//! from the input sequence: alignment, pattern matching, inframe
//! position matching and full read counting.
use anyhow::{self};
use bio::alignment::AlignmentOperation;
use bio::alignment::pairwise::{Aligner, MatchFunc, Scoring};
use bio::alphabets::dna::revcomp;
use bio::bio_types::sequence::Sequence;
use clap::ValueEnum;
use itertools::{Itertools, izip};
use log::{debug, info};
use std::collections::HashMap;
use std::iter::zip;
use std::sync::Once;

use crate::containers::{
    CombinationKey, ObservedCombination, ObservedCombinations, ObservedRegion, RegionCompleteness,
    RegionKey,
};
use crate::errors::{AlignmentInfo, LibSpecError, ReadCountError, seq_to_string_or_log};
use crate::filters::{FilterConfig, FilterReason, mean_quality};
use crate::lib_spec::LibrarySpec;
use crate::logging::{Progress, ProgressStyle};
use crate::parsing::{ReadKey, ReadPair, ReadPairParser};

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
            Some(rev) => {
                ObservedRegion::key("rev".to_string(), rev.seq(), RegionCompleteness::Complete)
            }
            None => {
                return Err(ReadCountError::Error {
                    desc: "No reverse read found in paired raw mode".to_string(),
                }
                .into());
            }
        };

        let comb_key: CombinationKey = ObservedCombination::key(vec![fwd_key, rev_key]);

        counts.add_or_increment_combination(&comb_key, record.group)?;
        progress.inc(1);
    }
    progress.finish();

    Ok(counts)
}
