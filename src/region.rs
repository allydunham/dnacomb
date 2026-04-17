//! Observed variable-region sequences extracted from sequencing reads.
//!
//! This module defines the per-region objects used after read parsing and region
//! extraction but before full combination-level summarisation. It captures both
//! the observed sequence itself and how completely that region was observed,
//! along with optional comparison to an expected library.
use bio::bio_types::sequence::Sequence;
use itertools::Itertools;
use std::sync::Arc;

use crate::interning::{RegionID, SeqHandle, seq_from_bytes, seq_to_bytes};
use crate::library::{DistanceMetric, Library, LibraryRegion, PartialMatching, merge_matches};
use crate::seq_diff::SequenceDiff;

/// Key identifying one distinct observed region state.
///
/// Two observations are considered the same region key if they share:
/// - the same region ID,
/// - the same observed sequence,
/// - and the same completeness status.
///
/// This is used to deduplicate repeated region observations across reads.
#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub struct RegionKey {
    pub id: RegionID,
    pub sequence: Sequence,
    pub completeness: RegionCompleteness,
}

impl RegionKey {
    pub fn new(id: RegionID, sequence: Sequence, completeness: RegionCompleteness) -> Self {
        Self {
            id,
            sequence,
            completeness,
        }
    }
}

/// Observed sequence for one variable region.
///
/// An `ObservedRegion` stores:
/// - the region ID,
/// - the observed sequence,
/// - how completely that region was observed in the read(s),
/// - and, once performed, the result of comparison to the expected library.
///
/// Identical region observations can be shared across multiple observed
/// combinations to avoid repeating library-comparison work.
#[derive(Debug)]
pub struct ObservedRegion {
    /// Name of the region from LibSpec
    pub id: RegionID,

    /// Observed Sequence
    pub seq: SeqHandle,

    /// Completeness status of the region
    pub completeness: RegionCompleteness,

    /// Library indeces of nearest matches (Vec of Vecs containing indeces from each matching seq)
    pub nearest_matches: RegionMatch,
}

impl ObservedRegion {
    /// Create a new observed region in the uncompared state.
    pub fn new(id: RegionID, seq: &[u8], complete: RegionCompleteness) -> Self {
        Self {
            id,
            seq: seq_from_bytes(seq),
            completeness: complete,
            nearest_matches: RegionMatch::Uncompared,
        }
    }

    /// Length of the observed sequence.
    pub fn len(&self) -> usize {
        seq_to_bytes(self.seq).len()
    }

    /// Return `true` if the observed sequence is empty.
    pub fn is_empty(&self) -> bool {
        seq_to_bytes(self.seq).is_empty()
    }

    /// Return `true` if library comparison has already been performed for this region.
    pub fn is_compared_to_library(&self) -> bool {
        !matches!(self.nearest_matches, RegionMatch::Uncompared)
    }

    /// Compare this observed region to the expected library.
    ///
    /// The lookup strategy depends on `completeness`:
    /// - `Complete` regions are matched against the full library sequence,
    /// - `Partial5Prime` regions are matched against the library 3' end,
    /// - `Partial3Prime` regions are matched against the library 5' end,
    /// - `MissingCenter` and `Overlapping` regions are split into left/right
    ///   pieces and matched independently before intersecting compatible results.
    ///
    /// The result is returned as a `RegionMatch`, which may represent a unique
    /// match, multiple equally good matches, too many matches, no match, or the
    /// absence of this region from the supplied library.
    pub fn compare_to_library(
        &self,
        library: &Library,
        distance_metric: DistanceMetric,
        max_matches: usize,
    ) -> RegionMatch {
        let lib_match = match self.completeness {
            RegionCompleteness::Complete => {
                match library.lookup(&self.id, self.seq, distance_metric, PartialMatching::Full) {
                    // lookup can only return Err(LibraryError::MissingRegion) which implies the
                    // region isn't in the library. Needs changing if more errors are added to it.
                    Err(_) => return RegionMatch::NoLibrary { seq: None },
                    Ok(x) => x,
                }
            }
            RegionCompleteness::Partial5Prime => {
                match library.lookup(
                    &self.id,
                    self.seq,
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
                    self.seq,
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
                let seq = seq_to_bytes(self.seq);
                let left_seq = seq_from_bytes(&seq[0..split_ind]);
                let right_seq = seq_from_bytes(&seq[split_ind + 1..seq.len()]);

                let left_match = library.lookup(
                    &self.id,
                    left_seq,
                    distance_metric,
                    PartialMatching::FivePrimeOnly,
                );

                let right_match = library.lookup(
                    &self.id,
                    right_seq,
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
                        diff: SequenceDiff::compute_ids(self.seq, x.matches[0].sequence),
                    }
                } else if x.matches.len() > max_matches {
                    RegionMatch::Overmatched {
                        distance: x.distance,
                        matches: x.matches.len(),
                    }
                } else {
                    RegionMatch::MultiMatch {
                        distance: x.distance,
                        diffs: x
                            .matches
                            .iter()
                            .map(|m| SequenceDiff::compute_ids(self.seq, m.sequence))
                            .collect(),
                        seq_matches: x.matches,
                    }
                }
            }
        }
    }

    /// Convert this region into display strings for TSV/output generation.
    ///
    /// The returned tuple contains:
    /// * Sequence (with ^ at start or end to represent a partially truncated end)
    /// * nearest match(s)
    /// * difference to match(s)
    /// * distance from match
    /// * number of matches
    pub fn to_strings(&self) -> (String, String, String, String, String) {
        let seq = match self.completeness {
            RegionCompleteness::Complete
            | RegionCompleteness::MissingCenter { .. }
            | RegionCompleteness::Overlapping { .. } => self.seq.to_str_or_log(),
            RegionCompleteness::Partial5Prime => format!("^{}", self.seq.to_str_or_log()),
            RegionCompleteness::Partial3Prime => format!("{}^", self.seq.to_str_or_log()),
        };

        match &self.nearest_matches {
            RegionMatch::Uncompared | RegionMatch::Unmatched | RegionMatch::NoLibrary { .. } => (
                seq,
                "".to_string(),
                "".to_string(),
                "".to_string(),
                "0".to_string(),
            ),
            RegionMatch::Overmatched { distance, matches } => (
                seq,
                "".to_string(),
                "".to_string(),
                distance.to_string(),
                matches.to_string(),
            ),
            RegionMatch::Match {
                seq_match,
                distance,
                diff,
            } => (
                seq,
                seq_match.sequence.to_str_or_log(),
                diff.to_string(),
                distance.to_string(),
                "1".to_string(),
            ),
            RegionMatch::MultiMatch {
                seq_matches,
                distance,
                diffs,
            } => {
                let seqs: String = seq_matches
                    .iter()
                    .map(|x| x.sequence.to_str_or_log())
                    .join(",");

                let diff_str: String = diffs.iter().map(|x| x.to_string()).join(",");

                let count = seq_matches.len().to_string();

                (seq, seqs, diff_str, distance.to_string(), count)
            }
        }
    }
}

/// How completely a region was observed in the read data.
///
/// Completeness captures whether the full region sequence is known, whether one
/// end may be truncated by read boundaries, or whether paired-read evidence
/// produced gapped or overlapping partial observations.
///
/// For `MissingCenter` and `Overlapping`, the stored sequence is expected to
/// contain a separator byte at `split_ind`, with left and right observed
/// fragments stored on either side.
#[derive(Debug, Eq, Hash, PartialEq, Clone, Copy)]
pub enum RegionCompleteness {
    /// The full region sequence is believed to be observed.
    Complete,

    /// The 5' end may be truncated; the observed sequence begins inside the region.
    Partial5Prime,

    /// The 3' end may be truncated; the observed sequence ends inside the region.
    Partial3Prime,

    /// Forward/reverse evidence covers both ends of the region but leaves a gap
    /// in the middle. The stored sequence contains a `/` separator at `split_ind`.
    MissingCenter { split_ind: usize },

    /// Forward/reverse evidence covers both ends of the region and overlaps.
    /// The stored sequence contains a `/` separator at `split_ind`.
    Overlapping { split_ind: usize },
}

/// Status of the match between an ObservedRegion and a Library
///
/// This stores both the qualitative match status and any associated library
/// sequence(s), distance(s), and sequence-difference information.
#[derive(Debug, Eq, Hash, PartialEq, Clone)]
pub enum RegionMatch {
    /// Library comparison has not yet been performed.
    Uncompared,

    /// A unique best library-region match was found.
    Match {
        seq_match: Arc<LibraryRegion>,
        distance: u64,
        diff: SequenceDiff,
    },

    /// Multiple equally good best library-region matches were found.
    MultiMatch {
        seq_matches: Vec<Arc<LibraryRegion>>,
        distance: u64,
        diffs: Vec<SequenceDiff>,
    },

    /// More than `max_matches` equally good matches were found.
    Overmatched { distance: u64, matches: usize },

    /// No acceptable library-region match was found.
    Unmatched,

    /// This region is not represented in the supplied library.
    ///
    /// The observed sequence may optionally be retained for output.
    NoLibrary { seq: Option<SeqHandle> },
}

impl RegionMatch {
    /// Convert the library-match portion of this region comparison into strings
    /// for TSV/output generation.
    ///
    /// Returns:
    /// - matching sequence(s),
    /// - sequence difference(s),
    /// - match distance,
    /// - number of matches.
    pub fn to_strings(&self) -> (String, String, String, String) {
        match self {
            RegionMatch::Uncompared | RegionMatch::Unmatched | RegionMatch::NoLibrary { .. } => (
                "".to_string(),
                "".to_string(),
                "".to_string(),
                "0".to_string(),
            ),
            RegionMatch::Overmatched { distance, matches } => (
                "".to_string(),
                "".to_string(),
                distance.to_string(),
                matches.to_string(),
            ),
            RegionMatch::Match {
                seq_match,
                distance,
                diff,
            } => (
                seq_match.sequence.to_str_or_log(),
                diff.to_string(),
                distance.to_string(),
                "1".to_string(),
            ),
            RegionMatch::MultiMatch {
                seq_matches,
                distance,
                diffs,
            } => {
                let seqs: String = seq_matches
                    .iter()
                    .map(|x| x.sequence.to_str_or_log())
                    .join(",");

                let diff_str: String = diffs.iter().map(|x| x.to_string()).join(",");

                let count = seq_matches.len().to_string();

                (seqs, diff_str, distance.to_string(), count)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::interning::{region_id_from_str, region_id_to_str};

    use super::*;

    fn make_region(id: &str, seq: &[u8], c: RegionCompleteness) -> ObservedRegion {
        ObservedRegion::new(region_id_from_str(id), seq, c)
    }

    /// Creating a Complete region should preserve id, bytes, length, and completeness.
    #[test]
    fn new_complete_region_holds_data() {
        let r = make_region("barcode", b"ACGTACGT", RegionCompleteness::Complete);

        assert_eq!(region_id_to_str(r.id).to_string(), "barcode");
        assert_eq!(seq_to_bytes(r.seq).as_ref(), b"ACGTACGT");
        assert_eq!(r.len(), 8);
        assert!(matches!(r.completeness, RegionCompleteness::Complete));
    }

    /// Empty sequences are allowed and correctly reported.
    #[test]
    fn empty_sequence_is_valid() {
        let r = make_region("empty", b"", RegionCompleteness::Complete);
        assert_eq!(region_id_to_str(r.id).to_string(), "empty");
        assert_eq!(r.len(), 0);
        assert!(r.is_empty());
        assert_eq!(seq_to_bytes(r.seq).as_ref(), b"");
    }

    /// Byte content must be preserved exactly; no implicit normalisation should occur.
    #[test]
    fn preserves_bytes_verbatim() {
        let weird = b"ACGTNN--acgt\x00\xff";
        let r = make_region("weird", weird, RegionCompleteness::Complete);
        assert_eq!(
            seq_to_bytes(r.seq).as_ref(),
            weird,
            "region bytes changed unexpectedly"
        );
        assert_eq!(r.len(), weird.len());
    }

    /// Regions should not alias their input slice; subsequent mutations to the source
    /// buffer (if any) must not affect the stored sequence.
    #[test]
    fn owns_its_sequence() {
        let mut buf = b"AAAA".to_vec();
        let r = make_region("id", &buf, RegionCompleteness::Complete);
        buf[0] = b'T'; // mutate the source buffer
        assert_eq!(
            seq_to_bytes(r.seq).as_ref(),
            b"AAAA",
            "region leaked aliasing to input slice"
        );
    }

    /// Extremely short and extremely long sequences should not panic.
    #[test]
    fn size_extremes_smoke() {
        // Short (including 1-base)
        let r1 = make_region("s", b"A", RegionCompleteness::Complete);
        assert_eq!(r1.len(), 1);

        // Long (1e5 bytes) – keep modest to avoid slow tests; still catches realloc issues.
        let big = vec![b'G'; 100_000];
        let r2 = make_region("big", &big, RegionCompleteness::Complete);
        assert_eq!(r2.len(), 100_000);
        // spot-check end bytes to ensure contiguous storage
        assert_eq!(seq_to_bytes(r2.seq).as_ref()[0], b'G');
        assert_eq!(seq_to_bytes(r2.seq).as_ref()[99_999], b'G');
    }
}
