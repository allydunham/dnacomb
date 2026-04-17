//! Observed sequence regions from a sequencing experiment
//!
//! Structures and functions to store and manipulate library reqions extracted
//! from sequencing data, as defined by a LibSpec.
use bio::bio_types::sequence::Sequence;
use itertools::Itertools;
use std::sync::Arc;

use crate::interning::{RegionID, SeqHandle, seq_from_bytes, seq_to_bytes};
use crate::library::{DistanceMetric, Library, LibraryRegion, PartialMatching, merge_matches};
use crate::seq_diff::SequenceDiff;

/// Key identifying an observed Region
///
/// Contains a subset of the region information to be used as a hash key
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

/// Observed sequence for a given library region
///
/// Each time this sequence is observed for this region it will be linked back to this object to
/// save library comparison overhead. The struct includes the observed sequence, a marker for
/// whether the whole region has necessarily been observed (i.e. if it is at the the start or end
///  of a read) and a RegionMatch object detailing the library match
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
    /// Create a new ObservedRegion
    pub fn new(id: RegionID, seq: &[u8], complete: RegionCompleteness) -> Self {
        Self {
            id,
            seq: seq_from_bytes(seq),
            completeness: complete,
            nearest_matches: RegionMatch::Uncompared,
        }
    }

    pub fn len(&self) -> usize {
        seq_to_bytes(self.seq).len()
    }

    pub fn is_empty(&self) -> bool {
        seq_to_bytes(self.seq).is_empty()
    }

    /// Check if library comparison has been performed
    pub fn is_compared_to_library(&self) -> bool {
        !matches!(self.nearest_matches, RegionMatch::Uncompared)
    }

    /// Compare the region sequence to a library of oligo options, dispatching to the appropriate implementation
    ///
    /// The distance to each is calculated and all nearest match indeces below the threshold
    /// supplied to library are returned.
    /// If no applicable matches are found nearest_matches will be an empty.
    /// If over max_matches hits are returned the region will be considered overmatched and
    /// indeterminant.
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

    /// Generate display string representing the region
    ///
    /// Output a tuple of values that can be used to represent the region (may be replaced
    /// a struct in future). The values are:
    ///
    /// * Sequence (with ^ at start or end to represent incompleteness)
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
        seq_match: Arc<LibraryRegion>,
        distance: u64,
        diff: SequenceDiff,
    },

    /// Multiple equidistant matches and the distance
    MultiMatch {
        seq_matches: Vec<Arc<LibraryRegion>>,
        distance: u64,
        diffs: Vec<SequenceDiff>,
    },

    /// Too many matches
    Overmatched { distance: u64, matches: usize },

    /// No match found
    Unmatched,

    /// Region not in library, with option to store the observed sequence
    NoLibrary { seq: Option<SeqHandle> },
}

impl RegionMatch {
    /// Get library Sequence(s), difference(s), distance, number of matches as strings for output.
    ///
    /// NoLibrary matches are considered not to have a matching sequence
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
