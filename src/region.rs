//! Observed sequence regions from a sequencing experiment
//!
//! Structures and functions to store and manipulate library reqions extracted
//! from sequencing data, as defined by a LibSpec.
use bio::bio_types::sequence::Sequence;
use std::sync::Arc;

use crate::errors::seq_to_string_or_log;
use crate::library::{DistanceMetric, Library, LibraryRegion, PartialMatching, merge_matches};

/// Key identifying an observed Region
///
/// Contains a subset of the region information to be used as a hash key
#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub struct RegionKey {
    pub id: String,
    pub sequence: Sequence,
    pub completeness: RegionCompleteness,
}

impl RegionKey {
    pub fn new(id: String, sequence: Sequence, completeness: RegionCompleteness) -> Self {
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
    pub id: String,

    /// Observed Sequence
    pub seq: Sequence,

    /// Completeness status of the region
    pub completeness: RegionCompleteness,

    /// Library indeces of nearest matches (Vec of Vecs containing indeces from each matching seq)
    pub nearest_matches: RegionMatch,
}

impl ObservedRegion {
    /// Create a new ObservedRegion
    pub fn new(id: String, seq: &[u8], complete: RegionCompleteness) -> Self {
        Self {
            id,
            seq: seq.to_vec(),
            completeness: complete,
            nearest_matches: RegionMatch::Uncompared,
        }
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
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
    pub fn to_tsv_chunk(&self) -> String {
        // Generate "{region} {region}_nearest {region}_distance {region}_n_matches" String
        let mut seq = match String::from_utf8(self.seq.clone()) {
            Ok(x) => x,
            Err(_) => panic!("Error converting Vec<u8> sequence to string"),
        };

        // Add appropriate marker to incomplete regions
        match self.completeness {
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
    },

    /// Multiple equidistant matches and the distance
    MultiMatch {
        seq_matches: Vec<Arc<LibraryRegion>>,
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
            } => {
                let mut seqs = seq_matches
                    .iter()
                    .map(|x| seq_to_string_or_log(&x.sequence))
                    .collect::<Vec<_>>();
                seqs.sort();
                (seqs.join(","), distance.to_string(), seqs.len().to_string())
            }
        }
    }

    /// Get matching Sequence as a string, including the passed through raw
    /// Sequence for NoLibrary matches.
    pub fn str_sequence(&self) -> String {
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

#[cfg(test)]
mod tests {
    use super::*;

    fn make_region(id: &str, seq: &[u8], c: RegionCompleteness) -> ObservedRegion {
        ObservedRegion::new(id.to_string(), seq, c)
    }

    /// Creating a Complete region should preserve id, bytes, length, and completeness.
    #[test]
    fn new_complete_region_holds_data() {
        let r = make_region("barcode", b"ACGTACGT", RegionCompleteness::Complete);

        assert_eq!(r.id, "barcode");
        assert_eq!(r.seq, b"ACGTACGT");
        assert_eq!(r.len(), 8);
        assert!(matches!(r.completeness, RegionCompleteness::Complete));
    }

    /// Empty sequences are allowed and correctly reported.
    #[test]
    fn empty_sequence_is_valid() {
        let r = make_region("empty", b"", RegionCompleteness::Complete);
        assert_eq!(r.id, "empty");
        assert_eq!(r.len(), 0);
        assert!(r.is_empty());
        assert_eq!(r.seq, b"");
    }

    /// Byte content must be preserved exactly; no implicit normalisation should occur.
    #[test]
    fn preserves_bytes_verbatim() {
        let weird = b"ACGTNN--acgt\x00\xff";
        let r = make_region("weird", weird, RegionCompleteness::Complete);
        assert_eq!(r.seq, weird, "region bytes changed unexpectedly");
        assert_eq!(r.len(), weird.len());
    }

    /// Regions should not alias their input slice; subsequent mutations to the source
    /// buffer (if any) must not affect the stored sequence.
    #[test]
    fn owns_its_sequence() {
        let mut buf = b"AAAA".to_vec();
        let r = make_region("id", &buf, RegionCompleteness::Complete);
        buf[0] = b'T'; // mutate the source buffer
        assert_eq!(r.seq, b"AAAA", "region leaked aliasing to input slice");
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
        assert_eq!(r2.seq[0], b'G');
        assert_eq!(r2.seq[99_999], b'G');
    }
}
