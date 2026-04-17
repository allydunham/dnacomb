//! Representation of individual observed region combinations.
//!
//! This module defines [`ObservedCombination`] and related types describing a
//! single distinct combination of observed variable regions extracted from reads.
//! It also includes the logic for combining per-region library matches into an
//! overall combination-level assignment relative to an expected library.
use std::collections::{HashMap, HashSet};
use std::sync::{Arc, Mutex};

use crate::errors::LibraryError;
use crate::groups::ReadGroup;
use crate::interning::{LibraryID, RegionID, library_id_to_str};
use crate::library::{DistanceMetric, Library};
use crate::region::{ObservedRegion, RegionKey, RegionMatch};
use crate::seqs::SeqPair;

/// Key identifying an observed combination
///
/// Contains a subset of the combination information to use as a hash key
#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub struct CombinationKey {
    pub sequence: Option<SeqPair>,
    pub regions: Vec<RegionKey>,
}

impl CombinationKey {
    pub fn new(sequence: Option<SeqPair>, regions: Vec<RegionKey>) -> Self {
        Self { sequence, regions }
    }
}

/// Combination of observed variable regions seen in one or more reads.
///
/// An `ObservedCombination` represents one distinct combination of ObservedRegions.
/// It stores:
/// - the observed region for each variable region,
/// - optional full read sequence
/// - counts split by read group,
/// - the inferred relationship between this combination and
///   the expected library (after comparison is run)
#[derive(Debug, Clone)]
pub struct ObservedCombination {
    /// Count of observations for each read group. Ungrouped reads are stored in None
    pub counts: HashMap<ReadGroup, u32>,

    /// Full sequence
    pub sequence: Option<SeqPair>,

    /// ObservedRegions defining the sequence form. References to ObservedRegion which
    /// should be stored in the parent ObservedCombinations object.
    pub regions: HashMap<RegionID, Arc<Mutex<ObservedRegion>>>,

    /// Status and result of comparison with the expected library of sequences
    pub library_matches: CombinationMatch,
}

impl ObservedCombination {
    pub fn new(
        regions: HashMap<RegionID, Arc<Mutex<ObservedRegion>>>,
        sequence: Option<SeqPair>,
    ) -> Self {
        Self {
            counts: HashMap::new(),
            sequence,
            regions,
            library_matches: CombinationMatch::Uncompared,
        }
    }

    /// Total count across all read groups
    pub fn total_count(&self) -> u32 {
        self.counts.values().sum()
    }

    /// Incremment count for a read group
    pub fn increment_count(&mut self, group: ReadGroup) {
        match self.counts.get_mut(&group) {
            Some(x) => *x += 1,
            None => {
                self.counts.insert(group, 1);
            }
        }
    }

    /// Compare this observed combination to the expected library.
    ///
    /// Each observed region is first compared to its corresponding library region
    /// if that comparison has not already been performed. The per-region matches
    /// are then combined across all expected regions to determine whether the
    /// overall observed combination is:
    ///
    /// - a unique library match,
    /// - a multimatch to several equally good library combinations,
    /// - a recombination of valid library elements in an unexpected combination,
    /// - a mismatch because at least one observed region could not be assigned,
    /// - or a nonmatch because one or more expected regions are missing.
    ///
    /// Regions that are not present in the supplied library are ignored when
    /// determining the overall combination match.
    pub fn compare_to_library(
        &self,
        region_ids: &Vec<RegionID>,
        library: &Library,
        distance_metric: DistanceMetric,
        max_matches: usize,
    ) -> CombinationMatch {
        let mut comb_dist: u64 = 0;
        let mut candidate_matches: Vec<Option<HashSet<LibraryID>>> =
            vec![None; library.sublibraries.len()];

        for reg_id in region_ids {
            let reg = match self.regions.get(reg_id) {
                None => return CombinationMatch::Nonmatch,
                Some(x) => {
                    // Should never need this with the implementation in ObservedCombinations, but here as a back-up as otherwise could panic later. Do all regions first as slightly more efficient and easier to follow in log
                    let mut reg = x.lock().unwrap();
                    if !reg.is_compared_to_library() {
                        let val = reg.compare_to_library(library, distance_metric, max_matches);
                        reg.nearest_matches = val;
                    }
                    reg
                }
            };

            let sublib = match library.get_sublibrary_index(reg_id) {
                Ok(x) => x,
                Err(_) => continue, // If region doesn't match a sublib can skip it,
            };

            match &reg.nearest_matches {
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
                    ..
                } => {
                    comb_dist += distance;
                    match candidate_matches[sublib] {
                        // If this is the first region, candidates are all it's inds
                        None => candidate_matches[sublib] = Some(seq_match.ids.clone()),

                        // Otherwise remove any inds it doesn't overlap to narrow down
                        Some(ref mut x) => {
                            x.retain(|i| seq_match.ids.contains(i));
                        }
                    };
                }
                RegionMatch::MultiMatch {
                    seq_matches,
                    distance,
                    ..
                } => {
                    comb_dist += distance;

                    // Identify the union of inds the multimatch covers
                    let mut match_ind_union: HashSet<LibraryID> = HashSet::new();
                    for mat in seq_matches {
                        match_ind_union.extend(mat.ids.iter());
                    }

                    // Set as the search space or remove anything not overlapping it
                    match candidate_matches[sublib] {
                        None => {
                            candidate_matches[sublib] = Some(match_ind_union);
                        }
                        Some(ref mut x) => {
                            x.retain(|i| match_ind_union.contains(i));
                        }
                    }
                }
            }
        }

        // If all sublibs haven't matched (e.g. empty read) is a nonmatch
        if candidate_matches.iter().all(|x| x.is_none()) {
            return CombinationMatch::Nonmatch;
        }

        // If any of the matched libraries have no remaining indeces must be a recombination
        if candidate_matches.iter().any(|x| match x {
            Some(x) => x.is_empty(),
            None => false,
        }) {
            return CombinationMatch::Recombination {
                distance: comb_dist,
            };
        }

        // If all None (e.g. no regions from that lib) or length 1 then a full match
        if candidate_matches.iter().all(|x| match x {
            Some(x) => x.len() == 1,
            None => true,
        }) {
            return CombinationMatch::Match {
                // Can unwrap because we know x.len() == 1
                inds: candidate_matches
                    .iter()
                    .map(|x| x.as_ref().map(|x| *x.iter().next().unwrap()))
                    .collect(),
                distance: comb_dist,
            };
        }

        // Only remaining case is a multi-match with a mix of match lengths
        CombinationMatch::MultiMatch {
            inds: candidate_matches,
            distance: comb_dist,
        }
    }
}

/// Status of the match between ObservedCombination and a Library
///
/// This enum summarises the result of combining per-region library matches
/// across all expected variable regions.
/// It includes a status plus distance and the match(s) as potential payloads.
#[derive(Debug, Clone)]
pub enum CombinationMatch {
    /// Library comparison has not been performed.
    Uncompared,

    /// All relevant regions resolve to a single consistent library assignment.
    ///
    /// For multi-sublibrary designs, `inds` contains one entry per sublibrary.
    /// `None` means that this combination did not include any region from that
    /// sublibrary.
    Match {
        inds: Vec<Option<LibraryID>>,
        distance: u64,
    },

    /// The observed combination is consistent with multiple equally good library
    /// assignments in at least one sublibrary.
    ///
    /// For multi-sublibrary designs, `inds` contains one entry per sublibrary.
    /// `None` means that this combination did not include any region from that
    /// sublibrary.
    MultiMatch {
        inds: Vec<Option<HashSet<LibraryID>>>,
        distance: u64,
    },

    /// The individual regions match library elements, but not in a valid expected
    /// combination.
    Recombination { distance: u64 },

    /// All regions exist but at least one could not be assigned to the library
    Mismatch,

    /// One or more expected regions were missing from the observed combination.
    Nonmatch,
}

impl CombinationMatch {
    /// Format the matched library IDs for display or TSV output.
    ///
    /// For a unique match this returns one ID per sublibrary, joined by `/`.
    /// For a multimatch, IDs within a sublibrary are joined by `,` and different
    /// sublibraries are joined by `/`.
    /// Returns an empty string for statuses without a meaningful match.
    pub fn id_string(&self) -> Result<String, LibraryError> {
        Ok(match self {
            CombinationMatch::Match { inds, .. } => inds
                .iter()
                .map(|id| match id {
                    Some(id) => Ok(library_id_to_str(*id).to_string()),
                    None => Ok("".to_string()),
                })
                .collect::<Result<Vec<String>, LibraryError>>()?
                .join("/"),

            CombinationMatch::MultiMatch { inds, .. } => {
                let mut names = Vec::with_capacity(inds.len());

                for ids in inds.iter() {
                    names.push(match ids {
                        Some(ids) => ids
                            .iter()
                            .map(|id| library_id_to_str(*id))
                            .collect::<Vec<_>>()
                            .join(","),
                        None => "".to_string(),
                    });
                }

                names.join("/")
            }
            _ => "".to_string(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use std::sync::{Arc, Mutex};

    use crate::SubLibrary;
    use crate::interning::{RegionID, region_id_from_str};
    use crate::region::RegionCompleteness;

    // Library members
    // seq1 - ATAT GGGG
    // seq2 - AAAA CCCC
    // seq3 - AAAT CCCC

    #[test]
    fn compare_to_library_table() {
        use RegionCompleteness::Complete;

        let cases: Vec<TestCase> = vec![
            // Exact matching
            TestCase {
                name: "exact_match",
                regions: vec![("r1", b"AAAA", Complete), ("r2", b"CCCC", Complete)],
                metric: DistanceMetric::Exact,
                expected: Expected::Match {
                    name: "seq2",
                    distance: 0,
                },
            },
            TestCase {
                name: "exact_mismatch",
                regions: vec![("r1", b"AAAA", Complete), ("r2", b"AAAA", Complete)],
                metric: DistanceMetric::Exact,
                expected: Expected::Mismatch,
            },
            TestCase {
                name: "exact_nonmatch",
                regions: vec![("r1", b"AAAA", Complete)],
                metric: DistanceMetric::Exact,
                expected: Expected::Nonmatch,
            },
            TestCase {
                name: "exact_recombination",
                regions: vec![("r1", b"AAAA", Complete), ("r2", b"GGGG", Complete)],
                metric: DistanceMetric::Exact,
                expected: Expected::Recombination { distance: 0 },
            },
            // Hamming lookup
            TestCase {
                name: "hamming_match",
                regions: vec![("r1", b"ATAT", Complete), ("r2", b"GGGC", Complete)],
                metric: DistanceMetric::Hamming,
                expected: Expected::Match {
                    name: "seq1",
                    distance: 1,
                },
            },
            TestCase {
                name: "hamming_multimatch",
                regions: vec![("r1", b"AAAG", Complete), ("r2", b"CCCC", Complete)],
                metric: DistanceMetric::Hamming,
                expected: Expected::MultiMatch {
                    inds_len: 2,
                    distance: 1,
                },
            },
            TestCase {
                name: "hamming_mismatch",
                regions: vec![("r1", b"TT", Complete), ("r2", b"CCCC", Complete)],
                metric: DistanceMetric::Hamming,
                expected: Expected::Mismatch,
            },
            TestCase {
                name: "hamming_nonmatch",
                regions: vec![("r1", b"AAAA", Complete)],
                metric: DistanceMetric::Hamming,
                expected: Expected::Nonmatch,
            },
            TestCase {
                name: "hamming_recombination",
                regions: vec![("r1", b"AAAA", Complete), ("r2", b"GGGC", Complete)],
                metric: DistanceMetric::Hamming,
                expected: Expected::Recombination { distance: 1 },
            },
            // Bounded-levenshtein lookup
            TestCase {
                name: "bounded_levenshtein_match",
                regions: vec![("r1", b"ATT", Complete), ("r2", b"GGGC", Complete)],
                metric: DistanceMetric::BoundedLevenshtein,
                expected: Expected::Match {
                    name: "seq1",
                    distance: 2,
                },
            },
            TestCase {
                name: "bounded_levenshtein_multimatch",
                regions: vec![("r1", b"AAAG", Complete), ("r2", b"CCCC", Complete)],
                metric: DistanceMetric::BoundedLevenshtein,
                expected: Expected::MultiMatch {
                    inds_len: 2,
                    distance: 1,
                },
            },
            TestCase {
                name: "bounded_levenshtein_mismatch",
                regions: vec![("r1", b"CG", Complete), ("r2", b"CCCC", Complete)],
                metric: DistanceMetric::BoundedLevenshtein,
                expected: Expected::Mismatch,
            },
            TestCase {
                name: "bounded_levenshtein_nonmatch",
                regions: vec![("r2", b"AAAA", Complete)],
                metric: DistanceMetric::BoundedLevenshtein,
                expected: Expected::Nonmatch,
            },
            TestCase {
                name: "bounded_levenshtein_recombination",
                regions: vec![("r1", b"AAAA", Complete), ("r2", b"GGG", Complete)],
                metric: DistanceMetric::BoundedLevenshtein,
                expected: Expected::Recombination { distance: 1 },
            },
            // Levenshtein lookup
            TestCase {
                name: "levenshtein_match",
                regions: vec![("r1", b"AAT", Complete), ("r2", b"CCTC", Complete)],
                metric: DistanceMetric::Levenshtein,
                expected: Expected::Match {
                    name: "seq3",
                    distance: 2,
                },
            },
            TestCase {
                name: "levenshtein_multimatch",
                regions: vec![("r1", b"AAAG", Complete), ("r2", b"CCCC", Complete)],
                metric: DistanceMetric::Levenshtein,
                expected: Expected::MultiMatch {
                    inds_len: 2,
                    distance: 1,
                },
            },
            TestCase {
                name: "levenshtein_mismatch",
                regions: vec![("r1", b"CG", Complete), ("r2", b"CCCC", Complete)],
                metric: DistanceMetric::Levenshtein,
                expected: Expected::Mismatch,
            },
            TestCase {
                name: "levenshtein_nonmatch",
                regions: vec![("r2", b"AAAA", Complete)],
                metric: DistanceMetric::Levenshtein,
                expected: Expected::Nonmatch,
            },
            TestCase {
                name: "levenshtein_recombination",
                regions: vec![("r1", b"AAAA", Complete), ("r2", b"GGG", Complete)],
                metric: DistanceMetric::Levenshtein,
                expected: Expected::Recombination { distance: 1 },
            },
        ];

        for tc in &cases {
            run_row(tc);
        }
    }

    /// Expected outcome for a test row
    #[derive(Debug)]
    enum Expected {
        Match { name: &'static str, distance: u64 },
        MultiMatch { inds_len: usize, distance: u64 },
        Recombination { distance: u64 },
        Mismatch,
        Nonmatch,
    }

    /// Case to test
    struct TestCase {
        name: &'static str,
        regions: Vec<(&'static str, &'static [u8], RegionCompleteness)>,
        metric: DistanceMetric,
        expected: Expected,
    }

    /// Construct a library to compare to
    fn make_library() -> Library {
        use std::collections::HashMap;

        let mut map: HashMap<RegionID, Vec<Vec<u8>>> = HashMap::new();
        let ids = Some(vec![
            "seq1".to_string(),
            "seq2".to_string(),
            "seq3".to_string(),
        ]);
        let region_max: HashMap<RegionID, u64> = HashMap::new();

        map.insert(
            region_id_from_str("r1"),
            vec![b"ATAT".to_vec(), b"AAAA".to_vec(), b"AAAT".to_vec()],
        );
        map.insert(
            region_id_from_str("r2"),
            vec![b"GGGG".to_vec(), b"CCCC".to_vec(), b"CCCC".to_vec()],
        );

        // default_max_distance is only relevant for edit-distance modes; 2 is fine.
        let sublib = SubLibrary::new(map, ids, region_max, 2, None);

        Library::new(vec![sublib.expect("Sublib failed")]).expect("library builds")
    }

    /// Generate observed combination
    fn make_combination(
        regs: &[(&'static str, &'static [u8], RegionCompleteness)],
    ) -> ObservedCombination {
        let mut map: HashMap<RegionID, Arc<Mutex<ObservedRegion>>> = HashMap::new();
        for (id, seq, comp) in regs.iter().copied() {
            map.insert(
                region_id_from_str(id),
                Arc::new(Mutex::new(ObservedRegion::new(
                    region_id_from_str(id),
                    seq,
                    comp,
                ))),
            );
        }
        ObservedCombination::new(map, None)
    }

    /// Run tests on a row
    fn run_row(tc: &TestCase) {
        let lib = make_library();
        let comb = make_combination(&tc.regions);
        let region_ids: Vec<RegionID> = vec![region_id_from_str("r1"), region_id_from_str("r2")];

        let got = comb.compare_to_library(&region_ids, &lib, tc.metric, 3);

        match (&tc.expected, got.clone()) {
            (Expected::Match { name, distance }, CombinationMatch::Match { distance: d, .. }) => {
                assert_eq!(d, *distance, "[{}] distance mismatch", tc.name);
                let got_name = got.id_string().expect("name exists");
                assert_eq!(got_name, *name, "[{}] matched name mismatch", tc.name);
            }
            (
                Expected::MultiMatch { inds_len, distance },
                CombinationMatch::MultiMatch { inds, distance: d },
            ) => {
                assert_eq!(d, *distance, "[{}] distance mismatch", tc.name);
                assert_eq!(
                    inds.iter()
                        .map(|x| match x {
                            Some(x) => x.len(),
                            None => 1,
                        })
                        .reduce(|x, y| x * y)
                        .unwrap_or(0),
                    *inds_len,
                    "[{}] candidate set size mismatch",
                    tc.name
                );
            }
            (
                Expected::Recombination { distance },
                CombinationMatch::Recombination { distance: d },
            ) => {
                assert_eq!(
                    d, *distance,
                    "[{}] recombination distance mismatch",
                    tc.name
                );
            }
            (Expected::Mismatch, CombinationMatch::Mismatch) => {}
            (Expected::Nonmatch, CombinationMatch::Nonmatch) => {}
            (exp, got_actual) => {
                panic!(
                    "[{}] unexpected result.\n  expected: {:?}\n  got: {:?}",
                    tc.name, exp, got_actual
                );
            }
        }
    }
}
