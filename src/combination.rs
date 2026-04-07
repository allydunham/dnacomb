//! Observed combination of sequence regions from a sequencing experiment
//!
//! Structures and functions to store and manipulate combinations of reqions extracted from sequencing data
use std::collections::{HashMap, HashSet};
use std::sync::{Arc, Mutex};

use crate::errors::{LibraryError, seq_to_string_or_log};
use crate::library::{DistanceMetric, Library};
use crate::region::{ObservedRegion, RegionKey, RegionMatch};
use crate::seqs::ReadGroup;
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

/// Combination of ObservedRegions seen in sequence reads
///
/// A set of observed regions determining the "type" of read, as defined in the LibSpec.
/// Also includes a count, the read grouping (for instance for different cells in single
/// cell studies) and whether it matches an expected library member.
#[derive(Debug, Clone)]
pub struct ObservedCombination {
    /// Count of observations for each read group. Ungrouped reads are stored in None
    pub counts: HashMap<ReadGroup, u32>,

    /// Full sequence
    pub sequence: Option<SeqPair>,

    /// ObservedRegions defining the sequence form. References to ObservedRegion which
    /// should be stored in the parent ObservedCombinations object.
    pub regions: HashMap<String, Arc<Mutex<ObservedRegion>>>,

    /// Status and result of comparison with the expected library of sequences
    pub library_matches: CombinationMatch,
}

impl ObservedCombination {
    pub fn new(
        regions: HashMap<String, Arc<Mutex<ObservedRegion>>>,
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

    /// Compare the combination to expected combinations in the library
    ///
    /// Returns a CombinationMatch object which can also be added to the ObservedCombination
    /// library matches field. Looks at each region in turn and identifies which library
    /// combinations are possible overall matches.
    pub fn compare_to_library(
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
                    let mut reg = x.lock().unwrap();
                    if !reg.is_compared_to_library() {
                        let val = reg.compare_to_library(library, distance_metric, max_matches);
                        reg.nearest_matches = val;
                    }
                    reg
                }
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
                } => {
                    comb_dist += distance;
                    match candidate_matches {
                        // If this is the first region, candidates are all it's inds
                        None => candidate_matches = Some(seq_match.inds.clone()),

                        // Otherwise remove any inds it doesn't overlap to narrow down
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

                    // Identify the union of inds the multimatch covers
                    let mut match_ind_union: HashSet<usize> = HashSet::new();
                    for mat in seq_matches {
                        match_ind_union.extend(mat.inds.iter());
                    }

                    // Set as the search space or remove anything not overlapping it
                    match candidate_matches {
                        None => {
                            candidate_matches = Some(match_ind_union);
                        }
                        Some(ref mut x) => {
                            x.retain(|i| match_ind_union.contains(i));
                        }
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
    pub fn to_tsv(
        &self,
        region_ids: &Vec<String>,
        library: Option<&Library>,
    ) -> Result<String, LibraryError> {
        // Line has \t separated format:
        // group forward reverse[{region} {region}_nearest {region}_distance {region}_n_matches for each region] status combination_distance combinations_in_library combination_indexes count

        let mut output = String::with_capacity(100 * self.counts.len());

        for (group, count) in self.counts.iter() {
            // Read group
            output.push_str(&group.to_string());
            output.push('\t');

            match &self.sequence {
                Some(seq) => {
                    output.push_str(&seq_to_string_or_log(&seq.forward));
                    output.push('\t');
                    match &seq.reverse {
                        Some(rev) => {
                            output.push_str(&seq_to_string_or_log(rev));
                            output.push('\t');
                        }
                        None => output.push('\t'),
                    }
                }
                None => output.push_str("\t\t"),
            }

            // Region seq/nearest match(s)/distance per region
            for reg_id in region_ids {
                let region = self.regions.get(reg_id);

                match region {
                    None => output.push_str("\t\t\t\t"), // Missing regions 4 blanks
                    Some(r) => {
                        output.push_str(&r.lock().unwrap().to_tsv_chunk());
                        output.push('\t');
                    }
                }
            }

            output.push_str(&self.library_matches.to_tsv_chunk(library)?);
            output.push_str(&count.to_string());
            output.push('\n');
        }

        Ok(output)
    }
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
    fn to_tsv_chunk(&self, library: Option<&Library>) -> Result<String, LibraryError> {
        Ok(match self {
            CombinationMatch::Uncompared => "uncompared\t\t\t\t".to_string(),
            CombinationMatch::Match { ind, distance } => {
                let name = match library {
                    None => ind.to_string(),
                    Some(l) => l.get_name(*ind)?,
                };
                format!("match\t{distance}\t1\t{name}\t",)
            }
            CombinationMatch::MultiMatch { inds, distance } => {
                let names = match library {
                    None => inds
                        .iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>()
                        .join(","),
                    Some(l) => inds
                        .iter()
                        .map(|x| l.get_name(*x))
                        .collect::<Result<Vec<_>, _>>()?
                        .join(","),
                };

                format!("match\t{}\t{}\t{}\t", distance, inds.len(), names)
            }
            CombinationMatch::Recombination { distance } => {
                format!("recombination\t{distance}\t0\t\t",)
            }
            CombinationMatch::Mismatch => "mismatch\t\t0\t\t".to_string(),
            CombinationMatch::Nonmatch => "nonmatch\t\t0\t\t".to_string(),
        })
    }

    /// Output a summary TSV chunk for the combination match status
    ///
    /// Has the \t separated format:
    /// status combinations_in_library combination_indexes
    pub fn to_summary_tsv_chunk(&self, library: Option<&Library>) -> Result<String, LibraryError> {
        Ok(match self {
            CombinationMatch::Uncompared => "uncompared\t\t\t".to_string(),
            CombinationMatch::Match { ind, .. } => {
                let name = match library {
                    None => ind.to_string(),
                    Some(l) => l.get_name(*ind)?,
                };

                format!("match\t1\t{name}\t",)
            }
            CombinationMatch::MultiMatch { inds, .. } => {
                let names = match library {
                    None => inds
                        .iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>()
                        .join(","),
                    Some(l) => inds
                        .iter()
                        .map(|x| l.get_name(*x))
                        .collect::<Result<Vec<_>, _>>()?
                        .join(","),
                };

                format!("match\t{}\t{}\t", inds.len(), names)
            }
            CombinationMatch::Recombination { .. } => "recombination\t0\t\t".to_string(),
            CombinationMatch::Mismatch => "mismatch\t0\t\t".to_string(),
            CombinationMatch::Nonmatch => "nonmatch\t0\t\t".to_string(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use std::sync::{Arc, Mutex};

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
    fn make_library() -> crate::lib_spec::Library {
        use std::collections::HashMap;

        let mut map: HashMap<String, Vec<Vec<u8>>> = HashMap::new();
        let ids = Some(vec![
            "seq1".to_string(),
            "seq2".to_string(),
            "seq3".to_string(),
        ]);
        let region_max: HashMap<String, u64> = HashMap::new();

        map.insert(
            "r1".into(),
            vec![b"ATAT".to_vec(), b"AAAA".to_vec(), b"AAAT".to_vec()],
        );
        map.insert(
            "r2".into(),
            vec![b"GGGG".to_vec(), b"CCCC".to_vec(), b"CCCC".to_vec()],
        );

        // default_max_distance is only relevant for edit-distance modes; 2 is fine.
        crate::lib_spec::Library::new(map, ids, region_max, 2).expect("library builds")
    }

    /// Generate observed combination
    fn make_combination(
        regs: &[(&'static str, &'static [u8], RegionCompleteness)],
    ) -> ObservedCombination {
        let mut map: HashMap<String, Arc<Mutex<ObservedRegion>>> = HashMap::new();
        for (id, seq, comp) in regs.iter().copied() {
            map.insert(
                id.to_string(),
                Arc::new(Mutex::new(ObservedRegion::new(id.to_string(), seq, comp))),
            );
        }
        ObservedCombination::new(map, None)
    }

    /// Run tests on a row
    fn run_row(tc: &TestCase) {
        let lib = make_library();
        let comb = make_combination(&tc.regions);
        let region_ids: Vec<String> = vec!["r1".to_string(), "r2".to_string()];

        let got = comb.compare_to_library(&region_ids, &lib, tc.metric, 3);

        match (&tc.expected, got) {
            (Expected::Match { name, distance }, CombinationMatch::Match { ind, distance: d }) => {
                assert_eq!(d, *distance, "[{}] distance mismatch", tc.name);
                let got_name = lib.get_name(ind).expect("name exists");
                assert_eq!(got_name, *name, "[{}] matched name mismatch", tc.name);
            }
            (
                Expected::MultiMatch { inds_len, distance },
                CombinationMatch::MultiMatch { inds, distance: d },
            ) => {
                assert_eq!(d, *distance, "[{}] distance mismatch", tc.name);
                assert_eq!(
                    inds.len(),
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
