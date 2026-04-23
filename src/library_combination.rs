//! Summary representations of library-assigned region combinations.
//!
//! These types are used after library comparison to collapse full observed
//! combinations into a library-centric summary table. In particular, they remove
//! per-observation detail such as sequence diffs while retaining enough
//! information to group counts by inferred library assignment.
use std::collections::HashMap;
use std::sync::Arc;

use crate::combination::CombinationMatch;
use crate::groups::ReadGroup;
use crate::interning::{RegionID, SeqHandle};
use crate::library::LibraryRegion;
use crate::region::RegionMatch;

/// Key identifying one distinct summarised library-assignment pattern.
///
/// This key is used to collapse multiple `ObservedCombination`s that differ in
/// observed sequence detail but share the same per-region library-assignment
/// summary.
#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub struct LibraryCombinationKey {
    pub regions: Vec<(RegionID, LibraryRegionMatch)>,
}

impl LibraryCombinationKey {
    pub fn new(regions: Vec<(RegionID, LibraryRegionMatch)>) -> Self {
        Self { regions }
    }
}

/// Summary form of a region-to-library match.
///
/// This is a reduced version of [`RegionMatch`] used for library-summary output.
/// Unlike `RegionMatch`, it does not retain additional per-match information (SequenceDiff,
/// distance or number of matches), so multiple observed regions with the same library
/// assignment can be grouped together.
#[derive(Debug, Eq, Hash, PartialEq, Clone)]
pub enum LibraryRegionMatch {
    /// Library comparison has not been performed.
    Uncompared,

    /// A unique library-region assignment was found.
    Match { seq_match: Arc<LibraryRegion> },

    /// Multiple equally good library-region assignments were found.
    MultiMatch {
        seq_matches: Vec<Arc<LibraryRegion>>,
    },

    /// Too many equally good matches were found to report individually.
    Overmatched,

    /// No library-region assignment was found.
    Unmatched,

    /// This region is not represented in the library summary.
    ///
    /// The observed sequence may optionally be retained so outputs can still
    /// report the raw sequence for regions intentionally absent from the library.
    NoLibrary { seq: Option<SeqHandle> },
}

impl LibraryRegionMatch {
    /// Convert a full per-region match into its summary representation.
    ///
    /// This drops sequence-difference detail so that library-summary rows group
    /// by library assignment rather than by exact observed variant.
    pub fn from_region_match(region_match: &RegionMatch) -> Self {
        match region_match {
            RegionMatch::Uncompared => Self::Uncompared,
            RegionMatch::Unmatched => Self::Unmatched,
            RegionMatch::NoLibrary { seq } => Self::NoLibrary { seq: seq.clone() },
            RegionMatch::Overmatched { .. } => Self::Overmatched,
            RegionMatch::Match { seq_match, .. } => Self::Match {
                seq_match: seq_match.clone(),
            },
            RegionMatch::MultiMatch { seq_matches, .. } => Self::MultiMatch {
                seq_matches: seq_matches.clone(),
            },
        }
    }

    /// Return the library-associated sequence(s) for display or TSV output.
    ///
    /// - unique matches return one sequence,
    /// - multimatches return comma-separated sequences,
    /// - `NoLibrary` returns the stored observed sequence if available,
    /// - unmatched/uncompared/overmatched states return an empty string.
    pub fn str_sequence(&self) -> String {
        match self {
            LibraryRegionMatch::Uncompared
            | LibraryRegionMatch::Unmatched
            | LibraryRegionMatch::Overmatched { .. } => "".to_string(),
            LibraryRegionMatch::NoLibrary { seq } => match seq {
                Some(s) => s.to_str_or_log(),
                None => "".to_string(),
            },
            LibraryRegionMatch::Match { seq_match, .. } => seq_match.sequence.to_str_or_log(),
            LibraryRegionMatch::MultiMatch { seq_matches, .. } => seq_matches
                .iter()
                .map(|x| x.sequence.to_str_or_log())
                .collect::<Vec<_>>()
                .join(","),
        }
    }
}

/// Summarised library-assignment counts across observed combinations.
///
/// This groups together observed combinations that share the same per-region
/// summary library assignments and the same overall combination-match status.
/// Counts are accumulated per read group.
#[derive(Debug, Clone)]
pub struct LibraryCombination {
    /// Count of observations for each read group. Ungrouped reads are stored in None
    pub counts: HashMap<ReadGroup, u32>,

    /// RegionMatches determine the connection to the library
    pub regions: HashMap<RegionID, LibraryRegionMatch>,

    /// Status and result of comparison with the expected library of sequences
    pub library_matches: CombinationMatch,
}

impl LibraryCombination {
    /// Create an empty library-summary combination with no counts.
    pub fn new(
        regions: HashMap<RegionID, LibraryRegionMatch>,
        library_matches: CombinationMatch,
    ) -> Self {
        Self {
            counts: HashMap::new(),
            regions,
            library_matches,
        }
    }

    /// Get the total count across all read groups
    pub fn total_count(&self) -> u32 {
        self.counts.values().sum()
    }

    /// Increment the count for the desired read group by n
    pub fn increment_count(&mut self, group: &ReadGroup, n: u32) {
        match self.counts.get_mut(group) {
            Some(x) => *x += n,
            None => {
                self.counts.insert(group.clone(), n);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::collections::HashMap;
    use std::sync::Arc;

    use crate::combination::CombinationMatch;
    use crate::groups::ReadGroup;
    use crate::interning::{library_id_from_str, region_id_from_str, seq_from_bytes};
    use crate::library::LibraryRegion;
    use crate::region::RegionMatch;
    use crate::seq_diff::{EditOperation, SequenceDiff};

    fn make_library_region(seq: &[u8], inds: &[usize], ids: &[&str]) -> Arc<LibraryRegion> {
        Arc::new(LibraryRegion {
            sequence: seq_from_bytes(seq),
            inds: inds.iter().copied().collect(),
            ids: ids.iter().map(|x| library_id_from_str(x)).collect(),
        })
    }

    #[test]
    fn library_combination_key_new_preserves_regions() {
        let key = LibraryCombinationKey::new(vec![
            (
                region_id_from_str("r1"),
                LibraryRegionMatch::NoLibrary {
                    seq: Some(seq_from_bytes(b"ACGT")),
                },
            ),
            (region_id_from_str("r2"), LibraryRegionMatch::Unmatched),
        ]);

        assert_eq!(key.regions.len(), 2);
        assert_eq!(key.regions[0].0, region_id_from_str("r1"));
        assert_eq!(key.regions[1].0, region_id_from_str("r2"));
        assert!(matches!(key.regions[1].1, LibraryRegionMatch::Unmatched));
    }

    #[test]
    fn from_region_match_preserves_summary_state_for_simple_variants() {
        let seq = seq_from_bytes(b"ACGT");
        let unique = make_library_region(b"AAAA", &[0], &["seq1"]);
        let alt1 = make_library_region(b"CCCC", &[1], &["seq2"]);
        let alt2 = make_library_region(b"GGGG", &[2], &["seq3"]);

        let cases = vec![
            (RegionMatch::Uncompared, LibraryRegionMatch::Uncompared),
            (RegionMatch::Unmatched, LibraryRegionMatch::Unmatched),
            (
                RegionMatch::NoLibrary {
                    seq: Some(seq.clone()),
                },
                LibraryRegionMatch::NoLibrary { seq: Some(seq) },
            ),
            (
                RegionMatch::Overmatched {
                    distance: 2,
                    matches: 99,
                },
                LibraryRegionMatch::Overmatched,
            ),
            (
                RegionMatch::Match {
                    seq_match: unique.clone(),
                    distance: 1,
                    diff: SequenceDiff::new(vec![EditOperation::Sub(0, b'A', b'T')]),
                },
                LibraryRegionMatch::Match {
                    seq_match: unique.clone(),
                },
            ),
            (
                RegionMatch::MultiMatch {
                    seq_matches: vec![alt1.clone(), alt2.clone()],
                    distance: 1,
                    diffs: vec![
                        SequenceDiff::new(vec![EditOperation::Sub(0, b'C', b'A')]),
                        SequenceDiff::new(vec![EditOperation::Sub(0, b'G', b'A')]),
                    ],
                },
                LibraryRegionMatch::MultiMatch {
                    seq_matches: vec![alt1.clone(), alt2.clone()],
                },
            ),
        ];

        for (full, expected) in cases {
            let got = LibraryRegionMatch::from_region_match(&full);
            assert_eq!(got, expected);
        }
    }

    #[test]
    fn str_sequence_formats_each_variant_correctly() {
        let unique = make_library_region(b"AAAA", &[0], &["seq1"]);
        let alt1 = make_library_region(b"CCCC", &[1], &["seq2"]);
        let alt2 = make_library_region(b"GGGG", &[2], &["seq3"]);

        let cases = vec![
            (LibraryRegionMatch::Uncompared, ""),
            (LibraryRegionMatch::Unmatched, ""),
            (LibraryRegionMatch::Overmatched, ""),
            (
                LibraryRegionMatch::NoLibrary {
                    seq: Some(seq_from_bytes(b"TTTT")),
                },
                "TTTT",
            ),
            (LibraryRegionMatch::NoLibrary { seq: None }, ""),
            (
                LibraryRegionMatch::Match {
                    seq_match: unique.clone(),
                },
                "AAAA",
            ),
            (
                LibraryRegionMatch::MultiMatch {
                    seq_matches: vec![alt1.clone(), alt2.clone()],
                },
                "CCCC,GGGG",
            ),
        ];

        for (input, expected) in cases {
            assert_eq!(input.str_sequence(), expected);
        }
    }

    #[test]
    fn library_combination_new_starts_empty() {
        let mut regions = HashMap::new();
        regions.insert(
            region_id_from_str("r1"),
            LibraryRegionMatch::NoLibrary {
                seq: Some(seq_from_bytes(b"ACGT")),
            },
        );

        let comb = LibraryCombination::new(regions.clone(), CombinationMatch::Mismatch);

        assert!(comb.counts.is_empty());
        assert_eq!(comb.regions.len(), 1);
        assert_eq!(comb.regions.keys().next(), Some(&region_id_from_str("r1")));
        assert!(matches!(comb.library_matches, CombinationMatch::Mismatch));
    }

    #[test]
    fn increment_count_accumulates_by_group() {
        let comb_match = CombinationMatch::Match {
            inds: vec![Some(library_id_from_str("seq1"))],
            distance: 0,
        };

        let mut comb = LibraryCombination::new(HashMap::new(), comb_match);

        comb.increment_count(&ReadGroup::ungrouped(), 1);
        comb.increment_count(&ReadGroup::ungrouped(), 2);
        comb.increment_count(&ReadGroup::grouped("g1"), 3);

        assert_eq!(comb.counts.get(&ReadGroup::ungrouped()), Some(&3));
        assert_eq!(comb.counts.get(&ReadGroup::grouped("g1")), Some(&3));
        assert_eq!(comb.total_count(), 6);
    }

    #[test]
    fn total_count_sums_all_groups() {
        let mut comb = LibraryCombination::new(HashMap::new(), CombinationMatch::Nonmatch);

        comb.increment_count(&ReadGroup::ungrouped(), 5);
        comb.increment_count(&ReadGroup::grouped("a"), 7);
        comb.increment_count(&ReadGroup::grouped("b"), 11);

        assert_eq!(comb.total_count(), 23);
    }

    #[test]
    fn library_region_match_equality_is_summary_based() {
        let reg = make_library_region(b"AAAA", &[0, 1], &["seq1", "seq2"]);

        let a = LibraryRegionMatch::Match {
            seq_match: reg.clone(),
        };
        let b = LibraryRegionMatch::Match {
            seq_match: reg.clone(),
        };

        assert_eq!(a, b);

        let mm1 = LibraryRegionMatch::MultiMatch {
            seq_matches: vec![
                make_library_region(b"CCCC", &[0], &["seqA"]),
                make_library_region(b"GGGG", &[1], &["seqB"]),
            ],
        };
        let mm2 = LibraryRegionMatch::MultiMatch {
            seq_matches: vec![
                make_library_region(b"CCCC", &[0], &["seqA"]),
                make_library_region(b"GGGG", &[1], &["seqB"]),
            ],
        };

        assert_eq!(mm1, mm2);
    }

    #[test]
    fn library_combination_key_hash_and_eq_are_stable() {
        use std::collections::HashSet;

        let key1 = LibraryCombinationKey::new(vec![
            (
                region_id_from_str("r1"),
                LibraryRegionMatch::NoLibrary {
                    seq: Some(seq_from_bytes(b"ACGT")),
                },
            ),
            (
                region_id_from_str("r2"),
                LibraryRegionMatch::Match {
                    seq_match: make_library_region(b"AAAA", &[0], &["seq1"]),
                },
            ),
        ]);

        let key2 = LibraryCombinationKey::new(vec![
            (
                region_id_from_str("r1"),
                LibraryRegionMatch::NoLibrary {
                    seq: Some(seq_from_bytes(b"ACGT")),
                },
            ),
            (
                region_id_from_str("r2"),
                LibraryRegionMatch::Match {
                    seq_match: make_library_region(b"AAAA", &[0], &["seq1"]),
                },
            ),
        ]);

        let key3 = LibraryCombinationKey::new(vec![
            (
                region_id_from_str("r1"),
                LibraryRegionMatch::NoLibrary {
                    seq: Some(seq_from_bytes(b"TGCA")),
                },
            ),
            (
                region_id_from_str("r2"),
                LibraryRegionMatch::Match {
                    seq_match: make_library_region(b"AAAA", &[0], &["seq1"]),
                },
            ),
        ]);

        let mut set = HashSet::new();
        assert!(set.insert(key1.clone()));
        assert!(!set.insert(key2), "equivalent key should not insert twice");
        assert!(set.insert(key3), "different summary key should insert");
        assert_eq!(set.len(), 2);
    }
}
