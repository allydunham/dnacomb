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
            RegionMatch::NoLibrary { seq } => Self::NoLibrary { seq: *seq },
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
#[derive(Debug)]
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
                self.counts.insert(*group, n);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    // use super::*;
}
