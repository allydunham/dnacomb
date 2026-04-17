//! Library matches from sequencing data
//!
//! Structures and functions to store and manipulate library matches extracted from sequencing data
use std::collections::HashMap;
use std::sync::Arc;

use crate::combination::CombinationMatch;
use crate::groups::ReadGroup;
use crate::interning::{RegionID, SeqHandle};
use crate::library::LibraryRegion;
use crate::region::RegionMatch;

/// Key identifying a particular library match
///
/// Contains a subset of library match information for use as a hash key
#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub struct LibraryCombinationKey {
    pub regions: Vec<(RegionID, LibraryRegionMatch)>,
}

impl LibraryCombinationKey {
    pub fn new(regions: Vec<(RegionID, LibraryRegionMatch)>) -> Self {
        Self { regions }
    }
}

/// Status of the match between a region and a Library
///
/// Includes the status plus reference(s) to the matched sequence in
/// the library and how distant it is. Different from RegionMatch as
/// it doesn't include SeqDiffs, meaning regions can be combined over
/// sequences.
#[derive(Debug, Eq, Hash, PartialEq, Clone)]
pub enum LibraryRegionMatch {
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
    NoLibrary { seq: Option<SeqHandle> },
}

impl LibraryRegionMatch {
    pub fn from_region_match(region_match: &RegionMatch) -> Self {
        match region_match {
            RegionMatch::Uncompared => Self::Uncompared,
            RegionMatch::Unmatched => Self::Unmatched,
            RegionMatch::NoLibrary { seq } => Self::NoLibrary { seq: *seq },
            RegionMatch::Overmatched { .. } => Self::Overmatched {
                distance: 0,
                matches: 0,
            },
            RegionMatch::Match { seq_match, .. } => Self::Match {
                seq_match: seq_match.clone(),
                distance: 0,
            },
            RegionMatch::MultiMatch { seq_matches, .. } => Self::MultiMatch {
                seq_matches: seq_matches.clone(),
                distance: 0,
            },
        }
    }

    /// Get matching Sequence as a string, including the passed through raw
    /// Sequence for NoLibrary matches.
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

/// Summary version of ObservedCombination
///
/// This counts a particular form of match with the expected library instead of a particular
/// combination of observed reads and so uses a hash map of RegionMatches instead of ObservedRegions
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
