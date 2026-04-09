//! Library matches from sequencing data
//!
//! Structures and functions to store and manipulate library matches extracted from sequencing data
use std::collections::HashMap;

use crate::combination::CombinationMatch;
use crate::errors::LibraryError;
use crate::library::Library;
use crate::region::RegionMatch;
use crate::seqs::ReadGroup;

/// Key identifying a particular library match
///
/// Contains a subset of library match information for use as a hash key
#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub struct LibraryCombinationKey {
    pub regions: Vec<(String, RegionMatch)>,
}

impl LibraryCombinationKey {
    pub fn new(regions: Vec<(String, RegionMatch)>) -> Self {
        Self { regions }
    }
}

/// Summary version of ObservedCombination
///
/// This counts a particular form of match with the expected library instead of a particular
/// combination of observed reads and so uses a hash map of RegionMatches instead of ObservedRegions
#[derive(Debug)]
pub struct LibraryCombination {
    /// Count of observations for each read group. Ungrouped reads are stored in None
    counts: HashMap<ReadGroup, u32>,

    /// RegionMatches determine the connection to the library
    regions: HashMap<String, RegionMatch>,

    /// Status and result of comparison with the expected library of sequences
    library_matches: CombinationMatch,
}

impl LibraryCombination {
    pub fn new(regions: HashMap<String, RegionMatch>, library_matches: CombinationMatch) -> Self {
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

    /// Generate tsv line(s) corresponding to this combination. Each read group
    /// the combination is observed is given a separate line
    pub fn to_tsv(
        &self,
        region_ids: &Vec<String>,
        library: Option<&Library>,
    ) -> Result<String, LibraryError> {
        // Line has \t separated format:
        // group [{region} for each region] status combinations_in_library combination_indexes count

        let mut output = String::with_capacity(100 * self.counts.len());

        for (group, count) in self.counts.iter() {
            // Read group
            match group {
                ReadGroup::Ungrouped => output.push('\t'),
                ReadGroup::Unmatched => output.push_str("_unmatched_\t"),
                ReadGroup::Match(x) => {
                    output.push_str(x);
                    output.push('\t');
                }
            };

            // Region seq
            for reg_id in region_ids {
                let region = self.regions.get(reg_id);

                match region {
                    None => output.push('\t'), // Missing regions 1 blanks
                    Some(r) => {
                        output.push_str(&r.str_sequence());
                        output.push('\t');
                    }
                }
            }

            output.push_str(&self.library_matches.to_summary_tsv_chunk(library)?);
            output.push_str(&count.to_string());
            output.push('\n');
        }

        Ok(output)
    }
}

#[cfg(test)]
mod tests {
    // use super::*;
}
