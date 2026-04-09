//! Specification for DNA constructs and libraries
//!
//! Provides methdods for importing JSON based DNA construct specifications
//! and manipulating them. Additionally supports TSV libraries corresponding
//! to these constructs, with lookup capabilities.
use bio::alphabets::dna::revcomp;
use bio::bio_types::sequence::Sequence;
use serde::{Deserialize, Serialize};
use std::cmp;
use std::collections::{HashMap, HashSet};
use std::fmt::Display;
use std::fs::read_to_string;
use std::str::FromStr;

use crate::errors::{LibSpecError, seq_to_string_or_log};

/// LibSpec region types
///
/// Specification for serde json to parse LibSpec regions
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(tag = "seq_type")]
pub enum Region {
    /// A variable region with a list of possible values/combinations in a library. For
    /// example a CRISPR spacer or barcode.
    ///
    /// id: region id
    /// min_length: minimum expected length
    /// max_length: maximum expected length
    /// max_distance: maximum number of mismatches for an observed sequence to be considered
    /// a library match
    Library {
        id: String,
        min_length: usize,
        max_length: usize,
        max_distance: Option<u64>,
    },

    /// A fixed region such as a primer or scaffold sequence
    ///
    /// id: region id
    /// seq: expected sequence
    Fixed { id: String, seq: String },
}

impl Region {
    /// Get the region ID
    pub fn id(&self) -> &String {
        match self {
            Region::Library { id, .. } => id,
            Region::Fixed { id, .. } => id,
        }
    }

    /// Get the region length
    pub fn len(&self) -> usize {
        match self {
            Region::Fixed { seq, .. } => seq.len(),
            Region::Library { max_length, .. } => *max_length,
        }
    }

    // Is the region "empty" (i.e. of 0 length)
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// If the region is variable (i.e. to be extracted during counting)
    pub fn is_variable(&self) -> bool {
        matches!(self, Region::Library { .. })
    }

    /// Check the region is valid
    pub fn validate(&self) -> Result<(), LibSpecError> {
        match self {
            Region::Library {
                id,
                min_length,
                max_length,
                ..
            } => {
                if min_length > max_length {
                    return Err(LibSpecError::MinGreaterThanMax {
                        id: id.clone(),
                        min: *min_length,
                        max: *max_length,
                    });
                }
            }
            Region::Fixed { .. } => {}
        }

        Ok(())
    }
}

/// Sequence of flanking regions around a sequence of interest
#[derive(Debug)]
pub enum FlankingSequences {
    Unflanked,
    OpenStart(Sequence),
    Internal(Sequence, Sequence),
    OpenEnd(Sequence),
}

impl Display for FlankingSequences {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FlankingSequences::Unflanked => write!(f, "Unflanked"),
            FlankingSequences::OpenStart(end) => write!(f, "(Open, {})", seq_to_string_or_log(end)),
            FlankingSequences::Internal(start, end) => write!(
                f,
                "({}, {})",
                seq_to_string_or_log(start),
                seq_to_string_or_log(end)
            ),
            FlankingSequences::OpenEnd(start) => {
                write!(f, "({}, Open)", seq_to_string_or_log(start))
            }
        }
    }
}

/// LibSpec definition
///
/// Specification for serde json to parse LibSpec JSON files
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct LibrarySpec {
    /// The name of the library/sequence type
    pub id: String,

    /// Id of the region forward reads start from
    pub forward_start_region: String,

    /// Length of forward reads
    pub forward_read_length: u32,

    /// Region reverse reads start in
    pub reverse_start_region: String,

    /// Reverse read length
    pub reverse_read_length: u32,

    /// Array of Region objects
    pub regions: Vec<Region>,
}

impl LibrarySpec {
    /// Read a LibrarySpec from a JSON file
    pub fn from_file(
        path: &str,
        forward_start: Option<String>,
        forward_length: Option<u32>,
        reverse_start: Option<String>,
        reverse_length: Option<u32>,
    ) -> Result<LibrarySpec, LibSpecError> {
        let json_str: String = read_to_string(path)?;
        let mut lib_spec: LibrarySpec = LibrarySpec::from_str(&json_str)?;

        // Optionally override read structure
        if let Some(x) = forward_start {
            lib_spec.forward_start_region = x
        }

        if let Some(x) = forward_length {
            lib_spec.forward_read_length = x
        }

        if let Some(x) = reverse_start {
            lib_spec.reverse_start_region = x
        }

        if let Some(x) = reverse_length {
            lib_spec.reverse_read_length = x
        }

        lib_spec.validate()?;

        Ok(lib_spec)
    }

    /// Fetch a specified region
    pub fn get_region(&self, id: &str) -> Result<&Region, LibSpecError> {
        for r in &self.regions {
            if r.id() == id {
                return Ok(r);
            }
        }

        Err(LibSpecError::MissingRegion { id: id.to_string() })
    }

    /// Get a HashMap of max_distances per region
    pub fn get_max_distances(&self) -> HashMap<String, u64> {
        let mut max_dists = HashMap::new();

        for r in &self.regions {
            match r {
                Region::Fixed { .. } => (),
                Region::Library {
                    id, max_distance, ..
                } => match max_distance {
                    None => (),
                    Some(x) => {
                        max_dists.insert(id.clone(), *x);
                    }
                },
            }
        }

        max_dists
    }

    /// Validate the integrity of a LibrarySpec, raising an error if any annomolies are identified
    ///
    /// Much of the LibrarySpec format is checked by Serde as part of it's definition and
    /// deserialisation process but some properties are not amenable to this. These are validated
    /// here instead.
    pub fn validate(&self) -> Result<(), LibSpecError> {
        let mut errors: Vec<String> = Vec::new();

        // Check indicated read start regions are present
        match &self.get_region(&self.forward_start_region) {
            Ok(_) => {}
            Err(err) => errors.push(format!("{}", err)),
        }

        match &self.get_region(&self.reverse_start_region) {
            Ok(_) => {}
            Err(err) => errors.push(format!("{}", err)),
        }

        let mut observed_regions: HashSet<String> = HashSet::new();

        // Validate each region
        for region in &self.regions {
            if observed_regions.contains(region.id()) {
                errors.push(format!(
                    "{}",
                    LibSpecError::DuplicateRegion {
                        id: region.id().to_string()
                    }
                ))
            }
            observed_regions.insert(region.id().clone());

            match region.validate() {
                Ok(_) => {}
                Err(err) => errors.push(format!("{}", err)),
            }
        }

        // Check no variable regions next to each other
        let mut last_variable = false;
        for region in &self.regions {
            if region.is_variable() {
                if last_variable {
                    errors.push(format!(
                        "{}",
                        LibSpecError::NeighbouringVariable {
                            id: region.id().to_string()
                        }
                    ))
                }
                last_variable = true;
            } else {
                last_variable = false;
            }
        }

        if !errors.is_empty() {
            return Err(LibSpecError::InvalidLibSpec { errs: errors });
        }

        Ok(())
    }

    /// Generate a template sequence from the sequence specification
    ///
    /// This function returns a template sequence with the expected structure of sequences
    /// coming from this library, for example to align reads against. It is the concatenation
    /// of each region in the library with max_length Ns included for variable regions.
    pub fn template_sequence(&self) -> Sequence {
        let mut len: usize = 0;
        for region in &self.regions {
            match region {
                Region::Library { max_length, .. } => len += max_length,
                Region::Fixed { seq, .. } => len += seq.len(),
            }
        }

        let mut template: Sequence = Vec::with_capacity(len);

        for region in &self.regions {
            match region {
                Region::Library { max_length, .. } => {
                    for _ in 0..*max_length {
                        template.push(b'N')
                    }
                }
                Region::Fixed { seq, .. } => template.extend_from_slice(&seq.clone().into_bytes()),
            }
        }

        template
    }

    /// Expected forward read
    ///
    /// Generate a minimum length expected forward read as a lower bound for
    /// alignment to the template
    pub fn expected_forward_read(&self) -> Sequence {
        let mut template: Sequence = Vec::with_capacity(self.forward_read_length as usize);
        let mut read_started = false;

        for region in &self.regions {
            // Ignore regions until the start region is found
            if !read_started && *region.id() == self.forward_start_region {
                read_started = true;
            } else if !read_started {
                continue;
            }

            match region {
                Region::Library { min_length, .. } => {
                    for _ in 0..*min_length {
                        template.push(b'N')
                    }
                }
                Region::Fixed { seq, .. } => template.extend_from_slice(&seq.clone().into_bytes()),
            }

            // Add regions until exhausted or longer than the expected read length
            if template.len() >= self.forward_read_length as usize {
                break;
            }
        }

        template[0..cmp::min(self.forward_read_length as usize, template.len())].to_vec()
    }

    /// Expected reverse read
    ///
    /// Generate a minimum length expected reverse read as a lower bound for
    /// alignment to the template
    pub fn expected_reverse_read(&self) -> Sequence {
        let mut template: Sequence = Vec::with_capacity(self.forward_read_length as usize);
        let mut read_started = false;

        for region in self.regions.iter().rev() {
            // Ignore regions until the start region is found
            if !read_started && *region.id() == self.reverse_start_region {
                read_started = true;
            } else if !read_started {
                continue;
            }

            match region {
                Region::Library { min_length, .. } => {
                    for _ in 0..*min_length {
                        template.push(b'N')
                    }
                }
                Region::Fixed { seq, .. } => {
                    template.extend_from_slice(&revcomp(seq.clone().into_bytes()))
                }
            }

            // Add regions until exhausted or longer than the expected read length
            if template.len() >= self.reverse_read_length as usize {
                break;
            }
        }

        template[0..cmp::min(self.reverse_read_length as usize, template.len())].to_vec()
    }

    /// Identify the position of a region in the library template sequence
    ///
    /// This function returns a tuple of the start/end indeces of the passed region, as a half
    /// open interval [a, b) as used for rust vector slices.
    pub fn template_position(&self, region: &str) -> Result<(usize, usize), LibSpecError> {
        let mut start: usize = 0;

        for r in &self.regions {
            if r.id() == region {
                return Ok((start, start + r.len()));
            }
            start += r.len()
        }

        Err(LibSpecError::MissingRegion {
            id: region.to_string(),
        })
    }

    /// Identify the variable regions in the library
    pub fn variable_regions(&self) -> Vec<String> {
        self.regions
            .iter()
            .filter(|x| x.is_variable())
            .map(|x| x.id().clone())
            .collect()
    }

    /// Get flanking sequences for all variable regions
    pub fn get_all_flanking_regions(
        &self,
        len: usize,
    ) -> Result<Vec<FlankingSequences>, LibSpecError> {
        let regions = self.variable_regions();
        let flanks = regions
            .iter()
            .map(|x| self.flanking_regions(x, len))
            .collect::<Result<Vec<FlankingSequences>, LibSpecError>>()?;

        Self::validate_flank_seqs(&flanks)?;

        Ok(flanks)
    }

    /// Validate flanking regions
    ///
    /// Currently check that they form a valid and findable sequence of region types
    pub fn validate_flank_seqs(flanks: &[FlankingSequences]) -> Result<(), LibSpecError> {
        for (i, r) in flanks.iter().enumerate() {
            match r {
                FlankingSequences::Unflanked => {
                    return Err(LibSpecError::LibSpec {
                        desc: "Unflanked region after all flank patterns found".to_string(),
                    });
                }
                FlankingSequences::OpenStart(..) => {
                    if i == 0 {
                        continue;
                    }
                    return Err(LibSpecError::LibSpec {
                        desc: "Region with an open start found after first region in flanking patterns".to_string()
                    });
                }
                FlankingSequences::Internal(..) => continue,
                FlankingSequences::OpenEnd(..) => {
                    if i == flanks.len() - 1 {
                        continue;
                    }
                    return Err(LibSpecError::LibSpec {
                        desc:
                            "Region with an open end found before final region in flanking patterns"
                                .to_string(),
                    });
                }
            }
        }

        Ok(())
    }

    /// Identify the sequences flanking a region of interest
    ///
    /// If the region is first/last the corresponding flanking region is None, otherwise
    /// it is `Some<vec<u8>>` up to len long (less if a variable region or the end is reached).
    pub fn flanking_regions(
        &self,
        region: &str,
        len: usize,
    ) -> Result<FlankingSequences, LibSpecError> {
        let reg_ind = self
            .regions
            .iter()
            .enumerate()
            .find_map(|(i, x)| if x.id() == region { Some(i) } else { None })
            .unwrap();

        // Identify before
        let mut before: Vec<u8> = Vec::new();
        if reg_ind > 0 {
            // reg_ind == 0 means first region
            let mut i = reg_ind - 1;
            loop {
                let len_needed = len - before.len();
                if len_needed == 0 {
                    break;
                }

                let seq = match &self.regions[i] {
                    Region::Library { .. } => break,
                    Region::Fixed { seq, .. } => seq.as_bytes(),
                };

                // Extend before with up to len_needed in reverse
                before.extend(seq.iter().rev().take(len_needed));

                if i == 0 {
                    break;
                }
                i -= 1;
            }
            // Return to the cannonical order
            before.reverse();
        }

        // Identify after
        let mut after: Vec<u8> = Vec::new();
        if reg_ind < self.regions.len() - 1 {
            // Similarly for non-end regions
            let mut i = reg_ind + 1;
            while i < self.regions.len() {
                let len_needed = len - after.len();
                if len_needed == 0 {
                    break;
                }

                let seq = match &self.regions[i] {
                    Region::Library { .. } => break,
                    Region::Fixed { seq, .. } => seq.as_bytes(),
                };

                // Extend with up to len_needed
                after.extend(seq.iter().take(len_needed));

                i += 1;
            }
        }

        Ok(match (before.is_empty(), after.is_empty()) {
            (true, true) => FlankingSequences::Unflanked,
            (true, false) => FlankingSequences::OpenStart(after),
            (false, true) => FlankingSequences::OpenEnd(before),
            (false, false) => FlankingSequences::Internal(before, after),
        })
    }

    /// Determine number of variable length region
    ///
    /// Returns the count of variable regions with min_length != max_length
    pub fn variable_length_regions(&self) -> usize {
        let mut count: usize = 0;

        for region in &self.regions {
            match region {
                Region::Library {
                    min_length,
                    max_length,
                    ..
                } => {
                    if min_length != max_length {
                        count += 1
                    }
                }
                Region::Fixed { .. } => continue,
            }
        }

        count
    }
}

impl FromStr for LibrarySpec {
    type Err = LibSpecError;

    /// Parse a LibrarySpec from a JSON string
    fn from_str(spec: &str) -> Result<Self, Self::Err> {
        let lib_spec: LibrarySpec = serde_json::from_str::<LibrarySpec>(spec)?;
        lib_spec.validate()?;
        Ok(lib_spec)
    }
}

#[cfg(test)]
mod tests {
    // use super::*;
}
