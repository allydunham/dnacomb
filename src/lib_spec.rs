//! Definition and validation of DNA construct specifications (`LibSpec`).
//!
//! A `LibrarySpec` describes the expected structure of a sequenced construct:
//! the ordered regions that make up the construct, which regions are fixed or
//! variable, where forward and reverse reads are expected to start, and optional
//! per-region matching tolerances for later library comparison.
//!
//! This module also provides helpers for deriving template sequences, expected
//! read sequences, variable-region flanks, and other information needed by the
//! counting algorithms.
use bio::alphabets::dna::revcomp;
use bio::bio_types::sequence::Sequence;
use serde::ser::Error as SerError;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::cmp;
use std::collections::{HashMap, HashSet};
use std::fmt::Display;
use std::fs::read_to_string;
use std::str::FromStr;

use crate::errors::{LibSpecError, seq_to_string_or_log};
use crate::interning::{RegionID, region_id_from_str, region_id_to_str};

fn serialize_region_id<S>(id: &RegionID, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let s = region_id_to_str(id);
    s.serialize(serializer)
}

fn deserialize_region_id<'de, D>(deserializer: D) -> Result<RegionID, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    Ok(region_id_from_str(&s))
}

fn serialize_sequence<S>(seq: &Sequence, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    match std::str::from_utf8(seq) {
        Ok(i) => i.serialize(serializer),
        Err(e) => Err(S::Error::custom(e)),
    }
}

fn deserialize_sequence<'de, D>(deserializer: D) -> Result<Sequence, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    Ok(s.into_bytes())
}

/// Region in a `LibrarySpec`.
///
/// A construct is described as an ordered list of regions. Regions are either:
/// - `Fixed`: constant sequence used as known scaffold/anchor sequence,
/// - `Library`: variable sequence to be extracted from reads and optionally
///   compared to an expected library.
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
        #[serde(
            deserialize_with = "deserialize_region_id",
            serialize_with = "serialize_region_id"
        )]
        id: RegionID,
        min_length: usize,
        max_length: usize,
        max_distance: Option<u64>,
    },

    /// A fixed region such as a primer or scaffold sequence
    ///
    /// id: region id
    /// seq: expected sequence
    Fixed {
        #[serde(
            deserialize_with = "deserialize_region_id",
            serialize_with = "serialize_region_id"
        )]
        id: RegionID,
        #[serde(
            deserialize_with = "deserialize_sequence",
            serialize_with = "serialize_sequence"
        )]
        seq: Sequence,
    },
}

impl Region {
    /// Get the region ID
    pub fn id(&self) -> &RegionID {
        match self {
            Region::Library { id, .. } => id,
            Region::Fixed { id, .. } => id,
        }
    }

    /// Return the region length
    ///
    /// For fixed regions this is the exact sequence length. For library regions
    /// this is the maximum allowed length.
    pub fn len(&self) -> usize {
        match self {
            Region::Fixed { seq, .. } => seq.len(),
            Region::Library { max_length, .. } => *max_length,
        }
    }

    // Return true if the region is empty (i.e. of 0 length)
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Return true if the region is variable (i.e. to be extracted during counting)
    pub fn is_variable(&self) -> bool {
        matches!(self, Region::Library { .. })
    }

    /// Validate region-level constraints.
    ///
    /// Currently this just checks that `min_length <= max_length` for library regions.
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

/// Fixed-sequence context flanking a variable region.
///
/// This describes how a variable region is bounded for pattern-based extraction:
/// - `Unflanked`: no usable fixed sequence on either side,
/// - `OpenStart`: only a downstream flank exists,
/// - `Internal`: fixed flanks exist on both sides,
/// - `OpenEnd`: only an upstream flank exists.
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

/// Specification of a sequenced construct and expected read layout.
///
/// A `LibrarySpec` defines:
/// - the ordered regions composing the construct,
/// - where forward and reverse reads are expected to start,
/// - the expected read lengths,
/// - and optional matching tolerances for variable regions.
///
/// It is the main source of structural information used by all structured
/// counting modes. The struct is used by `serde-json` to parse LibSpec
/// JSON files.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct LibrarySpec {
    /// The name of the library/sequence type
    pub id: String,

    /// Id of the region forward reads start from
    #[serde(
        deserialize_with = "deserialize_region_id",
        serialize_with = "serialize_region_id"
    )]
    pub forward_start_region: RegionID,

    /// Length of forward reads
    pub forward_read_length: u32,

    /// Region reverse reads start in
    #[serde(
        deserialize_with = "deserialize_region_id",
        serialize_with = "serialize_region_id"
    )]
    pub reverse_start_region: RegionID,

    /// Reverse read length
    pub reverse_read_length: u32,

    /// Array of Region objects
    pub regions: Vec<Region>,
}

impl LibrarySpec {
    /// Read a `LibrarySpec` from a JSON file and optionally override read-layout fields.
    ///
    /// The JSON file is parsed and validated, then any provided CLI-style
    /// overrides are applied to:
    /// - forward start region,
    /// - forward read length,
    /// - reverse start region,
    /// - reverse read length.
    ///
    /// The resulting modified spec is validated again before being returned.
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
            lib_spec.forward_start_region = region_id_from_str(&x)
        }

        if let Some(x) = forward_length {
            lib_spec.forward_read_length = x
        }

        if let Some(x) = reverse_start {
            lib_spec.reverse_start_region = region_id_from_str(&x)
        }

        if let Some(x) = reverse_length {
            lib_spec.reverse_read_length = x
        }

        lib_spec.validate()?;

        Ok(lib_spec)
    }

    /// Fetch a region by ID
    pub fn get_region(&self, id: &RegionID) -> Result<&Region, LibSpecError> {
        for r in &self.regions {
            if r.id() == id {
                return Ok(r);
            }
        }

        Err(LibSpecError::MissingRegion { id: id.clone() })
    }

    /// Return per-region maximum library-match distances defined in the spec.
    ///
    /// Only variable regions with an explicit `max_distance` are included.
    pub fn get_max_distances(&self) -> HashMap<RegionID, u64> {
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

    /// Validate the logical integrity of the `LibrarySpec`.
    ///
    /// In addition to serde-level structural validation, this checks:
    /// - that the forward and reverse start regions exist,
    /// - that region IDs are unique,
    /// - that each region is individually valid,
    /// - and that variable regions are not adjacent.
    ///
    /// Adjacent variable regions are rejected because region boundaries cannot be
    /// inferred reliably without fixed anchor sequence between them.
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

        let mut observed_regions: HashSet<RegionID> = HashSet::new();

        // Validate each region
        for region in &self.regions {
            if observed_regions.contains(region.id()) {
                errors.push(format!(
                    "{}",
                    LibSpecError::DuplicateRegion {
                        id: region.id().clone()
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
                            id: region.id().clone()
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

    /// Build the full template sequence implied by the `LibrarySpec`.
    ///
    /// Fixed regions contribute their literal sequence. Variable regions are
    /// represented by `N` repeated to their maximum allowed length.
    ///
    /// This template is primarily used for alignment-based extraction.
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
                Region::Fixed { seq, .. } => template.extend_from_slice(seq),
            }
        }

        template
    }

    /// Construct the minimum expected forward read sequence.
    ///
    /// Starting at `forward_start_region`, this walks forward through the construct
    /// and appends:
    /// - fixed-region sequence verbatim,
    /// - `N` repeated to each variable region's minimum length.
    ///
    /// The result is truncated to `forward_read_length` and acts as a lower-bound
    /// expected read for alignment-threshold calculations.
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
                Region::Fixed { seq, .. } => template.extend_from_slice(seq),
            }

            // Add regions until exhausted or longer than the expected read length
            if template.len() >= self.forward_read_length as usize {
                break;
            }
        }

        template[0..cmp::min(self.forward_read_length as usize, template.len())].to_vec()
    }

    /// Construct the minimum expected reverse read sequence.
    ///
    /// Starting at `reverse_start_region`, this walks backward through the construct
    /// and appends:
    /// - reverse-complemented fixed-region sequence,
    /// - `N` repeated to each variable region's minimum length.
    ///
    /// The result is truncated to `reverse_read_length` and acts as a lower-bound
    /// expected reverse read for alignment-threshold calculations.
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
                Region::Fixed { seq, .. } => template.extend_from_slice(&revcomp(seq)),
            }

            // Add regions until exhausted or longer than the expected read length
            if template.len() >= self.reverse_read_length as usize {
                break;
            }
        }

        template[0..cmp::min(self.reverse_read_length as usize, template.len())].to_vec()
    }

    /// Return the half-open template interval `[start, end)` for a region.
    ///
    /// Coordinates are relative to the template sequence returned by
    /// `template_sequence()`.
    pub fn template_position(&self, region: &RegionID) -> Result<(usize, usize), LibSpecError> {
        let mut start: usize = 0;

        for r in &self.regions {
            if r.id() == region {
                return Ok((start, start + r.len()));
            }
            start += r.len()
        }

        Err(LibSpecError::MissingRegion { id: region.clone() })
    }

    /// Return the variable-region IDs in construct order.
    pub fn variable_regions(&self) -> Vec<RegionID> {
        self.regions
            .iter()
            .filter(|x| x.is_variable())
            .map(|x| x.id().clone())
            .collect()
    }

    /// Compute flanking fixed-sequence patterns for all variable regions.
    ///
    /// Each variable region is assigned up to `len` bases of fixed sequence from
    /// its upstream and/or downstream context, stopping early at construct ends or
    /// at neighbouring variable regions.
    ///
    /// The resulting flank descriptors are validated to ensure they form a
    /// consistent sequence for pattern matching.
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

    /// Validate a list of flank descriptors for use in pattern matching.
    ///
    /// This checks that the flank descriptors form a consistent ordered pattern:
    /// - at most the first region may have an open start,
    /// - at most the last region may have an open end,
    /// - and no unflanked region appears once usable flank patterns are expected.
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

    /// Compute the fixed-sequence context flanking a single variable region.
    ///
    /// Up to `len` bases are collected from the nearest upstream and downstream
    /// fixed regions, stopping if another variable region or a construct boundary
    /// is encountered first.
    ///
    /// This is used by pattern-based region extraction to identify variable
    /// regions from surrounding constant sequence.
    pub fn flanking_regions(
        &self,
        region: &RegionID,
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
                    Region::Fixed { seq, .. } => seq,
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
                    Region::Fixed { seq, .. } => seq,
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

    /// Parse and validate a `LibrarySpec` from a JSON string.
    fn from_str(spec: &str) -> Result<Self, Self::Err> {
        let lib_spec: LibrarySpec = serde_json::from_str::<LibrarySpec>(spec)?;
        lib_spec.validate()?;
        Ok(lib_spec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_region_fixed(id: &str, seq: &[u8]) -> Region {
        Region::Fixed {
            id: region_id_from_str(id),
            seq: seq.to_vec(),
        }
    }

    fn create_test_region_library(
        id: &str,
        min_len: usize,
        max_len: usize,
        max_dist: Option<u64>,
    ) -> Region {
        Region::Library {
            id: region_id_from_str(id),
            min_length: min_len,
            max_length: max_len,
            max_distance: max_dist,
        }
    }

    fn create_test_spec(
        regions: Vec<Region>,
        forward_start: &str,
        forward_len: u32,
        reverse_start: &str,
        reverse_len: u32,
    ) -> LibrarySpec {
        LibrarySpec {
            id: "test_spec".to_string(),
            forward_start_region: region_id_from_str(forward_start),
            forward_read_length: forward_len,
            reverse_start_region: region_id_from_str(reverse_start),
            reverse_read_length: reverse_len,
            regions,
        }
    }

    #[test]
    fn test_region_fixed_id() {
        let region = create_test_region_fixed("r1", b"ATCG");
        assert_eq!(*region.id(), region_id_from_str("r1"));
    }

    #[test]
    fn test_region_fixed_len() {
        let region = create_test_region_fixed("r1", b"ATCG");
        assert_eq!(region.len(), 4);
    }

    #[test]
    fn test_region_fixed_is_empty() {
        let region_empty = create_test_region_fixed("r1", b"");
        let region_full = create_test_region_fixed("r2", b"ATCG");
        assert!(region_empty.is_empty());
        assert!(!region_full.is_empty());
    }

    #[test]
    fn test_region_fixed_is_not_variable() {
        let region = create_test_region_fixed("r1", b"ATCG");
        assert!(!region.is_variable());
    }

    #[test]
    fn test_region_library_id() {
        let region = create_test_region_library("lib1", 10, 20, None);
        assert_eq!(*region.id(), region_id_from_str("lib1"));
    }

    #[test]
    fn test_region_library_len_returns_max() {
        let region = create_test_region_library("lib1", 10, 20, None);
        assert_eq!(region.len(), 20);
    }

    #[test]
    fn test_region_library_is_variable() {
        let region = create_test_region_library("lib1", 10, 20, None);
        assert!(region.is_variable());
    }

    #[test]
    fn test_region_library_validate_valid() {
        let region = create_test_region_library("lib1", 10, 20, None);
        assert!(region.validate().is_ok());
    }

    #[test]
    fn test_region_library_validate_min_greater_than_max() {
        let region = create_test_region_library("lib1", 30, 20, None);
        assert!(region.validate().is_err());
    }

    #[test]
    fn test_region_library_validate_equal_min_max() {
        let region = create_test_region_library("lib1", 20, 20, None);
        assert!(region.validate().is_ok());
    }

    #[test]
    fn test_template_sequence_all_fixed() {
        let regions = vec![
            create_test_region_fixed("r1", b"ATG"),
            create_test_region_fixed("r2", b"TAG"),
        ];
        let spec = create_test_spec(regions, "r1", 10, "r2", 10);
        let template = spec.template_sequence();
        assert_eq!(template, b"ATGTAG");
    }

    #[test]
    fn test_template_sequence_mixed() {
        let regions = vec![
            create_test_region_fixed("r1", b"ATG"),
            create_test_region_library("lib1", 5, 10, None),
            create_test_region_fixed("r2", b"TAG"),
        ];
        let spec = create_test_spec(regions, "r1", 20, "r2", 20);
        let template = spec.template_sequence();
        assert_eq!(template, b"ATGNNNNNNNNNNTAG");
    }

    #[test]
    fn test_expected_forward_read() {
        let regions = vec![
            create_test_region_fixed("r1", b"ATGC"),
            create_test_region_library("lib1", 3, 10, None),
            create_test_region_fixed("r2", b"TAGA"),
        ];
        let spec = create_test_spec(regions, "r1", 10, "r2", 10);
        let read = spec.expected_forward_read();
        assert_eq!(read.len(), 10);
        assert!(read.starts_with(b"ATGC"));
    }

    #[test]
    fn test_expected_forward_read_from_middle_region() {
        let regions = vec![
            create_test_region_fixed("r1", b"ATGC"),
            create_test_region_fixed("r2", b"GGGG"),
            create_test_region_library("lib1", 2, 5, None),
        ];
        let spec = create_test_spec(regions, "r2", 8, "r1", 8);
        let read = spec.expected_forward_read();
        assert!(read.starts_with(b"GGGG"));
    }

    #[test]
    fn test_expected_reverse_read() {
        let regions = vec![
            create_test_region_fixed("r1", b"ATGC"),
            create_test_region_library("lib1", 3, 10, None),
            create_test_region_fixed("r2", b"TAGA"),
        ];
        let spec = create_test_spec(regions, "r1", 10, "r2", 10);
        let read = spec.expected_reverse_read();
        assert_eq!(read.len(), 10);
    }

    #[test]
    fn test_template_position_single_region() {
        let regions = vec![create_test_region_fixed("r1", b"ATGC")];
        let spec = create_test_spec(regions, "r1", 10, "r1", 10);
        let (start, end) = spec.template_position(&region_id_from_str("r1")).unwrap();
        assert_eq!(start, 0);
        assert_eq!(end, 4);
    }

    #[test]
    fn test_template_position_multiple_regions() {
        let regions = vec![
            create_test_region_fixed("r1", b"ATGC"),
            create_test_region_fixed("r2", b"GGGG"),
            create_test_region_library("lib1", 5, 10, None),
        ];
        let spec = create_test_spec(regions, "r1", 20, "lib1", 20);
        let (start, end) = spec.template_position(&region_id_from_str("r2")).unwrap();
        assert_eq!(start, 4);
        assert_eq!(end, 8);

        let (start, end) = spec.template_position(&region_id_from_str("lib1")).unwrap();
        assert_eq!(start, 8);
        assert_eq!(end, 18);
    }

    #[test]
    fn test_template_position_missing_region() {
        let regions = vec![create_test_region_fixed("r1", b"ATGC")];
        let spec = create_test_spec(regions, "r1", 10, "r1", 10);
        assert!(
            spec.template_position(&region_id_from_str("missing"))
                .is_err()
        );
    }

    #[test]
    fn test_variable_regions() {
        let regions = vec![
            create_test_region_fixed("r1", b"ATGC"),
            create_test_region_library("lib1", 5, 10, None),
            create_test_region_fixed("r2", b"GGGG"),
            create_test_region_library("lib2", 3, 8, None),
        ];
        let spec = create_test_spec(regions, "r1", 20, "r2", 20);
        let var_regions = spec.variable_regions();
        assert_eq!(var_regions.len(), 2);
        assert_eq!(var_regions[0], region_id_from_str("lib1"));
        assert_eq!(var_regions[1], region_id_from_str("lib2"));
    }

    #[test]
    fn test_get_max_distances() {
        let regions = vec![
            create_test_region_library("lib1", 5, 10, Some(2)),
            create_test_region_library("lib2", 3, 8, None),
            create_test_region_library("lib3", 5, 10, Some(5)),
        ];
        let spec = create_test_spec(regions, "lib1", 20, "lib3", 20);
        let distances = spec.get_max_distances();
        assert_eq!(distances.len(), 2);
        assert_eq!(distances.get(&region_id_from_str("lib1")), Some(&2));
        assert_eq!(distances.get(&region_id_from_str("lib3")), Some(&5));
        assert_eq!(distances.get(&region_id_from_str("lib2")), None);
    }

    #[test]
    fn test_validate_valid_spec() {
        let regions = vec![
            create_test_region_fixed("r1", b"ATGC"),
            create_test_region_library("lib1", 5, 10, None),
            create_test_region_fixed("r2", b"GGGG"),
        ];
        let spec = create_test_spec(regions, "r1", 20, "r2", 20);
        assert!(spec.validate().is_ok());
    }

    #[test]
    fn test_validate_missing_forward_start() {
        let regions = vec![create_test_region_fixed("r1", b"ATGC")];
        let spec = create_test_spec(regions, "missing", 10, "r1", 10);
        assert!(spec.validate().is_err());
    }

    #[test]
    fn test_validate_missing_reverse_start() {
        let regions = vec![create_test_region_fixed("r1", b"ATGC")];
        let spec = create_test_spec(regions, "r1", 10, "missing", 10);
        assert!(spec.validate().is_err());
    }

    #[test]
    fn test_validate_duplicate_regions() {
        let regions = vec![
            create_test_region_fixed("r1", b"ATGC"),
            create_test_region_fixed("r1", b"GGGG"),
        ];
        let spec = create_test_spec(regions, "r1", 10, "r1", 10);
        assert!(spec.validate().is_err());
    }

    #[test]
    fn test_validate_adjacent_variable_regions() {
        let regions = vec![
            create_test_region_library("lib1", 5, 10, None),
            create_test_region_library("lib2", 3, 8, None),
        ];
        let spec = create_test_spec(regions, "lib1", 20, "lib2", 20);
        assert!(spec.validate().is_err());
    }

    #[test]
    fn test_validate_invalid_library_region() {
        let regions = vec![create_test_region_library("lib1", 20, 10, None)];
        let spec = create_test_spec(regions, "lib1", 20, "lib1", 20);
        assert!(spec.validate().is_err());
    }

    #[test]
    fn test_flanking_regions_internal() {
        let regions = vec![
            create_test_region_fixed("r1", b"ATGC"),
            create_test_region_library("lib1", 5, 10, None),
            create_test_region_fixed("r2", b"GGGG"),
        ];
        let spec = create_test_spec(regions, "r1", 20, "r2", 20);
        let flanks = spec
            .flanking_regions(&region_id_from_str("lib1"), 3)
            .unwrap();
        match flanks {
            FlankingSequences::Internal(before, after) => {
                assert_eq!(before, b"TGC");
                assert_eq!(after, b"GGG");
            }
            _ => panic!("Expected Internal flanks"),
        }
    }

    #[test]
    fn test_flanking_regions_open_start() {
        let regions = vec![
            create_test_region_library("lib1", 5, 10, None),
            create_test_region_fixed("r2", b"GGGG"),
        ];
        let spec = create_test_spec(regions, "lib1", 20, "r2", 20);
        let flanks = spec
            .flanking_regions(&region_id_from_str("lib1"), 3)
            .unwrap();
        match flanks {
            FlankingSequences::OpenStart(after) => {
                assert_eq!(after, b"GGG");
            }
            _ => panic!("Expected OpenStart flanks"),
        }
    }

    #[test]
    fn test_flanking_regions_open_end() {
        let regions = vec![
            create_test_region_fixed("r1", b"ATGC"),
            create_test_region_library("lib1", 5, 10, None),
        ];
        let spec = create_test_spec(regions, "r1", 20, "lib1", 20);
        let flanks = spec
            .flanking_regions(&region_id_from_str("lib1"), 3)
            .unwrap();
        match flanks {
            FlankingSequences::OpenEnd(before) => {
                assert_eq!(before, b"TGC");
            }
            _ => panic!("Expected OpenEnd flanks"),
        }
    }

    #[test]
    fn test_flanking_regions_unflanked() {
        let regions = vec![create_test_region_library("lib1", 5, 10, None)];
        let spec = create_test_spec(regions, "lib1", 20, "lib1", 20);
        let flanks = spec
            .flanking_regions(&region_id_from_str("lib1"), 3)
            .unwrap();
        match flanks {
            FlankingSequences::Unflanked => {}
            _ => panic!("Expected Unflanked"),
        }
    }

    #[test]
    fn test_variable_length_regions() {
        let regions = vec![
            create_test_region_library("lib1", 5, 10, None),
            create_test_region_library("lib2", 8, 8, None),
            create_test_region_library("lib3", 3, 7, None),
            create_test_region_fixed("r1", b"ATGC"),
        ];
        let spec = create_test_spec(regions, "lib1", 20, "r1", 20);
        assert_eq!(spec.variable_length_regions(), 2);
    }

    #[test]
    fn test_get_region_exists() {
        let regions = vec![
            create_test_region_fixed("r1", b"ATGC"),
            create_test_region_library("lib1", 5, 10, None),
        ];
        let spec = create_test_spec(regions, "r1", 20, "lib1", 20);
        let region = spec.get_region(&region_id_from_str("lib1")).unwrap();
        assert!(region.is_variable());
    }

    #[test]
    fn test_get_region_missing() {
        let regions = vec![create_test_region_fixed("r1", b"ATGC")];
        let spec = create_test_spec(regions, "r1", 20, "r1", 20);
        assert!(spec.get_region(&region_id_from_str("missing")).is_err());
    }

    #[test]
    fn test_validate_flank_seqs_valid_internal() {
        let flanks = vec![
            FlankingSequences::OpenStart(b"ATG".to_vec()),
            FlankingSequences::Internal(b"GC".to_vec(), b"TA".to_vec()),
            FlankingSequences::OpenEnd(b"G".to_vec()),
        ];
        assert!(LibrarySpec::validate_flank_seqs(&flanks).is_ok());
    }

    #[test]
    fn test_validate_flank_seqs_invalid_open_start_middle() {
        let flanks = vec![
            FlankingSequences::Internal(b"ATG".to_vec(), b"CG".to_vec()),
            FlankingSequences::OpenStart(b"TA".to_vec()),
        ];
        assert!(LibrarySpec::validate_flank_seqs(&flanks).is_err());
    }

    #[test]
    fn test_validate_flank_seqs_invalid_open_end_middle() {
        let flanks = vec![
            FlankingSequences::OpenEnd(b"ATG".to_vec()),
            FlankingSequences::Internal(b"CG".to_vec(), b"TA".to_vec()),
        ];
        assert!(LibrarySpec::validate_flank_seqs(&flanks).is_err());
    }

    #[test]
    fn test_validate_flank_seqs_invalid_unflanked() {
        let flanks = vec![
            FlankingSequences::Internal(b"ATG".to_vec(), b"CG".to_vec()),
            FlankingSequences::Unflanked,
        ];
        assert!(LibrarySpec::validate_flank_seqs(&flanks).is_err());
    }

    #[test]
    fn test_flanking_regions_display() {
        let unflanked = FlankingSequences::Unflanked;
        assert_eq!(format!("{}", unflanked), "Unflanked");

        let open_start = FlankingSequences::OpenStart(b"ATG".to_vec());
        let open_start_str = format!("{}", open_start);
        assert!(open_start_str.contains("Open"));

        let open_end = FlankingSequences::OpenEnd(b"TAG".to_vec());
        let open_end_str = format!("{}", open_end);
        assert!(open_end_str.contains("Open"));

        let internal = FlankingSequences::Internal(b"ATG".to_vec(), b"TAG".to_vec());
        let internal_str = format!("{}", internal);
        assert!(!internal_str.contains("Open"));
    }
}
