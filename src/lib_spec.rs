//! Specification for DNA constructs and libraries
//!
//! Provides methdods for importing JSON based DNA construct specifications
//! and manipulating them. Additionally supports TSV libraries corresponding
//! to these constructs, with lookup capabilities.
use bio::alignment::distance;
use bio::alphabets::dna::revcomp;
use bio::bio_types::sequence::Sequence;
use clap::ValueEnum;
use csv::ReaderBuilder;
use serde::{Deserialize, Serialize};
use std::cmp;
use std::collections::{HashMap, HashSet};
use std::fs::read_to_string;
use std::hash::{Hash, Hasher};
use std::sync::Arc;
use std::str::FromStr;

use crate::errors::{LibSpecError, LibraryError, seq_to_string_or_log};

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
    /// length: length (must match seq)
    Fixed {
        id: String,
        seq: String,
        length: usize,
    },
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
            Region::Fixed { length, .. } => *length,
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
    OpenEnd(Sequence)
}

impl ToString for FlankingSequences {
    fn to_string(&self) -> String {
        match self {
            FlankingSequences::Unflanked => "Unflanked".to_string(),
            FlankingSequences::OpenStart(end) => format!("(Open, {})", seq_to_string_or_log(end)),
            FlankingSequences::Internal(start, end) => format!(
                "({}, {})", seq_to_string_or_log(start), seq_to_string_or_log(end)
            ),
            FlankingSequences::OpenEnd(start) => format!("({}, Open)", seq_to_string_or_log(start))
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

    /// TSV file containing expected sequences
    ///
    /// This must be strictly tab separated and have columns
    /// for each region that will be compared to the library.
    /// Each row gives an expected combination. It will be used
    /// to initialise a Library object.
    pub library: Option<String>,
}

impl LibrarySpec {
    /// Read a LibrarySpec from a JSON file
    pub fn from_file(path: &str) -> Result<LibrarySpec, LibSpecError> {
        let json_str: String = read_to_string(path)?;
        let lib_spec: LibrarySpec = LibrarySpec::from_str(&json_str)?;
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
                Region::Fixed { length, .. } => len += length,
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
    pub fn get_all_flanking_regions(&self, len: usize) -> Result<Vec<FlankingSequences>, LibSpecError>{
        let regions = self.variable_regions();
        let flanks = regions
            .iter()
            .map(|x| self.flanking_regions(x, len))
            .collect::<Result<Vec<FlankingSequences>, LibSpecError>>()?;

        Self::validate_flank_seqs(&flanks)?;

        return Ok(flanks)
    }

    /// Validate flanking regions
    ///
    /// Currently check that they form a valid and findable sequence of region types
    pub fn validate_flank_seqs(flanks: &[FlankingSequences]) -> Result<(), LibSpecError>{
        for (i, r) in flanks.iter().enumerate() {
            match r {
                FlankingSequences::Unflanked => {
                    return Err(LibSpecError::LibSpec {
                        desc: "Unflanked region after all flank patterns found".to_string()
                    });
                },
                FlankingSequences::OpenStart( .. ) => {
                    if i == 0 {continue;}
                    return Err(LibSpecError::LibSpec {
                        desc: "Region with an open start found after first region in flanking patterns".to_string()
                    });
                },
                FlankingSequences::Internal( .. ) => continue,
                FlankingSequences::OpenEnd( .. ) => {
                    if i == flanks.len() - 1 {continue;}
                    return Err(LibSpecError::LibSpec {
                        desc: "Region with an open end found before final region in flanking patterns".to_string()
                    });
                },
            }
        }

        return Ok(());
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

/// A compiled sequence library, carrying the expected sequence combinations in each variable region
///
/// This allows efficient lookup of candidate sequences against the library
#[derive(Debug)]
pub struct Library {
    /// Full sequences for each member of the library, divided into region vectors. The full nth
    /// sequence contains the nth sequence from each region vector
    pub library: HashMap<String, Vec<Arc<LibraryRegion>>>,

    /// Unique sequences for each region, mapping back to which full combinations they are part
    /// of by index
    pub regions: HashMap<String, Vec<Arc<LibraryRegion>>>,

    /// HashMap of exact hits to Library regions for quick initial lookup and
    /// exact matching
    exact_matches: HashMap<String, HashMap<Sequence, Arc<LibraryRegion>>>,

    /// Max distance to consider for each region
    region_max_distance: HashMap<String, u64>,

    /// Default max distance to consider
    default_max_distance: u64,
}

/// A sequence region from a compiled library
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct LibraryRegion {
    /// The `Vec<u8>` sequence
    pub sequence: Sequence,

    /// The library members that contain this sequence
    pub inds: HashSet<usize>,
}

impl Hash for LibraryRegion {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.sequence.hash(state);
        self.inds
            .iter()
            .copied()
            .collect::<Vec<usize>>()
            .hash(state);
    }
}

/// Match with a LibraryRegion at a given distance
#[derive(Debug)]
pub struct LibraryMatch {
    pub matches: Vec<Arc<LibraryRegion>>,
    pub distance: u64,
}

/// Combine two library matches to matches consistent with both
/// Distance is summed, which makes sense for the desired case of partial
/// matches at both ends but may double count if overlapping
pub fn merge_matches(x: Option<LibraryMatch>, y: Option<LibraryMatch>) -> Option<LibraryMatch> {
    match (x, y) {
        (None, _) | (_, None) => None,
        (Some(x), Some(y)) => {
            let matches: Vec<Arc<LibraryRegion>> = x
                .matches
                .iter()
                .filter(|m| y.matches.contains(m))
                .cloned()
                .collect();

            if matches.is_empty() {
                None
            } else {
                Some(LibraryMatch {
                    matches,
                    distance: x.distance + y.distance,
                })
            }
        }
    }
}

impl Library {
    pub fn new(
        library: HashMap<String, Vec<Sequence>>,
        region_max_distance: HashMap<String, u64>,
        default_max_distance: u64,
    ) -> Result<Library, LibraryError> {
        let mut regions = HashMap::new();

        // This allows an empty library
        if library.len() > 1 {
            let exp_len: usize = library
                .values()
                .next()
                .expect("Just checked library has at least 2 elements")
                .len();
            if !library.values().all(|x| x.len() == exp_len) {
                return Err(LibraryError::Library {
                    desc: "Library must contain the same number of sequences for each region"
                        .to_string(),
                });
            }
        }

        for key in library.keys() {
            let mut reg_map: HashMap<Sequence, HashSet<usize>> = HashMap::new();
            let seqs = match library.get(key) {
                Some(x) => x,
                None => {
                    return Err(LibraryError::Library {
                        desc: "Library sequence vec missing unexpectedly during compilation"
                            .to_string(),
                    });
                }
            };

            for (ind, seq) in seqs.iter().enumerate() {
                match reg_map.get_mut(seq) {
                    Some(x) => {
                        x.insert(ind);
                    }
                    None => {
                        reg_map.insert(seq.clone(), HashSet::from([ind]));
                    }
                }
            }

            if regions
                .insert(
                    key.clone(),
                    Vec::from_iter(reg_map.into_iter().map(|x| {
                        Arc::new(LibraryRegion {
                            sequence: x.0,
                            inds: x.1,
                        })
                    })),
                )
                .is_some()
            {
                return Err(LibraryError::DuplicateRegion {
                    id: key.to_string(),
                });
            }
        }

        // Compile Exact matches HashMap
        let mut exact_matches = HashMap::new();
        for key in regions.keys() {
            exact_matches.insert(key.clone(), HashMap::new());
            for reg in regions
                .get(key)
                .expect("Key known to be in regions HashMap")
            {
                exact_matches
                    .get_mut(key)
                    .expect("Key just added to exact_matchs")
                    .insert(reg.sequence.clone(), reg.clone());
            }
        }

        // Reconstruct library with links to correct region Rcs
        let mut library_compiled = HashMap::new();

        for (key, seqs) in library {
            let mut rc_vec: Vec<Option<Arc<LibraryRegion>>> = vec![None; seqs.len()];

            for reg in regions
                .get(&key)
                .expect("regions should contain key as just inserted")
            {
                for ind in &reg.inds {
                    rc_vec[*ind] = Some(reg.clone());
                }
            }

            library_compiled.insert(
                key.clone(),
                rc_vec.into_iter().map(
                    |x| x.expect("All library members should have been assigned Some(Arc<LibraryRegion>) by construction")).collect()
            );
        }

        Ok(Library {
            library: library_compiled,
            regions,
            exact_matches,
            region_max_distance,
            default_max_distance,
        })
    }

    #[allow(dead_code)] // not used in count_reads but useful for users
    /// Number of elements in the library (i.e. length of one seq vector)
    pub fn len(&self) -> usize {
        if self.library.is_empty() {
            return 0;
        }

        self.library
            .values()
            .next()
            .expect("Returned previously if library empty")
            .len()
    }

    #[allow(dead_code)] // not used in count_reads but useful for users
    /// Check is the library is empty
    pub fn is_empty(&self) -> bool {
        self.library.is_empty()
    }

    /// Compare an observed sequence to the library
    ///
    /// Itentify the library members that most closely match a query sequence, with options
    /// for distance metric to use and whether to require matches to the full sequence or
    /// just one end. This is useful where you know your query is incomplete compared to the
    /// library regions.
    ///
    /// The implementation dispatches to the appropriate lookup funciton based on distance matric
    /// and partial match type. We use the SIMD optimised versions of Hamming and Levenshtein
    /// distance if AVX2 and SSE4.1 are available at compile time. The Rust Bio SIMD distance
    /// metrics fall back to standard versions if SIMD isn't available so this is just a minor
    /// optimisation and either version should be portable without undefined behaviour. Bounded
    /// levenshtein is only available as a SIMD implementation with fallback so we rely on
    /// the Rust Bio and editdistancek authors for the check.
    ///
    /// The max distance is used as the upper bound for bounded Levenshteinso this give identical
    /// results to Levenshtein in less time. Therefore generally bounded should be prefered to
    /// Levenshtein but both options are available in case of edge cases.
    pub fn lookup(
        &self,
        region: &str,
        seq: &[u8],
        metric: DistanceMetric,
        partial: PartialMatching,
    ) -> Result<Option<LibraryMatch>, LibraryError> {
        // Try exact matching first - short circuit if we find the region
        if let Some(exact) = self.exact_matches.get(region) {
            if let Some(hit) = exact.get(seq) {
                return Ok(Some(LibraryMatch {
                    matches: vec![hit.clone()],
                    distance: 0,
                }));
            }
        }

        // Else try lookup
        let regions: &Vec<Arc<LibraryRegion>> = match self.regions.get(region) {
            Some(x) => x,
            None => {
                return Err(LibraryError::MissingRegion {
                    id: region.to_string(),
                });
            }
        };

        let max_dist = match self.region_max_distance.get(region) {
            None => self.default_max_distance,
            Some(x) => *x,
        };

        let (hits, best_dist) = match (metric, partial) {
            (DistanceMetric::Exact, PartialMatching::Full) => {
                // Already checked, if reached here no exact match
                return Ok(None);
            }
            (DistanceMetric::Exact, PartialMatching::FivePrimeOnly) => {
                // Exact FivePrimeOnly is the same as Hamming lookup on 5' end with 0 dist
                Self::lookup_hamming_5prime(seq, regions, 0)
            }
            (DistanceMetric::Exact, PartialMatching::ThreePrimeOnly) => {
                // Exact ThreePrimeOnly is the same as Hamming lookup on 3' end with 0 dist
                Self::lookup_hamming_3prime(seq, regions, 0)
            }

            (DistanceMetric::Hamming, PartialMatching::Full) => {
                Self::lookup_hamming(seq, regions, max_dist)
            }
            (DistanceMetric::Hamming, PartialMatching::FivePrimeOnly) => {
                Self::lookup_hamming_5prime(seq, regions, max_dist)
            }
            (DistanceMetric::Hamming, PartialMatching::ThreePrimeOnly) => {
                Self::lookup_hamming_3prime(seq, regions, max_dist)
            }

            (DistanceMetric::Levenshtein, PartialMatching::Full) => {
                Self::lookup_levenshtein(seq, regions, max_dist as u32)
            }
            (DistanceMetric::Levenshtein, PartialMatching::FivePrimeOnly) => {
                Self::lookup_levenshtein_5prime(seq, regions, max_dist as u32)
            }
            (DistanceMetric::Levenshtein, PartialMatching::ThreePrimeOnly) => {
                Self::lookup_levenshtein_3prime(seq, regions, max_dist as u32)
            }

            (DistanceMetric::BoundedLevenshtein, PartialMatching::Full) => {
                Self::lookup_bounded_levenshtein(seq, regions, max_dist as u32)
            }
            (DistanceMetric::BoundedLevenshtein, PartialMatching::FivePrimeOnly) => {
                Self::lookup_bounded_levenshtein_5prime(seq, regions, max_dist as u32)
            }
            (DistanceMetric::BoundedLevenshtein, PartialMatching::ThreePrimeOnly) => {
                Self::lookup_bounded_levenshtein_3prime(seq, regions, max_dist as u32)
            }
        };

        if hits.is_empty() {
            return Ok(None);
        }

        Ok(Some(LibraryMatch {
            matches: hits,
            distance: best_dist,
        }))
    }

    /// Compare an observed sequence to the library via Hamming distance
    fn lookup_hamming(
        seq: &[u8],
        regions: &[Arc<LibraryRegion>],
        max_dist: u64,
    ) -> (Vec<Arc<LibraryRegion>>, u64) {
        let mut dist: u64;
        let mut best_dist: u64 = u64::MAX;
        let mut hits: Vec<Arc<LibraryRegion>> = Vec::new();

        for reg in regions.iter() {
            // Hamming distance only applicaple for matching length, ignore
            // non-matching lengths
            if seq.len() != reg.sequence.len() {
                continue;
            }

            // Use appropriate hamming implemntation (other branch should be
            // pruned at compile time)
            if cfg!(all(target_feature = "avx2", target_feature = "sse4.1")) {
                dist = distance::simd::hamming(seq, &reg.sequence);
            } else {
                dist = distance::hamming(seq, &reg.sequence);
            }

            // Ignore too distant seqs - could make custom dist functions that short
            // circuit sooner to squeeze extra performance potentially
            if (dist > max_dist) || (dist > best_dist) {
                continue;
            } else if dist < best_dist {
                hits.clear();
                hits.push(reg.clone());
                best_dist = dist;
            } else if dist == best_dist {
                hits.push(reg.clone());
            }

            if best_dist == 0 {
                break;
            }
        }

        (hits, best_dist)
    }

    /// Compare an observed sequence to the library via Hamming distance at the library sequences
    /// 5 prime end
    fn lookup_hamming_5prime(
        seq: &[u8],
        regions: &[Arc<LibraryRegion>],
        max_dist: u64,
    ) -> (Vec<Arc<LibraryRegion>>, u64) {
        let mut dist: u64;
        let mut best_dist: u64 = u64::MAX;
        let mut hits: Vec<Arc<LibraryRegion>> = Vec::new();
        let query_len = seq.len();

        for reg in regions.iter() {
            // Hamming distance only defined for equal length - if query is
            // longer than region discard
            if query_len > reg.sequence.len() {
                continue;
            }

            // Use appropriate hamming implemntation (other branch should be
            // pruned at compile time)
            if cfg!(all(target_feature = "avx2", target_feature = "sse4.1")) {
                dist = distance::simd::hamming(seq, &reg.sequence[0..query_len]);
            } else {
                dist = distance::hamming(seq, &reg.sequence[0..query_len]);
            }

            // Ignore too distant seqs - could make custom dist functions that short
            // circuit sooner to squeeze extra performance potentially
            if (dist > max_dist) || (dist > best_dist) {
                continue;
            } else if dist < best_dist {
                hits.clear();
                hits.push(reg.clone());
                best_dist = dist;
            } else if dist == best_dist {
                hits.push(reg.clone());
            }
        }

        (hits, best_dist)
    }

    /// Compare an observed sequence to the library via Hamming distance at the library sequences
    /// 3 prime end
    fn lookup_hamming_3prime(
        seq: &[u8],
        regions: &[Arc<LibraryRegion>],
        max_dist: u64,
    ) -> (Vec<Arc<LibraryRegion>>, u64) {
        let mut dist: u64;
        let mut best_dist: u64 = u64::MAX;
        let mut hits: Vec<Arc<LibraryRegion>> = Vec::new();
        let query_len = seq.len();

        for reg in regions.iter() {
            let end = reg.sequence.len();
            let start = end.saturating_sub(query_len);

            // Hamming distance only defined for equal length - if query is
            // different length than region discard
            if end - start != seq.len() {
                continue
            }

            // Use appropriate hamming implemntation (other branch should be
            // pruned at compile time)
            if cfg!(all(target_feature = "avx2", target_feature = "sse4.1")) {
                dist = distance::simd::hamming(seq, &reg.sequence[start..end]);
            } else {
                dist = distance::hamming(seq, &reg.sequence[start..end]);
            }

            // Ignore too distant seqs - could make custom dist functions that short
            // circuit sooner to squeeze extra performance potentially
            if (dist > max_dist) || (dist > best_dist) {
                continue;
            } else if dist < best_dist {
                hits.clear();
                hits.push(reg.clone());
                best_dist = dist;
            } else if dist == best_dist {
                hits.push(reg.clone());
            }
        }

        (hits, best_dist)
    }

    /// Compare an observed sequence to the library via Levenshtein distance
    fn lookup_levenshtein(
        seq: &[u8],
        regions: &[Arc<LibraryRegion>],
        max_dist: u32,
    ) -> (Vec<Arc<LibraryRegion>>, u64) {
        let mut dist: u32;
        let mut best_dist: u32 = u32::MAX;
        let mut hits: Vec<Arc<LibraryRegion>> = Vec::new();

        for reg in regions.iter() {
            // Use appropriate levenshtein implemntation (other branch should be
            // pruned at compile time)
            if cfg!(all(target_feature = "avx2", target_feature = "sse4.1")) {
                dist = distance::simd::levenshtein(seq, &reg.sequence);
            } else {
                dist = distance::levenshtein(seq, &reg.sequence);
            }

            // Ignore too distant seqs - could make custom dist functions that short
            // circuit sooner to squeeze extra performance potentially
            if (dist > max_dist) || (dist > best_dist) {
                continue;
            } else if dist < best_dist {
                hits.clear();
                hits.push(reg.clone());
                best_dist = dist;
            } else if dist == best_dist {
                hits.push(reg.clone());
            }

            if best_dist == 0 {
                break;
            }
        }

        (hits, best_dist as u64)
    }

    /// Compare an observed sequence to the library via Levenshtein distance at the library
    /// sequences 5 prime end
    fn lookup_levenshtein_5prime(
        seq: &[u8],
        regions: &[Arc<LibraryRegion>],
        max_dist: u32,
    ) -> (Vec<Arc<LibraryRegion>>, u64) {
        let mut dist: u32;
        let mut best_dist: u32 = u32::MAX;
        let mut hits: Vec<Arc<LibraryRegion>> = Vec::new();
        let query_len = seq.len();

        for reg in regions.iter() {
            let reg_end = cmp::min(query_len, reg.sequence.len());

            // Use appropriate levenshtein implemntation (other branch should be
            // pruned at compile time)
            if cfg!(all(target_feature = "avx2", target_feature = "sse4.1")) {
                dist = distance::simd::levenshtein(seq, &reg.sequence[0..reg_end]);
            } else {
                dist = distance::levenshtein(seq, &reg.sequence[0..reg_end]);
            }

            // Ignore too distant seqs - could make custom dist functions that short
            // circuit sooner to squeeze extra performance potentially
            if (dist > max_dist) || (dist > best_dist) {
                continue;
            } else if dist < best_dist {
                hits.clear();
                hits.push(reg.clone());
                best_dist = dist;
            } else if dist == best_dist {
                hits.push(reg.clone());
            }
        }

        (hits, best_dist as u64)
    }

    /// Compare an observed sequence to the library via Levenshtein distance at the library
    /// sequences 3 prime end
    fn lookup_levenshtein_3prime(
        seq: &[u8],
        regions: &[Arc<LibraryRegion>],
        max_dist: u32,
    ) -> (Vec<Arc<LibraryRegion>>, u64) {
        let mut dist: u32;
        let mut best_dist: u32 = u32::MAX;
        let mut hits: Vec<Arc<LibraryRegion>> = Vec::new();
        let query_len = seq.len();

        for reg in regions.iter() {
            let end = reg.sequence.len();
            let start = end.saturating_sub(query_len);

            // Use appropriate levenshtein implemntation (other branch should be
            // pruned at compile time)
            if cfg!(all(target_feature = "avx2", target_feature = "sse4.1")) {
                dist = distance::simd::levenshtein(seq, &reg.sequence[start..end]);
            } else {
                dist = distance::levenshtein(seq, &reg.sequence[start..end]);
            }

            // Ignore too distant seqs - could make custom dist functions that short
            // circuit sooner to squeeze extra performance potentially
            if (dist > max_dist) || (dist > best_dist) {
                continue;
            } else if dist < best_dist {
                hits.clear();
                hits.push(reg.clone());
                best_dist = dist;
            } else if dist == best_dist {
                hits.push(reg.clone());
            }
        }

        (hits, best_dist as u64)
    }

    /// Compare an observed sequence to the library via Bounded Levenshtein distance
    fn lookup_bounded_levenshtein(
        seq: &[u8],
        regions: &[Arc<LibraryRegion>],
        max_dist: u32,
    ) -> (Vec<Arc<LibraryRegion>>, u64) {
        let mut dist: Option<u32>;
        let mut best_dist: u32 = u32::MAX;
        let mut hits: Vec<Arc<LibraryRegion>> = Vec::new();

        for reg in regions.iter() {
            dist = distance::simd::bounded_levenshtein(
                seq,
                &reg.sequence,
                cmp::min(best_dist, max_dist),
            );

            match dist {
                // Ignore cases where d is over k/max_dist
                None => continue,
                Some(d) if d < best_dist => {
                    hits.clear();
                    hits.push(reg.clone());
                    best_dist = d;
                }
                Some(d) if d == best_dist => hits.push(reg.clone()),
                Some(_) => continue,
            }

            if best_dist == 0 {
                break;
            }
        }

        (hits, best_dist as u64)
    }

    /// Compare an observed sequence to the library via Bounded Levenshtein distance at the library
    /// sequences 5 prime end
    fn lookup_bounded_levenshtein_5prime(
        seq: &[u8],
        regions: &[Arc<LibraryRegion>],
        max_dist: u32,
    ) -> (Vec<Arc<LibraryRegion>>, u64) {
        let mut dist: Option<u32>;
        let mut best_dist: u32 = u32::MAX;
        let mut hits: Vec<Arc<LibraryRegion>> = Vec::new();
        let query_len = seq.len();

        for reg in regions.iter() {
            let reg_end = cmp::min(query_len, reg.sequence.len());

            dist = distance::simd::bounded_levenshtein(
                seq,
                &reg.sequence[0..reg_end],
                cmp::min(best_dist, max_dist),
            );

            match dist {
                // Ignore cases where d is over k/max_dist
                None => continue,
                Some(d) if d < best_dist => {
                    hits.clear();
                    hits.push(reg.clone());
                    best_dist = d;
                }
                Some(d) if d == best_dist => hits.push(reg.clone()),
                Some(_) => continue,
            }
        }

        (hits, best_dist as u64)
    }

    /// Compare an observed sequence to the library via Bounded Levenshtein distance at the library
    /// sequences 3 prime end
    fn lookup_bounded_levenshtein_3prime(
        seq: &[u8],
        regions: &[Arc<LibraryRegion>],
        max_dist: u32,
    ) -> (Vec<Arc<LibraryRegion>>, u64) {
        let mut dist: Option<u32>;
        let mut best_dist: u32 = u32::MAX;
        let mut hits: Vec<Arc<LibraryRegion>> = Vec::new();
        let query_len = seq.len();

        for reg in regions.iter() {
            let end = reg.sequence.len();
            let start = end.saturating_sub(query_len);

            dist = distance::simd::bounded_levenshtein(
                seq,
                &reg.sequence[start..end],
                cmp::min(best_dist, max_dist),
            );

            match dist {
                // Ignore cases where d is over k/max_dist
                None => continue,
                Some(d) if d < best_dist => {
                    hits.clear();
                    hits.push(reg.clone());
                    best_dist = d;
                }
                Some(d) if d == best_dist => hits.push(reg.clone()),
                Some(_) => continue,
            }
        }

        (hits, best_dist as u64)
    }

    /// Initialise a Library from a LibrarySpec object
    pub fn from_lib_spec(
        lib_spec: &LibrarySpec,
        default_max_distance: u64,
    ) -> Result<Option<Library>, LibraryError> {
        let spec_regions = lib_spec.variable_regions();
        let region_max_distance = lib_spec.get_max_distances();

        let lib: Library = match &lib_spec.library {
            None => return Ok(None),
            Some(x) => Library::from_file(x, region_max_distance, default_max_distance)?,
        };

        if !lib.library.keys().all(|x| spec_regions.contains(x)) {
            return Err(LibraryError::Library {
                desc: "Library region ids don't match variable LibSpec region ids".to_string(),
            });
        }

        Ok(Some(lib))
    }

    /// Import a Library from a TSV file
    pub fn from_file(
        path: &str,
        region_max_distance: HashMap<String, u64>,
        default_max_distance: u64,
    ) -> Result<Library, LibraryError> {
        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .comment(Some(b'#'))
            .trim(csv::Trim::All)
            .from_path(path)?;

        // Prepare a HashMap to store column name to values
        let names = reader.headers()?.clone();
        let mut regions: HashMap<String, Vec<Sequence>> = HashMap::new();

        // Initialize empty Vec<Sequence> for each column
        for name in names.iter() {
            regions.insert(name.to_string(), Vec::new());
        }

        // Iterate through records, appending values to the respective columns
        for result in reader.records() {
            let record = result?;
            for (name, seq) in names.iter().zip(record.iter()) {
                match regions.get_mut(name) {
                    None => {
                        return Err(LibraryError::MissingRegion {
                            id: name.to_string(),
                        });
                    }
                    Some(v) => v.push(seq.as_bytes().to_vec()),
                }
            }
        }

        Library::new(regions, region_max_distance, default_max_distance)
    }
}

/// Distance metric types
#[derive(Clone, ValueEnum, Debug, Copy)]
pub enum DistanceMetric {
    Exact,
    Hamming,
    Levenshtein,
    BoundedLevenshtein,
}

/// Partial matching options
///
/// Refers to the end of the library sequence to include - so a query
/// that is trunctated at the 3 prime end would use FivePrimeOnly.
#[derive(Clone, Debug, Copy)]
pub enum PartialMatching {
    Full,
    FivePrimeOnly,
    ThreePrimeOnly,
}

#[cfg(test)]
mod tests {
    // use super::*;

}
