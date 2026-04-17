//! Compiled expected-sequence libraries and lookup logic.
//!
//! This module defines:
//! - [`Library`], which coordinates lookup across one or more independent
//!   sublibraries,
//! - [`SubLibrary`], which stores expected sequences for a specific subset of
//!   variable regions,
//! - and the distance-based matching logic used to compare observed region
//!   sequences to those expected libraries.
//!
//! A single library TSV defines one `SubLibrary`. Multiple TSVs can be combined
//! into a top-level `Library` to represent combinatorial designs where each TSV
//! constrains one independent subset of regions.
use bio::alignment::distance;
use bio::bio_types::sequence::Sequence;
use clap::ValueEnum;
use csv::ReaderBuilder;
use std::cmp;
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use crate::errors::LibraryError;
use crate::interning::{
    LibraryID, RegionID, SeqHandle, library_id_from_str, region_id_from_str, seq_from_bytes,
    seq_to_bytes,
};
use crate::lib_spec::LibrarySpec;

/// Compiled expected-sequence library spanning one or more sublibraries.
///
/// A `Library` dispatches region lookups to the appropriate [`SubLibrary`] based
/// on region ID. This allows independent library TSVs to define separate parts of
/// a combinatorial design while presenting a single lookup interface.
#[derive(Debug)]
pub struct Library {
    pub regions: HashMap<RegionID, usize>,
    pub sublibraries: Vec<SubLibrary>,
}

impl Library {
    /// Construct a top-level library from compiled sublibraries.
    ///
    /// Each region may belong to at most one sublibrary. An error is returned if
    /// the same region appears in multiple sublibraries.
    pub fn new(sublibraries: Vec<SubLibrary>) -> Result<Self, LibraryError> {
        let mut regions = HashMap::new();

        for (i, lib) in sublibraries.iter().enumerate() {
            let new_regions = lib.regions();

            for r in new_regions {
                if regions.insert(*r, i).is_some() {
                    return Err(LibraryError::DuplicateSubLibraryRegion { id: *r });
                }
            }
        }

        Ok(Self {
            regions,
            sublibraries,
        })
    }

    /// Return true if the library contains no SubLibraries
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.sublibraries.is_empty()
    }

    /// Return the sublibrary index responsible for a given region ID.
    pub fn get_sublibrary_index(&self, region: &RegionID) -> Result<usize, LibraryError> {
        match self.regions.get(region) {
            Some(x) => Ok(*x),
            None => Err(LibraryError::MissingRegion { id: *region }),
        }
    }

    /// Look up an observed sequence against the expected library for one region.
    ///
    /// The region ID determines which sublibrary should be queried. Matching is
    /// then delegated to that sublibrary using the requested distance metric and
    /// partial-matching mode.
    ///
    /// Returns:
    /// - `Ok(Some(...))` for one or more best matches within the allowed distance,
    /// - `Ok(None)` if no acceptable match is found,
    /// - `Err(...)` if the region is not represented in the library.
    pub fn lookup(
        &self,
        region: &RegionID,
        seq: SeqHandle,
        metric: DistanceMetric,
        partial: PartialMatching,
    ) -> Result<Option<LibraryMatch>, LibraryError> {
        self.sublibraries[self.get_sublibrary_index(region)?].lookup(region, seq, metric, partial)
    }

    /// Build a top-level library from one or more library TSV files.
    ///
    /// Each TSV becomes one independent [`SubLibrary`]. Together these define a
    /// potentially combinatorial expected design, where each region must appear in
    /// at most one input file.
    pub fn from_files(
        paths: &[String],
        lib_spec: &LibrarySpec,
        default_max_distance: u64,
    ) -> Result<Library, LibraryError> {
        let libs: Result<Vec<SubLibrary>, LibraryError> = paths
            .iter()
            .enumerate()
            .map(|(i, x)| {
                SubLibrary::from_file_with_lib_spec(
                    x,
                    lib_spec,
                    default_max_distance,
                    Some(format!("lib{i}")),
                )
            })
            .collect();

        Library::new(libs?)
    }
}

/// Compiled expected-sequence library for a specific subset of variable regions.
///
/// A `SubLibrary` stores the expected sequence combinations from one library TSV
/// and supports efficient lookup of observed sequences against the unique
/// sequences present in each region.
#[derive(Debug)]
pub struct SubLibrary {
    /// Full sequences for each member of the library, divided into region vectors. The full nth
    /// sequence contains the nth sequence from each region vector
    pub library: HashMap<RegionID, Vec<Arc<LibraryRegion>>>,

    /// Unique sequences for each region, mapping back to which full combinations they are part
    /// of by index plus a copy of the real sequence to avoid hitting hte interning during the hot loop
    regions: HashMap<RegionID, Vec<LibrarySequence>>,

    /// Library member IDs
    pub ids: Vec<LibraryID>,

    /// HashMap of exact hits to Library regions for quick initial lookup and
    /// exact matching
    exact_matches: HashMap<RegionID, HashMap<SeqHandle, Arc<LibraryRegion>>>,

    /// Max distance to consider for each region
    region_max_distance: HashMap<RegionID, u64>,

    /// Default max distance to consider
    default_max_distance: u64,
}

/// A LibraryRegion and it's paired real sequence
///
/// Storing the sequence raw avoids hitting the interner during the hot loop
#[derive(Debug, Eq, PartialEq, Clone)]
struct LibrarySequence {
    region: Arc<LibraryRegion>,
    sequence: Sequence,
}

impl LibrarySequence {
    fn from_region(region: Arc<LibraryRegion>) -> Self {
        Self {
            sequence: seq_to_bytes(region.sequence).to_vec(),
            region,
        }
    }
}

/// One unique expected sequence for a region within a compiled sublibrary.
///
/// A `LibraryRegion` stores the canonical sequence plus the set of library-member
/// indices and IDs that contain that sequence. This allows repeated identical
/// region sequences across multiple library members to be represented once.
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct LibraryRegion {
    /// The sequence
    pub sequence: SeqHandle,

    /// The indeces of library members that contain this sequence
    pub inds: HashSet<usize>,

    /// The library members that contain this sequence
    pub ids: HashSet<LibraryID>,
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

/// Best-match result for a region lookup.
///
/// Contains all equally good best matches and the distance shared by those matches.
#[derive(Debug)]
pub struct LibraryMatch {
    pub matches: Vec<Arc<LibraryRegion>>,
    pub distance: u64,
}

/// Intersect two partial region-match results.
///
/// This is used when a region is observed only as two partial pieces, for example
/// from opposite ends of a read pair. Only library-region matches consistent with
/// both partial observations are retained.
///
/// Distances are summed across the two partial matches. This is appropriate when
/// the partial observations cover distinct parts of the region, but may
/// overcount if they overlap.
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

impl SubLibrary {
    /// Compile a sublibrary from per-region expected sequences.
    ///
    /// The input `library` maps each region ID to the full column of expected
    /// sequences from one library TSV. All region vectors must have the same
    /// length, representing the same ordered library members.
    ///
    /// During compilation this:
    /// - validates dimensions and IDs,
    /// - assigns library-member IDs,
    /// - deduplicates identical region sequences,
    /// - builds exact-match lookup maps,
    /// - and reconstructs full per-member region assignments.
    pub fn new(
        library: HashMap<RegionID, Vec<Sequence>>,
        ids: Option<Vec<String>>,
        region_max_distance: HashMap<RegionID, u64>,
        default_max_distance: u64,
        default_id: Option<String>,
    ) -> Result<SubLibrary, LibraryError> {
        let mut regions: HashMap<RegionID, Vec<LibrarySequence>> = HashMap::new();

        let mut exp_len: usize = 0;

        // This allows an empty library
        if !library.is_empty() {
            exp_len = library
                .values()
                .next()
                .expect("Just checked library has at least 1 elements")
                .len();

            if !library.values().all(|x| x.len() == exp_len) {
                return Err(LibraryError::Library {
                    desc: "Library must contain the same number of sequences for each region"
                        .to_string(),
                });
            }

            if let Some(i) = &ids {
                if i.len() != exp_len {
                    return Err(LibraryError::Library {
                        desc: "Library must have as many IDs as the number of elements".to_string(),
                    });
                }
            }
        }

        let baseid = default_id.unwrap_or("seq".to_string());
        let lib_ids: Vec<LibraryID> = match ids {
            Some(ids) => ids.iter().map(|x| library_id_from_str(x)).collect(),
            None => (0..exp_len)
                .map(|x| library_id_from_str(&format!("{baseid}_{x}")))
                .collect(),
        };

        // Check if Libary IDs are unique
        if lib_ids.iter().collect::<HashSet<_>>().len() != exp_len {
            return Err(LibraryError::Library {
                desc: "Library IDs are not unique".to_string(),
            });
        }

        for key in library.keys() {
            let mut reg_map: HashMap<SeqHandle, HashSet<usize>> = HashMap::new();
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
                match reg_map.get_mut(&seq_from_bytes(seq)) {
                    Some(x) => {
                        x.insert(ind);
                    }
                    None => {
                        reg_map.insert(seq_from_bytes(seq), HashSet::from([ind]));
                    }
                }
            }

            if regions
                .insert(
                    *key,
                    Vec::from_iter(reg_map.into_iter().map(|x| {
                        LibrarySequence::from_region(Arc::new(LibraryRegion {
                            sequence: x.0,
                            ids: x.1.iter().map(|x| lib_ids[*x]).collect(),
                            inds: x.1,
                        }))
                    })),
                )
                .is_some()
            {
                return Err(LibraryError::DuplicateRegion { id: *key });
            }
        }

        // Compile Exact matches HashMap
        let mut exact_matches: HashMap<RegionID, HashMap<SeqHandle, Arc<LibraryRegion>>> =
            HashMap::new();
        for key in regions.keys() {
            exact_matches.insert(*key, HashMap::new());
            for reg in regions
                .get(key)
                .expect("Key known to be in regions HashMap")
            {
                exact_matches
                    .get_mut(key)
                    .expect("Key just added to exact_matchs")
                    .insert(reg.region.sequence, reg.region.clone());
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
                for ind in &reg.region.inds {
                    rc_vec[*ind] = Some(reg.region.clone());
                }
            }

            library_compiled.insert(
                key,
                rc_vec.into_iter().map(
                    |x| x.expect("All library members should have been assigned Some(Arc<LibraryRegion>) by construction")).collect()
            );
        }

        Ok(SubLibrary {
            library: library_compiled,
            regions,
            ids: lib_ids,
            exact_matches,
            region_max_distance,
            default_max_distance,
        })
    }

    /// Number of library members in this sublibrary.
    #[allow(dead_code)]
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

    /// Return `true` if this sublibrary contains no library members.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.library.is_empty()
    }

    /// Return the region IDs represented in this sublibrary.
    pub fn regions(&self) -> Vec<&RegionID> {
        self.regions.keys().collect()
    }

    /// Look up an observed sequence against one region of this sublibrary.
    ///
    /// Matching proceeds in two stages:
    /// 1. exact-match lookup is attempted first,
    /// 2. if needed, the requested distance metric is used to find the best
    ///    sequence(s) within the configured maximum distance.
    ///
    /// `partial` controls whether the query must match the full library sequence
    /// or only one end of it. This is used for incompletely observed regions,
    /// such as truncation at the 5' or 3' end of a read.
    ///
    /// Returns the best match set and its distance, or `None` if no match lies
    /// within the allowed threshold.
    ///
    /// `DistanceMetric::BoundedLevenshtein` uses the configured maximum distance
    /// as an upper bound and is generally equivalent to Levenshtein for accepted
    /// matches while often being faster.
    ///
    /// The implementation dispatches to the appropriate lookup funciton based on distance matric
    /// and partial match type. We use the SIMD optimised versions of Hamming and Levenshtein
    /// distance if AVX2 and SSE4.1 are available at compile time. The Rust Bio SIMD distance
    /// metrics fall back to standard versions if SIMD isn't available so this is just a minor
    /// optimisation and either version should be portable without undefined behaviour. Bounded
    /// levenshtein is only available as a SIMD implementation with fallback so we rely on
    /// the Rust Bio and editdistancek authors for the check.
    pub fn lookup(
        &self,
        region: &RegionID,
        seq: SeqHandle,
        metric: DistanceMetric,
        partial: PartialMatching,
    ) -> Result<Option<LibraryMatch>, LibraryError> {
        // Try exact matching first - short circuit if we find the region
        if let Some(exact) = self.exact_matches.get(region) {
            if let Some(hit) = exact.get(&seq) {
                return Ok(Some(LibraryMatch {
                    matches: vec![hit.clone()],
                    distance: 0,
                }));
            }
        }

        // Else try lookup
        let regions: &Vec<LibrarySequence> = match self.regions.get(region) {
            Some(x) => x,
            None => {
                return Err(LibraryError::MissingRegion { id: *region });
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
                Self::lookup_hamming_5prime(seq_to_bytes(seq).as_ref(), regions, 0)
            }
            (DistanceMetric::Exact, PartialMatching::ThreePrimeOnly) => {
                // Exact ThreePrimeOnly is the same as Hamming lookup on 3' end with 0 dist
                Self::lookup_hamming_3prime(seq_to_bytes(seq).as_ref(), regions, 0)
            }

            (DistanceMetric::Hamming, PartialMatching::Full) => {
                Self::lookup_hamming(seq_to_bytes(seq).as_ref(), regions, max_dist)
            }
            (DistanceMetric::Hamming, PartialMatching::FivePrimeOnly) => {
                Self::lookup_hamming_5prime(seq_to_bytes(seq).as_ref(), regions, max_dist)
            }
            (DistanceMetric::Hamming, PartialMatching::ThreePrimeOnly) => {
                Self::lookup_hamming_3prime(seq_to_bytes(seq).as_ref(), regions, max_dist)
            }

            (DistanceMetric::Levenshtein, PartialMatching::Full) => {
                Self::lookup_levenshtein(seq_to_bytes(seq).as_ref(), regions, max_dist as u32)
            }
            (DistanceMetric::Levenshtein, PartialMatching::FivePrimeOnly) => {
                Self::lookup_levenshtein_5prime(
                    seq_to_bytes(seq).as_ref(),
                    regions,
                    max_dist as u32,
                )
            }
            (DistanceMetric::Levenshtein, PartialMatching::ThreePrimeOnly) => {
                Self::lookup_levenshtein_3prime(
                    seq_to_bytes(seq).as_ref(),
                    regions,
                    max_dist as u32,
                )
            }

            (DistanceMetric::BoundedLevenshtein, PartialMatching::Full) => {
                Self::lookup_bounded_levenshtein(
                    seq_to_bytes(seq).as_ref(),
                    regions,
                    max_dist as u32,
                )
            }
            (DistanceMetric::BoundedLevenshtein, PartialMatching::FivePrimeOnly) => {
                Self::lookup_bounded_levenshtein_5prime(
                    seq_to_bytes(seq).as_ref(),
                    regions,
                    max_dist as u32,
                )
            }
            (DistanceMetric::BoundedLevenshtein, PartialMatching::ThreePrimeOnly) => {
                Self::lookup_bounded_levenshtein_3prime(
                    seq_to_bytes(seq).as_ref(),
                    regions,
                    max_dist as u32,
                )
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
        regions: &[LibrarySequence],
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
                hits.push(reg.region.clone());
                best_dist = dist;
            } else if dist == best_dist {
                hits.push(reg.region.clone());
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
        regions: &[LibrarySequence],
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
                hits.push(reg.region.clone());
                best_dist = dist;
            } else if dist == best_dist {
                hits.push(reg.region.clone());
            }
        }

        (hits, best_dist)
    }

    /// Compare an observed sequence to the library via Hamming distance at the library sequences
    /// 3 prime end
    fn lookup_hamming_3prime(
        seq: &[u8],
        regions: &[LibrarySequence],
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
                continue;
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
                hits.push(reg.region.clone());
                best_dist = dist;
            } else if dist == best_dist {
                hits.push(reg.region.clone());
            }
        }

        (hits, best_dist)
    }

    /// Compare an observed sequence to the library via Levenshtein distance
    fn lookup_levenshtein(
        seq: &[u8],
        regions: &[LibrarySequence],
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
                hits.push(reg.region.clone());
                best_dist = dist;
            } else if dist == best_dist {
                hits.push(reg.region.clone());
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
        regions: &[LibrarySequence],
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
                hits.push(reg.region.clone());
                best_dist = dist;
            } else if dist == best_dist {
                hits.push(reg.region.clone());
            }
        }

        (hits, best_dist as u64)
    }

    /// Compare an observed sequence to the library via Levenshtein distance at the library
    /// sequences 3 prime end
    fn lookup_levenshtein_3prime(
        seq: &[u8],
        regions: &[LibrarySequence],
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
                hits.push(reg.region.clone());
                best_dist = dist;
            } else if dist == best_dist {
                hits.push(reg.region.clone());
            }
        }

        (hits, best_dist as u64)
    }

    /// Compare an observed sequence to the library via Bounded Levenshtein distance
    fn lookup_bounded_levenshtein(
        seq: &[u8],
        regions: &[LibrarySequence],
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
                    hits.push(reg.region.clone());
                    best_dist = d;
                }
                Some(d) if d == best_dist => hits.push(reg.region.clone()),
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
        regions: &[LibrarySequence],
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
                    hits.push(reg.region.clone());
                    best_dist = d;
                }
                Some(d) if d == best_dist => hits.push(reg.region.clone()),
                Some(_) => continue,
            }
        }

        (hits, best_dist as u64)
    }

    /// Compare an observed sequence to the library via Bounded Levenshtein distance at the library
    /// sequences 3 prime end
    fn lookup_bounded_levenshtein_3prime(
        seq: &[u8],
        regions: &[LibrarySequence],
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
                    hits.push(reg.region.clone());
                    best_dist = d;
                }
                Some(d) if d == best_dist => hits.push(reg.region.clone()),
                Some(_) => continue,
            }
        }

        (hits, best_dist as u64)
    }

    /// Build a sublibrary from a TSV file and validate it against a `LibrarySpec`.
    ///
    /// This checks that:
    /// - the reserved `_id` name is not used as a variable region in the spec,
    /// - all TSV region columns correspond to variable regions in the spec,
    /// - and per-region max-distance defaults are taken from the spec where present.
    pub fn from_file_with_lib_spec(
        path: &str,
        lib_spec: &LibrarySpec,
        default_max_distance: u64,
        default_id: Option<String>,
    ) -> Result<SubLibrary, LibraryError> {
        let spec_regions = lib_spec.variable_regions();
        let region_max_distance = lib_spec.get_max_distances();

        if lib_spec
            .variable_regions()
            .contains(&region_id_from_str("_id"))
        {
            return Err(LibraryError::Library {
                desc: "Region named '_id'. This is reserved for element names when doing library comparison".to_string(),
            });
        }

        let lib: SubLibrary =
            SubLibrary::from_file(path, region_max_distance, default_max_distance, default_id)?;

        if !lib.library.keys().all(|x| spec_regions.contains(x)) {
            return Err(LibraryError::Library {
                desc: "Library region ids don't match variable LibSpec region ids".to_string(),
            });
        }

        Ok(lib)
    }

    /// Build a sublibrary directly from a library TSV file.
    ///
    /// The TSV must contain one column per region and one row per library member.
    /// An optional `_id` column supplies library-member names; otherwise names are
    /// generated automatically from `default_id`.
    pub fn from_file(
        path: &str,
        region_max_distance: HashMap<RegionID, u64>,
        default_max_distance: u64,
        default_id: Option<String>,
    ) -> Result<SubLibrary, LibraryError> {
        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .comment(Some(b'#'))
            .trim(csv::Trim::All)
            .from_path(path)?;

        // Prepare a HashMap to store column name to values
        let names = reader.headers()?.clone();
        let mut id_vec: Vec<String> = Vec::new();
        let mut regions: HashMap<RegionID, Vec<Sequence>> = HashMap::new();

        // Initialize empty Vec<Sequence> for each column
        for name in names.iter() {
            if name != "_id" {
                regions.insert(region_id_from_str(name), Vec::new());
            }
        }

        // Iterate through records, appending values to the respective columns
        for result in reader.records() {
            let record = result?;
            for (name, val) in names.iter().zip(record.iter()) {
                if name == "_id" {
                    id_vec.push(val.to_string());
                    continue;
                }

                match regions.get_mut(&region_id_from_str(name)) {
                    None => {
                        return Err(LibraryError::MissingRegion {
                            id: region_id_from_str(name),
                        });
                    }
                    Some(v) => v.push(val.as_bytes().to_vec()),
                }
            }
        }

        let ids = if !id_vec.is_empty() {
            Some(id_vec)
        } else {
            None
        };

        SubLibrary::new(
            regions,
            ids,
            region_max_distance,
            default_max_distance,
            default_id,
        )
    }
}

/// Distance metric used for library lookup.
#[derive(Clone, ValueEnum, Debug, Copy)]
pub enum DistanceMetric {
    /// Require exact sequence equality.
    Exact,

    /// Count substitutions only; lengths must match for full matching.
    Hamming,

    /// Full edit distance allowing substitutions, insertions, and deletions.
    Levenshtein,

    /// Edit distance capped at the configured maximum threshold.
    ///
    /// Usually equivalent to Levenshtein for accepted matches, but often faster.
    BoundedLevenshtein,
}

/// How an observed query sequence should be aligned against a library region.
///
/// Partial matching is used for truncated observed regions. The enum names refer
/// to which end of the *library sequence* must be matched by the query.
#[derive(Clone, Debug, Copy)]
pub enum PartialMatching {
    /// Require the full query to match the full library region.
    Full,

    /// Match the query to the 5' end of the library sequence.
    FivePrimeOnly,

    /// Match the query to the 3' end of the library sequence.
    ThreePrimeOnly,
}

#[cfg(test)]
mod tests {
    // use super::*;
}
