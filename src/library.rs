//! Specification for DNA constructs and libraries
//!
//! Provides methdods for importing JSON based DNA construct specifications
//! and manipulating them. Additionally supports TSV libraries corresponding
//! to these constructs, with lookup capabilities.
use bio::alignment::distance;
use bio::bio_types::sequence::Sequence;
use clap::ValueEnum;
use csv::ReaderBuilder;
use std::cmp;
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use crate::errors::LibraryError;
use crate::interning::{LibraryID, RegionID, library_id_from_str, region_id_from_str};
use crate::lib_spec::LibrarySpec;

/// A compiled sequence library
///
/// Dispatces to sub-libraries for matching each subset of regions
#[derive(Debug)]
pub struct Library {
    pub regions: HashMap<RegionID, usize>,
    pub sublibraries: Vec<SubLibrary>,
}

impl Library {
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

    #[allow(dead_code)]
    /// Check if the library is empty
    pub fn is_empty(&self) -> bool {
        self.sublibraries.is_empty()
    }

    pub fn get_sublibrary_index(&self, region: &RegionID) -> Result<usize, LibraryError> {
        match self.regions.get(region) {
            Some(x) => Ok(*x),
            None => Err(LibraryError::MissingRegion { id: *region }),
        }
    }

    /// Compare an observed sequence to the library
    ///
    /// Itentify the subpool a candidate is from and dispatch to the appropriate SubLibrary
    pub fn lookup(
        &self,
        region: &RegionID,
        seq: &[u8],
        metric: DistanceMetric,
        partial: PartialMatching,
    ) -> Result<Option<LibraryMatch>, LibraryError> {
        self.sublibraries[self.get_sublibrary_index(region)?].lookup(region, seq, metric, partial)
    }

    /// Import a Library from a series of TSV files
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

/// A compiled sequence sub-library, carrying the expected sequence combinations in each variable region
///
/// This allows efficient lookup of candidate sequences against the library
#[derive(Debug)]
pub struct SubLibrary {
    /// Full sequences for each member of the library, divided into region vectors. The full nth
    /// sequence contains the nth sequence from each region vector
    pub library: HashMap<RegionID, Vec<Arc<LibraryRegion>>>,

    /// Unique sequences for each region, mapping back to which full combinations they are part
    /// of by index
    pub regions: HashMap<RegionID, Vec<Arc<LibraryRegion>>>,

    /// Library member IDs
    pub ids: Vec<LibraryID>,

    /// HashMap of exact hits to Library regions for quick initial lookup and
    /// exact matching
    exact_matches: HashMap<RegionID, HashMap<Sequence, Arc<LibraryRegion>>>,

    /// Max distance to consider for each region
    region_max_distance: HashMap<RegionID, u64>,

    /// Default max distance to consider
    default_max_distance: u64,
}

/// A sequence region from a compiled library
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct LibraryRegion {
    /// The `Vec<u8>` sequence
    pub sequence: Sequence,

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

impl SubLibrary {
    pub fn new(
        library: HashMap<RegionID, Vec<Sequence>>,
        ids: Option<Vec<String>>,
        region_max_distance: HashMap<RegionID, u64>,
        default_max_distance: u64,
        default_id: Option<String>,
    ) -> Result<SubLibrary, LibraryError> {
        let mut regions = HashMap::new();

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
                    *key,
                    Vec::from_iter(reg_map.into_iter().map(|x| {
                        Arc::new(LibraryRegion {
                            sequence: x.0,
                            ids: x.1.iter().map(|x| lib_ids[*x]).collect(),
                            inds: x.1,
                        })
                    })),
                )
                .is_some()
            {
                return Err(LibraryError::DuplicateRegion { id: *key });
            }
        }

        // Compile Exact matches HashMap
        let mut exact_matches = HashMap::new();
        for key in regions.keys() {
            exact_matches.insert(*key, HashMap::new());
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

    /// Get the names of regions in the sublibrary
    pub fn regions(&self) -> Vec<&RegionID> {
        self.regions.keys().collect()
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
        region: &RegionID,
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

    /// Import a Library from a TSV file
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
