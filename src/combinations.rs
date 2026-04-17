//! Container and aggregation logic for observed region combinations.
//!
//! This module defines [`ObservedCombinations`], the central data structure used
//! during counting. It accumulates distinct observed combinations, deduplicates
//! regions, tracks read counts and filtered reads, and provides methods for
//! merging results, performing library comparison, and generating summary views
//! for downstream output.
use anyhow::{self};
use bio::bio_types::alignment::Alignment;
use crossbeam::channel::{Receiver, unbounded};
use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use std::thread::scope;

use crate::combination::{CombinationKey, CombinationMatch, ObservedCombination};
use crate::errors::{LibraryError, ReadCountError};
use crate::filters::{FilterConfig, FilterReason, FilteredCounts, FilteredReads};
use crate::groups::ReadGroup;
use crate::interning::RegionID;
use crate::library::{DistanceMetric, Library};
use crate::library_combination::{LibraryCombination, LibraryCombinationKey, LibraryRegionMatch};
use crate::logging::{Progress, ProgressStyle};
use crate::region::{ObservedRegion, RegionKey};
use crate::seqs::ReadPair;
use crate::seqs::SeqPair;

/// Container for all distinct `ObservedCombination` objects in a dataset
///
/// This is the main accumulation structure used during counting. It stores:
/// - the ordered list of variable region IDs expected from the `LibSpec`,
/// - deduplicated `ObservedRegion` objects keyed by sequence and completeness,
/// - deduplicated `ObservedCombination` objects keyed by their constituent regions, with
///   the contained ObservedRegions linking back to the region container.
/// - optional library-comparison results,
/// - filtered read counts,
/// - and an optional read-level cache used to avoid repeating extraction work.
///
/// The container is designed for incremental counting followed by an optional
/// library-comparison step and final output summarisation.
#[derive(Debug)]
pub struct ObservedCombinations {
    pub region_ids: Vec<RegionID>,
    regions: HashMap<RegionKey, Arc<Mutex<ObservedRegion>>>,
    combinations: HashMap<CombinationKey, ObservedCombination>,
    library: Option<Library>,
    library_combinations: Option<HashMap<LibraryCombinationKey, LibraryCombination>>,
    filtered_reads: FilteredReads,
    cache: ObservedReads,
}

impl ObservedCombinations {
    /// Create an empty `ObservedCombinations` container for a specific ordered
    /// set of variable region IDs and filtering configuration.
    pub fn new(region_ids: Vec<RegionID>, filter_config: FilterConfig) -> Self {
        Self {
            region_ids,
            regions: HashMap::new(),
            combinations: HashMap::new(),
            library: None,
            library_combinations: None,
            filtered_reads: FilteredReads::new(filter_config),
            cache: HashMap::new(),
        }
    }

    /// Borrow the filtered-read summary accumulated during counting.
    pub fn filtered_reads(&self) -> &FilteredReads {
        &self.filtered_reads
    }

    /// Merge another `ObservedCombinations` into this one.
    ///
    /// This is primarily used to combine per-thread counting results. Regions and
    /// combinations are deduplicated into the receiving container, and counts are
    /// summed across matching read groups.
    ///
    /// Merge is only valid before library comparison has been performed, because
    /// library-comparison results depend on shared region state and are not
    /// currently merged.
    ///
    /// Returns an error if:
    /// - either container has already been compared to a library,
    /// - the ordered `region_ids` differ,
    /// - or the filtered-read configurations are incompatible.
    pub fn merge(&mut self, new_counts: ObservedCombinations) -> Result<(), ReadCountError> {
        if self.is_compared_to_library() || new_counts.is_compared_to_library() {
            return Err(ReadCountError::Error {
                desc: "Can't merge ObservedCombinations once library comparison has been run"
                    .to_string(),
            });
        }

        if self.region_ids != new_counts.region_ids {
            return Err(ReadCountError::Error {
                desc: "Can't merge ObservedCombinations with different region_ids".to_string(),
            });
        }

        self.filtered_reads.merge(new_counts.filtered_reads)?;

        for (k, v) in new_counts.regions.iter() {
            if !self.regions.contains_key(k) {
                self.regions.insert(k.clone(), v.to_owned());
            }
        }

        for (comb_key, mut new_comb) in new_counts.combinations.into_iter() {
            match self.combinations.get_mut(&comb_key) {
                Some(old_comb) => {
                    for (group, new_count) in &new_comb.counts {
                        match old_comb.counts.get_mut(group) {
                            Some(old_count) => *old_count += new_count,
                            None => {
                                old_comb.counts.insert(*group, *new_count);
                            }
                        }
                    }
                }
                None => {
                    // If combination isn't in the map we need to update it's Arc<Mutex<ObservedRegion>>
                    // references to point to the internal regions

                    // Clear old regions - the information is in the key and the region objects are merged already
                    new_comb.regions.clear();

                    for reg_key in &comb_key.regions {
                        let arc = match self.regions.get(reg_key) {
                            Some(x) => x,
                            None => return Err(ReadCountError::Error {
                                desc: "Region key missing during combination merge after merging regions".to_string(),
                            }),
                        };

                        new_comb.regions.insert(reg_key.id, arc.clone());
                    }

                    self.combinations.insert(comb_key, new_comb);
                }
            }
        }

        Ok(())
    }

    /// Number of distinct observed combination types
    pub fn len(&self) -> usize {
        self.combinations.len()
    }

    /// Check if there are no observed combinations
    pub fn is_empty(&self) -> bool {
        self.combinations.is_empty()
    }

    /// Total number of filtered reads
    pub fn total_filtered(&self) -> u64 {
        self.filtered_reads.total()
    }

    /// Add a new observed combination or increment an existing one.
    ///
    /// Region keys are resolved against the container-wide deduplicated region map,
    /// creating new `ObservedRegion` objects only when a region sequence/completeness
    /// combination has not been seen before.
    pub fn add_or_increment_combination(
        &mut self,
        comb_key: &CombinationKey,
        group: ReadGroup,
    ) -> Result<(), anyhow::Error> {
        match self.combinations.get_mut(comb_key) {
            Some(comb) => comb.increment_count(group),
            None => {
                let mut reg_map = HashMap::new();

                for reg_key in &comb_key.regions {
                    if !self.region_ids.contains(&reg_key.id) {
                        return Err(ReadCountError::UnexpectedRegion { region: reg_key.id }.into());
                    }

                    match self.regions.get(reg_key) {
                        None => {
                            let new_reg = Arc::new(Mutex::new(ObservedRegion::new(
                                reg_key.id,
                                &reg_key.sequence,
                                reg_key.completeness,
                            )));
                            self.regions.insert(reg_key.clone(), new_reg.clone());
                            reg_map.insert(reg_key.id, new_reg.clone());
                        }
                        Some(r) => {
                            reg_map.insert(reg_key.id, r.clone());
                        }
                    }
                }

                let mut comb = ObservedCombination::new(reg_map, comb_key.sequence.clone());
                comb.increment_count(group);
                self.combinations.insert(comb_key.clone(), comb);
            }
        }

        Ok(())
    }

    /// Update filter counts without checking the read
    ///
    /// Passes through to self.filtered_reads.update_count, useful when using
    /// cached FilterReasons to prevent needing to re-align.
    pub fn update_filter_count(&mut self, read: &ReadPair, reason: FilterReason) {
        self.filtered_reads.increment_count(read, reason)
    }

    /// Determine if a read should be filtered
    ///
    /// Pass through to self.filtered_reads.filter_read, which Checks whether the read should
    /// be filtered, adding it to the appropriate count if increment is true, and returns a the filter reason.
    pub fn filter_readpair(&mut self, record: &ReadPair, increment: bool) -> Option<FilterReason> {
        self.filtered_reads.filter_readpair(record, increment)
    }

    /// Determine if an alignment should be filtered
    ///
    /// Passes through to self.filtered_reads.filter_alignment, which checks if an alignment should
    /// be filtered, adding it to the appropriate count if increment is true, and returns a the filter reason.
    pub fn filter_alignment(
        &mut self,
        record: &ReadPair,
        f_alignment: &Alignment,
        r_alignment: Option<&Alignment>,
        increment: bool,
    ) -> Option<FilterReason> {
        self.filtered_reads
            .filter_alignment(record, f_alignment, r_alignment, increment)
    }

    /// Store a read-level cache entry.
    ///
    /// Cache entries are keyed by full read sequence(s), not by quality values or
    /// other metadata, so only sequence-derived outcomes should be cached.
    pub fn cache(&mut self, key: SeqPair, value: CacheHit) {
        self.cache.insert(key, value);
    }

    /// Check whether this read has a cached result and optionally apply it.
    ///
    /// If `increment` is true, the cached combination or filter result is replayed
    /// into the current counts before returning the cached value.
    pub fn check_cache(
        &mut self,
        record: &ReadPair,
        increment: bool,
    ) -> Result<Option<CacheHit>, anyhow::Error> {
        let key = SeqPair::from_readpair(record);

        let hit: CacheHit = match self.cache.get(&key) {
            Some(x) => x.clone(),
            None => return Ok(None),
        };

        if increment {
            match hit {
                CacheHit::Comb(ref k) => {
                    self.add_or_increment_combination(k, record.group)?;
                }
                CacheHit::Filter(r) => {
                    self.update_filter_count(record, r);
                }
            }
        }

        Ok(Some(hit))
    }

    /// Compare all observed regions and combinations to an expected library.
    ///
    /// This proceeds in three stages:
    /// 1. compare each distinct observed region to the appropriate library region,
    /// 2. combine those per-region matches into per-combination library assignments,
    /// 3. build a summarised set of `LibraryCombination` counts for output.
    ///
    /// This mutates the container in place by:
    /// - storing the compiled library,
    /// - updating each `ObservedRegion` with its nearest library match state,
    /// - updating each `ObservedCombination` with its overall combination match,
    /// - and constructing the library-level summary table.
    ///
    /// Once this has been run, the container is considered library-compared and
    /// can no longer be merged with uncompared containers.
    pub fn compare_to_library(
        &mut self,
        library: Library,
        progress_style: Option<&ProgressStyle>,
        distance_metric: DistanceMetric,
        max_matches: usize,
        threads: usize,
    ) -> Result<(), LibraryError> {
        let n_regs = self.regions.len() as u64;
        let n_combs = self.combinations.len() as u64;

        // Compare each region to the library
        let mut reg_progress: Progress = Progress::from_style(
            progress_style.unwrap_or(&ProgressStyle::new(None, false)),
            "Matching regions:",
            "Matched regions:",
            Some(n_regs),
            match distance_metric {
                DistanceMetric::Hamming | DistanceMetric::Exact => 250000,
                DistanceMetric::BoundedLevenshtein => std::cmp::max(n_regs / 10, 50000),
                DistanceMetric::Levenshtein => std::cmp::max(n_regs / 20, 10000),
            },
        );

        match threads.cmp(&1) {
            std::cmp::Ordering::Less => {
                return Err(LibraryError::Library {
                    desc: "Threads must be >0".to_string(),
                });
            }
            std::cmp::Ordering::Equal => {
                for r in self.regions.values() {
                    let mut reg = r.lock().unwrap();
                    let val = reg.compare_to_library(&library, distance_metric, max_matches);
                    reg.nearest_matches = val;
                    reg_progress.inc(1);
                }
                reg_progress.finish();
            }
            std::cmp::Ordering::Greater => {
                let (reg_tx, reg_rx) = unbounded();
                let (done_tx, done_rx) = unbounded();
                let lib_arc = Arc::new(&library);

                scope(|scope| {
                    // Spin up threads to do work
                    for _ in 0..threads {
                        let rx: Receiver<Arc<Mutex<ObservedRegion>>> = reg_rx.clone();
                        let tx = done_tx.clone();
                        let lib = lib_arc.clone();

                        scope.spawn(move || {
                            while let Ok(region) = rx.recv() {
                                let mut reg = region.lock().unwrap();
                                let val =
                                    reg.compare_to_library(&lib, distance_metric, max_matches);
                                reg.nearest_matches = val;
                                tx.send(()).expect("Main thread comparison reciever failed");
                            }
                        });
                    }

                    // Send regions to workers
                    for region in self.regions.values() {
                        reg_tx
                            .send(region.clone())
                            .expect("Library comparison thread send failed");
                    }
                    drop(reg_tx);
                    drop(reg_rx);
                    drop(done_tx);

                    // Collect results to check all regions processed
                    for _ in done_rx.iter() {
                        reg_progress.inc(1);
                    }
                    drop(done_rx);
                });
                reg_progress.finish();
            }
        }

        // Compare the combinations to the libary
        let mut comb_progress: Progress = Progress::from_style(
            progress_style.unwrap_or(&ProgressStyle::new(None, false)),
            "Comparing combinations:",
            "Compared combinations:",
            Some(n_combs),
            2000000,
        );

        for value in self.combinations.values_mut() {
            value.library_matches =
                value.compare_to_library(&self.region_ids, &library, distance_metric, max_matches);
            comb_progress.inc(1);
        }
        comb_progress.finish();

        // Make summary counts of observed library combinations
        let mut lib_summary_progress: Progress = Progress::from_style(
            progress_style.unwrap_or(&ProgressStyle::new(None, false)),
            "Summarising library matches:",
            "Summarised library matches:",
            Some(n_combs),
            std::cmp::max(n_combs / 4, 250000),
        );

        let mut lib_combs: HashMap<LibraryCombinationKey, LibraryCombination> = HashMap::new();

        for comb in self.combinations.values_mut() {
            let mut key: LibraryCombinationKey = LibraryCombinationKey::new(Vec::with_capacity(5));

            for reg in &self.region_ids {
                // Add the appropriate library match to the key. Simplify over distance,
                // match count, etc. to summarise the library and store NoLibrary
                // region seqs to capture e.g. barcodes.
                match comb.regions.get(reg) {
                    None => key.regions.push((*reg, LibraryRegionMatch::Unmatched)),
                    Some(x) => {
                        let or = x.lock().unwrap();
                        key.regions.push((
                            *reg,
                            LibraryRegionMatch::from_region_match(&or.nearest_matches),
                        ))
                    }
                }
            }

            match lib_combs.get_mut(&key) {
                Some(x) => {
                    for (group, count) in &comb.counts {
                        x.increment_count(group, *count);
                    }
                }
                None => {
                    let mut x: LibraryCombination = LibraryCombination::new(
                        HashMap::from_iter(key.regions.clone()),
                        match &comb.library_matches {
                            CombinationMatch::Uncompared => CombinationMatch::Uncompared,
                            CombinationMatch::Match { inds, .. } => CombinationMatch::Match {
                                inds: inds.clone(),
                                distance: 0,
                            },
                            CombinationMatch::MultiMatch { inds, .. } => {
                                CombinationMatch::MultiMatch {
                                    inds: inds.clone(),
                                    distance: 0,
                                }
                            }
                            CombinationMatch::Recombination { .. } => {
                                CombinationMatch::Recombination { distance: 0 }
                            }
                            CombinationMatch::Mismatch => CombinationMatch::Mismatch,
                            CombinationMatch::Nonmatch => CombinationMatch::Nonmatch,
                        },
                    );
                    for (group, count) in &comb.counts {
                        x.increment_count(group, *count);
                    }
                    lib_combs.insert(key, x);
                }
            }

            lib_summary_progress.inc(1);
        }
        lib_summary_progress.finish();

        self.library = Some(library);
        self.library_combinations = Some(lib_combs);

        Ok(())
    }

    /// Return `true` if library comparison has been run on this container.
    pub fn is_compared_to_library(&self) -> bool {
        self.library.is_some()
    }

    /// Summarise observed counts into high-level read categories.
    ///
    /// The returned `ReadSummary` includes both filtered-read totals and counts of
    /// each `CombinationMatch` category across unfiltered reads.
    ///
    /// Before library comparison, all unfiltered reads fall into the `uncompared`
    /// category.
    pub fn summarise(&self) -> ReadSummary {
        let mut read_summary = ReadSummary::empty();

        read_summary.filtered_reads = self.filtered_reads.totals.clone();

        for comb in self.combinations.values() {
            let count: u64 = comb.total_count() as u64;

            match comb.library_matches {
                CombinationMatch::Uncompared => read_summary.uncompared += count,
                CombinationMatch::Match { distance, .. } => {
                    if distance == 0 {
                        read_summary.exact_match += count
                    } else {
                        read_summary.nearest_match += count
                    }
                }
                CombinationMatch::MultiMatch { .. } => read_summary.multimatch += count,
                CombinationMatch::Recombination { distance } => {
                    if distance == 0 {
                        read_summary.exact_recombination += count
                    } else {
                        read_summary.nearest_recombination += count
                    }
                }
                CombinationMatch::Mismatch => read_summary.mismatch += count,
                CombinationMatch::Nonmatch => read_summary.nonmatch += count,
            }
        }

        read_summary
    }

    /// Output the library-combination summary as a vector for iteration/output.
    ///
    /// If `sort` is true, combinations are returned in descending total-count order.
    /// This allocates a new vector of references.
    ///
    /// Returns an error if library comparison has not yet been performed.
    pub fn to_library_vector(
        &self,
        sort: bool,
    ) -> Result<Vec<&LibraryCombination>, ReadCountError> {
        if let Some(combs) = &self.library_combinations {
            let mut vec: Vec<&LibraryCombination> = combs.iter().map(|(_, v)| v).collect();

            if sort {
                vec.sort_unstable_by_key(|c| std::cmp::Reverse(c.total_count()));
            }

            Ok(vec)
        } else {
            Err(ReadCountError::Error {
                desc: "Combinations uncompared, compare before summarising".to_string(),
            })
        }
    }

    /// Output observed combinations as a vector for iteration/output.
    ///
    /// If `sort` is true, combinations are returned in descending total-count order.
    /// This allocates a new vector of references.
    pub fn to_vector(&self, sort: bool) -> Vec<&ObservedCombination> {
        let mut vec: Vec<&ObservedCombination> = self.combinations.values().collect();

        if sort {
            vec.sort_unstable_by_key(|c| std::cmp::Reverse(c.total_count()));
        }

        vec
    }
}

/// Read-level cache mapping full read sequences to previously computed outcomes.
pub type ObservedReads = HashMap<SeqPair, CacheHit>;

/// Cached outcome for a previously seen read sequence.
///
/// The cache works on `SeqPair`s only, which only store sequences, so only
/// sequence-derived outcomes should be cached. In particular, filter results
/// that depend on qualities or alignment scoring should not be reused purely from
/// sequence identity unless that behaviour is known to be correct for the calling
/// context.
#[derive(Debug, Clone)]
pub enum CacheHit {
    Comb(CombinationKey),
    Filter(FilterReason),
}

/// High-level summary of read outcomes.
///
/// This struct collapses the full observed-combination table into broad categories
/// used for reporting, including:
/// - unmatched/uncompared reads,
/// - exact and inexact library matches,
/// - multimatches,
/// - recombinations,
/// - mismatches,
/// - nonmatches,
/// - and filtered-read totals by reason.
pub struct ReadSummary {
    /// Comparison hasn't occured
    pub uncompared: u64,

    /// Full match with a specific library member
    pub exact_match: u64,

    /// Nearest match with a specific library member
    pub nearest_match: u64,

    /// Fully matches multiple library members and the distance
    pub multimatch: u64,

    /// Exactly matches multiple library members but recombined
    pub exact_recombination: u64,

    /// Partially matches multiple library members but recombinaed
    pub nearest_recombination: u64,

    /// Regions exist but at least one cannot be assigned to the library
    pub mismatch: u64,

    /// Not all regions exist
    pub nonmatch: u64,

    /// Filtered Reads
    pub filtered_reads: FilteredCounts,
}

impl ReadSummary {
    #[allow(dead_code)]
    fn new(
        uncompared: u64,
        exact_match: u64,
        nearest_match: u64,
        multimatch: u64,
        exact_recombination: u64,
        nearest_recombination: u64,
        mismatch: u64,
        nonmatch: u64,
        filtered_reads: FilteredCounts,
    ) -> Self {
        Self {
            uncompared,
            exact_match,
            nearest_match,
            multimatch,
            exact_recombination,
            nearest_recombination,
            mismatch,
            nonmatch,
            filtered_reads,
        }
    }

    /// Inititalise an empty ReadSummary
    ///
    /// Useful shortcut for using it as a counter
    fn empty() -> Self {
        Self {
            uncompared: 0,
            exact_match: 0,
            nearest_match: 0,
            multimatch: 0,
            exact_recombination: 0,
            nearest_recombination: 0,
            mismatch: 0,
            nonmatch: 0,
            filtered_reads: FilteredCounts::new(),
        }
    }

    /// Total number of unfiltered reads represented in the summary.
    pub fn total_unfiltered(&self) -> u64 {
        self.uncompared
            + self.exact_match
            + self.nearest_match
            + self.multimatch
            + self.exact_recombination
            + self.nearest_recombination
            + self.mismatch
            + self.nonmatch
    }

    /// Total number of reads processed, including filtered reads.
    pub fn total(&self) -> u64 {
        self.total_unfiltered() + self.filtered_reads.total()
    }
}

#[cfg(test)]
mod tests {
    // use super::*;
}
