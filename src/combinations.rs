//! Store collections of observed combinations
//!
//! Structures and functions to store and manipulate groups of observed combinations extracted from sequencing data
use anyhow::{self};
use bio::bio_types::alignment::Alignment;
use crossbeam::channel::{Receiver, unbounded};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::{Arc, Mutex};
use std::thread::scope;

use crate::combination::{CombinationKey, CombinationMatch, ObservedCombination};
use crate::errors::{LibraryError, ReadCountError};
use crate::filters::{FilterConfig, FilterReason, FilteredCounts, FilteredReads};
use crate::lib_spec::{DistanceMetric, Library};
use crate::library_combination::{LibraryCombination, LibraryCombinationKey};
use crate::logging::{Progress, ProgressStyle};
use crate::region::{ObservedRegion, RegionKey, RegionMatch};
use crate::seqs::ReadGroup;
use crate::seqs::ReadPair;
use crate::utils::div_or_zero;

/// Container for ObservedCombination objects
///
/// The core is a HashMap of the ObservedRegions seen so far and a HashMap
/// of ObservedCombination objects that link to the contained regions.
/// An additional HashMap plus the library is added when library comparison
/// is run. It also carries the region ids to be considered in order.
#[derive(Debug)]
pub struct ObservedCombinations {
    region_ids: Vec<String>,
    regions: HashMap<RegionKey, Arc<Mutex<ObservedRegion>>>,
    combinations: HashMap<CombinationKey, ObservedCombination>,
    library: Option<Library>,
    library_combinations: Option<HashMap<LibraryCombinationKey, LibraryCombination>>,
    filtered_reads: FilteredReads,
}

impl ObservedCombinations {
    pub fn new(region_ids: Vec<String>, filter_config: FilterConfig) -> Self {
        Self {
            region_ids,
            regions: HashMap::new(),
            combinations: HashMap::new(),
            library: None,
            library_combinations: None,
            filtered_reads: FilteredReads::new(filter_config),
        }
    }

    /// Merge counts from another ObservedCombinations object into this one
    ///
    /// New counts and regions are added to this container with an empty return on
    /// successful operation and an error returned if the two ObservedCombinations objects
    /// are incompatible.
    ///
    /// Reasons for incompatibility:
    /// - Library comparison has occured
    /// - Region ids are not the same
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
                                old_comb.counts.insert(group.clone(), *new_count);
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

                        new_comb.regions.insert(reg_key.id.clone(), arc.clone());
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

    /// Increment a combination count or add a new combination if it hasn't been seen yet
    pub fn add_or_increment_combination(
        &mut self,
        key: &CombinationKey,
        group: ReadGroup,
    ) -> Result<(), anyhow::Error> {
        match self.combinations.get_mut(key) {
            Some(comb) => comb.increment_count(group),
            None => {
                let mut reg_map = HashMap::new();

                for reg_key in &key.regions {
                    if !self.region_ids.contains(&reg_key.id) {
                        return Err(ReadCountError::UnexpectedRegion {
                            region: reg_key.id.clone(),
                        }
                        .into());
                    }

                    match self.regions.get(reg_key) {
                        None => {
                            let new_reg = Arc::new(Mutex::new(ObservedRegion::new(
                                reg_key.id.clone(),
                                &reg_key.sequence,
                                reg_key.completeness,
                            )));
                            self.regions.insert(reg_key.clone(), new_reg.clone());
                            reg_map.insert(reg_key.id.clone(), new_reg.clone());
                        }
                        Some(r) => {
                            reg_map.insert(reg_key.id.clone(), r.clone());
                        }
                    }
                }

                let mut comb = ObservedCombination::new(reg_map, key.sequence.clone());
                comb.increment_count(group);
                self.combinations.insert(key.clone(), comb);
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
    /// be filtered, adding it to the appropriate count if is, and returns a bool determining
    /// if it was filtered.
    pub fn filter_readpair(&mut self, record: &ReadPair) -> Option<FilterReason> {
        self.filtered_reads.filter_readpair(record)
    }

    /// Determine if an alignment should be filtered
    ///
    /// Passes through to self.filtered_reads.filter_alignment, which checks if an alignment should
    /// be filtered, adding it to the appropriate count if so, and returns a bool determining if it
    /// was filtered
    pub fn filter_alignment(
        &mut self,
        record: &ReadPair,
        f_alignment: &Alignment,
        r_alignment: Option<&Alignment>,
    ) -> Option<FilterReason> {
        self.filtered_reads
            .filter_alignment(record, f_alignment, r_alignment)
    }

    /// Compare observed combinations to those expected in a Library
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
                    None => key.regions.push((reg.to_string(), RegionMatch::Unmatched)),
                    Some(x) => {
                        let or = x.lock().unwrap();
                        key.regions.push((
                            reg.to_string(),
                            match &or.nearest_matches {
                                RegionMatch::Unmatched => RegionMatch::Unmatched,
                                RegionMatch::Overmatched { .. } => RegionMatch::Overmatched {
                                    distance: 0,
                                    matches: 0,
                                },
                                RegionMatch::NoLibrary { .. } => RegionMatch::NoLibrary {
                                    seq: Some(or.seq.clone()),
                                },
                                RegionMatch::Uncompared => RegionMatch::Uncompared,
                                RegionMatch::Match { seq_match, .. } => RegionMatch::Match {
                                    seq_match: seq_match.clone(),
                                    distance: 0,
                                },
                                RegionMatch::MultiMatch { seq_matches, .. } => {
                                    RegionMatch::MultiMatch {
                                        seq_matches: seq_matches.to_vec(),
                                        distance: 0,
                                    }
                                }
                            },
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
                            CombinationMatch::Match { ind, .. } => CombinationMatch::Match {
                                ind: *ind,
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

    /// Check if library comparison has occured
    pub fn is_compared_to_library(&self) -> bool {
        self.library.is_some()
    }

    /// Summarise the observed  read count categories
    ///
    /// Counts the occurance of each CombinationMatch, so makes little sense if library comparison
    /// hasn't occured first as it will just sum the total reads
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

    /// Write all counts to file
    pub fn write_tsv(&self, file: File, sort: bool) -> Result<(), anyhow::Error> {
        let mut count_writer = BufWriter::new(file);
        let mut keys: Vec<(&CombinationKey, u32)> = self
            .combinations
            .iter()
            .map(|x| (x.0, x.1.total_count()))
            .collect();

        if sort {
            // Invert count to get desc order
            keys.sort_unstable_by_key(|x| 0 - i64::from(x.1));
        }

        // Write header
        write!(count_writer, "group\tforward\treverse\t")?;
        for r in &self.region_ids {
            write!(
                count_writer,
                "{}\t{}_nearest\t{}_distance\t{}_n_matches\t",
                r, r, r, r
            )?;
        }
        writeln!(
            count_writer,
            "combination_status\tcombination_distance\tcombinations_in_library\tcombination_id\tcount"
        )?;

        for (key, _) in keys {
            let combination = self.combinations.get(key).expect(
                "Combination key from extracted key list missing from ObservedCombinations",
            );
            write!(
                count_writer,
                "{}",
                combination.to_tsv(&self.region_ids, self.library.as_ref())?
            )?;
        }

        count_writer.flush()?;
        Ok(())
    }

    /// Write matched library combination counts to file
    pub fn write_summary_tsv(&self, file: File, sort: bool) -> Result<(), anyhow::Error> {
        let combs = match &self.library_combinations {
            None => {
                return Err(ReadCountError::Error {
                    desc: "Combinations uncompared, compare before summarising".to_string(),
                }
                .into());
            }
            Some(x) => x,
        };

        let mut writer = BufWriter::new(file);
        let mut keys: Vec<(&LibraryCombinationKey, u32)> =
            combs.iter().map(|x| (x.0, x.1.total_count())).collect();

        // Sort
        if sort {
            // Invert count to get desc order
            keys.sort_unstable_by_key(|x| 0 - i64::from(x.1));
        }

        // Write
        write!(writer, "group\t")?;
        for r in &self.region_ids {
            write!(writer, "{}\t", r)?;
        }
        writeln!(
            writer,
            "combination_status\tcombinations_in_library\tcombination_id\tcount"
        )?;

        for (key, _) in keys {
            let combination = combs.get(key).expect(
                "Combination key from extracted key list missing from ObservedCombinations",
            );
            write!(
                writer,
                "{}",
                combination.to_tsv(&self.region_ids, self.library.as_ref())?
            )?;
        }

        writer.flush()?;
        Ok(())
    }

    /// Write all counts to file
    pub fn write_filtered_tsv(&self, file: File, sort: bool) -> Result<(), anyhow::Error> {
        self.filtered_reads.write_filter_tsv(file, sort)?;
        Ok(())
    }
}

/// Summary counts of read types
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

    /// Sum of unfiltered reads
    fn total_unfiltered(&self) -> u64 {
        self.uncompared
            + self.exact_match
            + self.nearest_match
            + self.multimatch
            + self.exact_recombination
            + self.nearest_recombination
            + self.mismatch
            + self.nonmatch
    }

    /// Total reads observed across categories
    pub fn total(&self) -> u64 {
        self.total_unfiltered() + self.filtered_reads.total()
    }

    /// Write to file in TSV format
    pub fn write_tsv(self, file: File) -> Result<(), anyhow::Error> {
        let mut writer = BufWriter::new(file);
        let unfiltered_total = self.total_unfiltered();
        let total = self.total();

        // Write data
        writeln!(
            writer,
            "group\tmetric\tcount\toverall_proportion\tgroup_proportion"
        )?;
        writeln!(
            writer,
            "all\ttotal\t{}\t{:.4}\t{:.4}",
            total,
            if total > 0 { 1.0 } else { 0.0 },
            if total > 0 { 1.0 } else { 0.0 },
        )?;
        writeln!(
            writer,
            "unfiltered\ttotal\t{}\t{:.4}\t{:.4}",
            unfiltered_total,
            div_or_zero(unfiltered_total as f32, total as f32),
            div_or_zero(unfiltered_total as f32, unfiltered_total as f32),
        )?;
        writeln!(
            writer,
            "unfiltered\tuncompared\t{}\t{:.4}\t{:.4}",
            self.uncompared,
            div_or_zero(self.uncompared as f32, total as f32),
            div_or_zero(self.uncompared as f32, unfiltered_total as f32),
        )?;
        writeln!(
            writer,
            "unfiltered\texact_match\t{}\t{:.4}\t{:.4}",
            self.exact_match,
            div_or_zero(self.exact_match as f32, total as f32),
            div_or_zero(self.exact_match as f32, unfiltered_total as f32),
        )?;
        writeln!(
            writer,
            "unfiltered\tnearest_match\t{}\t{:.4}\t{:.4}",
            self.nearest_match,
            div_or_zero(self.nearest_match as f32, total as f32),
            div_or_zero(self.nearest_match as f32, unfiltered_total as f32),
        )?;
        writeln!(
            writer,
            "unfiltered\tmultimatch\t{}\t{:.4}\t{:.4}",
            self.multimatch,
            div_or_zero(self.multimatch as f32, total as f32),
            div_or_zero(self.multimatch as f32, unfiltered_total as f32),
        )?;
        writeln!(
            writer,
            "unfiltered\texact_recombination\t{}\t{:.4}\t{:.4}",
            self.exact_recombination,
            div_or_zero(self.exact_recombination as f32, total as f32),
            div_or_zero(self.exact_recombination as f32, unfiltered_total as f32),
        )?;
        writeln!(
            writer,
            "unfiltered\tnearest_recombination\t{}\t{:.4}\t{:.4}",
            self.nearest_recombination,
            div_or_zero(self.nearest_recombination as f32, total as f32),
            div_or_zero(self.nearest_recombination as f32, unfiltered_total as f32),
        )?;
        writeln!(
            writer,
            "unfiltered\tmismatch\t{}\t{:.4}\t{:.4}",
            self.mismatch,
            div_or_zero(self.mismatch as f32, total as f32),
            div_or_zero(self.mismatch as f32, unfiltered_total as f32),
        )?;
        writeln!(
            writer,
            "unfiltered\tnonmatch\t{}\t{:.4}\t{:.4}",
            self.nonmatch,
            div_or_zero(self.nonmatch as f32, total as f32),
            div_or_zero(self.nonmatch as f32, unfiltered_total as f32),
        )?;

        write!(writer, "{}", self.filtered_reads.to_long_tsv_lines(total))?;

        writer.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    // use super::*;
}
