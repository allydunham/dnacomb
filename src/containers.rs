//! Container for observed sequence reads, built around representing
//! them as combinations of regions
//!
//! Contains structs and supporting methods to count the type of reads observed
//! in your input data.
//! Reads are compressed to a combination of specific regions of interest,  as
//! defined by a LibSpec.
//! Supporting methods help use the containers, for instance outputing TSV files
//! and assigning each observed combination to a library of expected sequences defined
//! in a Library object.
use anyhow::{self};
use bio::bio_types::{alignment::Alignment, sequence::Sequence};
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::rc::Rc;

use crate::errors::{LibraryError, ReadCountError, seq_to_string_or_log};
use crate::filters::{FilterConfig, FilterReason, FilteredReads};
use crate::lib_spec::{DistanceMetric, Library, LibraryRegion, PartialMatching, merge_matches};
use crate::logging::{Progress, ProgressStyle};
use crate::parsing::{ReadGroup, ReadPair};

/// Region keys identify via the region name, the observed sequence and completeness status
pub type RegionKey = (String, Sequence, RegionCompleteness);

/// Combination keys identify via a vector of RegionKeys
pub type CombinationKey = Vec<RegionKey>;

/// Library keys identify a particular library match type
pub type LibraryKey = Vec<(String, RegionMatch)>;

/// Container for ObservedCombinations
///
/// The core is a HashMap of the ObservedRegions seen so far and a HashMap
/// of ObservedCombination objects that link to the contained regions.
/// An additional HashMap plus the library is added when library comparison
/// is run. It also carries the region ids to be considered in order.
#[derive(Debug)]
pub struct ObservedCombinations {
    region_ids: Vec<String>,
    regions: HashMap<RegionKey, Rc<RefCell<ObservedRegion>>>,
    combinations: HashMap<CombinationKey, ObservedCombination>,
    library: Option<Library>,
    library_combinations: Option<HashMap<LibraryKey, LibraryCombination>>,
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

    /// Number of distinct observed combination types
    pub fn len(&self) -> usize {
        self.combinations.len()
    }

    /// Check if there are no observed combinations
    pub fn is_empty(&self) -> bool {
        self.combinations.is_empty()
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

                for reg_key in key {
                    if !self.region_ids.contains(&reg_key.0) {
                        return Err(ReadCountError::UnexpectedRegion {
                            region: reg_key.0.clone(),
                        }
                        .into());
                    }

                    match self.regions.get(reg_key) {
                        None => {
                            let new_reg = Rc::new(RefCell::new(ObservedRegion::new(
                                reg_key.0.clone(),
                                &reg_key.1,
                                reg_key.2,
                            )));
                            self.regions.insert(reg_key.clone(), new_reg.clone());
                            reg_map.insert(reg_key.0.clone(), new_reg.clone());
                        }
                        Some(r) => {
                            reg_map.insert(r.borrow().id.clone(), r.clone());
                        }
                    }
                }

                let mut comb = ObservedCombination::new(reg_map);
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
    pub fn update_filter_count(&mut self, reason: &FilterReason) {
        self.filtered_reads.update_count(reason)
    }

    /// Determine if a read should be filtered
    ///
    /// Pass through to self.filtered_reads.filter_read, which Checks whether the read should
    /// be filtered, adding it to the appropriate count if is, and returns a bool determining
    /// if it was filtered.
    pub fn filter_readpair(&mut self, record: &ReadPair) -> FilterReason {
        self.filtered_reads.filter_readpair(record)
    }

    /// Determine if an alignment should be filtered
    ///
    /// Passes through to self.filtered_reads.filter_alignment, which checks if an alignment should
    /// be filtered, adding it to the appropriate count if so, and returns a bool determining if it
    /// was filtered
    pub fn filter_alignment(
        &mut self,
        f_alignment: &Alignment,
        r_alignment: Option<&Alignment>,
    ) -> FilterReason {
        self.filtered_reads
            .filter_alignment(f_alignment, r_alignment)
    }

    /// Compare observed combinations to those expected in a Library
    pub fn compare_to_library(
        &mut self,
        library: Library,
        progress_style: Option<&ProgressStyle>,
        distance_metric: DistanceMetric,
        max_matches: usize,
    ) -> Result<(), LibraryError> {
        let n_regs = self.regions.len() as u64;
        let n_combs = self.combinations.len() as u64;

        // Compare each region to the library
        let mut reg_progress: Progress = Progress::from_style(
            progress_style.unwrap_or(&ProgressStyle::new(None)),
            "Matching regions:",
            "Matched regions:",
            Some(n_regs),
            match distance_metric {
                DistanceMetric::Hamming | DistanceMetric::Exact => 250000,
                DistanceMetric::BoundedLevenshtein => std::cmp::max(n_regs / 10, 50000),
                DistanceMetric::Levenshtein => std::cmp::max(n_regs / 20, 10000),
            },
        );

        for value in self.regions.values() {
            let val = value
                .borrow()
                .compare_to_library(&library, distance_metric, max_matches);
            value.borrow_mut().nearest_matches = val;
            reg_progress.inc(1);
        }
        reg_progress.finish();

        // Compare the combinations to the libary
        let mut comb_progress: Progress = Progress::from_style(
            progress_style.unwrap_or(&ProgressStyle::new(None)),
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
            progress_style.unwrap_or(&ProgressStyle::new(None)),
            "Summarising library matches:",
            "Summarised library matches:",
            Some(n_combs),
            std::cmp::max(n_combs / 4, 250000),
        );

        let mut lib_combs: HashMap<Vec<(String, RegionMatch)>, LibraryCombination> = HashMap::new();

        for comb in self.combinations.values_mut() {
            let mut key: LibraryKey = Vec::with_capacity(5);

            for reg in &self.region_ids {
                // Add the appropriate library match to the key. Simplify over distance,
                // match count, etc. to summarise the library and store NoLibrary
                // region seqs to capture e.g. barcodes.
                match comb.regions.get(reg) {
                    None => key.push((reg.to_string(), RegionMatch::Unmatched)),
                    Some(x) => key.push((
                        reg.to_string(),
                        match &x.borrow().nearest_matches {
                            RegionMatch::Unmatched => RegionMatch::Unmatched,
                            RegionMatch::Overmatched { .. } => RegionMatch::Overmatched {
                                distance: 0,
                                matches: 0,
                            },
                            RegionMatch::NoLibrary { .. } => RegionMatch::NoLibrary {
                                seq: Some(x.borrow().seq.clone()),
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
                    )),
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
                        HashMap::from_iter(key.clone()),
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

        read_summary.filtered_reads = Some(self.filtered_reads.clone());

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
        write!(count_writer, "group\t")?;
        for r in &self.region_ids {
            write!(
                count_writer,
                "{}\t{}_nearest\t{}_distance\t{}_n_matches\t",
                r, r, r, r
            )?;
        }
        writeln!(
            count_writer,
            "combination_status\tcombination_distance\tcombinations_in_library\tcombination_indexes\tcount"
        )?;

        for (key, _) in keys {
            let combination = self.combinations.get(key).expect(
                "Combination key from extracted key list missing from ObservedCombinations",
            );
            write!(count_writer, "{}", combination.to_tsv(&self.region_ids))?;
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
        let mut keys: Vec<(&LibraryKey, u32)> =
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
            "combination_status\tcombinations_in_library\tcombination_indexes\tcount"
        )?;

        for (key, _) in keys {
            let combination = combs.get(key).expect(
                "Combination key from extracted key list missing from ObservedCombinations",
            );
            write!(writer, "{}", combination.to_tsv(&self.region_ids))?;
        }

        writer.flush()?;
        Ok(())
    }
}

/// Combination of ObservedRegions seen in sequence reads
///
/// A set of observed regions determining the "type" of read, as defined in the LibSpec.
/// Also includes a count, the read grouping (for instance for different cells in single
/// cell studies) and whether it matches an expected library member.
#[derive(Debug)]
pub struct ObservedCombination {
    /// Count of observations for each read group. Ungrouped reads are stored in None
    counts: HashMap<ReadGroup, u32>,

    /// ObservedRegions defining the sequence form. References to ObservedRegion which
    /// should be stored in the parent ObservedCombinations object.
    regions: HashMap<String, Rc<RefCell<ObservedRegion>>>,

    /// Status and result of comparison with the expected library of sequences
    library_matches: CombinationMatch,
}

/// Status of the match between ObservedCombination and a Library
///
/// Includes the match status and the indeces of matches in the libray
/// plus the distance to the library.
#[derive(Debug, Clone)]
pub enum CombinationMatch {
    /// Comparison hasn't occured
    Uncompared,

    /// Full match with a specific library member and the distance
    Match { ind: usize, distance: u64 },

    /// Fully matches multiple library members and the distance
    MultiMatch { inds: HashSet<usize>, distance: u64 },

    /// Partially matches multiple library members and the total distance
    Recombination { distance: u64 },

    /// Regions exist but at least one cannot be assigned to the library
    Mismatch,

    /// Not all regions exist
    Nonmatch,
}

impl CombinationMatch {
    /// Output a TSV chunk for the combination match status
    ///
    /// Has the \t separated format:
    /// status combination_distance combinations_in_library combination_indexes
    fn to_tsv_chunk(&self) -> String {
        match self {
            CombinationMatch::Uncompared => "uncompared\t\t\t\t".to_string(),
            CombinationMatch::Match { ind, distance } => format!("match\t{distance}\t1\t{ind}\t",),
            CombinationMatch::MultiMatch { inds, distance } => format!(
                "match\t{}\t{}\t{}\t",
                distance,
                inds.len(),
                inds.iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            ),
            CombinationMatch::Recombination { distance } => {
                format!("recombination\t{distance}\t0\t\t",)
            }
            CombinationMatch::Mismatch => "mismatch\t\t0\t\t".to_string(),
            CombinationMatch::Nonmatch => "nonmatch\t\t0\t\t".to_string(),
        }
    }

    /// Output a summary TSV chunk for the combination match status
    ///
    /// Has the \t separated format:
    /// status combinations_in_library combination_indexes
    fn to_summary_tsv_chunk(&self) -> String {
        match self {
            CombinationMatch::Uncompared => "uncompared\t\t\t".to_string(),
            CombinationMatch::Match { ind, .. } => format!("match\t1\t{ind}\t",),
            CombinationMatch::MultiMatch { inds, .. } => format!(
                "match\t{}\t{}\t",
                inds.len(),
                inds.iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            ),
            CombinationMatch::Recombination { .. } => "recombination\t0\t\t".to_string(),
            CombinationMatch::Mismatch => "mismatch\t0\t\t".to_string(),
            CombinationMatch::Nonmatch => "nonmatch\t0\t\t".to_string(),
        }
    }
}

impl ObservedCombination {
    fn new(regions: HashMap<String, Rc<RefCell<ObservedRegion>>>) -> Self {
        Self {
            counts: HashMap::new(),
            regions,
            library_matches: CombinationMatch::Uncompared,
        }
    }

    /// Generate a key uniquely identifying a combination
    pub fn key(regions: Vec<RegionKey>) -> CombinationKey {
        regions
    }

    /// Total count across all read groups
    fn total_count(&self) -> u32 {
        self.counts.values().sum()
    }

    /// Incremment count for a read group
    fn increment_count(&mut self, group: ReadGroup) {
        match self.counts.get_mut(&group) {
            Some(x) => *x += 1,
            None => {
                self.counts.insert(group, 1);
            }
        }
    }

    /// Compare the combination to expected combinations in the library
    ///
    /// Returns a CombinationMatch object which can also be added to the ObservedCombination
    /// library matches field. Looks at each region in turn and identifies which library
    /// combinations are possible overall matches.
    fn compare_to_library(
        &self,
        region_ids: &Vec<String>,
        library: &Library,
        distance_metric: DistanceMetric,
        max_matches: usize,
    ) -> CombinationMatch {
        let mut comb_dist: u64 = 0;
        let mut candidate_matches: Option<HashSet<usize>> = None;

        for reg_id in region_ids {
            let reg = match self.regions.get(reg_id) {
                None => return CombinationMatch::Nonmatch,
                Some(x) => {
                    // Should never need this with the implementation in ObservedCombinations, but here as a back-up as otherwise could panic later. Do all regions first as slightly more efficient and easier to follow in log
                    if !x.borrow().is_compared_to_library() {
                        let val =
                            x.borrow()
                                .compare_to_library(library, distance_metric, max_matches);
                        x.borrow_mut().nearest_matches = val;
                    }
                    x
                }
            };

            match &reg.borrow().nearest_matches {
                RegionMatch::Uncompared => panic!("Region uncompared despite just comparing"),
                RegionMatch::Unmatched | RegionMatch::Overmatched { .. } => {
                    // An unmatched region or an indeterminate one means the
                    // combination cannot be assigned
                    return CombinationMatch::Mismatch;
                }
                RegionMatch::NoLibrary { .. } => {}
                RegionMatch::Match {
                    seq_match,
                    distance,
                } => {
                    comb_dist += distance;
                    match candidate_matches {
                        None => candidate_matches = Some(seq_match.inds.clone()),
                        Some(ref mut x) => {
                            x.retain(|i| seq_match.inds.contains(i));
                        }
                    };
                }
                RegionMatch::MultiMatch {
                    seq_matches,
                    distance,
                } => {
                    comb_dist += distance;
                    for mat in seq_matches {
                        match candidate_matches {
                            None => candidate_matches = Some(mat.inds.clone()),
                            Some(ref mut x) => {
                                x.retain(|i| mat.inds.contains(i));
                            }
                        };
                    }
                }
            }
        }

        match candidate_matches {
            None => CombinationMatch::Nonmatch, // Only occurs if no regions (e.g. empty reads)
            Some(x) => {
                if x.len() == 1 {
                    CombinationMatch::Match {
                        // Can unwrap because we know x.len() == 1
                        ind: *x.iter().next().unwrap(),
                        distance: comb_dist,
                    }
                } else if x.is_empty() {
                    CombinationMatch::Recombination {
                        distance: comb_dist,
                    }
                } else {
                    CombinationMatch::MultiMatch {
                        inds: x,
                        distance: comb_dist,
                    }
                }
            }
        }
    }

    /// Generate tsv line(s) corresponding to this combination. Each read group
    /// the combination is observed is given a separate line
    fn to_tsv(&self, region_ids: &Vec<String>) -> String {
        // Line has \t separated format:
        // group [{region} {region}_nearest {region}_distance {region}_n_matches for each region] status combination_distance combinations_in_library combination_indexes count

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

            // Region seq/nearest match(s)/distance per region
            for reg_id in region_ids {
                let region = self.regions.get(reg_id);

                match region {
                    None => output.push_str("\t\t\t\t"), // Missing regions 4 blanks
                    Some(r) => {
                        output.push_str(&r.borrow().to_tsv_chunk());
                        output.push('\t');
                    }
                }
            }

            output.push_str(&self.library_matches.to_tsv_chunk());
            output.push_str(&count.to_string());
            output.push('\n');
        }

        output
    }
}

/// Observed sequence for a given library region
///
/// Each time this sequence is observed for this region it will be linked back to this object to
/// save library comparison overhead. The struct includes the observed sequence, a marker for
/// whether the whole region has necessarily been observed (i.e. if it is at the the start or end
///  of a read) and a RegionMatch object detailing the library match
#[derive(Debug)]
pub struct ObservedRegion {
    /// Name of the region from LibSpec
    id: String,

    /// Observed Sequence
    seq: Sequence,

    /// Completeness status of the region
    complete: RegionCompleteness,

    /// Library indeces of nearest matches (Vec of Vecs containing indeces from each matching seq)
    nearest_matches: RegionMatch,
}

/// Status of the match between an ObservedRegion and a Library
///
/// Includes the status plus reference(s) to the matched sequence in
/// the library and how distant it is.
#[derive(Debug, Eq, Hash, PartialEq, Clone)]
pub enum RegionMatch {
    /// No comparison has occured yet
    Uncompared,

    /// A single match (`Vec<u8>` sequence and library indeces) and associated distance
    Match {
        seq_match: Rc<LibraryRegion>,
        distance: u64,
    },

    /// Multiple equidistant matches and the distance
    MultiMatch {
        seq_matches: Vec<Rc<LibraryRegion>>,
        distance: u64,
    },

    /// Too many matches
    Overmatched { distance: u64, matches: usize },

    /// No match found
    Unmatched,

    /// Region not in library, with option to store the observed sequence
    NoLibrary { seq: Option<Sequence> },
}

impl RegionMatch {
    /// Get library Sequence, distance and number of matches as strings for output.
    ///
    /// NoLibrary matches are considered not to have a matching sequence
    fn to_str_fields(&self) -> (String, String, String) {
        match self {
            RegionMatch::Uncompared | RegionMatch::Unmatched | RegionMatch::NoLibrary { .. } => {
                ("".to_string(), "".to_string(), "0".to_string())
            }
            RegionMatch::Overmatched { distance, matches } => {
                ("".to_string(), distance.to_string(), matches.to_string())
            }
            RegionMatch::Match {
                seq_match,
                distance,
            } => (
                seq_to_string_or_log(&seq_match.sequence),
                distance.to_string(),
                "1".to_string(),
            ),
            RegionMatch::MultiMatch {
                seq_matches,
                distance,
            } => (
                seq_matches
                    .iter()
                    .map(|x| seq_to_string_or_log(&x.sequence))
                    .collect::<Vec<_>>()
                    .join(","),
                distance.to_string(),
                seq_matches.len().to_string(),
            ),
        }
    }

    /// Get matching Sequence as a string, including the passed through raw
    /// Sequence for NoLibrary matches.
    fn str_sequence(&self) -> String {
        match self {
            RegionMatch::Uncompared | RegionMatch::Unmatched | RegionMatch::Overmatched { .. } => {
                "".to_string()
            }
            RegionMatch::NoLibrary { seq } => match seq {
                Some(s) => seq_to_string_or_log(s),
                None => "".to_string(),
            },
            RegionMatch::Match { seq_match, .. } => seq_to_string_or_log(&seq_match.sequence),
            RegionMatch::MultiMatch { seq_matches, .. } => seq_matches
                .iter()
                .map(|x| seq_to_string_or_log(&x.sequence))
                .collect::<Vec<_>>()
                .join(","),
        }
    }
}

/// Whether an observed region is complete
///
/// Complete regions are known to be full, either being flanked by neighbours on
/// both side in the read or being the full expected length. Partial regions might
/// have more sequence before/after because this is not the case.
#[derive(Debug, Eq, Hash, PartialEq, Clone, Copy)]
pub enum RegionCompleteness {
    /// A full region
    Complete,

    /// 5 prime end is potentially incomplete
    Partial5Prime,

    /// 3 prime end is potentially incomplete
    Partial3Prime,

    /// Reads come from both ends without meeting. Stored as
    /// both seqs divided by a / at split_ind
    MissingCenter { split_ind: usize },

    /// Reads come from both ends and overlap. Stored as
    /// both seqs divided by a / at split_ind
    Overlapping { split_ind: usize },
}

impl ObservedRegion {
    /// Create a new ObservedRegion
    fn new(id: String, seq: &[u8], complete: RegionCompleteness) -> Self {
        Self {
            id,
            seq: seq.to_vec(),
            complete,
            nearest_matches: RegionMatch::Uncompared,
        }
    }

    /// Generate a key for a region for mapping it in an ObservedRegions HashMap
    pub fn key(name: String, seq: &[u8], complete: RegionCompleteness) -> RegionKey {
        (name, seq.to_vec(), complete)
    }

    /// Check if library comparison has been performed
    fn is_compared_to_library(&self) -> bool {
        !matches!(self.nearest_matches, RegionMatch::Uncompared)
    }

    /// Compare the region sequence to a library of oligo options, dispatching to the appropriate implementation
    ///
    /// The distance to each is calculated and all nearest match indeces below the threshold
    /// supplied to library are returned.
    /// If no applicable matches are found nearest_matches will be an empty.
    /// If over max_matches hits are returned the region will be considered overmatched and
    /// indeterminant.
    fn compare_to_library(
        &self,
        library: &Library,
        distance_metric: DistanceMetric,
        max_matches: usize,
    ) -> RegionMatch {
        let lib_match = match self.complete {
            RegionCompleteness::Complete => {
                match library.lookup(&self.id, &self.seq, distance_metric, PartialMatching::Full) {
                    // lookup can only return Err(LibraryError::MissingRegion) which implies the
                    // region isn't in the library. Needs changing if more errors are added to it.
                    Err(_) => return RegionMatch::NoLibrary { seq: None },
                    Ok(x) => x,
                }
            }
            RegionCompleteness::Partial5Prime => {
                match library.lookup(
                    &self.id,
                    &self.seq,
                    distance_metric,
                    PartialMatching::ThreePrimeOnly,
                ) {
                    // lookup can only return Err(LibraryError::MissingRegion) which implies the
                    // region isn't in the library. Needs changing if more errors are added to it.
                    Err(_) => return RegionMatch::NoLibrary { seq: None },
                    Ok(x) => x,
                }
            }
            RegionCompleteness::Partial3Prime => {
                match library.lookup(
                    &self.id,
                    &self.seq,
                    distance_metric,
                    PartialMatching::FivePrimeOnly,
                ) {
                    // lookup can only return Err(LibraryError::MissingRegion) which implies the
                    // region isn't in the library. Needs changing if more errors are added to it.
                    Err(_) => return RegionMatch::NoLibrary { seq: None },
                    Ok(x) => x,
                }
            }
            RegionCompleteness::MissingCenter { split_ind }
            | RegionCompleteness::Overlapping { split_ind } => {
                let left_match = library.lookup(
                    &self.id,
                    &self.seq[0..split_ind],
                    distance_metric,
                    PartialMatching::FivePrimeOnly,
                );

                let right_match = library.lookup(
                    &self.id,
                    &self.seq[split_ind + 1..self.seq.len()],
                    distance_metric,
                    PartialMatching::ThreePrimeOnly,
                );

                match (left_match, right_match) {
                    (Ok(left), Ok(right)) => merge_matches(left, right),

                    // lookup can only return Err(LibraryError::MissingRegion) which implies the
                    // region isn't in the library. Needs changing if more errors are added to it.
                    (_, Err(_)) | (Err(_), _) => return RegionMatch::NoLibrary { seq: None },
                }
            }
        };

        match lib_match {
            None => RegionMatch::Unmatched,
            Some(x) => {
                if x.matches.len() == 1 {
                    RegionMatch::Match {
                        seq_match: x.matches[0].clone(),
                        distance: x.distance,
                    }
                } else if x.matches.len() > max_matches {
                    RegionMatch::Overmatched {
                        distance: x.distance,
                        matches: x.matches.len(),
                    }
                } else {
                    RegionMatch::MultiMatch {
                        seq_matches: x.matches,
                        distance: x.distance,
                    }
                }
            }
        }
    }

    /// Generate output TSV chunk describing the region
    fn to_tsv_chunk(&self) -> String {
        // Generate "{region} {region}_nearest {region}_distance {region}_n_matches" String
        let mut seq = match String::from_utf8(self.seq.clone()) {
            Ok(x) => x,
            Err(_) => panic!("Error converting Vec<u8> sequence to string"),
        };

        // Add appropriate marker to incomplete regions
        match self.complete {
            RegionCompleteness::Complete
            | RegionCompleteness::MissingCenter { .. }
            | RegionCompleteness::Overlapping { .. } => {}
            RegionCompleteness::Partial5Prime => seq.insert(0, '^'),
            RegionCompleteness::Partial3Prime => seq.push('^'),
        };

        let (nearest, distance, matches) = self.nearest_matches.to_str_fields();

        format!("{}\t{}\t{}\t{}", seq, nearest, distance, matches)
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
    fn new(regions: HashMap<String, RegionMatch>, library_matches: CombinationMatch) -> Self {
        Self {
            counts: HashMap::new(),
            regions,
            library_matches,
        }
    }

    /// Get the total count across all read groups
    fn total_count(&self) -> u32 {
        self.counts.values().sum()
    }

    /// Increment the count for the desired read group by n
    fn increment_count(&mut self, group: &ReadGroup, n: u32) {
        match self.counts.get_mut(group) {
            Some(x) => *x += n,
            None => {
                self.counts.insert(group.clone(), n);
            }
        }
    }

    /// Generate tsv line(s) corresponding to this combination. Each read group
    /// the combination is observed is given a separate line
    fn to_tsv(&self, region_ids: &Vec<String>) -> String {
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

            output.push_str(&self.library_matches.to_summary_tsv_chunk());
            output.push_str(&count.to_string());
            output.push('\n');
        }

        output
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
    pub filtered_reads: Option<FilteredReads>,
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
        filtered_reads: Option<FilteredReads>,
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
            filtered_reads: None,
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
        self.total_unfiltered()
            + match &self.filtered_reads {
                Some(f) => f.total(),
                None => 0,
            }
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
        writeln!(writer, "all\ttotal\t{}\t1.0000\t1.000", total)?;
        writeln!(
            writer,
            "unfiltered\ttotal\t{}\t{:.4}\t1.000",
            unfiltered_total,
            unfiltered_total as f32 / total as f32,
        )?;
        writeln!(
            writer,
            "unfiltered\tuncompared\t{}\t{:.4}\t{:.4}",
            self.uncompared,
            self.uncompared as f32 / total as f32,
            self.uncompared as f32 / unfiltered_total as f32,
        )?;
        writeln!(
            writer,
            "unfiltered\texact_match\t{}\t{:.4}\t{:.4}",
            self.exact_match,
            self.exact_match as f32 / total as f32,
            self.exact_match as f32 / unfiltered_total as f32,
        )?;
        writeln!(
            writer,
            "unfiltered\tnearest_match\t{}\t{:.4}\t{:.4}",
            self.nearest_match,
            self.nearest_match as f32 / total as f32,
            self.nearest_match as f32 / unfiltered_total as f32,
        )?;
        writeln!(
            writer,
            "unfiltered\tmultimatch\t{}\t{:.4}\t{:.4}",
            self.multimatch,
            self.multimatch as f32 / total as f32,
            self.multimatch as f32 / unfiltered_total as f32,
        )?;
        writeln!(
            writer,
            "unfiltered\texact_recombination\t{}\t{:.4}\t{:.4}",
            self.exact_recombination,
            self.exact_recombination as f32 / total as f32,
            self.exact_recombination as f32 / unfiltered_total as f32,
        )?;
        writeln!(
            writer,
            "unfiltered\tnearest_recombination\t{}\t{:.4}\t{:.4}",
            self.nearest_recombination,
            self.nearest_recombination as f32 / total as f32,
            self.nearest_recombination as f32 / unfiltered_total as f32,
        )?;
        writeln!(
            writer,
            "unfiltered\tmismatch\t{}\t{:.4}\t{:.4}",
            self.mismatch,
            self.mismatch as f32 / total as f32,
            self.mismatch as f32 / unfiltered_total as f32,
        )?;
        writeln!(
            writer,
            "unfiltered\tnonmatch\t{}\t{:.4}\t{:.4}",
            self.nonmatch,
            self.nonmatch as f32 / total as f32,
            self.nonmatch as f32 / unfiltered_total as f32,
        )?;

        match self.filtered_reads {
            Some(f) => {
                write!(writer, "{}", f.to_tsv(total))?;
            }
            None => {
                writeln!(writer, "filtered\ttotal\t0\t0.0000\t0.0000",)?;
            }
        }

        writer.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;


}

