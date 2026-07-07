//! TSV output generation for DNAComb results.
//!
//! This module contains functions for writing the main output tables produced by
//! DNAComb, including:
//! - full observed-combination counts,
//! - library-assignment summaries,
//! - read-level summary statistics,
//! - and filtered-read summaries.
//!
//! These functions format internal data structures into stable TSV schemas for
//! downstream analysis.
use std::fs::File;
use std::io::{BufWriter, Write};

use itertools::Itertools;

use crate::combination::CombinationMatch;
use crate::combinations::{ObservedCombinations, ReadSummary};
use crate::filters::{FilterReason, FilteredReads};
use crate::interning::region_id_to_str;
use crate::utils::div_or_zero;

/// Write all observed combinations to a TSV file
///
/// Write a full ouput TSV from an ObservedCombinations object with each
/// `ReadGroup`/`ObservedCombination` pair getting one output row.
///
/// This TSV has the format:
/// - group: ReadGroup
/// - forward: Forward read (if tracked)
/// - reverse: Reverse read (if tracked/paired)
/// - Per region:
///     - {region}: the observed sequence
///     - {region}_nearest: the matching library sequence(s)
///     - {region}_diff: difference with the matched library sequence(s)
///     - {region}_distance: distance to the library match(s)
///     - {region}_n_matches: number of matches of distance found
/// - combination_status: Status of the match with the library
/// - combination_distance: Overall distance of the combination to the library
/// - combinations_in_library: Number of matches in the library
/// - combination_id: The library ID of the matching combinations
/// - count: Number of times observed
///
/// Skipping variants removes the {region}_variants column, which simplifies output
/// slightly but is primarily useful when the variants haven't been calculated
/// during library comparison to avoid many alignment operations, which would otherwise
/// produce a column of "NA" strings.
pub fn write_counts(
    combinations: &ObservedCombinations,
    file: File,
    sort: bool,
    skip_variants: bool,
) -> Result<(), anyhow::Error> {
    let mut writer = BufWriter::new(file);

    // Write header
    write!(writer, "group\tforward\treverse\t")?;
    for r in &combinations.region_ids {
        let s = region_id_to_str(r);
        if skip_variants {
            write!(writer, "{s}\t{s}_nearest\t{s}_distance\t{s}_n_matches\t")?;
        } else {
            write!(
                writer,
                "{s}\t{s}_nearest\t{s}_variants\t{s}_distance\t{s}_n_matches\t"
            )?;
        }
    }
    writeln!(
        writer,
        "combination_status\tcombination_distance\tcombinations_in_library\tcombination_id\tcount"
    )?;

    // Write lines for each combination
    for comb in combinations.to_vector(sort) {
        // Fwd/Rev sequences if tracked
        let (fwd, rev) = match &comb.sequence {
            Some(seq) => match &seq.reverse {
                Some(rev) => (seq.forward.to_str_or_log(), rev.to_str_or_log()),
                None => (seq.forward.to_str_or_log(), "".to_string()),
            },
            None => ("".to_string(), "".to_string()),
        };

        // String describing the regions matches
        // Each region is described by:
        // {region} {region}_nearest {region}_variants {region}_distance {region}_n_matches
        let region_str = &combinations
            .region_ids
            .iter()
            .map(|k| match comb.regions.get(k) {
                Some(r) => {
                    let (seq, nearest, diff, dist, n) = r.lock().unwrap().to_strings();

                    if skip_variants {
                        format!("{seq}\t{nearest}\t{dist}\t{n}")
                    } else {
                        format!("{seq}\t{nearest}\t{diff}\t{dist}\t{n}")
                    }
                }
                None => {
                    if skip_variants {
                        "\t\t\t".to_string()
                    } else {
                        "\t\t\t\t".to_string()
                    }
                }
            })
            .join("\t");

        // String describing the combination match with the library
        let name = comb.library_matches.id_string()?;
        let comb_str = match &comb.library_matches {
            CombinationMatch::Uncompared => "uncompared\t\t\t".to_string(),
            CombinationMatch::Match { distance, .. } => format!("match\t{distance}\t1\t{name}"),
            CombinationMatch::MultiMatch { inds, distance } => {
                // Number of matches is the product of sub-library matches as all combinations possible
                let n_matches = inds
                    .iter()
                    .map(|x| match x {
                        Some(x) => x.len(),
                        None => 1,
                    })
                    .reduce(|x, y| x * y)
                    .unwrap_or(0);

                format!("match\t{}\t{}\t{}", distance, n_matches, name)
            }
            CombinationMatch::Recombination { distance } => {
                format!("recombination\t{distance}\t0\t",)
            }
            CombinationMatch::Mismatch => "mismatch\t\t0\t".to_string(),
            CombinationMatch::Nonmatch => "nonmatch\t\t0\t".to_string(),
        };

        // Write line for each group, adding group and count to const strings for the combination
        for (group, count) in comb.counts.iter() {
            writeln!(
                writer,
                "{group}\t{fwd}\t{rev}\t{region_str}\t{comb_str}\t{count}"
            )?;
        }
    }

    writer.flush()?;
    Ok(())
}

/// Write a compressed count TSV file with library matches
///
/// Write a compressed count TSV from an ObservedCombinations object covering only library matches
/// and merging hits at different distances for each match.
///
/// This TSV has the format:
/// - group: ReadGroup
/// - Per region:
///     - {region}: the matching library sequence sequence
/// - combination_status: Status of the match with the library
/// - combinations_in_library: Number of matches in the library
/// - combination_id: The library ID of the matching combinations
/// - count: Number of times observed
pub fn write_library_counts(
    combinations: &ObservedCombinations,
    file: File,
    sort: bool,
) -> Result<(), anyhow::Error> {
    let mut writer = BufWriter::new(file);

    // Write headers
    write!(writer, "group\t")?;
    for r in combinations.region_ids.iter() {
        write!(writer, "{}\t", region_id_to_str(r))?;
    }
    writeln!(
        writer,
        "combination_status\tcombinations_in_library\tcombination_id\tcount"
    )?;

    // Write combinations
    for comb in combinations.to_library_vector(sort)? {
        // Region and library strins constant for all readgroups
        let region_str = combinations
            .region_ids
            .iter()
            .map(|r| match comb.regions.get(r) {
                Some(r) => r.str_sequence(),
                None => "".to_string(),
            })
            .join("\t");

        let name = comb.library_matches.id_string()?;

        let lib_str = match &comb.library_matches {
            CombinationMatch::Uncompared => "uncompared\t\t".to_string(),
            CombinationMatch::Match { .. } => format!("match\t1\t{name}"),
            CombinationMatch::MultiMatch { inds, .. } => {
                // Number of matches is the product of sub-library matches as all combinations possible
                let n_matches = inds
                    .iter()
                    .map(|x| match x {
                        Some(x) => x.len(),
                        None => 1,
                    })
                    .reduce(|x, y| x * y)
                    .unwrap_or(0);

                format!("match\t{}\t{}", n_matches, name)
            }
            CombinationMatch::Recombination { .. } => "recombination\t0\t".to_string(),
            CombinationMatch::Mismatch => "mismatch\t0\t".to_string(),
            CombinationMatch::Nonmatch => "nonmatch\t0\t".to_string(),
        };

        // Write each readgroup line
        for (group, count) in comb.counts.iter() {
            writeln!(writer, "{group}\t{region_str}\t{lib_str}\t{count}")?;
        }
    }

    writer.flush()?;
    Ok(())
}

/// Write an overall summary TSV file
///
/// Write a summary TSV from a ReadSummary object, detailing how many reads were filtered
/// how many had each match status
///
/// This TSV has the format:
/// - group: Subsection of summary for group_proportion
/// - metric: The thing being counted (e.g. total, uncompared, exact_match etc.)
/// - count: Read count
/// - overall_proportion: proportion of all reads
/// - group_proportion: proportion of reads in that group
pub fn write_summary(summary: &ReadSummary, file: File) -> Result<(), anyhow::Error> {
    let mut writer = BufWriter::new(file);

    let filtered_reads = &summary.filtered_reads;

    let total = summary.total();
    let unfiltered_total = summary.total_unfiltered();
    let filtered_total = filtered_reads.total();

    // Headers
    writeln!(
        writer,
        "group\tmetric\tcount\toverall_proportion\tgroup_proportion"
    )?;

    // Overall total - use total_denom so proportions are 0.0 if no reads recieved
    write_summary_row(&mut writer, "all", "total", total, total, total)?;

    // Unfiltered reads
    let rows = [
        ("total", unfiltered_total),
        ("uncompared", summary.uncompared),
        ("exact_match", summary.exact_match),
        ("nearest_match", summary.nearest_match),
        ("multimatch", summary.multimatch),
        ("exact_recombination", summary.exact_recombination),
        ("nearest_recombination", summary.nearest_recombination),
        ("mismatch", summary.mismatch),
        ("nonmatch", summary.nonmatch),
    ];

    for (metric, count) in rows {
        write_summary_row(
            &mut writer,
            "unfiltered",
            metric,
            count,
            total,
            unfiltered_total,
        )?;
    }

    // Filtered total
    write_summary_row(
        &mut writer,
        "filtered",
        "total",
        filtered_total,
        total,
        filtered_total,
    )?;

    // Filtered reads
    for reason in FilterReason::ALL_FILTERS {
        write_summary_row(
            &mut writer,
            "filtered",
            reason.meta().id,
            filtered_reads.get(reason),
            total,
            filtered_total,
        )?;
    }

    writer.flush()?;
    Ok(())
}

#[inline]
fn write_summary_row<W: Write>(
    writer: &mut W,
    group: &str,
    metric: &str,
    count: u64,
    overall_total: u64,
    group_total: u64,
) -> std::io::Result<()> {
    writeln!(
        writer,
        "{group}\t{metric}\t{count}\t{:.4}\t{:.4}",
        div_or_zero(count as f32, overall_total as f32),
        div_or_zero(count as f32, group_total as f32),
    )
}

/// Write a TSV file summarising filtered reads
///
/// Write a summary TSV from a FilteredReads object, detailing which reads were filtered and
/// for what reason.
///
/// This TSV has the format:
/// - group: ReadGroup the read came from, if any
/// - forward: Forward read
/// - reverse: Reverse read (if paired)
/// - count: Read count
/// - proportion: Proportion of all filtered reads
/// - One column per filter in FilterReason::ALL_FILTERS, giving the count of reads filtered for that reason.
pub fn write_filter_summary(
    filtered_reads: &FilteredReads,
    file: File,
    sort: bool,
) -> Result<(), anyhow::Error> {
    let total = filtered_reads.total() as f32;

    let mut writer = BufWriter::new(file);

    // Write headers
    writeln!(
        writer,
        "group\tforward\treverse\tcount\tproportion\t{}",
        FilterReason::ALL_FILTERS
            .iter()
            .map(|r| r.meta().id)
            .join("\t")
    )?;

    // Write filtered reads
    for (read, group, counts) in filtered_reads.to_vector(sort) {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{:.4}\t{}",
            group,
            read.forward.to_str_or_log(),
            read.reverse
                .clone()
                .map_or("".to_string(), |x| x.to_str_or_log()),
            counts.total(),
            div_or_zero(counts.total() as f32, total),
            counts.iter().join("\t")
        )?;
    }

    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::NamedTempFile;

    use crate::combinations::{ObservedCombinations, ReadSummary};
    use crate::filters::FilteredReads;
    use crate::interning::{RegionID, region_id_from_str};
    use crate::library::Library;
    use crate::{DistanceMetric, FilterConfig};

    fn empty_filter_config() -> FilterConfig {
        FilterConfig::new(None, None, None, None, true)
    }

    fn empty_combs(r: Vec<RegionID>) -> ObservedCombinations {
        ObservedCombinations::new(r, empty_filter_config())
    }

    // Write summary row
    #[test]
    fn write_summary_row_normal() -> std::io::Result<()> {
        let mut buf = Vec::new();
        write_summary_row(&mut buf, "test_group", "metric_x", 100, 1000, 500)?;

        let output = String::from_utf8(buf).unwrap();
        assert_eq!(output, "test_group\tmetric_x\t100\t0.1000\t0.2000\n");
        Ok(())
    }

    #[test]
    fn write_summary_row_zero_count() -> std::io::Result<()> {
        let mut buf = Vec::new();
        write_summary_row(&mut buf, "group", "metric", 0, 100, 100)?;

        let output = String::from_utf8(buf).unwrap();
        assert_eq!(output, "group\tmetric\t0\t0.0000\t0.0000\n");
        Ok(())
    }

    #[test]
    fn write_summary_row_count_equals_total() -> std::io::Result<()> {
        let mut buf = Vec::new();
        write_summary_row(&mut buf, "group", "metric", 100, 100, 100)?;

        let output = String::from_utf8(buf).unwrap();
        assert_eq!(output, "group\tmetric\t100\t1.0000\t1.0000\n");
        Ok(())
    }

    #[test]
    fn write_summary_row_zero_denominators() -> std::io::Result<()> {
        let mut buf = Vec::new();
        write_summary_row(&mut buf, "group", "metric", 0, 0, 0)?;

        let output = String::from_utf8(buf).unwrap();
        assert_eq!(output, "group\tmetric\t0\t0.0000\t0.0000\n");
        Ok(())
    }

    #[test]
    fn write_summary_row_large_numbers() -> std::io::Result<()> {
        let mut buf = Vec::new();
        write_summary_row(
            &mut buf, "group", "metric", 1_000_000, 10_000_000, 5_000_000,
        )?;

        let output = String::from_utf8(buf).unwrap();
        assert_eq!(output, "group\tmetric\t1000000\t0.1000\t0.2000\n");
        Ok(())
    }

    #[test]
    fn write_summary_row_fractional_results() -> std::io::Result<()> {
        let mut buf = Vec::new();
        write_summary_row(&mut buf, "group", "metric", 1, 3, 7)?;

        let output = String::from_utf8(buf).unwrap();
        assert_eq!(output, "group\tmetric\t1\t0.3333\t0.1429\n");
        Ok(())
    }

    // Write counts header
    #[test]
    fn write_counts_header_structure() -> anyhow::Result<()> {
        let file = NamedTempFile::new()?;

        let combinations = empty_combs(vec![
            region_id_from_str("region1"),
            region_id_from_str("region2"),
        ]);

        write_counts(&combinations, file.reopen()?, false, false)?;

        let mut content = String::new();
        file.reopen()?.read_to_string(&mut content)?;

        let header = content.lines().next().unwrap();
        assert_eq!(
            header,
            "group\tforward\treverse\tregion1\tregion1_nearest\tregion1_variants\tregion1_distance\tregion1_n_matches\tregion2\tregion2_nearest\tregion2_variants\tregion2_distance\tregion2_n_matches\tcombination_status\tcombination_distance\tcombinations_in_library\tcombination_id\tcount"
        );

        Ok(())
    }

    // Write counts header
    #[test]
    fn write_counts_header_structure_no_variants() -> anyhow::Result<()> {
        let file = NamedTempFile::new()?;

        let combinations = empty_combs(vec![
            region_id_from_str("region1"),
            region_id_from_str("region2"),
        ]);

        write_counts(&combinations, file.reopen()?, false, true)?;

        let mut content = String::new();
        file.reopen()?.read_to_string(&mut content)?;

        let header = content.lines().next().unwrap();
        assert_eq!(
            header,
            "group\tforward\treverse\tregion1\tregion1_nearest\tregion1_distance\tregion1_n_matches\tregion2\tregion2_nearest\tregion2_distance\tregion2_n_matches\tcombination_status\tcombination_distance\tcombinations_in_library\tcombination_id\tcount"
        );

        Ok(())
    }

    // Library counts
    #[test]
    fn write_library_counts_header_structure() -> anyhow::Result<()> {
        let file = NamedTempFile::new()?;

        let mut combinations = empty_combs(vec![region_id_from_str("region1")]);

        combinations.compare_to_library(
            Library::new(Vec::new()).unwrap(),
            None,
            DistanceMetric::Hamming,
            1,
            false,
            1,
        )?;

        write_library_counts(&combinations, file.reopen()?, false)?;

        let mut content = String::new();
        file.reopen()?.read_to_string(&mut content)?;

        let header = content.lines().next().unwrap();
        assert_eq!(
            header,
            "group\tregion1\tcombination_status\tcombinations_in_library\tcombination_id\tcount"
        );

        Ok(())
    }

    // Summary
    #[test]
    fn write_summary_basic_structure() -> anyhow::Result<()> {
        let file = NamedTempFile::new()?;
        let summary = ReadSummary::empty();

        write_summary(&summary, file.reopen()?)?;

        let mut content = String::new();
        file.reopen()?.read_to_string(&mut content)?;

        let lines: Vec<&str> = content.lines().collect();
        // All lines should always be present
        assert_eq!(lines.len(), 17);

        let header = lines[0];
        assert_eq!(
            header,
            "group\tmetric\tcount\toverall_proportion\tgroup_proportion"
        );

        Ok(())
    }

    #[test]
    fn write_summary_all_rows_present() -> anyhow::Result<()> {
        let file = NamedTempFile::new()?;
        let mut summary = ReadSummary::empty();
        summary.exact_match = 50;
        summary.nearest_match = 30;
        summary.mismatch = 20;

        write_summary(&summary, file.reopen()?)?;

        let mut content = String::new();
        file.reopen()?.read_to_string(&mut content)?;

        let lines: Vec<&str> = content.lines().collect();

        assert_eq!(lines.len(), 17);

        assert!(content.contains("50"));
        assert!(content.contains("30"));
        assert!(content.contains("20"));
        assert!(content.contains("100"));

        Ok(())
    }

    // Filter summary
    #[test]
    fn write_filter_summary_header() -> anyhow::Result<()> {
        let file = NamedTempFile::new()?;
        let filtered = FilteredReads::new(empty_filter_config());

        write_filter_summary(&filtered, file.reopen()?, false)?;

        let mut content = String::new();
        file.reopen()?.read_to_string(&mut content)?;

        let header = content.lines().next().unwrap();
        assert_eq!(
            header,
            "group\tforward\treverse\tcount\tproportion\tempty_read\tshort_read\tlong_read	low_mean_quality\tbad_alignment"
        );

        Ok(())
    }

    #[test]
    fn write_filter_summary_empty() -> anyhow::Result<()> {
        let file = NamedTempFile::new()?;
        let filtered = FilteredReads::new(empty_filter_config());

        write_filter_summary(&filtered, file.reopen()?, false)?;

        let mut content = String::new();
        file.reopen()?.read_to_string(&mut content)?;

        let lines: Vec<&str> = content.lines().collect();
        // Should have header + no data lines if empty
        assert_eq!(lines.len(), 1);

        Ok(())
    }
}
