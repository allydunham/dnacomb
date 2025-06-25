//! Filters to identify anb count failing reads
//!
//! Provides a range of filtering critera that a function
//! processing reads can tap into via a simple config and
//! interface, easily filtering reads for a range of reasons.
//! Also makes it easy to keep track of filtered reads and
//! print a summary to file.
use crate::{errors::ReadCountError, parsing::ReadPair};
use bio::bio_types::alignment::Alignment;

/// Calculate the mean of a fastq quality vector
pub fn mean_quality(qual: &[u8]) -> f32 {
    let total: u32 = qual.iter().fold(0, |a, e| a + *e as u32);
    total as f32 / qual.len() as f32
}

/// Container for filtered reads
///
/// Keeps track of counts for filtered reads during counting
#[derive(Debug, Clone)]
pub struct FilteredReads {
    config: FilterConfig,
    low_mean_quality: u64,
    bad_alignment: u64,
}

pub enum FilterReason {
    None,
    LowMeanQuality,
    BadAlignment,
}

impl FilteredReads {
    pub fn new(config: FilterConfig) -> Self {
        Self {
            config,
            low_mean_quality: 0,
            bad_alignment: 0,
        }
    }

    /// Update the count without checking a read
    ///
    /// This lets you update counts if manually checking or when
    /// using a cached result
    pub fn update_count(&mut self, reason: &FilterReason) {
        match reason {
            FilterReason::None => {}
            FilterReason::LowMeanQuality => self.low_mean_quality += 1,
            FilterReason::BadAlignment => self.bad_alignment += 1,
        };
    }

    /// Determine if a readpair should be filtered based on the supplied config
    ///
    /// Checks whether the read should be filtered, adding it to the appropriate count if so, and
    /// returns a bool determining if it was filtered.
    pub fn filter_readpair(&mut self, record: &ReadPair) -> FilterReason {
        let f_read = &record.forward;
        let r_read = match &record.reverse {
            Some(r) => Some(r),
            None => None,
        };

        // Check mean quality is high enough
        match (self.config.mean_quality_threshold, r_read) {
            (None, _) => {}
            (Some(q), None) => {
                if mean_quality(f_read.qual()) < q {
                    self.low_mean_quality += 1;
                    return FilterReason::LowMeanQuality;
                }
            }
            (Some(q), Some(r)) => {
                if mean_quality(f_read.qual()) < q || mean_quality(r.qual()) < q {
                    self.low_mean_quality += 1;
                    return FilterReason::LowMeanQuality;
                }
            }
        }

        FilterReason::None
    }

    /// Determine if an alignment should be filtered based on the supplied config
    ///
    /// Checks if an alignment should be filtered, adding it to the appropriate count if so, and
    /// returns a bool determining if it was filtered
    pub fn filter_alignment(
        &mut self,
        f_alignment: &Alignment,
        r_alignment: Option<&Alignment>,
    ) -> FilterReason {
        match &self.config.alignment_tolerance {
            None => {}
            Some(t) => {
                if f_alignment.score < t.minimum_f_score {
                    self.bad_alignment += 1;
                    return FilterReason::BadAlignment;
                }

                if let Some(r_al) = r_alignment {
                    if r_al.score < t.minimum_r_score {
                        self.bad_alignment += 1;
                        return FilterReason::BadAlignment;
                    }
                }
            }
        }

        FilterReason::None
    }

    /// Total count of filtered reads
    pub fn total(&self) -> u64 {
        self.low_mean_quality + self.bad_alignment
    }

    /// Generate a string of TSV lines representing the filter counts
    ///
    /// Creates a TSV string with columns for filter reason, count and
    /// proportion of the supplied total read count. Mainly for use
    /// wihen outputing from ReadSummary
    pub fn to_tsv(&self, total: u64) -> String {
        let filtered_total = self.total();

        let mut out = String::with_capacity(300);

        out.push_str(&format!(
            "filtered\ttotal\t{}\t{:.4}\t1.000\n",
            filtered_total,
            filtered_total as f32 / total as f32,
        ));

        out.push_str(&format!(
            "filtered\tlow_mean_quality\t{}\t{:.4}\t{:.4}\n",
            self.low_mean_quality,
            self.low_mean_quality as f32 / total as f32,
            self.low_mean_quality as f32 / filtered_total as f32,
        ));

        out.push_str(&format!(
            "filtered\tbad_alignment\t{}\t{:.4}\t{:.4}\n",
            self.bad_alignment,
            self.bad_alignment as f32 / total as f32,
            self.bad_alignment as f32 / filtered_total as f32,
        ));

        out
    }
}

/// Configuration for filtering
///
/// Instructions for how to filter reads
#[derive(Debug, Clone)]
pub struct FilterConfig {
    mean_quality_threshold: Option<f32>,
    alignment_tolerance: Option<AlignmentTolerance>,
}

impl FilterConfig {
    pub fn new(
        mean_quality_threshold: Option<f32>,
        alignment_tolerance: Option<AlignmentTolerance>,
    ) -> Self {
        Self {
            mean_quality_threshold,
            alignment_tolerance,
        }
    }
}

#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct AlignmentTolerance {
    tolerance: f32,
    expected_f_score: i32,
    expected_r_score: i32,
    minimum_f_score: i32,
    minimum_r_score: i32,
}

impl AlignmentTolerance {
    pub fn new(
        tolerance: f32,
        expected_f_score: i32,
        expected_r_score: i32,
    ) -> Result<Self, ReadCountError> {
        if !(0.0..=1.0).contains(&tolerance) {
            return Err(ReadCountError::FilterConfigError {
                desc: "Alignment tolerance must be between 0 and 1".to_string(),
            });
        };

        Ok(Self {
            tolerance,
            expected_f_score,
            expected_r_score,
            minimum_f_score: (expected_f_score as f32 * tolerance) as i32,
            minimum_r_score: (expected_r_score as f32 * tolerance) as i32,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

}

