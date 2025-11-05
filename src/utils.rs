//! Utility fucntions for use throughout the library
//!
//! Provides a range of utility functions that are needed across modules

/// Calculate the mean of a fastq quality vector
pub fn mean_quality(qual: &[u8]) -> f32 {
    if qual.is_empty() {
        return 0.0;
    };

    let total: u32 = qual.iter().fold(0, |a, e| a + *e as u32);
    total as f32 / qual.len() as f32 - 33.0 // Subtract 33 as Phred scores are shifted 33 in byte codepoints
}

/// Divide x / y or return 0.0 if y == 0
pub fn div_or_zero(x: f32, y: f32) -> f32 {
    if y == 0.0 {
        return 0.0;
    }

    x / y
}
