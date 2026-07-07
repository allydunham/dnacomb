//! Small utility helpers used across DNAComb.

/// Compute the mean Phred quality score from FASTQ quality bytes.
///
/// Input qualities are assumed to use Phred+33 encoding.
pub fn mean_quality(qual: &[u8]) -> f32 {
    if qual.is_empty() {
        return 0.0;
    };

    let total: u32 = qual.iter().fold(0, |a, e| a + *e as u32);
    total as f32 / qual.len() as f32 - 33.0 // Subtract 33 as Phred scores are shifted 33 in byte codepoints
}

/// Divide `x / y`, returning `0.0` when `y == 0.0`.
///
/// This is mainly used when reporting proportions in summary/output tables.
pub fn div_or_zero(x: f32, y: f32) -> f32 {
    if y == 0.0 {
        return 0.0;
    }

    x / y
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mean_quality() {
        assert_eq!(mean_quality(b"FFFF"), 37.0);
        assert_eq!(mean_quality(b"AAAA"), 32.0);
        assert_eq!(mean_quality(b"!!!!"), 0.0);
        assert_eq!(mean_quality(b"0101"), 15.5);
        assert_eq!(mean_quality(b"F"), 37.0);
        assert_eq!(mean_quality(b"!"), 0.0);
        assert_eq!(mean_quality(b"~"), 93.0);
        assert_eq!(mean_quality(b"!~"), 46.5);
        assert_eq!(mean_quality(b""), 0.0);
    }

    #[test]
    fn test_mean_quality_very_long_string() {
        let long_qual = vec![b'A'; 100_000];
        assert_eq!(mean_quality(&long_qual), 32.0);
    }

    #[test]
    fn test_mean_quality_precision() {
        // Test that floating-point precision is maintained
        let qual = b"ABC";
        let expected = ((b'A' as f32 + b'B' as f32 + b'C' as f32) / 3.0) - 33.0;
        assert!((mean_quality(qual) - expected).abs() < 1e-6);
    }

    #[test]
    fn test_division_by_nonzero() {
        assert_eq!(div_or_zero(0.0, 5.0), 0.0);
        assert_eq!(div_or_zero(0.0, 1.0), 0.0);
        assert_eq!(div_or_zero(5.0, 1.0), 5.0);
        assert_eq!(div_or_zero(10.0, 2.0), 5.0);
        assert_eq!(div_or_zero(24.0, 3.0), 24.0 / 3.0);
        assert_eq!(div_or_zero(13.0, 2.5), 13.0 / 2.5);
    }

    #[test]
    fn test_division_by_zero() {
        assert_eq!(div_or_zero(1.0, 0.0), 0.0);
        assert_eq!(div_or_zero(1.0, 0.0), 0.0);
        assert_eq!(div_or_zero(100.0, 0.0), 0.0);
        assert_eq!(div_or_zero(0.0, 0.0), 0.0);
    }

    #[test]
    fn test_div_or_zero_very_small_denominator() {
        // Non-zero but very small
        let result = div_or_zero(1.0, 1e-10);
        assert!(result > 1e9);
    }

    #[test]
    fn test_div_or_zero_very_large_numbers() {
        assert_eq!(div_or_zero(1e10, 1e5), 1e5);
    }

    #[test]
    fn test_div_or_zero_infinity() {
        let result = div_or_zero(f32::INFINITY, 2.0);
        assert!(result.is_infinite());
    }

    #[test]
    fn test_div_or_zero_nan() {
        let result = div_or_zero(f32::NAN, 2.0);
        assert!(result.is_nan());
    }

    #[test]
    fn test_div_or_zero_result_very_close_to_zero() {
        // Non-zero result that's very small
        let result = div_or_zero(1e-10, 1.0);
        assert!(result > 0.0);
        assert!(result < 1e-9);
    }
}
