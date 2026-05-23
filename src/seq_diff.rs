//! Sequence-difference tracking for observed versus expected DNA sequences.
//!
//! This module computes edit operations between an observed sequence and an
//! expected/library sequence and formats those differences in a compact,
//! HGVS-inspired string representation for output.
use std::cell::RefCell;
use std::fmt::{self, Display};

use bio::alignment::pairwise::{Aligner, MatchFunc, Scoring};
use bio::alignment::{Alignment, AlignmentOperation};

use crate::interning::{SeqHandle, seq_to_bytes};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TerminalFilter {
    None,
    Leading,
    Trailing,
    Both,
}

#[derive(Clone, Copy, Debug)]
struct DiffScoring;

impl MatchFunc for DiffScoring {
    #[inline]
    fn score(&self, a: u8, b: u8) -> i32 {
        if a == b { 1 } else { -1 }
    }
}

thread_local! {
    static ALIGNER: RefCell<Aligner<DiffScoring>> = RefCell::new(
        Aligner::with_capacity_and_scoring(
            64,
            64,
            Scoring::new(-5, -1, DiffScoring),
        )
    );
}

/// One edit operation transforming an expected sequence into an observed sequence.
#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub enum EditOperation {
    /// Substitution at the given zero-based expected-sequence position:
    /// `(position, expected_base, observed_base)`.
    Sub(usize, u8, u8),

    /// Insertion in the observed sequence after the given expected-sequence position.
    Ins(usize, Vec<u8>),

    /// Deletion from the expected sequence starting at the given zero-based position.
    Del(usize, Vec<u8>),
}

impl EditOperation {
    /// Format this edit operation in a compact HGVS-like string form.
    ///
    /// Substitutions are reported using 1-based positions. Insertions and
    /// deletions are reported as interval-style events relative to the expected
    /// sequence coordinates.
    pub fn to_hgvs_string(&self) -> String {
        match self {
            EditOperation::Sub(pos, exp, obs) => {
                format!("{}{}>{}", pos + 1, *exp as char, *obs as char)
            }
            EditOperation::Ins(pos, seq) => {
                format!("{}_{}_ins{}", pos, pos + 1, String::from_utf8_lossy(seq))
            }
            EditOperation::Del(pos, seq) => {
                format!(
                    "{}_{}_del{}",
                    pos,
                    pos + seq.len() - 1,
                    String::from_utf8_lossy(seq)
                )
            }
        }
    }
}

/// Collection of edit operations describing the difference between an observed
/// sequence and an expected sequence.
#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub struct SequenceDiff {
    pub operations: Vec<EditOperation>,
}

impl SequenceDiff {
    /// Construct a sequence diff from a precomputed list of edit operations.
    pub fn new(operations: Vec<EditOperation>) -> Self {
        Self { operations }
    }

    /// Compute diff between observed and expected sequences from SeqHandles
    #[inline]
    pub fn compute_ids(
        observed: &SeqHandle,
        expected: &SeqHandle,
        terminal_filter: TerminalFilter,
    ) -> Self {
        Self::compute(
            &seq_to_bytes(observed),
            &seq_to_bytes(expected),
            terminal_filter,
        )
    }

    /// Compute the edit operations needed to describe an observed sequence
    /// relative to an expected sequence.
    ///
    /// A global alignment is used to derive substitutions, insertions, and
    /// deletions. Consecutive insertion or deletion operations are merged into
    /// single multi-base events where possible. Positions are reported relative to the expected sequence.
    pub fn compute(observed: &[u8], expected: &[u8], terminal_filter: TerminalFilter) -> Self {
        if observed.is_empty() && expected.is_empty() || observed == expected {
            return Self::new(Vec::new());
        }

        let alignment: Alignment =
            ALIGNER.with(|aligner_cell| aligner_cell.borrow_mut().global(observed, expected));

        let mut edits: Vec<EditOperation> = Vec::new();

        let mut obs_pos = alignment.xstart;
        let mut exp_pos = alignment.ystart;

        // Currently accumalating indel operation
        let mut pending: Option<EditOperation> = None;

        for op in alignment.operations {
            match op {
                AlignmentOperation::Match => {
                    if let Some(op) = pending.take() {
                        edits.push(op);
                    }
                    obs_pos += 1;
                    exp_pos += 1;
                }

                AlignmentOperation::Subst => {
                    if let Some(op) = pending.take() {
                        edits.push(op);
                    }
                    edits.push(EditOperation::Sub(
                        exp_pos,
                        expected[exp_pos],
                        observed[obs_pos],
                    ));
                    obs_pos += 1;
                    exp_pos += 1;
                }

                AlignmentOperation::Del => {
                    match pending.as_mut() {
                        Some(EditOperation::Sub(..)) => {
                            panic!("pending never set to Sub")
                        }
                        Some(EditOperation::Ins(..)) => {
                            let old = pending
                                .replace(EditOperation::Del(exp_pos, vec![expected[exp_pos]]));
                            edits.push(old.expect("pending was Some"));
                        }
                        Some(EditOperation::Del(_, items)) => {
                            items.push(expected[exp_pos]);
                        }
                        None => {
                            pending = Some(EditOperation::Del(exp_pos, vec![expected[exp_pos]]));
                        }
                    }

                    exp_pos += 1;
                }

                // Gap in observed / base(s) present in expected => insertion into observed
                AlignmentOperation::Ins => {
                    match pending.as_mut() {
                        Some(EditOperation::Sub(..)) => {
                            panic!("pending never set to Sub")
                        }
                        Some(EditOperation::Del(..)) => {
                            let old = pending
                                .replace(EditOperation::Ins(exp_pos, vec![observed[obs_pos]]));
                            edits.push(old.expect("pending was Some"));
                        }
                        Some(EditOperation::Ins(_, items)) => {
                            items.push(observed[obs_pos]);
                        }
                        None => {
                            pending = Some(EditOperation::Ins(exp_pos, vec![observed[obs_pos]]));
                        }
                    }

                    obs_pos += 1;
                }

                AlignmentOperation::Xclip(len) => {
                    if let Some(edit) = pending.take() {
                        edits.push(edit);
                    };
                    obs_pos += len;
                }

                AlignmentOperation::Yclip(len) => {
                    if let Some(edit) = pending.take() {
                        edits.push(edit);
                    };
                    exp_pos += len;
                }
            }
        }

        // Flush final op
        if let Some(edit) = pending.take() {
            edits.push(edit);
        };

        if !matches!(terminal_filter, TerminalFilter::None) {
            // Remove leading/trailing deletions if appropriate
            edits.retain(|op| {
                let leading = matches!(op, EditOperation::Del(0, _));

                let trailing = match op {
                    EditOperation::Del(pos, seq) => *pos + seq.len() == expected.len(),
                    EditOperation::Sub(..) | EditOperation::Ins(..) => false,
                };

                match terminal_filter {
                    TerminalFilter::None => true,
                    TerminalFilter::Leading => !leading,
                    TerminalFilter::Trailing => !trailing,
                    TerminalFilter::Both => !leading && !trailing,
                }
            });
        }

        Self::new(edits)
    }

    /// Offset all positions
    ///
    /// Useful if a diff is calculated against a subset of a full sequence, for
    /// instance for alignment anchored at one end.
    pub fn offset_expected_positions(mut self, offset: usize) -> Self {
        for op in &mut self.operations {
            match op {
                EditOperation::Sub(pos, ..)
                | EditOperation::Ins(pos, _)
                | EditOperation::Del(pos, _) => *pos += offset,
            }
        }
        self
    }
}

/// Display as `;`-separated HGVS-like edit operations.
impl Display for SequenceDiff {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            self.operations
                .iter()
                .map(|op| op.to_hgvs_string())
                .collect::<Vec<_>>()
                .join(";")
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn edit_operation_to_hgvs_string() {
        #[derive(Debug)]
        struct Case {
            name: &'static str,
            op: EditOperation,
            expected: String,
        }

        let cases = vec![
            Case {
                name: "substitution at position 0",
                op: EditOperation::Sub(0, b'A', b'T'),
                expected: "1A>T".to_string(),
            },
            Case {
                name: "substitution at position 5",
                op: EditOperation::Sub(5, b'G', b'C'),
                expected: "6G>C".to_string(),
            },
            Case {
                name: "insertion at start",
                op: EditOperation::Ins(0, vec![b'A']),
                expected: "0_1_insA".to_string(),
            },
            Case {
                name: "insertion of multiple bases",
                op: EditOperation::Ins(3, vec![b'G', b'T']),
                expected: "3_4_insGT".to_string(),
            },
            Case {
                name: "deletion of single base",
                op: EditOperation::Del(2, vec![b'C']),
                expected: "2_2_delC".to_string(),
            },
            Case {
                name: "deletion of multiple bases",
                op: EditOperation::Del(4, vec![b'A', b'T', b'G']),
                expected: "4_6_delATG".to_string(),
            },
        ];

        for c in cases {
            let got = c.op.to_hgvs_string();
            assert_eq!(got, c.expected, "Unexpected HGVS string (case: {})", c.name);
        }
    }

    #[test]
    fn sequence_diff_exact_match() {
        #[derive(Debug)]
        struct Case {
            name: &'static str,
            observed: &'static [u8],
            expected: &'static [u8],
            expected_ops_len: usize,
            expected_str: &'static str,
        }

        let cases = vec![
            Case {
                name: "identical sequences",
                observed: b"ACGT",
                expected: b"ACGT",
                expected_ops_len: 0,
                expected_str: "",
            },
            Case {
                name: "empty sequences",
                observed: b"",
                expected: b"",
                expected_ops_len: 0,
                expected_str: "",
            },
            Case {
                name: "long identical sequence",
                observed: b"ACGTACGTACGTACGT",
                expected: b"ACGTACGTACGTACGT",
                expected_ops_len: 0,
                expected_str: "",
            },
        ];

        for c in cases {
            let diff = SequenceDiff::compute(c.observed, c.expected, TerminalFilter::None);
            assert_eq!(
                diff.operations.len(),
                c.expected_ops_len,
                "Unexpected operation count (case: {})",
                c.name
            );
            assert_eq!(
                diff.to_string(),
                c.expected_str,
                "Unexpected diff string (case: {})",
                c.name
            );
        }
    }

    #[test]
    fn sequence_diff_substitutions() {
        #[derive(Debug)]
        struct Case {
            name: &'static str,
            observed: &'static [u8],
            expected: &'static [u8],
            expected_str: &'static str,
        }

        let cases = vec![
            Case {
                name: "single substitution at start",
                observed: b"ACGT",
                expected: b"TCGT",
                expected_str: "1T>A",
            },
            Case {
                name: "single substitution at end",
                observed: b"ACGT",
                expected: b"ACGA",
                expected_str: "4A>T",
            },
            Case {
                name: "single substitution in middle",
                observed: b"ACGT",
                expected: b"ACTT",
                expected_str: "3T>G",
            },
            Case {
                name: "multiple substitutions",
                observed: b"TATT",
                expected: b"AAGT",
                expected_str: "1A>T;3G>T",
            },
            Case {
                name: "all positions different",
                observed: b"AAAA",
                expected: b"TTTT",
                expected_str: "1T>A;2T>A;3T>A;4T>A",
            },
        ];

        for c in cases {
            let diff = SequenceDiff::compute(c.observed, c.expected, TerminalFilter::None);
            assert_eq!(
                diff.to_string(),
                c.expected_str,
                "Unexpected diff string (case: {})",
                c.name
            );
        }
    }

    #[test]
    fn sequence_diff_insertions() {
        #[derive(Debug)]
        struct Case {
            name: &'static str,
            observed: &'static [u8],
            expected: &'static [u8],
            expected_str: &'static str,
        }

        let cases = vec![
            Case {
                name: "insertion at start",
                observed: b"ACGT",
                expected: b"CGT",
                expected_str: "0_1_insA",
            },
            Case {
                name: "insertion at end",
                observed: b"ACGT",
                expected: b"ACG",
                expected_str: "3_4_insT",
            },
            Case {
                name: "insertion in middle",
                observed: b"ACGT",
                expected: b"AGT",
                expected_str: "1_2_insC",
            },
            Case {
                name: "multiple insertions",
                observed: b"ACGT",
                expected: b"AT",
                expected_str: "1_2_insCG",
            },
            Case {
                name: "insertion of multiple bases at once",
                observed: b"ACGTG",
                expected: b"A",
                expected_str: "1_2_insCGTG",
            },
        ];

        for c in cases {
            let diff = SequenceDiff::compute(c.observed, c.expected, TerminalFilter::None);
            assert_eq!(
                diff.to_string(),
                c.expected_str,
                "Unexpected diff string (case: {})",
                c.name
            );
        }
    }

    #[test]
    fn sequence_diff_deletions() {
        #[derive(Debug)]
        struct Case {
            name: &'static str,
            observed: &'static [u8],
            expected: &'static [u8],
            expected_str: &'static str,
        }

        let cases = vec![
            Case {
                name: "deletion at start",
                observed: b"ACGT",
                expected: b"TACGT",
                expected_str: "0_0_delT",
            },
            Case {
                name: "deletion at end",
                observed: b"ACGTAGT",
                expected: b"ACGTTAGT",
                expected_str: "3_3_delT",
            },
            Case {
                name: "deletion in middle",
                observed: b"ACGT",
                expected: b"ACTGT",
                expected_str: "2_2_delT",
            },
            Case {
                name: "multiple deletions",
                observed: b"TCGTTGGCCTAG",
                expected: b"TACGTTGGCCTGGAG",
                expected_str: "1_1_delA;11_12_delGG",
            },
            Case {
                name: "deletion of multiple bases at once",
                observed: b"ACGT",
                expected: b"ACGTTAGG",
                expected_str: "4_7_delTAGG",
            },
        ];

        for c in cases {
            let diff = SequenceDiff::compute(c.observed, c.expected, TerminalFilter::None);
            assert_eq!(
                diff.to_string(),
                c.expected_str,
                "Unexpected diff string (case: {})",
                c.name
            );
        }
    }

    #[test]
    fn sequence_diff_mixed_operations() {
        #[derive(Debug)]
        struct Case {
            name: &'static str,
            observed: &'static [u8],
            expected: &'static [u8],
            expected_str: &'static str,
        }

        let cases = vec![
            Case {
                name: "substitution and insertion",
                observed: b"ACGTA",
                expected: b"TCGT",
                expected_str: "1T>A;4_5_insA",
            },
            Case {
                name: "substitution and deletion",
                observed: b"TTTCGAGGCAGCA",
                expected: b"TACGTTCGATGCAGCA",
                expected_str: "1_3_delACG;10T>G",
            },
            Case {
                name: "insertion and deletion",
                observed: b"ACGTGCGCGACTAGAGTCCCTAG",
                expected: b"ACGCAATGCGCGACTAGCCCTAG",
                expected_str: "3_5_delCAA;17_18_insAGT",
            },
        ];

        for c in cases {
            let diff = SequenceDiff::compute(c.observed, c.expected, TerminalFilter::None);
            assert_eq!(
                diff.to_string(),
                c.expected_str,
                "Unexpected diff string (case: {})",
                c.name
            );
        }
    }

    #[test]
    fn sequence_diff_edge_cases() {
        #[derive(Debug)]
        struct Case {
            name: &'static str,
            observed: &'static [u8],
            expected: &'static [u8],
            expected_ops_len: usize,
        }

        let cases = vec![
            Case {
                name: "observed empty, expected non-empty",
                observed: b"",
                expected: b"ACGT",
                expected_ops_len: 1,
            },
            Case {
                name: "observed non-empty, expected empty",
                observed: b"ACGT",
                expected: b"",
                expected_ops_len: 1,
            },
            Case {
                name: "single base identical",
                observed: b"A",
                expected: b"A",
                expected_ops_len: 0,
            },
            Case {
                name: "single base different",
                observed: b"A",
                expected: b"T",
                expected_ops_len: 1,
            },
            Case {
                name: "long sequence identical",
                observed: b"ACGTACGTACGTACGTACGTACGTACGT",
                expected: b"ACGTACGTACGTACGTACGTACGTACGT",
                expected_ops_len: 0,
            },
        ];

        for c in cases {
            let diff = SequenceDiff::compute(c.observed, c.expected, TerminalFilter::None);
            assert_eq!(
                diff.operations.len(),
                c.expected_ops_len,
                "Unexpected operation count (case: {})",
                c.name
            );
        }
    }

    #[test]
    fn sequence_diff_display() {
        let ops = vec![
            EditOperation::Sub(0, b'A', b'T'),
            EditOperation::Ins(2, vec![b'G']),
            EditOperation::Del(5, vec![b'C']),
        ];
        let diff = SequenceDiff::new(ops);

        assert_eq!(
            diff.to_string(),
            "1A>T;2_3_insG;5_5_delC",
            "Display format incorrect"
        );
    }

    #[test]
    fn edit_operation_equality_sub() {
        let op1 = EditOperation::Sub(0, b'A', b'T');
        let op2 = EditOperation::Sub(0, b'A', b'T');
        let op3 = EditOperation::Sub(1, b'A', b'T');
        assert_eq!(op1, op2);
        assert_ne!(op1, op3);
    }

    #[test]
    fn edit_operation_equality_ins() {
        let op1 = EditOperation::Ins(2, vec![b'G', b'T']);
        let op2 = EditOperation::Ins(2, vec![b'G', b'T']);
        let op3 = EditOperation::Ins(2, vec![b'G']);
        assert_eq!(op1, op2);
        assert_ne!(op1, op3);
    }

    #[test]
    fn edit_operation_equality_del() {
        let op1 = EditOperation::Del(5, vec![b'C', b'A']);
        let op2 = EditOperation::Del(5, vec![b'C', b'A']);
        let op3 = EditOperation::Del(5, vec![b'C']);
        assert_eq!(op1, op2);
        assert_ne!(op1, op3);
    }

    #[test]
    fn edit_operation_hash_consistency() {
        use std::collections::HashSet;

        let op1 = EditOperation::Sub(0, b'A', b'T');
        let op2 = EditOperation::Sub(0, b'A', b'T');

        let mut set = HashSet::new();
        set.insert(op1);
        assert!(set.contains(&op2));
    }

    #[test]
    fn edit_operation_clone() {
        let op1 = EditOperation::Del(5, vec![b'A', b'T', b'G']);
        let op2 = op1.clone();
        assert_eq!(op1, op2);
    }

    #[test]
    fn sequence_diff_equality_empty() {
        let diff1 = SequenceDiff::new(vec![]);
        let diff2 = SequenceDiff::new(vec![]);
        assert_eq!(diff1, diff2);
    }

    #[test]
    fn sequence_diff_equality_same_ops() {
        let ops1 = vec![
            EditOperation::Sub(0, b'A', b'T'),
            EditOperation::Ins(2, vec![b'G']),
        ];
        let ops2 = vec![
            EditOperation::Sub(0, b'A', b'T'),
            EditOperation::Ins(2, vec![b'G']),
        ];
        let diff1 = SequenceDiff::new(ops1);
        let diff2 = SequenceDiff::new(ops2);
        assert_eq!(diff1, diff2);
    }

    #[test]
    fn sequence_diff_inequality_different_ops() {
        let ops1 = vec![EditOperation::Sub(0, b'A', b'T')];
        let ops2 = vec![EditOperation::Sub(1, b'A', b'T')];
        let diff1 = SequenceDiff::new(ops1);
        let diff2 = SequenceDiff::new(ops2);
        assert_ne!(diff1, diff2);
    }

    #[test]
    fn sequence_diff_inequality_different_order() {
        let ops1 = vec![
            EditOperation::Sub(0, b'A', b'T'),
            EditOperation::Ins(2, vec![b'G']),
        ];
        let ops2 = vec![
            EditOperation::Ins(2, vec![b'G']),
            EditOperation::Sub(0, b'A', b'T'),
        ];
        let diff1 = SequenceDiff::new(ops1);
        let diff2 = SequenceDiff::new(ops2);
        assert_ne!(diff1, diff2);
    }

    #[test]
    fn sequence_diff_hash_consistency() {
        use std::collections::HashSet;

        let ops = vec![EditOperation::Sub(0, b'A', b'T')];
        let diff1 = SequenceDiff::new(ops.clone());
        let diff2 = SequenceDiff::new(ops);

        let mut set = HashSet::new();
        set.insert(diff1);
        assert!(set.contains(&diff2));
    }

    #[test]
    fn sequence_diff_clone() {
        let ops = vec![
            EditOperation::Sub(0, b'A', b'T'),
            EditOperation::Del(5, vec![b'C']),
        ];
        let diff1 = SequenceDiff::new(ops);
        let diff2 = diff1.clone();
        assert_eq!(diff1, diff2);
    }

    #[test]
    fn sequence_diff_compute_ids() {
        use crate::interning::seq_from_bytes;

        let obs_handle = seq_from_bytes(b"ACGT");
        let exp_handle = seq_from_bytes(b"TCGT");
        let diff = SequenceDiff::compute_ids(&obs_handle, &exp_handle, TerminalFilter::None);

        assert_eq!(diff.to_string(), "1T>A");
    }

    #[test]
    fn sequence_diff_compute_ids_identical() {
        use crate::interning::seq_from_bytes;

        let obs_handle = seq_from_bytes(b"ACGTACGT");
        let exp_handle = seq_from_bytes(b"ACGTACGT");
        let diff = SequenceDiff::compute_ids(&obs_handle, &exp_handle, TerminalFilter::None);

        assert!(diff.operations.is_empty());
        assert_eq!(diff.to_string(), "");
    }

    #[test]
    fn edit_operation_hgvs_sub_with_non_dna_bases() {
        // Should handle non-standard bases gracefully
        let op = EditOperation::Sub(10, b'X', b'Y');
        assert_eq!(op.to_hgvs_string(), "11X>Y");
    }

    #[test]
    fn edit_operation_hgvs_del_single_vs_multiple() {
        let op_single = EditOperation::Del(2, vec![b'C']);
        let op_multi = EditOperation::Del(2, vec![b'C', b'A', b'T']);

        assert_eq!(op_single.to_hgvs_string(), "2_2_delC");
        assert_eq!(op_multi.to_hgvs_string(), "2_4_delCAT");
    }

    #[test]
    fn edit_operation_hgvs_ins_position_zero() {
        let op = EditOperation::Ins(0, vec![b'A', b'T']);
        assert_eq!(op.to_hgvs_string(), "0_1_insAT");
    }

    #[test]
    fn sequence_diff_large_identical() {
        let large = vec![b'A'; 500];
        let diff = SequenceDiff::compute(&large, &large, TerminalFilter::None);
        assert!(diff.operations.is_empty());
    }

    #[test]
    fn sequence_diff_large_single_sub_at_start() {
        let mut large_obs = vec![b'A'; 500];
        let large_exp = vec![b'A'; 500];
        large_obs[0] = b'T';

        let diff = SequenceDiff::compute(&large_obs, &large_exp, TerminalFilter::None);
        assert_eq!(diff.operations.len(), 1);
        assert_eq!(diff.to_string(), "1A>T");
    }

    #[test]
    fn sequence_diff_large_single_sub_at_end() {
        let mut large_obs = vec![b'A'; 500];
        let large_exp = vec![b'A'; 500];
        large_obs[499] = b'T';

        let diff = SequenceDiff::compute(&large_obs, &large_exp, TerminalFilter::None);
        assert_eq!(diff.operations.len(), 1);
        assert_eq!(diff.to_string(), "500A>T");
    }

    #[test]
    fn sequence_diff_large_with_multiple_edits() {
        let mut large_obs = vec![b'A'; 100];
        let large_exp = vec![b'A'; 100];
        large_obs[0] = b'T';
        large_obs[50] = b'C';
        large_obs[99] = b'G';

        let diff = SequenceDiff::compute(&large_obs, &large_exp, TerminalFilter::None);
        assert_eq!(diff.operations.len(), 3);
        assert!(diff.to_string().contains("1A>T"));
        assert!(diff.to_string().contains("51A>C"));
        assert!(diff.to_string().contains("100A>G"));
    }

    #[test]
    fn sequence_diff_many_consecutive_substitutions() {
        let obs = b"TTTTTTTTTT";
        let exp = b"AAAAAAAAAA";
        let diff = SequenceDiff::compute(obs, exp, TerminalFilter::None);

        // Should have 10 separate substitutions
        assert_eq!(diff.operations.len(), 10);
        for (i, op) in diff.operations.iter().enumerate() {
            assert!(matches!(op, EditOperation::Sub(pos, b'A', b'T') if *pos == i));
        }
    }

    #[test]
    fn sequence_diff_alternating_pattern() {
        let obs = b"ACACAC";
        let exp = b"AGAGAG";
        let diff = SequenceDiff::compute(obs, exp, TerminalFilter::None);

        assert_eq!(diff.operations.len(), 3);
        assert_eq!(diff.to_string(), "2G>C;4G>C;6G>C");
    }

    #[test]
    fn sequence_diff_with_ambiguous_bases() {
        let obs = b"ACNGT";
        let exp = b"ACXGT";
        let diff = SequenceDiff::compute(obs, exp, TerminalFilter::None);

        assert_eq!(diff.operations.len(), 1);
        assert_eq!(diff.to_string(), "3X>N");
    }

    #[test]
    fn sequence_diff_insertion_multiple_same_base() {
        let obs = b"AAAAA";
        let exp = b"A";
        let diff = SequenceDiff::compute(obs, exp, TerminalFilter::None);

        assert_eq!(diff.operations.len(), 1);
        assert_eq!(diff.to_string(), "0_1_insAAAA");
    }

    #[test]
    fn sequence_diff_deletion_multiple_same_base() {
        let obs = b"A";
        let exp = b"AAAAA";
        let diff = SequenceDiff::compute(obs, exp, TerminalFilter::None);

        assert_eq!(diff.operations.len(), 1);
        assert_eq!(diff.to_string(), "0_3_delAAAA");
    }

    #[test]
    fn sequence_diff_display_empty() {
        let diff = SequenceDiff::new(vec![]);
        assert_eq!(format!("{}", diff), "");
    }

    #[test]
    fn sequence_diff_display_single_op() {
        let ops = vec![EditOperation::Sub(5, b'G', b'C')];
        let diff = SequenceDiff::new(ops);
        assert_eq!(format!("{}", diff), "6G>C");
    }

    #[test]
    fn sequence_diff_display_many_ops() {
        let ops = (0..100)
            .map(|i| EditOperation::Sub(i, b'A', b'T'))
            .collect();
        let diff = SequenceDiff::new(ops);
        let output = format!("{}", diff);

        assert!(output.contains("1A>T"));
        assert!(output.contains("100A>T"));
        assert_eq!(output.matches(';').count(), 99);
    }

    #[test]
    fn sequence_diff_display_semicolon_separation() {
        let ops = vec![
            EditOperation::Sub(0, b'A', b'T'),
            EditOperation::Ins(2, vec![b'G']),
            EditOperation::Del(5, vec![b'C']),
        ];
        let diff = SequenceDiff::new(ops);
        let output = format!("{}", diff);

        let parts: Vec<&str> = output.split(';').collect();
        assert_eq!(parts.len(), 3);
        assert_eq!(parts[0], "1A>T");
        assert_eq!(parts[1], "2_3_insG");
        assert_eq!(parts[2], "5_5_delC");
    }

    #[test]
    fn terminal_filter_none_keeps_leading_deletion() {
        let diff = SequenceDiff::compute(b"CGT", b"ACGT", TerminalFilter::None);

        assert_eq!(diff.to_string(), "0_0_delA");
    }

    #[test]
    fn terminal_filter_leading_removes_leading_deletion() {
        let diff = SequenceDiff::compute(b"CGT", b"ACGT", TerminalFilter::Leading);

        assert_eq!(diff.to_string(), "");
    }

    #[test]
    fn terminal_filter_trailing_removes_trailing_deletion() {
        let diff = SequenceDiff::compute(b"ACG", b"ACGT", TerminalFilter::Trailing);

        assert_eq!(diff.to_string(), "");
    }

    #[test]
    fn terminal_filter_leading_does_not_remove_trailing_deletion() {
        let diff = SequenceDiff::compute(b"ACG", b"ACGT", TerminalFilter::Leading);

        assert_eq!(diff.to_string(), "3_3_delT");
    }

    #[test]
    fn terminal_filter_trailing_does_not_remove_leading_deletion() {
        let diff = SequenceDiff::compute(b"CGT", b"ACGT", TerminalFilter::Trailing);

        assert_eq!(diff.to_string(), "0_0_delA");
    }

    #[test]
    fn terminal_filter_both_removes_leading_and_trailing_deletions() {
        let diff = SequenceDiff::compute(b"CGCG", b"ACGCGT", TerminalFilter::Both);

        assert_eq!(diff.to_string(), "");
    }

    #[test]
    fn terminal_filter_preserves_internal_substitution() {
        let diff = SequenceDiff::compute(b"GCGCTT", b"AGCGCGT", TerminalFilter::Both);

        assert_eq!(diff.to_string(), "6G>T");
    }

    #[test]
    fn terminal_filter_preserves_internal_deletion() {
        let diff = SequenceDiff::compute(b"AGCGT", b"ACGCGT", TerminalFilter::Both);

        assert_eq!(diff.to_string(), "1_1_delC");
    }

    #[test]
    fn terminal_filter_preserves_internal_insertion() {
        let diff = SequenceDiff::compute(b"ACGCGT", b"ACCGT", TerminalFilter::Both);

        assert_eq!(diff.to_string(), "2_3_insG");
    }
}
