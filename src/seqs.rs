//! Core sequence container types used throughout DNAComb.
//!
//! This module provides lightweight containers for:
//! - full read sequences (`SeqPair`), used as stable keys and cached values,
//! - parsed sequencing records plus grouping metadata (`ReadPair`), used during
//!   parsing and counting.
use crate::interning::SeqHandle;
use crate::{groups::ReadGroup, interning::seq_from_bytes};
use bio::{bio_types::sequence::Sequence, io::fastq};

/// Forward/reverse sequence pair used as a stable content-based key.
///
/// `SeqPair` stores only sequence content, not qualities, IDs, or grouping
/// metadata. It is therefore suitable for hashing, deduplication, caching,
/// and output.
#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub struct SeqPair {
    pub forward: SeqHandle,
    pub reverse: Option<SeqHandle>,
}

impl SeqPair {
    /// Construct a `SeqPair` from owned forward and optional reverse sequences.
    pub fn new(forward: Sequence, reverse: Option<Sequence>) -> Self {
        Self {
            forward: seq_from_bytes(&forward),
            reverse: reverse.map(|x| seq_from_bytes(&x)),
        }
    }

    /// Build a sequence-only key from a parsed `ReadPair`.
    ///
    /// This drops read names, qualities, and grouping metadata while preserving
    /// forward/reverse sequence content.
    pub fn from_readpair(rp: &ReadPair) -> Self {
        Self::new(
            rp.forward.seq().to_vec(),
            rp.reverse.as_ref().map(|x| x.seq().to_vec()),
        )
    }
}

/// Parsed sequencing read pair plus grouping metadata.
///
/// This is the main record type yielded by parsers during counting. It stores
/// FASTQ-style forward and optional reverse reads together with the assigned
/// `ReadGroup`.
#[derive(Debug)]
pub struct ReadPair {
    pub forward: fastq::Record,
    pub reverse: Option<fastq::Record>,
    pub group: ReadGroup,
}

impl ReadPair {
    /// Generate the corresponding sequence-only `SeqPair` key.
    ///
    /// This is useful for caching and deduplication where only sequence content
    /// matters.
    pub fn key(&self) -> SeqPair {
        SeqPair::from_readpair(self)
    }

    /// Consume this `ReadPair` and convert it into a sequence-only `SeqPair`.
    pub fn into_seqpair(self) -> SeqPair {
        SeqPair::new(
            self.forward.seq().to_vec(),
            self.reverse.map(|x| x.seq().to_vec()),
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::interning::seq_to_bytes;

    use super::*;
    use bio::io::fastq;
    use std::collections::HashSet;

    fn make_record(id: &str, seq: &[u8]) -> fastq::Record {
        let qual = vec![b'I'; seq.len()]; // high, valid Phred+33
        fastq::Record::with_attrs(id, None, seq, &qual)
    }

    fn make_readpair(f_seq: &[u8], r_seq: Option<&[u8]>, group: ReadGroup) -> ReadPair {
        ReadPair {
            forward: make_record("F", f_seq),
            reverse: r_seq.map(|s| make_record("R", s)),
            group,
        }
    }

    /// Single-end: ReadPair::key() should match SeqPair::new with the same bytes.
    #[test]
    fn single_end_key_matches_constructor() {
        let rp = make_readpair(b"ACGTACGT", None, ReadGroup::ungrouped());
        let key = rp.key();

        let constructed = SeqPair::new(b"ACGTACGT".to_vec(), None);

        assert_eq!(
            key, constructed,
            "single-end key must equal constructor-derived SeqPair"
        );
        assert_eq!(seq_to_bytes(&key.forward).as_ref(), b"ACGTACGT");
        assert!(key.reverse.is_none());
    }

    /// Paired-end: reverse read must be part of the key and preserved byte-for-byte.
    #[test]
    fn paired_end_key_includes_reverse() {
        let rp = make_readpair(b"AAAA", Some(b"TTTT"), ReadGroup::grouped("g1".into()));
        let key = rp.key();

        assert_eq!(seq_to_bytes(&key.forward).as_ref(), b"AAAA");
        assert_eq!(
            seq_to_bytes(key.reverse.as_ref().unwrap()).as_ref(),
            b"TTTT"
        );

        // Cross-check with constructor:
        let constructed = SeqPair::new(b"AAAA".to_vec(), Some(b"TTTT".to_vec()));
        assert_eq!(
            key, constructed,
            "paired-end key should equal constructed SeqPair"
        );
    }

    /// Hash/Eq: identical SeqPairs must collapse in a HashSet; different ones must not.
    #[test]
    fn hash_equality_semantics() {
        let mut set: HashSet<SeqPair> = HashSet::new();
        let a = SeqPair::new(b"AAAA".to_vec(), Some(b"TTTT".to_vec()));
        let b = SeqPair::new(b"AAAA".to_vec(), Some(b"TTTT".to_vec())); // identical
        let c = SeqPair::new(b"AAAA".to_vec(), Some(b"TTTA".to_vec())); // reverse differs
        let d = SeqPair::new(b"AAAA".to_vec(), None); // reverse missing

        assert!(set.insert(a.clone()));
        assert!(!set.insert(b), "identical key should not insert twice");
        assert!(
            set.insert(c),
            "different reverse must produce a distinct key"
        );
        assert!(
            set.insert(d),
            "absence of reverse must produce a distinct key"
        );
        assert_eq!(set.len(), 3);
    }

    /// Reverse presence alone must change the key (single-end vs paired-end with same forward).
    #[test]
    fn reverse_presence_changes_key() {
        let rp_single = make_readpair(b"GGGG", None, ReadGroup::ungrouped());
        let rp_paired = make_readpair(b"GGGG", Some(b"A"), ReadGroup::ungrouped());

        let k_single = rp_single.key();
        let k_paired = rp_paired.key();

        assert_ne!(
            k_single, k_paired,
            "adding a reverse read must change the key"
        );
        assert_eq!(seq_to_bytes(&k_single.forward).as_ref(), b"GGGG");
        assert!(k_single.reverse.is_none());
        assert_eq!(seq_to_bytes(&k_paired.forward).as_ref(), b"GGGG");
        assert_eq!(
            seq_to_bytes(k_paired.reverse.as_ref().unwrap()).as_ref(),
            b"A"
        );
    }

    /// Determinism: repeated calls to ReadPair::key() must be stable.
    #[test]
    fn readpair_key_is_deterministic() {
        let rp = make_readpair(b"TACT", Some(b"AGGA"), ReadGroup::grouped("g2".into()));
        let k1 = rp.key();
        let k2 = rp.key();
        assert_eq!(k1, k2, "key generation must be deterministic");
    }

    /// Smoke test: very short reads (including empty reverse) should still produce valid keys.
    #[test]
    fn tiny_reads_smoke() {
        let rp1 = make_readpair(b"A", None, ReadGroup::ungrouped());
        let rp2 = make_readpair(b"A", Some(b""), ReadGroup::ungrouped());

        let k1 = rp1.key();
        let k2 = rp2.key();

        assert_eq!(seq_to_bytes(&k1.forward).as_ref(), b"A");
        assert!(k1.reverse.is_none());

        assert_eq!(seq_to_bytes(&k2.forward).as_ref(), b"A");
        assert_eq!(seq_to_bytes(k2.reverse.as_ref().unwrap()).as_ref(), b"");

        assert_ne!(
            k1, k2,
            "empty reverse is still a distinct key from no reverse"
        );
    }

    #[test]
    fn seqpair_new_single_end() {
        let sp = SeqPair::new(b"ACGT".to_vec(), None);
        assert_eq!(seq_to_bytes(&sp.forward).as_ref(), b"ACGT");
        assert!(sp.reverse.is_none());
    }

    #[test]
    fn seqpair_new_paired_end() {
        let sp = SeqPair::new(b"AAAA".to_vec(), Some(b"TTTT".to_vec()));
        assert_eq!(seq_to_bytes(&sp.forward).as_ref(), b"AAAA");
        assert_eq!(seq_to_bytes(sp.reverse.as_ref().unwrap()).as_ref(), b"TTTT");
    }

    #[test]
    fn seqpair_new_empty_forward() {
        let sp = SeqPair::new(b"".to_vec(), None);
        assert_eq!(seq_to_bytes(&sp.forward).as_ref(), b"");
        assert!(sp.reverse.is_none());
    }

    #[test]
    fn seqpair_new_empty_reverse() {
        let sp = SeqPair::new(b"ACGT".to_vec(), Some(b"".to_vec()));
        assert_eq!(seq_to_bytes(&sp.forward).as_ref(), b"ACGT");
        assert_eq!(seq_to_bytes(sp.reverse.as_ref().unwrap()).as_ref(), b"");
    }

    #[test]
    fn seqpair_new_both_empty() {
        let sp = SeqPair::new(b"".to_vec(), Some(b"".to_vec()));
        assert_eq!(seq_to_bytes(&sp.forward).as_ref(), b"");
        assert_eq!(seq_to_bytes(sp.reverse.as_ref().unwrap()).as_ref(), b"");
    }

    #[test]
    fn seqpair_new_long_sequences() {
        let long_f = vec![b'A'; 10_000];
        let long_r = vec![b'T'; 10_000];
        let sp = SeqPair::new(long_f.clone(), Some(long_r.clone()));
        assert_eq!(seq_to_bytes(&sp.forward).as_ref(), &long_f[..]);
        assert_eq!(
            seq_to_bytes(sp.reverse.as_ref().unwrap()).as_ref(),
            &long_r[..]
        );
    }

    #[test]
    fn seqpair_from_readpair_single_end() {
        let rp = make_readpair(b"ACGTACGT", None, ReadGroup::ungrouped());
        let sp = SeqPair::from_readpair(&rp);
        assert_eq!(seq_to_bytes(&sp.forward).as_ref(), b"ACGTACGT");
        assert!(sp.reverse.is_none());
    }

    #[test]
    fn seqpair_from_readpair_paired_end() {
        let rp = make_readpair(b"AAAA", Some(b"TTTT"), ReadGroup::grouped("g1"));
        let sp = SeqPair::from_readpair(&rp);
        assert_eq!(seq_to_bytes(&sp.forward).as_ref(), b"AAAA");
        assert_eq!(seq_to_bytes(sp.reverse.as_ref().unwrap()).as_ref(), b"TTTT");
    }

    #[test]
    fn seqpair_from_readpair_ignores_group() {
        let rp1 = make_readpair(b"ACGT", Some(b"TGCA"), ReadGroup::ungrouped());
        let rp2 = make_readpair(
            b"ACGT",
            Some(b"TGCA"),
            ReadGroup::grouped("different_group"),
        );
        let sp1 = SeqPair::from_readpair(&rp1);
        let sp2 = SeqPair::from_readpair(&rp2);
        assert_eq!(sp1, sp2, "group should not affect SeqPair");
    }

    #[test]
    fn seqpair_from_readpair_ignores_qualities() {
        let mut rp1 = make_readpair(b"ACGT", Some(b"TGCA"), ReadGroup::ungrouped());
        let mut rp2 = make_readpair(b"ACGT", Some(b"TGCA"), ReadGroup::ungrouped());

        // Manually set different qualities
        rp1.forward = fastq::Record::with_attrs("f1", None, b"ACGT", b"IIII");
        rp2.forward = fastq::Record::with_attrs("f2", None, b"ACGT", b"!!!!");

        let sp1 = SeqPair::from_readpair(&rp1);
        let sp2 = SeqPair::from_readpair(&rp2);
        assert_eq!(sp1, sp2, "qualities should not affect SeqPair");
    }

    #[test]
    fn seqpair_from_readpair_ignores_read_names() {
        let mut rp1 = make_readpair(b"ACGT", Some(b"TGCA"), ReadGroup::ungrouped());
        let mut rp2 = make_readpair(b"ACGT", Some(b"TGCA"), ReadGroup::ungrouped());

        rp1.forward = fastq::Record::with_attrs("name1", None, b"ACGT", b"IIII");
        rp2.forward = fastq::Record::with_attrs("name2", None, b"ACGT", b"IIII");

        let sp1 = SeqPair::from_readpair(&rp1);
        let sp2 = SeqPair::from_readpair(&rp2);
        assert_eq!(sp1, sp2, "read names should not affect SeqPair");
    }

    #[test]
    fn readpair_key_single_end() {
        let rp = make_readpair(b"ACGTACGT", None, ReadGroup::ungrouped());
        let key = rp.key();
        assert_eq!(seq_to_bytes(&key.forward).as_ref(), b"ACGTACGT");
        assert!(key.reverse.is_none());
    }

    #[test]
    fn readpair_key_paired_end() {
        let rp = make_readpair(b"AAAA", Some(b"TTTT"), ReadGroup::grouped("g1"));
        let key = rp.key();
        assert_eq!(seq_to_bytes(&key.forward).as_ref(), b"AAAA");
        assert_eq!(
            seq_to_bytes(key.reverse.as_ref().unwrap()).as_ref(),
            b"TTTT"
        );
    }

    #[test]
    fn readpair_key_matches_seqpair_new() {
        let rp = make_readpair(b"ACGTACGT", Some(b"TGCATGCA"), ReadGroup::ungrouped());
        let key = rp.key();
        let constructed = SeqPair::new(b"ACGTACGT".to_vec(), Some(b"TGCATGCA".to_vec()));
        assert_eq!(key, constructed);
    }

    #[test]
    fn readpair_key_is_stable() {
        let rp = make_readpair(b"TACT", Some(b"AGGA"), ReadGroup::grouped("g2"));
        let k1 = rp.key();
        let k2 = rp.key();
        assert_eq!(k1, k2);
    }

    #[test]
    fn readpair_into_seqpair_single_end() {
        let rp = make_readpair(b"ACGTACGT", None, ReadGroup::ungrouped());
        let sp = rp.into_seqpair();
        assert_eq!(seq_to_bytes(&sp.forward).as_ref(), b"ACGTACGT");
        assert!(sp.reverse.is_none());
    }

    #[test]
    fn readpair_into_seqpair_paired_end() {
        let rp = make_readpair(b"AAAA", Some(b"TTTT"), ReadGroup::grouped("g1"));
        let sp = rp.into_seqpair();
        assert_eq!(seq_to_bytes(&sp.forward).as_ref(), b"AAAA");
        assert_eq!(seq_to_bytes(sp.reverse.as_ref().unwrap()).as_ref(), b"TTTT");
    }

    #[test]
    fn seqpair_equality_identical() {
        let sp1 = SeqPair::new(b"ACGT".to_vec(), Some(b"TGCA".to_vec()));
        let sp2 = SeqPair::new(b"ACGT".to_vec(), Some(b"TGCA".to_vec()));
        assert_eq!(sp1, sp2);
    }

    #[test]
    fn seqpair_inequality_forward_differs() {
        let sp1 = SeqPair::new(b"ACGT".to_vec(), Some(b"TGCA".to_vec()));
        let sp2 = SeqPair::new(b"AAAA".to_vec(), Some(b"TGCA".to_vec()));
        assert_ne!(sp1, sp2);
    }

    #[test]
    fn seqpair_inequality_reverse_differs() {
        let sp1 = SeqPair::new(b"ACGT".to_vec(), Some(b"TGCA".to_vec()));
        let sp2 = SeqPair::new(b"ACGT".to_vec(), Some(b"AAAA".to_vec()));
        assert_ne!(sp1, sp2);
    }

    #[test]
    fn seqpair_inequality_reverse_presence() {
        let sp1 = SeqPair::new(b"ACGT".to_vec(), None);
        let sp2 = SeqPair::new(b"ACGT".to_vec(), Some(b"TGCA".to_vec()));
        assert_ne!(sp1, sp2);
    }

    #[test]
    fn seqpair_in_hashset() {
        let mut set: HashSet<SeqPair> = HashSet::new();
        let sp1 = SeqPair::new(b"AAAA".to_vec(), Some(b"TTTT".to_vec()));
        let sp2 = SeqPair::new(b"AAAA".to_vec(), Some(b"TTTT".to_vec()));
        let sp3 = SeqPair::new(b"AAAA".to_vec(), Some(b"TTTA".to_vec()));

        assert!(set.insert(sp1.clone()));
        assert!(!set.insert(sp2), "duplicate should not insert");
        assert!(set.insert(sp3), "different sequence should insert");
        assert_eq!(set.len(), 2);
    }

    #[test]
    fn seqpair_clone() {
        let sp1 = SeqPair::new(b"ACGT".to_vec(), Some(b"TGCA".to_vec()));
        let sp2 = sp1.clone();
        assert_eq!(sp1, sp2);
    }

    #[test]
    fn empty_reverse_vs_none_reverse() {
        let sp_none = SeqPair::new(b"ACGT".to_vec(), None);
        let sp_empty = SeqPair::new(b"ACGT".to_vec(), Some(b"".to_vec()));
        assert_ne!(
            sp_none, sp_empty,
            "empty reverse is distinct from no reverse"
        );
    }

    #[test]
    fn very_short_reads() {
        let sp1 = SeqPair::new(b"A".to_vec(), None);
        let sp2 = SeqPair::new(b"A".to_vec(), Some(b"T".to_vec()));
        assert_ne!(sp1, sp2);
        assert_eq!(seq_to_bytes(&sp1.forward).as_ref(), b"A");
        assert_eq!(seq_to_bytes(&sp2.forward).as_ref(), b"A");
    }

    #[test]
    fn readpair_with_different_groups_same_sequences() {
        let rp1 = make_readpair(b"ACGT", Some(b"TGCA"), ReadGroup::ungrouped());
        let rp2 = make_readpair(b"ACGT", Some(b"TGCA"), ReadGroup::grouped("g1"));
        let rp3 = make_readpair(b"ACGT", Some(b"TGCA"), ReadGroup::grouped("g2"));

        let k1 = rp1.key();
        let k2 = rp2.key();
        let k3 = rp3.key();

        assert_eq!(k1, k2, "group should not affect key");
        assert_eq!(k2, k3, "group should not affect key");
    }

    #[test]
    fn readpair_into_seqpair_equivalence_with_key() {
        let rp = make_readpair(b"AAAA", Some(b"TTTT"), ReadGroup::grouped("g1"));
        let sp_key = rp.key();

        let rp2 = make_readpair(b"AAAA", Some(b"TTTT"), ReadGroup::grouped("g1"));
        let sp_into = rp2.into_seqpair();

        assert_eq!(
            sp_key, sp_into,
            "key() and into_seqpair() should produce same result"
        );
    }

    #[test]
    fn various_dna_bases() {
        let sequences: Vec<&[u8]> = vec![b"ACGT", b"NNNNN", b"AAAAAA", b"GCGCGC"];
        let mut set: HashSet<SeqPair> = HashSet::new();

        for seq in sequences {
            let sp = SeqPair::new(seq.to_vec(), None);
            set.insert(sp);
        }

        assert_eq!(set.len(), 4, "all different sequences should be distinct");
    }

    #[test]
    fn asymmetric_forward_reverse_lengths() {
        let sp1 = SeqPair::new(b"A".to_vec(), Some(b"TTTTTTTTTT".to_vec()));
        let sp2 = SeqPair::new(b"AAAAAAAAAA".to_vec(), Some(b"T".to_vec()));

        assert_ne!(sp1, sp2);
        assert_eq!(seq_to_bytes(&sp1.forward).len(), 1);
        assert_eq!(seq_to_bytes(sp1.reverse.as_ref().unwrap()).len(), 10);
        assert_eq!(seq_to_bytes(&sp2.forward).len(), 10);
        assert_eq!(seq_to_bytes(sp2.reverse.as_ref().unwrap()).len(), 1);
    }
}
