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
        assert_eq!(seq_to_bytes(key.forward).as_ref(), b"ACGTACGT");
        assert!(key.reverse.is_none());
    }

    /// Paired-end: reverse read must be part of the key and preserved byte-for-byte.
    #[test]
    fn paired_end_key_includes_reverse() {
        let rp = make_readpair(b"AAAA", Some(b"TTTT"), ReadGroup::grouped("g1".into()));
        let key = rp.key();

        assert_eq!(seq_to_bytes(key.forward).as_ref(), b"AAAA");
        assert_eq!(seq_to_bytes(key.reverse.unwrap()).as_ref(), b"TTTT");

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
        assert_eq!(seq_to_bytes(k_single.forward).as_ref(), b"GGGG");
        assert!(k_single.reverse.is_none());
        assert_eq!(seq_to_bytes(k_paired.forward).as_ref(), b"GGGG");
        assert_eq!(seq_to_bytes(k_paired.reverse.unwrap()).as_ref(), b"A");
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

        assert_eq!(seq_to_bytes(k1.forward).as_ref(), b"A");
        assert!(k1.reverse.is_none());

        assert_eq!(seq_to_bytes(k2.forward).as_ref(), b"A");
        assert_eq!(seq_to_bytes(k2.reverse.unwrap()).as_ref(), b"");

        assert_ne!(
            k1, k2,
            "empty reverse is still a distinct key from no reverse"
        );
    }
}
