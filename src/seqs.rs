//! Core data structures for DNA sequence objects
//!
//! Provides various core data objects needed across the library for DNA sequences.
use bio::{bio_types::sequence::Sequence, io::fastq};
use std::fmt;

pub type ReadKey = (Vec<u8>, Option<Vec<u8>>);

/// Pair of sequences
#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub struct SeqPair {
    pub forward: Sequence,
    pub reverse: Option<Sequence>,
}

impl SeqPair {
    pub fn new(forward: Sequence, reverse: Option<Sequence>) -> Self {
        Self { forward, reverse }
    }

    pub fn from_readpair(rp: &ReadPair) -> Self {
        Self::new(
            rp.forward.seq().to_vec(),
            rp.reverse.as_ref().map(|x| x.seq().to_vec()),
        )
    }
}

/// Pair of linked Fastq reads
#[derive(Debug)]
pub struct ReadPair {
    pub forward: fastq::Record,
    pub reverse: Option<fastq::Record>,
    pub group: ReadGroup,
}

impl ReadPair {
    /// Generate a key to identify unique read types
    pub fn key(&self) -> ReadKey {
        if self.reverse.is_some() {
            (
                self.forward.seq().to_vec(),
                Some(self.reverse.as_ref().unwrap().seq().to_vec()),
            )
        } else {
            (self.forward.seq().to_vec(), None)
        }
    }

    pub fn into_seqpair(self) -> SeqPair {
        SeqPair::new(
            self.forward.seq().to_vec(),
            self.reverse.map(|x| x.seq().to_vec()),
        )
    }
}

/// Group status of a read
#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub enum ReadGroup {
    Ungrouped,
    Unmatched,
    Match(String),
}

impl fmt::Display for ReadGroup {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ReadGroup::Ungrouped => write!(f, ""),
            ReadGroup::Unmatched => write!(f, "_unmatched_"),
            ReadGroup::Match(x) => write!(f, "{}", x),
        }
    }
}

#[cfg(test)]
mod tests {
    // use super::*;
}
