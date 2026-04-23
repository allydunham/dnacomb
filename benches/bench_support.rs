//! Support functions for benchmarking
#![allow(dead_code)]
use dnacomb::{Compression, FilterConfig, LibrarySpec, ReadPairParser, SeqFormat, SeqPath};
use std::fs::File;
use std::io::Write;

use dnacomb::interning::seq_from_bytes;

/// Generate an interned SeqPair
pub fn seq_handle_pair(
    a: &[u8],
    b: &[u8],
) -> (dnacomb::interning::SeqHandle, dnacomb::interning::SeqHandle) {
    (seq_from_bytes(a), seq_from_bytes(b))
}

/// Write a temporary FASTQ file with repeated sequences
pub fn make_fastq(path: &str, seq: &[u8], n: usize) {
    let mut f = File::create(path).unwrap();
    for i in 0..n {
        writeln!(f, "@read{}", i).unwrap();
        writeln!(f, "{}", std::str::from_utf8(seq).unwrap()).unwrap();
        writeln!(f, "+").unwrap();
        writeln!(f, "{}", "I".repeat(seq.len())).unwrap();
    }
}

/// Simple LibSpec generator (fixed + one variable region)
pub fn simple_libspec() -> LibrarySpec {
    let json = r#"
    {
        "id": "test",
        "forward_start_region": "fixed",
        "forward_read_length": 50,
        "reverse_start_region": "fixed",
        "reverse_read_length": 50,
        "regions": [
            { "id": "fixed", "seq_type": "Fixed", "seq": "AAAA" },
            { "id": "var", "seq_type": "Library", "min_length": 4, "max_length": 4 }
        ]
    }
    "#;

    json.parse().unwrap()
}

/// Create parser from temp FASTQ
pub fn make_parser(path: &str, n_reads: u64) -> ReadPairParser {
    let forward = SeqPath::new(path.to_string(), SeqFormat::Fastq, Compression::None);
    ReadPairParser::from_paths(forward, None, None, n_reads, b'I').unwrap()
}

/// Default filter config (no filtering)
pub fn no_filter() -> FilterConfig {
    FilterConfig::new(None, None, None, None, true)
}

/// Deterministic pseudo-DNA sequence generator.
pub fn make_base_seq(i: usize, len: usize) -> Vec<u8> {
    let alphabet = [b'A', b'C', b'G', b'T'];
    let mut x = i ^ 0x9E37_79B9usize;
    let mut out = Vec::with_capacity(len);
    for _ in 0..len {
        out.push(alphabet[x & 0b11]);
        x = x.rotate_left(5) ^ 0x85EB_CA6Busize;
    }
    out
}

/// Add deterministic mutations to a sequence
pub fn mutate_sub(mut seq: Vec<u8>, pos: usize) -> Vec<u8> {
    seq[pos] = match seq[pos] {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    };
    seq
}

/// Apply a set of subs to a sequence
pub fn apply_n_subs(seq: &[u8], n_subs: usize) -> Vec<u8> {
    let mut out = seq.to_vec();
    if out.is_empty() {
        return out;
    }

    for k in 0..n_subs {
        let pos = (k * 7 + out.len() / 3) % out.len();
        out = mutate_sub(out, pos);
    }

    out
}

/// Insert a base into a sequence
pub fn insert_base(seq: &[u8], pos: usize, base: u8) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len() + 1);
    out.extend_from_slice(&seq[..pos]);
    out.push(base);
    out.extend_from_slice(&seq[pos..]);
    out
}

/// Remove a base from a sequence
pub fn delete_base(seq: &[u8], pos: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len().saturating_sub(1));
    out.extend_from_slice(&seq[..pos]);
    out.extend_from_slice(&seq[pos + 1..]);
    out
}

/// Make varied mutations through a sequence
pub fn apply_clustered_subs(seq: &[u8], start: usize, n_subs: usize) -> Vec<u8> {
    let mut out = seq.to_vec();
    if out.is_empty() {
        return out;
    }

    for k in 0..n_subs {
        let pos = (start + k).min(out.len() - 1);
        out = mutate_sub(out, pos);
    }

    out
}

/// Apply some mixed edits to a sequence
pub fn apply_mixed_edit(seq: &[u8]) -> Vec<u8> {
    if seq.len() < 4 {
        return insert_base(seq, seq.len() / 2, b'A');
    }

    let out = apply_n_subs(seq, 1);
    let out = insert_base(&out, out.len() / 2, b'T');
    delete_base(&out, out.len() / 3)
}
