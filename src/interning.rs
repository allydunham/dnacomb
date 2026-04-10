//! Facilities to intern and deduplicate frequently used objects
//!
//! Interning sequences (Vec<u8>) and strings (e.g. group and region names) can
//! drastically decrease memory usage as they are included in many objects.
//! It can therefore increase speed as copying load is a significant overhead and
//! moves initial comparisons to an O(1) opperation.
//!
//! Interning has no garbage collection because the standard workflow on the CLI accumalates
//! then outputs and everything is always in use. For use as a more general library where you
//! might process then drop sequences/regions/groups a feature flag is included to disable
//! interning in favour of a pass-through owned sequence/string.
//!
//! Public API (same in both modes):
//!   - SeqHandle,   seq_from_bytes(),   seq_bytes()
//!   - GroupHandle, group_from_str(),   group_str()
//!   - RegionHandle, region_from_str(), region_str()
//!
//! Since Groups, IDs and Regions are both strings they use ThreadedRodeo but are separated to allow
//! different tuning since they have different access patterns

use std::sync::Arc;

#[cfg(feature = "interning")]
mod enabled {
    use super::Arc;
    use std::num::NonZeroU32;
    use std::sync::atomic::{AtomicU32, Ordering};

    use ahash::RandomState;
    use dashmap::DashMap;
    use once_cell::sync::Lazy;
    use parking_lot::RwLock;

    // Sequences (Vec<u8>)
    /// Shared unique SeqInterner instance storing interned sequences
    static SEQ_INTERNER: Lazy<SeqInterner> = Lazy::new(SeqInterner::default);

    /// Get a SeqHandle for a Sequence vector, interning it if it's new
    #[inline]
    pub fn seq_from_bytes(bytes: &[u8]) -> SeqHandle {
        SEQ_INTERNER.intern(bytes)
    }

    /// Resolve seq handle to owned Sequence (Arc<[u8]>)
    #[inline]
    pub fn seq_to_bytes(h: SeqHandle) -> Arc<[u8]> {
        SEQ_INTERNER.resolve(h)
    }

    // Reserve more space in the Sequence interner
    pub fn reserve_seq_interner(estimated_unique: usize) {
        SEQ_INTERNER.reserve(estimated_unique);
    }

    /// Handle to access a Sequence object (Vec<u8>)
    ///
    /// This access point interns the underlying vector to save duplication. It is
    /// implemented as a NonZeroU32 ID under the hood to allow Option<SeqHandle> to
    /// fit in 4 bytes.
    #[repr(transparent)]
    #[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
    pub struct SeqHandle(NonZeroU32);

    impl SeqHandle {
        #[inline]
        pub fn get(self) -> u32 {
            self.0.get()
        }

        #[inline]
        fn from_index0(i: u32) -> Self {
            // Stored as 1-based so Option<SeqHandle> is niche-optimized
            SeqHandle(NonZeroU32::new(i + 1).expect("nonzero"))
        }

        #[inline]
        fn index0(self) -> usize {
            (self.0.get() - 1) as usize
        }
    }

    pub struct SeqInterner {
        /// Forward map: canonical bytes -> id
        forward: DashMap<Arc<[u8]>, SeqHandle, RandomState>,

        /// Reverse vec: id -> canonical bytes (index = id.get() - 1)
        reverse: RwLock<Vec<Arc<[u8]>>>,

        /// Next 0-based index to allocate
        next: AtomicU32,
    }

    impl Default for SeqInterner {
        fn default() -> Self {
            Self {
                forward: DashMap::with_hasher(RandomState::new()),
                reverse: RwLock::new(Vec::new()),
                next: AtomicU32::new(0),
            }
        }
    }

    impl SeqInterner {
        #[inline]
        pub fn intern(&self, bytes: &[u8]) -> SeqHandle {
            // Fast path
            if let Some(found) = self.forward.get(bytes) {
                return *found;
            }

            // Allocate once on miss
            let arc = Arc::<[u8]>::from(bytes);

            // Second check in case of race
            if let Some(found) = self.forward.get(&*arc) {
                return *found;
            }

            // Allocate an id
            let idx0 = self.next.fetch_add(1, Ordering::Relaxed);
            let id = SeqHandle::from_index0(idx0);

            // Publish reverse entry
            {
                let mut rev = self.reverse.write();
                let pos = idx0 as usize;

                match rev.len().cmp(&pos) {
                    std::cmp::Ordering::Less => {
                        // Extremely rare; keep safe
                        rev.resize_with(pos, || Arc::<[u8]>::from(&b""[..]));
                        rev.push(arc.clone())
                    }
                    std::cmp::Ordering::Equal => rev.push(arc.clone()),
                    std::cmp::Ordering::Greater => rev[pos] = arc.clone(),
                }
            }

            // Insert into forward map; if someone raced and inserted first, use theirs
            if let Some(prev_id) = self.forward.insert(arc, id) {
                return prev_id;
            }

            id
        }

        #[inline]
        pub fn resolve(&self, id: SeqHandle) -> Arc<[u8]> {
            self.reverse.read()[id.index0()].clone()
        }

        pub fn reserve(&self, additional_unique: usize) {
            self.reverse.write().reserve(additional_unique);
        }
    }

    // Group IDs (Str)
    /// Shared unique interner instance storing Group names
    static GROUPS: Lazy<lasso::ThreadedRodeo> = Lazy::new(lasso::ThreadedRodeo::default);

    /// Global interned group ID
    #[repr(transparent)]
    #[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
    pub struct Group(lasso::Spur);

    /// Get a GroupHandle for a Group ID, interning it if it's new
    #[inline]
    pub fn group_from_str(s: &str) -> Group {
        Group(GROUPS.get_or_intern(s))
    }

    /// Resolve group handle to owned string (Arc<str>)
    #[inline]
    pub fn group_to_str(id: Group) -> Arc<str> {
        Arc::<str>::from(GROUPS.resolve(&id.0))
    }

    // Group IDs (Str)
    /// Shared unique interner instance storing Group names
    static LIB_IDS: Lazy<lasso::ThreadedRodeo> = Lazy::new(lasso::ThreadedRodeo::default);

    /// Global interned Library ID
    #[repr(transparent)]
    #[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
    pub struct LibraryID(lasso::Spur);

    /// Get a LibraryID for a Group ID, interning it if it's new
    #[inline]
    pub fn library_id_from_str(s: &str) -> LibraryID {
        LibraryID(LIB_IDS.get_or_intern(s))
    }

    /// Resolve Library ID to owned string (Arc<str>)
    #[inline]
    pub fn library_id_to_str(id: LibraryID) -> Arc<str> {
        Arc::<str>::from(LIB_IDS.resolve(&id.0))
    }

    // Region IDs (Str)
    /// Shared unique interner instance storing Region IDs
    static REGIONS: Lazy<lasso::ThreadedRodeo> = Lazy::new(lasso::ThreadedRodeo::default);

    /// Global interned region ID
    #[repr(transparent)]
    #[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
    pub struct RegionID(lasso::Spur);

    /// Get a RegionHandle for a region ID, interning it if it's new
    #[inline]
    pub fn region_id_from_str(s: &str) -> RegionID {
        RegionID(REGIONS.get_or_intern(s))
    }

    /// Resolve region handle to owned string (Arc<str>).
    #[inline]
    pub fn region_id_to_str(id: RegionID) -> Arc<str> {
        Arc::<str>::from(REGIONS.resolve(&id.0))
    }
}

#[cfg(not(feature = "interning"))]
mod disabled {
    use super::Arc;

    // Sequences (Vec<u8>)
    /// Non-interning handle to access a Sequence object (Vec<u8>)
    #[derive(Clone, Eq, PartialEq, Hash, Debug)]
    pub struct SeqHandle(pub Arc<[u8]>);

    /// Get a SeqHandle for a Sequence vector (non-interning)
    #[inline]
    pub fn seq_from_bytes(bytes: &[u8]) -> SeqHandle {
        SeqHandle(Arc::<[u8]>::from(bytes))
    }

    /// Convert a non-interning SeqHandle to Sequence bytes
    #[inline]
    pub fn seq_to_bytes(h: SeqHandle) -> Arc<[u8]> {
        h.0
    }

    /// Non-interning reservation NoOp as no underlying interner)
    pub fn reserve_seq_interner(_estimated_unique: usize) {}

    // Group IDs (Str)
    /// Non-interning handle to access Group IDs
    #[derive(Clone, Eq, PartialEq, Hash, Debug)]
    pub struct Group(pub Arc<str>);

    /// Get a GroupHandle for a Group ID (non-interning)
    #[inline]
    pub fn group_from_str(s: &str) -> Group {
        Group(Arc::<str>::from(s))
    }

    /// Resolve non-interning group handle to owned string
    #[inline]
    pub fn group_to_str(h: Group) -> Arc<str> {
        h.0
    }

    // Library IDs (Str)
    /// Non-interning handle to access Group IDs
    #[derive(Clone, Eq, PartialEq, Hash, Debug)]
    pub struct LibraryID(pub Arc<str>);

    /// Get a GroupHandle for a Group ID (non-interning)
    #[inline]
    pub fn library_id_from_str(s: &str) -> LibraryID {
        LibraryID(Arc::<str>::from(s))
    }

    /// Resolve non-interning group handle to owned string
    #[inline]
    pub fn library_id_to_str(h: LibraryID) -> Arc<str> {
        h.0
    }

    // Region IDs (Str)
    /// Non-interning handle to access Region IDs
    #[derive(Clone, Eq, PartialEq, Hash, Debug)]
    pub struct RegionID(pub Arc<str>);

    /// Get a RegionHandle for a region ID (non-interning)
    #[inline]
    pub fn region_id_from_str(s: &str) -> RegionID {
        RegionID(Arc::<str>::from(s))
    }

    /// Resolve non-interning region handles to owned string (Arc<str>)
    #[inline]
    pub fn region_id_to_str(h: RegionID) -> Arc<str> {
        h.0
    }
}

// Re-export the chosen backend.
#[cfg(feature = "interning")]
pub use enabled::*;

#[cfg(not(feature = "interning"))]
pub use disabled::*;

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    // Common tests for interning and non-interning backends

    /// Test Sequence round trip
    #[test]
    fn seq_round_trip_bytes() {
        let h = seq_from_bytes(b"ACGTN");
        let out = seq_to_bytes(h);
        assert_eq!(&*out, b"ACGTN");
    }

    /// Test Group round trip
    #[test]
    fn group_round_trip_str() {
        let g = group_from_str("groupA");
        let s = group_to_str(g);
        assert_eq!(&*s, "groupA");
    }

    /// Test Library round trip
    #[test]
    fn library_round_trip_str() {
        let l = library_id_from_str("seqA");
        let s = library_id_to_str(l);
        assert_eq!(&*s, "seqA");
    }

    /// Test Region round trip
    #[test]
    fn region_round_trip_str() {
        let r = region_id_from_str("Reg1");
        let s = region_id_to_str(r);
        assert_eq!(&*s, "Reg1");
    }

    /// Handles need to work as Keys
    #[test]
    fn handles_work_as_hashmap_keys() {
        let a = seq_from_bytes(b"AAAA");
        let b = seq_from_bytes(b"GGGG");
        let mut m: HashMap<SeqHandle, u32> = HashMap::new();
        m.insert(a.clone(), 1);
        m.insert(b.clone(), 2);
        assert_eq!(m.get(&a), Some(&1));
        assert_eq!(m.get(&b), Some(&2));
    }

    // Test Interning module
    #[cfg(feature = "interning")]
    mod enabled_tests {
        use super::*;
        use std::collections::HashSet;
        use std::thread;

        /// Test same Sequence gets the same ID
        #[test]
        fn seq_dedup_same_bytes_get_same_handle() {
            let a1 = seq_from_bytes(b"ACGT");
            let a2 = seq_from_bytes(b"ACGT");
            assert_eq!(
                a1, a2,
                "interning enabled: identical bytes should dedupe to same handle"
            );

            // Optional stronger check: resolve pointers are identical (same Arc allocation)
            let p1 = Arc::as_ptr(&seq_to_bytes(a1));
            let p2 = Arc::as_ptr(&seq_to_bytes(a2));
            assert_eq!(
                p1, p2,
                "interning enabled: resolved Arc should be same allocation"
            );
        }

        /// Test Groups get the same handle
        #[test]
        fn group_dedup_same_str_get_same_handle() {
            let g1 = group_from_str("x");
            let g2 = group_from_str("x");
            assert_eq!(
                g1, g2,
                "interning enabled: identical group string should dedupe"
            );
        }

        /// Test regions get the same handle
        #[test]
        fn region_dedup_same_str_get_same_handle() {
            let r1 = region_id_from_str("r");
            let r2 = region_id_from_str("r");
            assert_eq!(
                r1, r2,
                "interning enabled: identical region string should dedupe"
            );
        }

        /// Test SeqHandle niche optimisation is working
        #[test]
        fn option_seqhandle_is_niche_optimized() {
            // If SeqHandle is NonZero-backed, Option<SeqHandle> should be same size as SeqHandle.
            assert_eq!(
                std::mem::size_of::<Option<SeqHandle>>(),
                std::mem::size_of::<SeqHandle>(),
                "enabled mode: SeqHandle should be niche-optimized so Option<SeqHandle> doesn't grow"
            );
        }

        /// Test concurrency is consistent
        #[test]
        fn concurrent_seq_interning_is_consistent() {
            // Intern the same set of sequences from multiple threads and verify same handle per sequence.
            let inputs: Vec<Vec<u8>> = (0..1000)
                .map(|i| format!("SEQ{:04}", i).into_bytes())
                .collect();

            let n_threads = 8;
            let handles_per_thread: Vec<Vec<SeqHandle>> = (0..n_threads)
                .map(|_| {
                    let inputs = inputs.clone();
                    thread::spawn(move || {
                        inputs.iter().map(|s| seq_from_bytes(s)).collect::<Vec<_>>()
                    })
                })
                .map(|j| j.join().unwrap())
                .collect();

            // Compare thread 0 to all others
            for t in 1..n_threads {
                for i in 0..inputs.len() {
                    assert_eq!(
                        handles_per_thread[0][i], handles_per_thread[t][i],
                        "enabled mode: same bytes should yield same handle across threads"
                    );
                }
            }
        }

        /// Smoke test checking that we can insert lots of sequences without clashes
        #[test]
        fn handles_have_reasonable_hash_distribution_smoke() {
            // Not a statistical test; just ensure we can insert lots of distinct keys without weird behavior.
            let mut set: HashSet<SeqHandle> = HashSet::new();
            for i in 0..10_000 {
                let s = format!("X{:05}", i).into_bytes();
                set.insert(seq_from_bytes(&s));
            }
            assert_eq!(set.len(), 10_000);
        }
    }

    // Test non-interning module
    #[cfg(not(feature = "interning"))]
    mod disabled_tests {
        use super::*;

        /// Test SeqHandles of the same bytes are equivalent
        #[test]
        fn seq_not_required_to_dedup() {
            // In disabled mode, we do not promise dedup. Handles may or may not compare equal depending on design.
            // BUT: our current disabled design uses Arc<[u8]> so equality is by content, thus they SHOULD be equal.
            // This test checks the semantic we rely on: identical sequences compare equal (important for HashMap keys).
            let a1 = seq_from_bytes(b"ACGT");
            let a2 = seq_from_bytes(b"ACGT");
            assert_eq!(
                a1, a2,
                "disabled mode: identical bytes should still compare equal"
            );
        }

        /// Check ownership over round trip
        #[test]
        fn are_owned_and_round_trip() {
            let g = group_from_str("hello");
            let r = region_id_from_str("regionX");
            let l = library_id_from_str("seqX");

            assert_eq!(&*group_to_str(g), "hello");
            assert_eq!(&*region_id_to_str(r), "regionX");
            assert_eq!(&*library_id_to_str(l), "seqX");
        }
    }
}
