//! Facilities for interning and deduplicating frequently reused sequences and IDs.
//!
//! DNAComb stores many repeated sequence values and short identifiers such as
//! region names, library IDs, and read-group labels. In the default `interning`
//! configuration, these values are deduplicated into shared global interners and
//! referred to through compact handle types.
//!
//! Benefits of interning include:
//! - reduced memory usage when the same values occur many times,
//! - cheaper cloning and copying of handles,
//! - and fast equality/hash behaviour driven by compact identifiers.
//!
//! When the `interning` feature is disabled, the same public API is preserved
//! but handles wrap owned shared data directly instead of global interned IDs.
//! This can be preferable in library-style workflows where values are created
//! and dropped over time and global accumulation is undesirable.
//!
//! Public API available in both modes:
//! - `SeqHandle`, `seq_from_bytes()`, `seq_to_bytes()`
//! - `GroupID`, `group_id_from_str()`, `group_id_to_str()`
//! - `LibraryID`, `library_id_from_str()`, `library_id_to_str()`
//! - `RegionID`, `region_id_from_str()`, `region_id_to_str()`
//!
//! SeqHandle uses a custom implementation to deal with the specifics of sequence data
//! and Groups, IDs and Regions use ThreadedRodeo but are separated to allow
//! different tuning since they have different access patterns
//!
//! In interning mode, interned values are stored in global process-wide tables
//! and are not garbage collected during execution. This matches DNAComb's main
//! CLI workflow, where values accumulate until final output. The sequence interner
//! can be cleared to reset it between analyses, but this should only be done carefully as
//! it invalidates all previously issued SeqHandles and using them will either produce
//! junk data or panic. The ThreadedRodeo's used for strings can't be cleared so will
//! accumalate until the process ends, but they are expected to be reasonably small for
//!  most use cases so that shouldn't matter. If it does the non-interning backend is
//! the only fall back.
use std::sync::Arc;

#[cfg(feature = "interning")]
mod enabled {
    use crate::errors::seq_to_string_or_log;

    use super::Arc;
    use std::num::NonZeroU32;
    use std::sync::atomic::{AtomicBool, AtomicU32, Ordering};

    use ahash::RandomState;
    use dashmap::DashMap;
    use lasso::{Key, Spur, ThreadedRodeo};
    use once_cell::sync::Lazy;
    use parking_lot::RwLock;

    // Sequences (Vec<u8>)
    /// Shared unique SeqInterner instance storing interned sequences
    static SEQ_INTERNER: Lazy<SeqInterner> = Lazy::new(SeqInterner::default);

    /// Return a handle for the given sequence bytes.
    ///
    /// In interning mode, identical byte sequences resolve to the same canonical
    /// handle across the process. In non-interning mode, a new owned shared value
    /// may be created, but equality and hashing remain content-based.
    #[inline]
    pub fn seq_from_bytes(bytes: &[u8]) -> SeqHandle {
        SEQ_INTERNER.intern(bytes)
    }

    /// Resolve a sequence handle to shared sequence bytes.
    ///
    /// The returned value is cheap to clone and may share storage with other
    /// equal handles.
    #[inline]
    pub fn seq_to_bytes(h: &SeqHandle) -> Arc<[u8]> {
        SEQ_INTERNER.resolve(h)
    }

    /// Reserve capacity for additional unique sequences in the sequence interner.
    ///
    /// This is only meaningful in interning mode and can reduce allocation churn
    /// when the approximate number of unique sequences is known in advance.
    pub fn reserve_seq_interner(estimated_unique: usize) {
        SEQ_INTERNER.reserve(estimated_unique);
    }

    /// Return the number of sequence entries currently stored in the reverse map.
    #[inline]
    pub fn num_interned_reverse() -> usize {
        SEQ_INTERNER.num_interned_reverse()
    }

    /// Return the number of distinct canonical sequences currently stored in the forward map.
    #[inline]
    pub fn num_interned_forward() -> usize {
        SEQ_INTERNER.num_interned_forward()
    }

    /// Clear all interned sequences.
    ///
    /// WARNING: This only is intended for advanced library users running multiple
    /// separate analyses in a long-lived process. It should only be called when
    /// all previously returned SeqHandle`s and all data structures containing
    /// them have been dropped. After clearing, resolving or reusing old sequence
    /// handles is invalid and will either return incorrect data or panic.
    pub fn clear_seq_interner() {
        SEQ_INTERNER.clear();
    }

    /// Handle identifying an observed sequence.
    ///
    /// `SeqHandle` is copyable, hashable, and equality-comparable. In interning
    /// mode it is backed by a compact non-zero integer ID, allowing
    /// `Option<SeqHandle>` to use niche optimisation and remain the same size as
    /// `SeqHandle` itself. In non-interning mode it is a thin wrapper around an owned
    /// Vec<u8>
    #[repr(transparent)]
    #[derive(Clone, Eq, PartialEq, Hash, Debug)]
    pub struct SeqHandle(NonZeroU32);

    impl SeqHandle {
        /// Return the raw underlying handle value.
        ///
        /// This is mainly useful for diagnostics and internal plumbing rather
        /// than normal library use.
        #[inline]
        pub fn get(&self) -> u32 {
            self.0.get()
        }

        #[inline]
        fn from_index0(i: u32) -> Self {
            // Stored as 1-based so Option<SeqHandle> is niche-optimized
            Self::from_raw(i + 1)
        }

        #[inline]
        fn from_raw(i: u32) -> Self {
            SeqHandle(NonZeroU32::new(i).expect("nonzero"))
        }

        #[inline]
        fn index0(&self) -> usize {
            (self.0.get() - 1) as usize
        }

        /// Length of the referenced sequence in bases/bytes.
        #[inline]
        pub fn len(&self) -> usize {
            seq_to_bytes(self).len()
        }

        /// Return `true` if the referenced sequence is empty.
        #[inline]
        pub fn is_empty(&self) -> bool {
            seq_to_bytes(self).is_empty()
        }

        /// Convert the sequence to a UTF-8 string for output, logging failures.
        ///
        /// This is intended for diagnostics and table writing; non-UTF-8 content
        /// yields an empty string after logging a warning.
        #[inline]
        pub fn to_str_or_log(&self) -> String {
            seq_to_string_or_log(&seq_to_bytes(self).to_vec())
        }
    }

    /// Forward-map entry for one canonical sequence.
    ///
    /// handle_raw stores the raw NonZeroU32 payload used by SeqHandle:
    ///   0 => not yet published
    ///   >0 => published and safe to use
    ///
    /// init_claimed elects a single thread to perform the first-time publication.
    struct ForwardEntry {
        bytes: Arc<[u8]>,
        handle_raw: AtomicU32,
        init_claimed: AtomicBool,
    }

    impl ForwardEntry {
        #[inline]
        pub fn new_unclaimed(bytes: Arc<[u8]>) -> Self {
            ForwardEntry {
                bytes,
                handle_raw: AtomicU32::new(0),
                init_claimed: AtomicBool::new(false),
            }
        }
    }

    /// Global concurrent sequence interner.
    ///
    /// Maintains:
    /// - a forward map from canonical byte content to handle metadata,
    /// - a reverse table from handle to canonical bytes,
    /// - and a monotonic ID allocator for new unique sequences.
    pub struct SeqInterner {
        /// Forward map: canonical bytes -> id
        forward: DashMap<Arc<[u8]>, Arc<ForwardEntry>, RandomState>,

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
        /// Intern a sequence and return its canonical handle.
        ///
        /// Concurrent callers interning the same byte sequence will converge on
        /// the same published handle.
        #[inline]
        pub fn intern(&self, bytes: &[u8]) -> SeqHandle {
            // First, get or set reference to the forward interner
            let entry = match self.forward.get(bytes) {
                Some(found) => Arc::clone(found.value()),
                None => {
                    let arc = Arc::<[u8]>::from(bytes);
                    let entry_ref = self
                        .forward
                        .entry(arc.clone())
                        .or_insert_with(|| Arc::new(ForwardEntry::new_unclaimed(arc)));

                    Arc::clone(entry_ref.value())
                }
            };

            // Loop until handle is initialised
            loop {
                let raw = entry.handle_raw.load(Ordering::Acquire);

                // Fast path: if handle is already published, return it
                if raw != 0 {
                    return SeqHandle::from_raw(raw);
                }

                // Otherwise try to claim the right to initialize
                if entry
                    .init_claimed
                    .compare_exchange(false, true, Ordering::AcqRel, Ordering::Acquire)
                    .is_ok()
                {
                    // We won the race. Allocate and publish the handle.
                    let idx0 = self.next.fetch_add(1, Ordering::Relaxed);
                    if idx0 == u32::MAX {
                        panic!("SeqInterner exhausted SeqHandle space");
                    }

                    let handle = SeqHandle::from_index0(idx0);

                    {
                        let mut rev = self.reverse.write();

                        let pos = idx0 as usize;
                        match rev.len().cmp(&pos) {
                            std::cmp::Ordering::Less => {
                                // Extremely rare; keep safe
                                rev.resize_with(pos, || Arc::<[u8]>::from(&b""[..]));
                                rev.push(entry.bytes.clone())
                            }
                            std::cmp::Ordering::Equal => rev.push(entry.bytes.clone()),
                            std::cmp::Ordering::Greater => rev[pos] = entry.bytes.clone(),
                        }
                    }

                    // Publish handle with Release so all threads see it
                    entry.handle_raw.store(handle.get(), Ordering::Release);
                    return handle;
                }

                // Lost race - yield and loop to wait for Release.
                std::thread::yield_now();
            }
        }

        /// Resolve a handle back to its canonical shared byte sequence.
        #[inline]
        pub fn resolve(&self, id: &SeqHandle) -> Arc<[u8]> {
            self.reverse.read()[id.index0()].clone()
        }

        /// Reserve space in the interner
        pub fn reserve(&self, additional_unique: usize) {
            self.reverse.write().reserve(additional_unique);
        }

        /// Return the number of sequence entries currently stored in the reverse map.
        #[inline]
        pub fn num_interned_reverse(&self) -> usize {
            self.reverse.read().len()
        }

        /// Return the number of distinct canonical sequences currently stored in the forward map.
        #[inline]
        pub fn num_interned_forward(&self) -> usize {
            self.forward.len()
        }

        /// Remove all interned sequence state.
        ///
        /// WARNING: This removes all sequences from the interner and
        /// invalidates all previously issued `SeqHandle`s. Use of previous
        /// `SeqHandle`s will produce nonsense results or panic. Only for
        /// advanced use on long running processes between fully distinct
        /// analyses.
        pub fn clear(&self) {
            self.forward.clear();
            self.reverse.write().clear();
            self.next.store(0, Ordering::Release);
        }
    }

    // GroupKey reserving 2 bits for None, Unmatched and Ungrouped
    #[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
    struct GroupKey(NonZeroU32);

    // Implement Lasso::Key for group Key
    unsafe impl Key for GroupKey {
        fn into_usize(self) -> usize {
            // Can do this as we add 0b100 previously
            (self.0.get() - 0b100) as usize
        }

        fn try_from_usize(int: usize) -> Option<Self> {
            // Reject values outside the niche range - Lasso starts at 0 so at 0b100 e.g.
            // the first bit we accept as outside the niche
            Some(Self(
                NonZeroU32::new(int.checked_add(0b100)? as u32).expect("Checked > 0b11"),
            ))
        }
    }

    // Group IDs (Str)
    /// Shared unique interner instance storing Group names
    static GROUP_IDS: Lazy<ThreadedRodeo<GroupKey>> = Lazy::new(ThreadedRodeo::<GroupKey>::new);

    /// Interned identifier for a read-group label string.
    #[repr(transparent)]
    #[derive(Clone, Eq, PartialEq, Hash, Debug)]
    pub struct GroupID(GroupKey);

    /// Get a GroupID for a read group string, interning it if it's new
    #[inline]
    pub fn group_id_from_str(s: &str) -> GroupID {
        GroupID(GROUP_IDS.get_or_intern(s))
    }

    /// Resolve GroupID to owned string (Arc<str>)
    #[inline]
    pub fn group_id_to_str(id: &GroupID) -> Arc<str> {
        Arc::<str>::from(GROUP_IDS.resolve(&id.0))
    }

    /// Get a GroupID from a raw int
    #[inline]
    pub fn group_id_from_raw(i: NonZeroU32) -> GroupID {
        GroupID(GroupKey(i))
    }

    /// Resolve GroupID to raw int
    #[inline]
    pub fn group_id_to_raw(id: &GroupID) -> NonZeroU32 {
        id.0.0
    }

    // Library IDs (Str)
    /// Shared unique interner instance storing Group names
    static LIB_IDS: Lazy<ThreadedRodeo> = Lazy::new(ThreadedRodeo::default);

    /// Global interned Library ID
    #[repr(transparent)]
    #[derive(Clone, Eq, PartialEq, Hash, Debug)]
    pub struct LibraryID(Spur);

    /// Get a LibraryID for a string name, interning it if it's new
    #[inline]
    pub fn library_id_from_str(s: &str) -> LibraryID {
        LibraryID(LIB_IDS.get_or_intern(s))
    }

    /// Resolve Library ID to owned string (Arc<str>)
    #[inline]
    pub fn library_id_to_str(id: &LibraryID) -> Arc<str> {
        Arc::<str>::from(LIB_IDS.resolve(&id.0))
    }

    // Region IDs (Str)
    /// Shared unique interner instance storing Region IDs
    static REGIONS: Lazy<ThreadedRodeo> = Lazy::new(ThreadedRodeo::default);

    /// Global interned region ID
    #[repr(transparent)]
    #[derive(Clone, Eq, PartialEq, Hash, Debug)]
    pub struct RegionID(Spur);

    /// Get a RegionHandle for a region ID, interning it if it's new
    #[inline]
    pub fn region_id_from_str(s: &str) -> RegionID {
        RegionID(REGIONS.get_or_intern(s))
    }

    /// Resolve region handle to owned string (Arc<str>).
    #[inline]
    pub fn region_id_to_str(id: &RegionID) -> Arc<str> {
        Arc::<str>::from(REGIONS.resolve(&id.0))
    }
}

#[cfg(not(feature = "interning"))]
mod disabled {
    use crate::errors::seq_to_string_or_log;

    use super::Arc;

    /// Shared owned sequence handle used when interning is disabled.
    ///
    /// Unlike the interned backend, identical sequences are not guaranteed to
    /// share a global canonical ID, but equality and hashing still behave by
    /// sequence content so the public API remains semantically compatible.
    #[derive(Eq, PartialEq, Hash, Debug)]
    pub struct SeqHandle(pub Arc<[u8]>);

    impl SeqHandle {
        /// Length of the referenced sequence in bases/bytes.
        #[inline]
        pub fn len(&self) -> usize {
            seq_to_bytes(self).len()
        }

        /// Return `true` if the referenced sequence is empty.
        #[inline]
        pub fn is_empty(&self) -> bool {
            seq_to_bytes(self).is_empty()
        }

        /// Convert the sequence to a UTF-8 string for output, logging failures.
        ///
        /// This is intended for diagnostics and table writing; non-UTF-8 content
        /// yields an empty string after logging a warning.
        #[inline]
        pub fn to_str_or_log(&self) -> String {
            seq_to_string_or_log(&seq_to_bytes(self).to_vec())
        }
    }

    impl Clone for SeqHandle {
        #[inline(always)]
        fn clone(&self) -> Self {
            SeqHandle(Arc::clone(&self.0))
        }
    }

    /// Return a handle for the given sequence bytes.
    ///
    /// In interning mode, identical byte sequences resolve to the same canonical
    /// handle across the process. In non-interning mode, a new owned shared value
    /// may be created, but equality and hashing remain content-based.
    #[inline]
    pub fn seq_from_bytes(bytes: &[u8]) -> SeqHandle {
        SeqHandle(Arc::<[u8]>::from(bytes))
    }

    /// Resolve a sequence handle to shared sequence bytes.
    ///
    /// The returned value is cheap to clone and may share storage with other
    /// equal handles.
    #[inline]
    pub fn seq_to_bytes(h: &SeqHandle) -> Arc<[u8]> {
        h.0.clone()
    }

    /// Non-interning reservation NoOp as no underlying interner)
    pub fn reserve_seq_interner(_estimated_unique: usize) {}

    /// Return the number of sequence entries currently stored in the reverse map.
    ///
    /// In the non-interning backend this is always 0
    #[inline]
    pub fn num_interned_reverse() -> usize {
        0
    }

    /// Return the number of distinct canonical sequences currently stored in the forward map.
    ///
    /// In the non-interning backend this is always 0
    #[inline]
    pub fn num_interned_forward() -> usize {
        0
    }

    /// Clear the sequence interning
    ///
    /// In the non-interning backend this is just a no-op, while in the interning
    /// backend it clears all currently interned sequences for resetting between
    /// distinct analyses.
    #[inline]
    pub fn clear_seq_interner() {}

    // Group IDs (Str)
    /// Shared owned group identifier used when interning is disabled.
    #[derive(Eq, PartialEq, Hash, Debug)]
    pub struct GroupID(pub Arc<str>);

    impl Clone for GroupID {
        #[inline(always)]
        fn clone(&self) -> Self {
            GroupID(Arc::clone(&self.0))
        }
    }

    /// Get a GroupHandle for a Group ID (non-interning)
    #[inline]
    pub fn group_id_from_str(s: &str) -> GroupID {
        GroupID(Arc::<str>::from(s))
    }

    /// Resolve non-interning group handle to owned string
    #[inline]
    pub fn group_id_to_str(h: &GroupID) -> Arc<str> {
        h.0.clone()
    }

    // Library IDs (Str)
    /// Shared owned library identifier used when interning is disabled.
    #[derive(Eq, PartialEq, Hash, Debug)]
    pub struct LibraryID(pub Arc<str>);

    impl Clone for LibraryID {
        #[inline(always)]
        fn clone(&self) -> Self {
            LibraryID(Arc::clone(&self.0))
        }
    }

    /// Get a GroupHandle for a Group ID (non-interning)
    #[inline]
    pub fn library_id_from_str(s: &str) -> LibraryID {
        LibraryID(Arc::<str>::from(s))
    }

    /// Resolve non-interning group handle to owned string
    #[inline]
    pub fn library_id_to_str(h: &LibraryID) -> Arc<str> {
        h.0.clone()
    }

    // Region IDs (Str)
    /// Shared owned region identifier used when interning is disabled.
    #[derive(Eq, PartialEq, Hash, Debug)]
    pub struct RegionID(pub Arc<str>);

    impl Clone for RegionID {
        #[inline(always)]
        fn clone(&self) -> Self {
            RegionID(Arc::clone(&self.0))
        }
    }

    /// Get a RegionHandle for a region ID (non-interning)
    #[inline]
    pub fn region_id_from_str(s: &str) -> RegionID {
        RegionID(Arc::<str>::from(s))
    }

    /// Resolve non-interning region handles to owned string (Arc<str>)
    #[inline]
    pub fn region_id_to_str(h: &RegionID) -> Arc<str> {
        h.0.clone()
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
        let out = seq_to_bytes(&h);
        assert_eq!(&*out, b"ACGTN");
    }

    /// Test Group round trip
    #[test]
    fn group_round_trip_str() {
        let g = group_id_from_str("groupA");
        let s = group_id_to_str(&g);
        assert_eq!(&*s, "groupA");
    }

    /// Test Library round trip
    #[test]
    fn library_round_trip_str() {
        let l = library_id_from_str("seqA");
        let s = library_id_to_str(&l);
        assert_eq!(&*s, "seqA");
    }

    /// Test Region round trip
    #[test]
    fn region_round_trip_str() {
        let r = region_id_from_str("Reg1");
        let s = region_id_to_str(&r);
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

    #[test]
    fn seq_round_trip_empty() {
        let h = seq_from_bytes(b"");
        let out = seq_to_bytes(&h);
        assert_eq!(&*out, b"");
    }

    #[test]
    fn seq_len_method() {
        let h = seq_from_bytes(b"ACGTACGT");
        assert_eq!(h.len(), 8);
        let h_empty = seq_from_bytes(b"");
        assert_eq!(h_empty.len(), 0);
    }

    #[test]
    fn seq_is_empty_method() {
        let h = seq_from_bytes(b"A");
        assert!(!h.is_empty());
        let h_empty = seq_from_bytes(b"");
        assert!(h_empty.is_empty());
    }

    #[test]
    fn seq_to_str_or_log_valid_utf8() {
        let h = seq_from_bytes(b"ACGT");
        let s = h.to_str_or_log();
        assert_eq!(s, "ACGT");
    }

    #[test]
    fn seq_to_str_or_log_empty() {
        let h = seq_from_bytes(b"");
        let s = h.to_str_or_log();
        assert_eq!(s, "");
    }

    #[test]
    fn group_round_trip_empty_string() {
        let g = group_id_from_str("");
        let s = group_id_to_str(&g);
        assert_eq!(&*s, "");
    }

    #[test]
    fn group_round_trip_special_chars() {
        let g = group_id_from_str("pool-A_123.xyz");
        let s = group_id_to_str(&g);
        assert_eq!(&*s, "pool-A_123.xyz");
    }

    #[test]
    fn library_round_trip_empty_string() {
        let l = library_id_from_str("");
        let s = library_id_to_str(&l);
        assert_eq!(&*s, "");
    }

    #[test]
    fn library_round_trip_special_chars() {
        let l = library_id_from_str("seq-001_variant.v2");
        let s = library_id_to_str(&l);
        assert_eq!(&*s, "seq-001_variant.v2");
    }

    #[test]
    fn region_round_trip_empty_string() {
        let r = region_id_from_str("");
        let s = region_id_to_str(&r);
        assert_eq!(&*s, "");
    }

    #[test]
    fn region_round_trip_special_chars() {
        let r = region_id_from_str("Region_A-1.upstream");
        let s = region_id_to_str(&r);
        assert_eq!(&*s, "Region_A-1.upstream");
    }

    #[test]
    fn seq_handles_clone() {
        let h1 = seq_from_bytes(b"ACGT");
        let h2 = h1.clone();
        assert_eq!(h1, h2);
    }

    #[test]
    fn group_handles_clone() {
        let g1 = group_id_from_str("group1");
        let g2 = g1.clone();
        assert_eq!(g1, g2);
    }

    #[test]
    fn library_handles_copy_and_clone() {
        let l1 = library_id_from_str("lib1");
        let l2 = l1.clone();
        assert_eq!(l1, l2);
    }

    #[test]
    fn region_handles_copy_and_clone() {
        let r1 = region_id_from_str("region1");
        let r2 = r1.clone();
        assert_eq!(r1, r2);
    }

    #[test]
    fn seq_handles_in_hashmap() {
        let h1 = seq_from_bytes(b"ACGT");
        let h2 = seq_from_bytes(b"GGGG");
        let h3 = seq_from_bytes(b"ACGT");
        let mut m: HashMap<SeqHandle, u32> = HashMap::new();
        m.insert(h1, 1);
        m.insert(h2, 2);
        assert_eq!(m.get(&h3), Some(&1));
        assert_eq!(m.len(), 2);
    }

    #[test]
    fn group_handles_in_hashmap() {
        let g1 = group_id_from_str("a");
        let g2 = group_id_from_str("b");
        let g3 = group_id_from_str("a");
        let mut m: HashMap<GroupID, u32> = HashMap::new();
        m.insert(g1, 1);
        m.insert(g2, 2);
        assert_eq!(m.get(&g3), Some(&1));
    }

    #[test]
    fn library_handles_in_hashmap() {
        let l1 = library_id_from_str("x");
        let l2 = library_id_from_str("y");
        let l3 = library_id_from_str("x");
        let mut m: HashMap<LibraryID, u32> = HashMap::new();
        m.insert(l1, 1);
        m.insert(l2, 2);
        assert_eq!(m.get(&l3), Some(&1));
    }

    #[test]
    fn region_handles_in_hashmap() {
        let r1 = region_id_from_str("alpha");
        let r2 = region_id_from_str("beta");
        let r3 = region_id_from_str("alpha");
        let mut m: HashMap<RegionID, u32> = HashMap::new();
        m.insert(r1, 1);
        m.insert(r2, 2);
        assert_eq!(m.get(&r3), Some(&1));
    }

    #[test]
    fn reserve_seq_interner_smoke_test() {
        // Should not panic regardless of feature
        reserve_seq_interner(1000);
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
            let p1 = Arc::as_ptr(&seq_to_bytes(&a1));
            let p2 = Arc::as_ptr(&seq_to_bytes(&a2));
            assert_eq!(
                p1, p2,
                "interning enabled: resolved Arc should be same allocation"
            );
        }

        /// Test Groups get the same handle
        #[test]
        fn group_dedup_same_str_get_same_handle() {
            let g1 = group_id_from_str("x");
            let g2 = group_id_from_str("x");
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

        /// Test multiple insertions of the same value don't add extra IDs to reverse
        #[test]
        fn concurrent_duplicate_inserts_do_not_allocate_extra_ids() {
            let interner = Arc::new(SeqInterner::default());
            let before_rev = interner.num_interned_reverse();
            let before_fwd = interner.num_interned_forward();

            let n_threads = 32;
            let n_iters = 1000;

            std::thread::scope(|s| {
                for _ in 0..n_threads {
                    s.spawn(|| {
                        for _ in 0..n_iters {
                            let _ = interner.intern(b"ACGT");
                        }
                    });
                }
            });

            let after_rev = interner.num_interned_reverse();
            let after_fwd = interner.num_interned_forward();
            assert_eq!(
                after_rev,
                before_rev + 1,
                "Concurrent inserts should add exactly 1 reverse ID"
            );
            assert_eq!(
                after_fwd,
                before_fwd + 1,
                "Concurrent inserts should add exactly 1 forward ID"
            );
        }

        /// Test concurrent inserts add exactly 1 ID each
        #[test]
        fn concurrent_unique_inserts_allocate_exactly_one_id_each() {
            let interner = Arc::new(SeqInterner::default());
            let before_rev = interner.num_interned_reverse();
            let before_fwd = interner.num_interned_forward();

            let n = 5_000;
            let seqs: Vec<Vec<u8>> = (0..n)
                .map(|i| format!("SEQ{:05}", i).into_bytes())
                .collect();

            std::thread::scope(|scope| {
                let chunk = 500;
                for part in seqs.chunks(chunk) {
                    let thread_interner = interner.clone();
                    scope.spawn(move || {
                        for s in part {
                            let _ = thread_interner.intern(s);
                        }
                    });
                }
            });

            let after_rev = interner.num_interned_reverse();
            let after_fwd = interner.num_interned_forward();
            assert_eq!(
                after_rev,
                before_rev + n,
                "unique inserts should allocate exactly one reverse entry per unique sequence"
            );
            assert_eq!(
                after_fwd,
                before_fwd + n,
                "unique inserts should allocate exactly one forward entry per unique sequence"
            );
        }

        /// Test concurrent interning and resolving
        #[test]
        fn concurrent_intern_and_resolve_is_safe() {
            let n_threads = 16;
            let n_iters = 5_000;

            std::thread::scope(|scope| {
                for _ in 0..n_threads {
                    scope.spawn(|| {
                        for _ in 0..n_iters {
                            let h = seq_from_bytes(b"GATTACA");
                            let seq = seq_to_bytes(&h);
                            assert_eq!(seq.as_ref(), b"GATTACA");
                        }
                    });
                }
            });
        }

        #[test]
        fn seq_get_raw_value() {
            let h = seq_from_bytes(b"ACGT");
            let raw = h.get();
            assert!(raw > 0, "SeqHandle should have non-zero raw value");
        }

        #[test]
        fn seq_handles_large_sequence() {
            let large_seq = vec![b'A'; 100_000];
            let h = seq_from_bytes(&large_seq);
            assert_eq!(h.len(), 100_000);
            let resolved = seq_to_bytes(&h);
            assert_eq!(resolved.len(), 100_000);
        }

        #[test]
        fn seq_handles_various_bytes() {
            // Not just valid DNA bases
            let h = seq_from_bytes(&[0, 255, 127, 128, 1, 2]);
            let resolved = seq_to_bytes(&h);
            assert_eq!(resolved.as_ref(), &[0, 255, 127, 128, 1, 2]);
        }

        #[test]
        fn multiple_different_sequences() {
            let h1 = seq_from_bytes(b"AAA");
            let h2 = seq_from_bytes(b"AAB");
            let h3 = seq_from_bytes(b"ABB");
            let h4 = seq_from_bytes(b"BBB");

            assert_ne!(h1, h2);
            assert_ne!(h2, h3);
            assert_ne!(h3, h4);
            assert_ne!(h1, h4);
        }

        #[test]
        fn group_concurrent_different_strings() {
            let n_threads = 10;
            let n_groups = 100;

            std::thread::scope(|scope| {
                for t in 0..n_threads {
                    scope.spawn(move || {
                        for i in 0..n_groups {
                            let name = format!("group_{}_{}", t, i);
                            let id = group_id_from_str(&name);
                            let resolved = group_id_to_str(&id);
                            assert_eq!(&*resolved, &name);
                        }
                    });
                }
            });
        }

        #[test]
        fn region_concurrent_different_strings() {
            let n_threads = 10;
            let n_regions = 100;

            std::thread::scope(|scope| {
                for t in 0..n_threads {
                    scope.spawn(move || {
                        for i in 0..n_regions {
                            let name = format!("region_{}_{}", t, i);
                            let id = region_id_from_str(&name);
                            let resolved = region_id_to_str(&id);
                            assert_eq!(&*resolved, &name);
                        }
                    });
                }
            });
        }

        #[test]
        fn library_concurrent_different_strings() {
            let n_threads = 10;
            let n_libs = 100;

            std::thread::scope(|scope| {
                for t in 0..n_threads {
                    scope.spawn(move || {
                        for i in 0..n_libs {
                            let name = format!("lib_{}_{}", t, i);
                            let id = library_id_from_str(&name);
                            let resolved = library_id_to_str(&id);
                            assert_eq!(&*resolved, &name);
                        }
                    });
                }
            });
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
            let g = group_id_from_str("hello");
            let r = region_id_from_str("regionX");
            let l = library_id_from_str("seqX");

            assert_eq!(&*group_id_to_str(&g), "hello");
            assert_eq!(&*region_id_to_str(&r), "regionX");
            assert_eq!(&*library_id_to_str(&l), "seqX");
        }

        #[test]
        fn seq_arc_equality_by_content() {
            let a1 = seq_from_bytes(b"ACGT");
            let a2 = seq_from_bytes(b"ACGT");
            let a3 = seq_from_bytes(b"GGGG");
            assert_eq!(a1, a2);
            assert_ne!(a1, a3);
        }

        #[test]
        fn group_arc_equality_by_content() {
            let g1 = group_id_from_str("test");
            let g2 = group_id_from_str("test");
            let g3 = group_id_from_str("other");
            assert_eq!(g1, g2);
            assert_ne!(g1, g3);
        }

        #[test]
        fn library_arc_equality_by_content() {
            let l1 = library_id_from_str("test");
            let l2 = library_id_from_str("test");
            let l3 = library_id_from_str("other");
            assert_eq!(l1, l2);
            assert_ne!(l1, l3);
        }

        #[test]
        fn region_arc_equality_by_content() {
            let r1 = region_id_from_str("test");
            let r2 = region_id_from_str("test");
            let r3 = region_id_from_str("other");
            assert_eq!(r1, r2);
            assert_ne!(r1, r3);
        }
    }
}
