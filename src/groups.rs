//! Read-group identifiers used to partition counts.
//!
//! A `ReadGroup` represents the grouping assigned to a read during parsing,
//! for example from a regex capture on the read name. Two sentinel values are
//! also supported:
//! - `ungrouped`, for workflows where no grouping is configured,
//! - `unmatched`, for reads where grouping was requested but no capture was found.
//!
//! This module provides both interned and non-interned implementations behind
//! the same public API.
#[cfg(feature = "interning")]
mod enabled {
    use std::fmt;
    use std::num::NonZeroU32;

    use crate::interning::{
        group_id_from_raw, group_id_from_str, group_id_to_raw, group_id_to_str,
    };

    /// Compact read-group identifier with special sentinel values for
    /// ungrouped and unmatched reads.
    ///
    /// In the interning backend, named groups are stored as interned IDs while
    /// `ungrouped` and `unmatched` are represented as reserved flag values.
    /// Internally this uses a NonZeroU32 for Option niche optimisation and the
    /// first 2 bits are reserved for the sentinel values.
    #[repr(transparent)]
    #[derive(Debug, Clone, Eq, PartialEq, Hash)]
    pub struct ReadGroup(NonZeroU32);

    impl ReadGroup {
        const UNGROUPED: u32 = 0b01;
        const UNMATCHED: u32 = 0b10;
        const FLAG_BITS: u32 = 0b11;
        // 4+ reserved for interned GroupIDs

        /// Construct the sentinel value used when no grouping is configured.
        #[inline]
        pub fn ungrouped() -> Self {
            ReadGroup(NonZeroU32::new(Self::UNGROUPED).expect("Know this is 1"))
        }

        /// Construct the sentinel value used when grouping is configured but
        /// no group could be extracted for a read.
        #[inline]
        pub fn unmatched() -> Self {
            ReadGroup(NonZeroU32::new(Self::UNMATCHED).expect("Know this is 2"))
        }

        /// Create a ReadGroup from a group name string (interning if needed)
        #[inline]
        pub fn grouped(s: &str) -> Self {
            // Send string to the interner, retrieving the ID and interning if necessary
            let group_id = group_id_from_str(s);
            ReadGroup(group_id_to_raw(&group_id))
        }

        /// Return `true` if this read belongs to the `ungrouped` sentinel class.
        #[inline]
        pub fn is_ungrouped(&self) -> bool {
            self.0.get() == Self::UNGROUPED
        }

        /// Return `true` if grouping was attempted but no group was matched.
        #[inline]
        pub fn is_unmatched(&self) -> bool {
            self.0.get() == Self::UNMATCHED
        }

        /// Return `true` if this value represents a concrete named group.
        #[inline]
        pub fn is_match(&self) -> bool {
            self.0.get() > Self::FLAG_BITS
        }
    }

    /// Display as the stable output form used in TSVs:
    /// - `ungrouped` renders as an empty string,
    /// - `unmatched` renders as `_unmatched_`,
    /// - named groups render as their group label.
    impl fmt::Display for ReadGroup {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            if self.is_match() {
                write!(f, "{}", group_id_to_str(&group_id_from_raw(self.0)))
            } else if self.is_unmatched() {
                write!(f, "_unmatched_")
            } else if self.is_ungrouped() {
                write!(f, "")
            } else {
                // Panic here as should never have values outside the allowed options - means a bug
                // or data corruption
                write!(f, "{:?}", self)
                //panic!("ReadGroup with illegal NonZeroU32 value")
            }
        }
    }
}

#[cfg(not(feature = "interning"))]
mod disabled {
    use std::fmt;

    /// Read-group assignment for a read.
    ///
    /// This enum distinguishes three cases:
    /// - `Ungrouped`: no grouping was requested,
    /// - `Unmatched`: grouping was requested but no capture was found,
    /// - `Match(String)`: a concrete extracted group label.
    #[derive(Debug, Clone, Eq, PartialEq, Hash)]
    pub enum ReadGroup {
        Ungrouped,
        Unmatched,
        Match(String),
    }

    impl ReadGroup {
        /// Construct the sentinel value used when no grouping is configured.
        pub fn ungrouped() -> Self {
            Self::Ungrouped
        }

        /// Construct the sentinel value used when grouping is configured but
        /// no group could be extracted for a read.
        pub fn unmatched() -> Self {
            Self::Unmatched
        }

        /// Construct a named read group from an extracted group label.
        pub fn grouped(s: &str) -> Self {
            ReadGroup::Match(s.to_string())
        }

        /// Return `true` if this read belongs to the `ungrouped` sentinel class.
        pub fn is_ungrouped(&self) -> bool {
            matches!(self, ReadGroup::Ungrouped)
        }

        /// Return `true` if grouping was attempted but no group was matched.
        pub fn is_unmatched(&self) -> bool {
            matches!(self, ReadGroup::Unmatched)
        }

        /// Return `true` if this value represents a concrete named group.
        pub fn is_match(&self) -> bool {
            matches!(self, ReadGroup::Match(_))
        }
    }

    /// Display as the stable output form used in TSVs:
    /// - `ungrouped` renders as an empty string,
    /// - `unmatched` renders as `_unmatched_`,
    /// - named groups render as their group label.
    impl fmt::Display for ReadGroup {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            match self {
                ReadGroup::Ungrouped => write!(f, ""),
                ReadGroup::Unmatched => write!(f, "_unmatched_"),
                ReadGroup::Match(x) => write!(f, "{}", x),
            }
        }
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

    /// ReadGroup Display: verify string forms are stable (UX-facing).
    #[test]
    fn readgroup_display_variants() {
        assert_eq!(ReadGroup::ungrouped().to_string(), "");
        assert_eq!(ReadGroup::unmatched().to_string(), "_unmatched_");

        let g = ReadGroup::grouped("poolA");
        assert_eq!(g.to_string(), "poolA");
    }

    #[test]
    fn readgroup_grouped_multiple_names() {
        let g1 = ReadGroup::grouped("pool_1");
        let g2 = ReadGroup::grouped("pool_2");
        let g3 = ReadGroup::grouped("pool_1");

        assert_eq!(g1.to_string(), "pool_1");
        assert_eq!(g2.to_string(), "pool_2");
        assert_eq!(g1.to_string(), g3.to_string());
    }

    #[test]
    fn readgroup_grouped_empty_string() {
        let g = ReadGroup::grouped("");
        assert!(g.is_match());
        assert_eq!(g.to_string(), "");
    }

    #[test]
    fn readgroup_grouped_special_chars() {
        let g = ReadGroup::grouped("pool-A_123.xyz");
        assert!(g.is_match());
        assert_eq!(g.to_string(), "pool-A_123.xyz");
    }

    #[test]
    fn readgroup_grouped_whitespace() {
        let g = ReadGroup::grouped("pool A");
        assert!(g.is_match());
        assert_eq!(g.to_string(), "pool A");
    }

    #[test]
    fn readgroup_ungrouped_predicates() {
        let ug = ReadGroup::ungrouped();
        assert!(ug.is_ungrouped());
        assert!(!ug.is_unmatched());
        assert!(!ug.is_match());
    }

    #[test]
    fn readgroup_unmatched_predicates() {
        let um = ReadGroup::unmatched();
        assert!(!um.is_ungrouped());
        assert!(um.is_unmatched());
        assert!(!um.is_match());
    }

    #[test]
    fn readgroup_grouped_predicates() {
        let m = ReadGroup::grouped("test");
        assert!(!m.is_ungrouped());
        assert!(!m.is_unmatched());
        assert!(m.is_match());
    }

    #[test]
    fn readgroup_equality_ungrouped() {
        let ug1 = ReadGroup::ungrouped();
        let ug2 = ReadGroup::ungrouped();
        assert_eq!(ug1, ug2);
    }

    #[test]
    fn readgroup_equality_unmatched() {
        let um1 = ReadGroup::unmatched();
        let um2 = ReadGroup::unmatched();
        assert_eq!(um1, um2);
    }

    #[test]
    fn readgroup_equality_grouped_same() {
        let g1 = ReadGroup::grouped("pool_A");
        let g2 = ReadGroup::grouped("pool_A");
        assert_eq!(g1, g2);
    }

    #[test]
    fn readgroup_equality_grouped_different() {
        let g1 = ReadGroup::grouped("pool_A");
        let g2 = ReadGroup::grouped("pool_B");
        assert_ne!(g1, g2);
    }

    #[test]
    fn readgroup_inequality_sentinels() {
        let ug = ReadGroup::ungrouped();
        let um = ReadGroup::unmatched();
        assert_ne!(ug, um);
    }

    #[test]
    fn readgroup_inequality_sentinel_vs_grouped() {
        let ug = ReadGroup::ungrouped();
        let m = ReadGroup::grouped("ungrouped");
        assert_ne!(ug, m);

        let um = ReadGroup::unmatched();
        let m2 = ReadGroup::grouped("unmatched");
        assert_ne!(um, m2);
    }

    #[test]
    fn readgroup_hash_consistency() {
        use std::collections::HashSet;

        let g1 = ReadGroup::grouped("pool_A");
        let g2 = ReadGroup::grouped("pool_A");

        let mut set = HashSet::new();
        set.insert(g1);
        assert!(set.contains(&g2));
    }

    #[test]
    fn readgroup_hash_different_values() {
        use std::collections::HashSet;

        let g1 = ReadGroup::grouped("pool_A");
        let g2 = ReadGroup::grouped("pool_B");

        let mut set = HashSet::new();
        set.insert(g1);
        set.insert(g2);
        assert_eq!(set.len(), 2);
    }

    #[test]
    fn readgroup_clone_equality() {
        let g1 = ReadGroup::grouped("pool_A");
        let g2 = g1.clone();
        assert_eq!(g1, g2);
    }

    #[cfg(feature = "interning")]
    #[test]
    fn niche_optimization_verified() {
        assert_eq!(std::mem::size_of::<ReadGroup>(), 4);
        assert_eq!(
            std::mem::size_of::<Option<ReadGroup>>(),
            4,
            "Option<ReadGroup> should be niche-optimized to 4 bytes"
        );
    }

    #[cfg(not(feature = "interning"))]
    #[test]
    fn readgroup_non_interning_backend_clone() {
        let g1 = ReadGroup::grouped("pool_A");
        let g2 = g1.clone();
        assert_eq!(g1, g2);
        match g1 {
            ReadGroup::Match(ref s) => assert_eq!(s, "pool_A"),
            _ => panic!("Expected Match variant"),
        }
    }

    #[cfg(not(feature = "interning"))]
    #[test]
    fn readgroup_non_interning_backend_size() {
        // Non-interning backend uses String, so will be larger
        let _g = ReadGroup::grouped("pool_A");
        assert!(std::mem::size_of::<ReadGroup>() > 4);
    }
}
