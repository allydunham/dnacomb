//! Read Group Data Structure
//!
//! ReadGroup struct with interning and non-interning forms

#[cfg(feature = "interning")]
mod enabled {
    use std::fmt;
    use std::num::NonZeroU32;

    use crate::interning::{
        group_id_from_raw, group_id_from_str, group_id_to_raw, group_id_to_str,
    };

    /// ReadGroup using interned GroupID with niche optimization.
    /// The first 2 bits are reserved for Option<ReadGroup> (via NonZeroU32) and the flag variants
    #[repr(transparent)]
    #[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
    pub struct ReadGroup(NonZeroU32);

    impl ReadGroup {
        const UNGROUPED: u32 = 0b01;
        const UNMATCHED: u32 = 0b10;
        const FLAG_BITS: u32 = 0b11;
        // 4+ reserved for interned GroupIDs

        /// Create an Ungrouped ReadGroup
        #[inline]
        pub fn ungrouped() -> Self {
            ReadGroup(NonZeroU32::new(Self::UNGROUPED).expect("Know this is 1"))
        }

        /// Create an Unmatched ReadGroup
        #[inline]
        pub fn unmatched() -> Self {
            ReadGroup(NonZeroU32::new(Self::UNMATCHED).expect("Know this is 2"))
        }

        /// Create a ReadGroup from a group name string (interns if needed)
        #[inline]
        pub fn grouped(s: &str) -> Self {
            // Send string to the interner, retrieving the ID and interning if necessary
            let group_id = group_id_from_str(s);
            ReadGroup(group_id_to_raw(group_id))
        }

        /// Check if this is Ungrouped
        #[inline]
        pub fn is_ungrouped(&self) -> bool {
            self.0.get() == Self::UNGROUPED
        }

        /// Check if this is Unmatched
        #[inline]
        pub fn is_unmatched(&self) -> bool {
            self.0.get() == Self::UNMATCHED
        }

        /// Check if this is a Match variant
        #[inline]
        pub fn is_match(&self) -> bool {
            self.0.get() > Self::FLAG_BITS
        }
    }

    impl fmt::Display for ReadGroup {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            if self.is_match() {
                write!(f, "{}", group_id_to_str(group_id_from_raw(self.0)))
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

    /// Group status of a read (non-interning)
    #[derive(Debug, Clone, Eq, PartialEq, Hash)]
    pub enum ReadGroup {
        Ungrouped,
        Unmatched,
        Match(String),
    }

    impl ReadGroup {
        /// Create an Unmatched ReadGroup
        pub fn ungrouped() -> Self {
            Self::Ungrouped
        }

        /// Create an Unmatched ReadGroup
        pub fn unmatched() -> Self {
            Self::Unmatched
        }

        /// Create a ReadGroup from a group name string (non-interning)
        pub fn grouped(s: &str) -> Self {
            ReadGroup::Match(s.to_string())
        }

        pub fn is_ungrouped(&self) -> bool {
            matches!(self, ReadGroup::Ungrouped)
        }

        pub fn is_unmatched(&self) -> bool {
            matches!(self, ReadGroup::Unmatched)
        }

        pub fn is_match(&self) -> bool {
            matches!(self, ReadGroup::Match(_))
        }
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

    #[cfg(feature = "interning")]
    #[test]
    fn readgroup_predicates() {
        let ug = ReadGroup::ungrouped();
        let um = ReadGroup::unmatched();
        let m = ReadGroup::grouped("test");

        assert!(ug.is_ungrouped());
        assert!(!ug.is_unmatched());
        assert!(!ug.is_match());

        assert!(!um.is_ungrouped());
        assert!(um.is_unmatched());
        assert!(!um.is_match());

        assert!(!m.is_ungrouped());
        assert!(!m.is_unmatched());
        assert!(m.is_match());
    }
}
