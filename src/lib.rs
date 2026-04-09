//! Library supporting the [DNAComb] structured sequence read processing tool.
//! The [GitHub] and [Crates.io] page give documentation for the CLI tool with
//! these pages documenting the library,
//!
//! Provides structs and functions for:
//! - Defining DNA molecule forms and libraries (libspec)
//! - Reading and filtering paired sequence reads (parsing, filters)
//! - Extracting the defined regions from reads (counting, containers)
//! - Comparing the observed region sequences to the library (libspec, containers)
//! - Supporting constructs for CLI usage (logging, errors)
//!
//! The main module brings this together to define the main CLI program.
//! The library can equally be used programatically to support other workflows although it
//! is designed primarily to support the CLI tool.
//!
//! [DNAComb]: https://github.com/allydunham/dnacomb
//! [GitHub]: https://github.com/allydunham/dnacomb
//! [Crates.io]: https://crates.io/crates/dnacomb
pub mod combination;
pub mod combinations;
pub mod counting;
pub mod errors;
pub mod filters;
pub mod interning;
pub mod lib_spec;
pub mod library;
pub mod library_combination;
pub mod logging;
pub mod parsing;
pub mod region;
pub mod seqs;
pub mod utils;

// Re-Exports
pub use combinations::ObservedCombinations;
pub use counting::{CountMode, count_reads};
pub use lib_spec::LibrarySpec;
pub use library::SubLibrary;
