# Changelog

This file contains a record of the major changes between versions.

## v0.6.0 Interning & combinatorial libraries

* Added interning of Sequences, Groups and Region IDs to reduce duplication, reduce memory usage and speed up performance

## v0.5.0 Deeper tracking, refined config and performance improvements

* Removed superfluous length field from fixed length regions in LibSpec - now the sequence is the core data.
* Cleaned up filters internally, making them more robust to changes and easier to add new ones, plus added a new filter
* Fixed bug preventing negative number input to alignment params
* Fixed bug in multi-match library lookup (they previously removed all possibilities rather than considering any that either match includes)
* Fixed bug in benchmark timing brake down (region matching main time not being tracked)

## v0.4.0 Tracking filtered reads

* Added a new output table detailing which reads are filtered and for what reason
* More comprehensive testing

## v0.3.0 Library IDs

* Added library IDs to the LibSpec library format and supported them through the workflow
* Logging handles threading better
* Comprehensive benchmark

## v0.2.0 Multithreading and redesigned pattern matching

* Added multithreading for region extraction and library mapping
* Redesigned the pattern matching algorithm to be more robust under mutations
* Added more comprehensive testing
* Added more comprehensive benchmarking

## v0.1.0 Initial release

Initial release of the package, containing basic functionality for the whole workflow.
