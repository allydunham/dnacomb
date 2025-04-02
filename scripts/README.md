# Scripts

A collection of scripts for testing and supporting the main tool:

* `utils.py` - Python utility to build simple shell wrappers used in testing/benchmarks
* `generate_test_data.py` - Scripts and a CLI for generating simulated sequence reads from a LibSpec JSON files, with a range of options for error and mutation models
* `benchmark.py` - Run a benchmark suite, generating data and measure performance across a range of parameter values
* `integration_test.py` - Run a suite of end-to-end tests, checking the tool gives the expected results for various inputs, again using generated data.
* `plots.R` - Plot the results of the above two scripts into some illustrative figures

The python scripts require `numpy` & `biopython` while the R script requires `tidyverse`, `ggpubr` & [plotlistr]().
These scripts are mainly useful for development and testing, and are not required for normal usage.
