#!/usr/bin/env Rscript
# Benchmark LUCA on DNAComb mutation_test.py and benchmark.py data
# Assumes those scripts have been run and the data is still available
library(tidyverse)
library(jsonlite)
library(fs)
library(glue)

manifest_path <- "data/accuracy_profile/accuracy_profile_manifest.tsv"
out_dir <- "data/luca_profile"
# luca_bin <- "luca"
luca_bin <- "/Users/ad44/software/LUCA/.venv/bin/luca"
samtools_bin <- "samtools"

libspec_paths <- c(
  grna = "config/grna_no_id.json",
  grna_sensor = "config/grna_sensor.json",
  pegrna = "config/pegrna.json"
)

library_paths <- c(
  grna = "config/grna_no_id.tsv",
  grna_sensor = "config/grna_sensor.tsv",
  pegrna = "config/pegrna.tsv"
)

dir_create(out_dir)

safe_unlink <- function(paths) {
  paths %>%
    discard(is.na) %>%
    discard(~.x == "") %>%
    walk(~unlink(.x, recursive = TRUE, force = TRUE))
}

run_checked <- function(cmd, args) {
  out <- tempfile()
  err <- tempfile()
  
  on.exit(unlink(c(out, err)), add = TRUE)
  
  start <- proc.time()[["elapsed"]]
  
  status <- system2(
    cmd,
    args,
    stdout = out,
    stderr = err
  )
  
  elapsed <- proc.time()[["elapsed"]] - start
  
  if (!identical(status, 0L)) {
    stop(glue(
      "Command failed with exit status {status}\n\n",
      "stdout:\n{read_file(out)}\n\n",
      "stderr:\n{read_file(err)}"
    ), call. = FALSE)
  }
  
  elapsed
}

write_luca_config <- function(design, spec_path, library_path, root) {
  spec <- read_json(spec_path, simplifyVector = FALSE)
  
  variable_ids <- spec$regions %>%
    keep(~.x$seq_type == "Library") %>%
    map_chr("id")
  
  lib_tbl <- read_tsv(library_path, show_col_types = FALSE) %>%
    mutate(.dnacomb_id = if ("_id" %in% names(.)) as.character(.data$`_id`) else as.character(row_number()))
  
  library_dir <- path(root, design, "libraries")
  dir_create(library_dir)
  
  seq_maps <- variable_ids %>%
    set_names() %>%
    map(~{
      seq_tbl <- lib_tbl %>%
        distinct(sequence = .data[[.x]]) %>%
        arrange(sequence) %>%
        mutate(id = row_number() - 1L) %>%
        select(id, sequence)
      
      write_tsv(seq_tbl, path(library_dir, .x))
      
      seq_tbl %>%
        mutate(region = .x)
    }) %>%
    bind_rows()
  
  fixed_libraries <- spec$regions %>%
    keep(~.x$seq_type == "Fixed") %>%
    map(~list(
      id = .x$id,
      values = list(.x$seq),
      reverse_on = "reverse_group"
    ))
  
  variable_libraries <- variable_ids %>%
    map(~list(
      id = .x,
      reverse_on = "reverse_group"
    ))
  
  filter_tbl <- lib_tbl %>%
    select(all_of(variable_ids))
  
  filter_path <- path(library_dir, "valid_combinations")
  write_tsv(filter_tbl, filter_path)
  
  regions <- spec$regions %>%
    map(~list(
      id = .x$id,
      libraries = list(.x$id),
      max_offset = 0
    ))
  
  experiment <- list(
    sequencing_type = "single_end",
    libraries = c(fixed_libraries, variable_libraries),
    read_templates = list(list(
      id = "default_template",
      anchor = "left",
      regions = regions
    )),
    read_group_templates = list(default = list("default_template")),
    combinations = list(list(
      id = "assignment",
      regions = variable_ids %>%
        map(~list(id = .x, filter = TRUE)),
      filters = list("valid_combinations")
    )),
    default_options = list(
      full_match_info = TRUE,
      count_mm_reads = FALSE,
      sort_mm_read_counts = FALSE,
      compress_mm_read_counts = FALSE,
      out_mm_match_info = FALSE,
      out_mm_reads = FALSE
    )
  )
  
  experiment_path <- path(root, design, "experiment.json")
  write_json(experiment, experiment_path, auto_unbox = TRUE, pretty = TRUE)
  
  tibble(
    design = design,
    experiment_json = experiment_path,
    library_dir = library_dir,
    variable_ids = list(variable_ids),
    seq_maps = list(seq_maps)
  )
}

fastq_to_unmapped_bam <- function(fq, bam) {
  run_checked(
    samtools_bin,
    c("import", "-o", shQuote(bam), shQuote(fq))
  )
  
  invisible(bam)
}

# Test accuracy profile results
analyse_luca_run <- function(path, true_counts, variable_ids, seq_maps) {
  # Load data from LUCA and true results
  exp_tbl <- read_tsv(true_counts, show_col_types = FALSE) %>%
    mutate(across(all_of(variable_ids), as.character),
           count = as.numeric(count))
  
  luca_tbl <- read_tsv(path, col_names = c(variable_ids, "count"), show_col_types = FALSE) %>%
    mutate(across(all_of(variable_ids), as.character),
           count = as.numeric(count))
  
  # Calculate fraction of reads assigned to source 
  all_regions <- full_join(
    exp_tbl %>%
      select(ends_with("_nearest"), count) %>%
      rename_with(.cols = ends_with("_nearest"), .fn = ~str_remove(., "_nearest")) %>%
      count(across(all_of(variable_ids)), wt = count, name = "expected"),
    luca_tbl %>%
      select(all_of(variable_ids), observed = count) %>%
      count(across(all_of(variable_ids)), wt = observed, name = "observed"),
    by = variable_ids
  ) %>%
    replace_na(list(expected = 0, observed = 0))
  
  total_reads <- sum(exp_tbl$count)
  assigned_reads <- sum(luca_tbl$count)
  unassigned_reads <- total_reads - assigned_reads
  
  correct_all <- sum(pmin(all_regions$expected, all_regions$observed))
  
  tibble(
    total_reads = total_reads,
    assigned_reads = assigned_reads,
    unassigned_reads = unassigned_reads,
    assignment_accuracy = correct_all / total_reads
  )
}

run_accuracy_dataset <- function(row, luca_configs) {
  # Prepare config
  cfg <- filter(luca_configs, design == row$design)
  
  run_id <- glue("{row$design}.level{row$level}.rep{row$replicate}")
  run_root <- path(out_dir, "tmp", run_id)
  result_root <- path(out_dir, "runs", run_id)
  
  dir_create(run_root)
  dir_create(result_root)
  
  on.exit(safe_unlink(c(run_root, result_root)), add = TRUE)
  
  # Generate BAM
  bam_path <- path(run_root, "reads.bam")
  fastq_to_unmapped_bam(row$input_f, bam_path)
  
  # Run LUCA
  elapsed <- run_checked(
    luca_bin,
    c(
      "count",
      "-o", shQuote(result_root),
      "-l", shQuote(cfg$library_dir[[1]]),
      "-s", "luca_test",
      shQuote(cfg$experiment_json[[1]]),
      shQuote(bam_path)
    )
  )
  
  combo_file <- dir_ls(result_root, regexp = "combination\\.[0-9]+\\.counts\\.tsv$")
  
  if (length(combo_file) != 1L) {
    stop(glue("Expected one LUCA combination counts file in {result_root}, found {length(combo_file)}"),
         call. = FALSE)
  }
  
  metrics <- analyse_luca_run(
    path = combo_file,
    true_counts = row$true_counts,
    variable_ids = cfg$variable_ids[[1]],
    seq_maps = cfg$seq_maps[[1]]
  )
  
  tibble(
    design = row$design,
    level = row$level,
    replicate = row$replicate,
    input_f = row$input_f,
    true_counts = row$true_counts,
    time = elapsed,
    
  ) %>%
    bind_cols(metrics) %>%
    mutate(reads_per_second = total_reads / time)
}

tmp_config_root <- path(out_dir, "tmp_config")
safe_unlink(tmp_config_root)
dir_create(tmp_config_root)

luca_configs <- imap_dfr(
  libspec_paths,
  ~write_luca_config(
    design = .y,
    spec_path = .x,
    library_path = library_paths[[.y]],
    root = tmp_config_root
  )
)

accuracy_manifest <- read_tsv(manifest_path, show_col_types = FALSE) %>%
  filter(end == "single") %>%
  distinct(design, level, replicate, sub_rate, indel_rate, recombination_rate,
           contamination_rate, mismatch_rate, input_f, true_counts) %>%
  arrange(design, replicate, level)

accuracy_results <- accuracy_manifest %>%
  group_split(row_number()) %>%
  imap(~{
    row <- .x[1, ]
    if (.y %% 10 == 1)
    message(glue("[{.y}/{nrow(manifest)}] {row$design} level {row$level} rep {row$replicate}"))
    run_accuracy_dataset(row, luca_configs)
  })  %>%
  bind_rows() %>%
  left_join(
    select(manifest, design, level, replicate, sub_rate, indel_rate,
           recombination_rate, contamination_rate, mismatch_rate),
    by = c("design", "level", "replicate")
  ) %>%
  relocate(sub_rate, indel_rate, recombination_rate, contamination_rate, mismatch_rate,
           .after = replicate)

write_tsv(accuracy_results, path(out_dir, "luca_accuracy_profile.tsv"))
safe_unlink(tmp_config_root)

# Run benchmark for speed
bench_reads <- 1000000

bench_inputs <- tibble(
  fq = dir("data/benchmark", pattern = "^bench[0-9]+_(grna|grna_sensor|pegrna)_1000\\.fq$", full.names = TRUE)
) %>%
  mutate(
    file = basename(fq),
    design = str_match(file, "^bench[0-9]+_(grna|grna_sensor|pegrna)_1000\\.fq$")[,2],
    bench = str_match(file, "^(bench[0-9]+)_")[,2],
    library = str_replace(fq, "\\.fq$", ".tsv")
  ) %>%
  arrange(design, bench)

run_benchmark <- function(row) {
  run_id <- glue("{row$bench}.{row$design}.1000")
  run_root <- path(out_dir, "tmp_benchmark", run_id)
  result_root <- path(out_dir, "benchmark_runs", run_id)
  
  dir_create(run_root)
  dir_create(result_root)
  
  on.exit(safe_unlink(c(run_root, result_root)), add = TRUE)
  
  cfg <- write_luca_config(
    design = row$design,
    spec_path = libspec_paths[[row$design]],
    library_path = row$library,
    root = run_root
  )
  
  bam_path <- path(run_root, "reads.bam")
  fastq_to_unmapped_bam(row$fq, bam_path)
  
  elapsed <- run_checked(
    luca_bin,
    c(
      "count",
      "-o", shQuote(result_root),
      "-l", shQuote(cfg$library_dir[[1]]),
      "-s", "luca_benchmark",
      shQuote(cfg$experiment_json[[1]]),
      shQuote(bam_path)
    )
  )
  
  message(glue("{run_id} completed in {elapsed}s"))
  
  tibble(
    bench = row$bench,
    design = row$design,
    library_size = 1000,
    input_f = row$fq,
    library = row$library,
    reads = bench_reads,
    time = elapsed,
    reads_per_second = bench_reads / elapsed
  )
}

benchmark_results <- bench_inputs %>%
  slice_head(n = 2) %>%
  group_split(row_number()) %>%
  map(~run_benchmark(.x[1, ])) %>%
  bind_rows()

write_tsv(benchmark_results, path(out_dir, "luca_benchmark_runs.tsv"))
