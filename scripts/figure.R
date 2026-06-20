#!/usr/bin/env Rscript
# Generate paper figure panels
library(tidyverse)
library(ggpubr)
library(ggh4x)
library(figpatch)
library(patchwork)
dir.create("plots/", showWarnings = FALSE, recursive = TRUE)

theme_set(theme_pubclean() + theme(legend.position = 'right',
                                   plot.title = element_text(hjust = 0.5),
                                   plot.subtitle = element_text(hjust = 0.5),
                                   strip.background = element_blank(),
                                   legend.key = element_blank()))

na_or_zero <- function(x) {
  if_else(is.na(x), 0, x)
}

time_label <- function(x) {
  x[is.na(x)] <- 0
  out <- str_c(x, "s")
  out[x >= 60] <- str_c(round(x[x >= 60]/60, 1), "m")
  out[x >= 3600] <- str_c(round(x[x >= 3600]/3600, 1), "h")
  return(out)
}

# Accuracy profile (assumes scripts/mutation_test.py has been run)
analyse_run <- function(expected_counts, observed_counts) {
  exp_tbl <- read_tsv(expected_counts)
  obs_tbl <- read_tsv(observed_counts)
  
  reg_names <- select(exp_tbl, -group, -ends_with("_nearest"), 
                      -combination_status, -combinations_in_library, -combination_indexes, -count) %>%
    colnames()
  
  total_reads <- sum(obs_tbl$count)
  
  # Matches found - proportion of proper reads with the real sequence extracted
  assignment <- left_join(
    filter(exp_tbl, combination_status %in% c("match", "recombination")) %>%
      select(one_of(reg_names), exp = count),
    select(obs_tbl, one_of(reg_names), obs = count)
  ) %>%
    replace_na(list(exp = 0, obs = 0))
  
  total_good <- sum(assignment$exp)
  match_accuracy <- (total_good - sum(abs(assignment$obs - assignment$exp))) / total_good
  
  # Region accuracy - proportion of all reads where the expected regions extracted
  regions <- full_join(
      select(exp_tbl, one_of(reg_names), exp = count),
      select(obs_tbl, one_of(reg_names), combination_status, obs = count) %>%
        mutate(across(one_of(reg_names), \(x) if_else(combination_status == "nonmatch", NA, x))) %>% 
        count(across(one_of(reg_names)), wt = obs, name = "obs")
    ) %>%
    replace_na(list(exp = 0, obs = 0))
  
  region_accuracy <- (total_reads - sum(abs(regions$obs - regions$exp)) / 2) / total_reads # Each mismatched read is counted twice
  
  # Match accuracy - of called matches, how many are the correct library member
  matches <- left_join(
      filter(exp_tbl, combination_status %in% c("match", "recombination")) %>%
        select(one_of(reg_names), ends_with("_nearest"), combination_status, exp = count),
      filter(obs_tbl, combination_status %in% c("match", "recombination")) %>%
        select(one_of(reg_names), ends_with("_nearest"), combination_status, obs = count)
    ) %>%
    replace_na(list(exp = 0, obs = 0))
  
  total_matches <- sum(matches$exp)
  library_accuracy <- (total_matches - sum(abs(matches$obs - matches$exp))) / total_matches
  
  # Output
  tibble(
    match_accuracy = match_accuracy, 
    region_accuracy = region_accuracy,
    library_accuracy = library_accuracy
  )
}
quiet_analysis <- purrr::quietly(analyse_run)

accuracy <- read_tsv("data/accuracy_profile/accuracy_profile_manifest.tsv") %>%
  {bind_cols(., bind_rows(map2(.$true_counts, .$counts_tsv, ~quiet_analysis(.x, .y)$result, .progress = TRUE)))}

# Check paired end and single end match - exact agreement
ends_correlation <- select(accuracy, design, level, replicate, mode, metric, end, match_accuracy, region_accuracy, library_accuracy) %>%
  pivot_longer(ends_with("accuracy")) %>%
  pivot_wider(names_from = end, values_from = value) %>%
  group_by(design, mode, metric, name) %>%
  group_modify(~broom::tidy(cor.test(.$single, .$paired)))

perturbation_labels <- c(
  "0 = Perfect",
  "1 = 0.5% Mutation Rate & 1% Mismatches",
  "2 = 1% Mutation Rate, 0.1% Indel Rate, 1% Mismatches & 1% Contaminants",
  "3 = 1% Mutation Rate, 0.5% Indel Rate, 1% Recombinations, 1% Mismatches & 1% Contaminants",
  "4 = 1% Mutation Rate, 1% Indel Rate, 2% Recombinations, 2% Mismatches & 2% Contaminants",
  "5 = 5% Mutation Rate, 1% Indel Rate, 5% Recombinations, 5% Mismatches & 5% Contaminants",
  "6 = 5% Mutation Rate, 2% Indel Rate, 10% Recombinations, 10% Mismatches & 5% Contaminants"
)

metric_labels <- c(
  region = "Regions (All)", match = "Regions (Expected Matches only)", 
  exact = "Matches (Exact)", hamming = "Matches (Hamming)", `bounded-levenshtein` = "Matches (Bounded-Levenshtein)"
)

p_accuracy <- filter(accuracy, end == "single") %>%
  select(design, level, replicate, mode, metric, match_accuracy, region_accuracy, library_accuracy) %>%
  {bind_rows(
    select(., design:metric, accuracy = library_accuracy), 
    filter(., metric == "exact") %>% select(design:metric, match_accuracy) %>% mutate(metric = "match") %>% rename(accuracy = match_accuracy),
    filter(., metric == "exact") %>% select(design:metric, region_accuracy) %>% mutate(metric = "region") %>% rename(accuracy = region_accuracy),
  )} %>%
  pivot_longer(ends_with("accuracy")) %>%
  group_by(design, level, mode, metric, name) %>%
  mutate(mean = mean(value),
         sd = sd(value)) %>%
  ungroup() %>%
  mutate(mode = factor(mode, levels = c("align", "pattern", "inframe")),
         design = factor(design, levels = c("grna_sensor", "pegrna"))) %>%
  {
    ggplot(., aes(x = level, colour = metric, linetype = design)) +
      facet_nested(cols = vars(mode), render_empty = FALSE, labeller = labeller(
        mode = c(align = "Alignment", pattern = "Pattern Matching", inframe = "Inframe"),
      )) +
      geom_point(aes(y = value), shape = 20) +
      geom_line(aes(y = mean)) +
      geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
      scale_y_continuous(name = "Fraction of Reads") +
      scale_x_continuous(name = "Perturbation Level", breaks = 0:6) +
      scale_colour_brewer(limits = names(metric_labels), labels = metric_labels, name = "", palette = "Set1", direction = -1) +
      scale_linetype_discrete(name = "", labels = c(grna_sensor = "gRNA + Sensor", pegrna = "pegRNA")) +
      theme(text = element_text(size = 12),
            legend.position = "bottom")
  }

# Benchmark
bench_cols <- c(
  "name", "fwd", "rev", "lib_spec", "mode", "metric", "no_cache", "sort", "group", "library", "library_size", 
  "read_length", "skip_variants", "additional_args", "total_time", "reads", "extraction_time",
  "extraction_rate", "unique_regions", "region_matching_time", "region_matching_rate", 
  "unique_combinations", "combination_time", "combination_rate", "summary_size",
  "summary_time", "summary_rate"
)
benchmark <- dir("data/benchmark", pattern = "bench_.*.tsv", full.names = TRUE) %>%
  set_names() %>%
  map(read_tsv, col_names = bench_cols, skip = 1) %>%
  bind_rows(.id = "rep") %>%
  separate_wider_regex(rep, c("data/benchmark/bench_", interning = "(?:no_intern_)?", rep = "[0-9]*", "_", threads = "[0-9]*", ".tsv")) %>%
  mutate(reads = as.integer(str_match(name, "reads:(10*)")[,2]), # Previously reads was mistakenly unique read count, this corrects for this
         extraction_rate = reads / extraction_time,
         interning = interning != "no_intern_") %>%
  drop_na(extraction_time)

p_mode_benchmark <- filter(benchmark, interning, sort, threads == 1) %>%
  select(mode, no_cache, reads, read_length, extraction_rate, rep) %>%
  mutate(mode = factor(mode, levels = c("full-read", "inframe", "pattern", "align"))) %>%
  {
    ggplot(., aes(x = as.factor(read_length), y = extraction_rate, fill = mode, linetype = no_cache)) +
      geom_boxplot(outlier.shape = 20, outlier.size = 0.5, linewidth = 0.5) +
      scale_fill_brewer(name = "", palette = "Reds") +
      scale_linetype_discrete(name = "", labels = c(`TRUE` = "Uncached", `FALSE` = "Cached")) +
      scale_y_continuous(name = "Reads/s", transform = "log10") +
      scale_x_discrete(name = "", labels = c("156bp\n(gRNA)", "244bp\n(gRNA + Sensor)", "313bp\n(pegRNA)")) + 
      theme(text = element_text(size = 12),
            legend.position = "bottom")
  }

p_metric_benchmark <- filter(benchmark, interning, sort, threads == 1, mode == "align") %>%
  select(metric, no_cache, reads, read_length, library_size, region_matching_rate, rep) %>%
  mutate(metric = factor(metric, levels = c("exact", "hamming", "bounded-levenshtein", "levenshtein"))) %>% 
  {
    ggplot(., aes(x = metric, y = region_matching_rate, fill = as.factor(library_size))) +
      geom_boxplot(outlier.shape = 20, outlier.size = 0.5, linewidth = 0.5) +
      scale_fill_brewer(name = "Library Size", palette = "Blues") +
      scale_y_continuous(name = "Regions/s", transform = "log10") +
      scale_x_discrete(name = "") + 
      theme(text = element_text(size = 12),
            legend.position = "bottom")
  }

# Assemble overall figure
dnacomb_schematic <- fig("plots/schematic.png", b_margin = margin())
pipeline_schematic <- fig("plots/pipeline.png", b_margin = margin())

figure <- dnacomb_schematic + pipeline_schematic + p_accuracy + p_mode_benchmark + p_metric_benchmark +
  plot_layout(design = "11\n22\n33\n45", heights = c(8, 8, 3, 3), widths = c(1, 1)) +
  plot_annotation(tag_levels = "A")
ggsave("plots/figure.pdf", figure, units = "cm", height = 16 * 1.8, width = 19 * 1.8)
ggsave("plots/figure.png", figure, units = "cm", height = 16 * 1.8, width = 19 * 1.8)

