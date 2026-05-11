#!/usr/bin/env Rscript
# Plot benchmark results
library(tidyverse)
library(ggpubr)
library(ggh4x)
dir.create("plots/bench", showWarnings = FALSE, recursive = TRUE)

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

# Captures any number of benchmark TSVs called bench_1_1, bench_2_1, ... for bench_rep_threads
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

# Replicate correlation
rep_cors <- mutate(benchmark, paired = rev != "None") %>%
  select(interning, rep, threads, paired, mode:group, library_size:skip_variants, reads, total_time, extraction_time, region_matching_time) %>%
  pivot_longer(ends_with("_time"), names_to = "part", values_to = "time") %>%
  pivot_wider(names_from = rep, names_prefix = "rep", values_from = time) %>%
  group_by(part) %>%
  group_modify(~as_tibble(cor(select(., starts_with("rep")), use = "pairwise.complete.obs"), rownames = "group1")) %>%
  ungroup() %>%
  pivot_longer(rep1:rep5, names_to = "group2", values_to = "cor")

p_rep_cors <- ggplot(rep_cors, aes(x = group1, y = group2, fill = cor, label = signif(cor, digits = 2))) +
  facet_grid(cols = vars(part)) +
  geom_raster() +
  geom_text(size = 2) +
  coord_fixed() +
  scale_fill_distiller(name = "Correlation", palette = "Reds", limits = c(0, 1), direction = 1) +
  labs(x = "", y = "") +
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size = 9))
ggsave("plots/bench/replicates.png", p_rep_cors, units = "cm", height = 10, width = 15)

# Reps correlate very well so can average metrics
averages <- mutate(benchmark, paired = rev != "None") %>%
  select(rep, interning:threads, paired, lib_spec:no_cache, library_size:skip_variants, reads, total_time, extraction_time, extraction_rate,
         region_matching_time, region_matching_rate, combination_time, combination_rate, summary_time, summary_rate) %>%
  group_by(interning, threads, paired, lib_spec, mode, metric, no_cache, library_size, read_length, skip_variants, reads) %>%
  summarise(across(total_time:summary_rate, sd, .names = "{.col}_sd"),
            across(total_time:summary_rate, mean),
            .groups = "drop") %>%
  mutate(target = str_c(str_match(lib_spec, "config/([a-z_]*)\\.json")[,2], " (", read_length, "bp)"),
         metric = factor(metric, levels = c("exact", "hamming", "bounded-levenshtein", "levenshtein")),
         mode = factor(mode, levels = c("full-read", "inframe", "pattern", "align")))

# Overall time variation
p_proportion <- select(averages, interning:reads, ends_with("_time")) %>%
  pivot_longer(extraction_time:summary_time, names_to = "type", values_to = "time") %>%
  mutate(f = time / total_time,
         type = factor(type, levels = c("extraction_time", "region_matching_time", "combination_time", "summary_time"))) %>%
  {
    ggplot(., aes(x = threads, y = f, fill = type)) +
      facet_nested(cols = vars(interning, no_cache), labeller = labeller(
        no_cache = c(`TRUE` = "Uncached", `FALSE` = "Cached"), interning = c(`TRUE` = "Interning", `FALSE` = "No Interning")
      )) +
      geom_boxplot(position = position_dodge(), outlier.size = 0.1, linewidth = 0.5) +
      scale_fill_brewer(name = "", palette = "Set1", labels = c(
        extraction_time = "Extraction", region_matching_time = "Matching", combination_time = "Combinations", summary_time = "Summary"
      )) +
      labs(x = "Threads", y = "Mean proportion of compute time")
  }
ggsave("plots/bench/fraction.png", p_proportion, units = "cm", height = 10, width = 15)
# Only extraction and matching take an appreciable time, so focus on them downstream.

# Interning
interning_ratios <- select(averages, interning:reads, total_time:region_matching_rate) %>%
  pivot_longer(total_time:region_matching_rate, names_to = "type", values_to = "value") %>%
  pivot_wider(names_from = interning, values_from = value) %>%
  rename(no_interning = `FALSE`, interning = `TRUE`) %>%
  mutate(diff = log2(interning/no_interning),
         target = str_c(str_match(lib_spec, "config/([a-z_]*)\\.json")[,2], " (", read_length, "bp)")) %>%
  pivot_wider(names_from = type, values_from = c(interning, no_interning, diff))

p_interning <- ggplot(interning_ratios, aes(x = no_interning_total_time, y = interning_total_time, colour = mode)) +
  facet_nested(cols = vars(threads), labeller = labeller(threads = ~str_c(., " thread", if_else(. > 1, "s", "")))) +
  geom_point(size = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_fixed() + 
  scale_x_continuous(breaks = c(0, 1, 30, 60, 300, 600, 1800, 3600), labels = time_label, transform = "pseudo_log") +
  scale_y_continuous(breaks = c(0, 1, 30, 60, 300, 600, 1800, 3600), labels = time_label, transform = "pseudo_log") +
  scale_colour_brewer(name = "", palette = "Set1") +
  labs(x = "Non-Interning Time", y = "Interning Time") +
  theme(text = element_text(size = 7))
ggsave("plots/bench/interning.png", p_interning, units = "cm", height = 7, width = 24)

p_interning_ratio <- mutate(interning_ratios, time_cat = case_when(
  interning_total_time < 10 ~ "< 10s",
  interning_total_time < 600 ~ "< 10 mins",
  TRUE ~ "> 10 mins"
) %>% factor(levels = c("< 10s", "< 10 mins", "> 10 mins"))) %>%
  select(threads:target, time_cat, starts_with("diff_")) %>%
  pivot_longer(starts_with("diff_"), names_to = "name", values_to = "ratio", names_prefix = "diff_") %>%
  filter(str_ends(name, "_time")) %>%
  {
    ggplot(., aes(x = name, y = ratio, fill = time_cat)) +
      facet_nested(cols = vars(threads), labeller = labeller(threads = ~str_c(., " thread", if_else(. > 1, "s", "")))) +
      geom_boxplot(outlier.size = 0.1) +
      geom_hline(yintercept = 0) +
      scale_x_discrete(limits = c("total_time", "extraction_time", "region_matching_time"),
                       labels = c("Total", "Extraction", "Matching", "Extraction")) +
      scale_fill_brewer(palette = "Reds", name = "") +
      labs(x = "", y = "log2(Interning / Non-interning)") +
      theme(text = element_text(size = 9),
            legend.position = "bottom")
  }
ggsave("plots/bench/interning_ratio.png", p_interning_ratio, units = "cm", height = 8, width = 18)

# Threading
p_threads <- filter(averages, paired, !skip_variants) %>%
  mutate(reads = factor(case_match(reads, 10000 ~ "10k", 100000 ~ "100k", 1000000 ~ "1M"), levels = c("10k", "100k", "1M"))) %>%
  {
    ggplot(., aes(x = as.integer(threads), y = total_time, colour = as.character(library_size), linetype = interning, shape = no_cache)) +
      facet_nested(rows = vars(mode, metric), cols = vars(target, reads), render_empty = FALSE, scales = "free_y") +
      geom_point() +
      geom_line() +
      scale_colour_brewer(name = "Library size", palette = "Set1") +
      scale_shape_discrete(name = "", labels = c(`TRUE` = "No Caching", `FALSE` = "Caching")) +
      scale_linetype_discrete(name = "", labels = c(`TRUE` = "Interning", `FALSE` = "No Interning")) +
      scale_x_continuous(breaks = c(1, 2, 4, 6)) +
      scale_y_continuous(breaks = c(1, 5, 60, 300, 600, 1800, 3600), labels = time_label, transform = "pseudo_log") +
      labs(x = "Number of Threads", y = "Total Time", caption = "Paired end only")
  }
ggsave("plots/bench/threads.png", p_threads, units = "cm", height = 30, width = 25)

# Extraction mode
p_modes_time <- filter(averages, interning, metric == "exact", threads == 1) %>%
  mutate(mode = factor(mode, levels = c("full-read", "inframe", "pattern", "align")),
         target = str_c(str_match(lib_spec, "config/([a-z_]*)\\.json")[,2], " (", read_length, "bp)")) %>%
  {
    ggplot(., aes(x = reads, y = extraction_time, colour = as.character(library_size), linetype = no_cache, shape = as.character(paired))) +
      facet_grid2(rows = vars(target), cols = vars(mode), render_empty = FALSE,
                  labeller = labeller(no_cache = c(`TRUE` = "Uncached", `FALSE` = "Cached"))) +
      geom_point() +
      geom_line() +
      scale_colour_brewer(name = "Library size", palette = "Set1") +
      scale_shape_discrete(name = "", labels = c(`TRUE` = "Paired-end", `FALSE` = "Single-end")) +
      scale_linetype_discrete(name = "", labels = c(`TRUE` = "No Caching", `FALSE` = "Caching")) +
      scale_x_log10(breaks = c(1e4, 1e5, 1e6), labels = c("10k", "100k", "1M")) +
      scale_y_continuous(breaks = c(0, 1, 30, 60, 300, 600, 1800, 3600), labels = time_label, transform = "pseudo_log") +
      labs(x = "Number of Reads", y = "Time")
  }
ggsave("plots/bench/modes_time.png", p_modes_time, units = "cm", height = 20, width = 30)

p_modes_rate <- filter(averages, interning, metric == "exact", threads == 1) %>%
  mutate(mode = factor(mode, levels = c("full-read", "inframe", "pattern", "align")),
         target = str_c(str_match(lib_spec, "config/([a-z_]*)\\.json")[,2], " (", read_length, "bp)")) %>%
  {
    ggplot(., aes(x = reads, y = extraction_rate, colour = as.character(library_size), linetype = no_cache, shape = as.character(paired))) +
      facet_grid2(rows = vars(target), cols = vars(mode), render_empty = FALSE) +
      geom_point() +
      geom_line() +
      scale_colour_brewer(name = "Library size", palette = "Set1") +
      scale_shape_discrete(name = "", labels = c(`TRUE` = "Paired-end", `FALSE` = "Single-end")) +
      scale_linetype_discrete(name = "", labels = c(`TRUE` = "No Caching", `FALSE` = "Caching")) +
      scale_x_log10(breaks = c(1e4, 1e5, 1e6), labels = c("10k", "100k", "1M")) +
      scale_y_continuous(breaks = c(1, 1e2, 1e4, 1e6), transform = "pseudo_log", limits = c(0, 2e6)) +
      labs(x = "Number of Reads", y = "Processing rate (reads/s)")
  }
ggsave("plots/bench/modes_rate.png", p_modes_rate, units = "cm", height = 20, width = 30)

# Comparison metric
p_metrics_time <- filter(averages, interning, threads == 1, mode == "align", !no_cache) %>%
  mutate(target = str_c(str_match(lib_spec, "config/([a-z_]*)\\.json")[,2], " (", read_length, "bp)"),
         metric = factor(metric, levels = c("exact", "hamming", "bounded-levenshtein", "levenshtein"))) %>%
  {
    ggplot(., aes(x = library_size, y = region_matching_time, colour = as.character(reads), linetype = skip_variants,
                  shape = as.character(paired))) +
      facet_grid2(rows = vars(target), cols = vars(metric), render_empty = FALSE) +
      geom_point() +
      geom_line() +
      scale_colour_brewer(name = "Number of reads", palette = "Set1") +
      scale_shape_discrete(name = "", labels = c(`TRUE` = "Paired-end", `FALSE` = "Single-end")) +
      scale_linetype_discrete(name = "", labels = c(`TRUE` = "No HGVS", `FALSE` = "Inc. HGVS")) +
      scale_x_log10(breaks = c(100, 1000, 10000), labels = c("100", "1k", "10k")) +
      scale_y_continuous(breaks = c(0, 1, 30, 60, 300, 600, 1800, 3600), labels = time_label, transform = "pseudo_log", limits = c(0, 1.5e4)) +
      labs(x = "Library Size", y = "Time")
  }
ggsave("plots/bench/metrics_time.png", p_metrics_time, units = "cm", height = 20, width = 30)

p_metrics_rate <- filter(averages, interning, threads == 1, mode == "align", !no_cache) %>%
  mutate(target = str_c(str_match(lib_spec, "config/([a-z_]*)\\.json")[,2], " (", read_length, "bp)"),
         metric = factor(metric, levels = c("exact", "hamming", "bounded-levenshtein", "levenshtein"))) %>%
  {
    ggplot(., aes(x = library_size, y = region_matching_rate, colour = as.character(reads), linetype = skip_variants,
                  shape = as.character(paired))) +
      facet_grid2(rows = vars(target), cols = vars(metric), render_empty = FALSE) +
      geom_point() +
      geom_line() +
      scale_colour_brewer(name = "Number of reads", palette = "Set1") +
      scale_shape_discrete(name = "", labels = c(`TRUE` = "Paired-end", `FALSE` = "Single-end")) +
      scale_linetype_discrete(name = "", labels = c(`TRUE` = "No HGVS", `FALSE` = "Inc. HGVS")) +
      scale_x_log10(breaks = c(100, 1000, 10000), labels = c("100", "1k", "10k")) +
      scale_y_continuous(breaks = c(0, 1e2, 1e4, 1e6), transform = "pseudo_log") +
      labs(x = "Library Size", y = "Processing rate (regions/s)")
  }
ggsave("plots/bench/metrics_rate.png", p_metrics_rate, units = "cm", height = 20, width = 30)
