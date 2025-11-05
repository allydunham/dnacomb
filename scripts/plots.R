#!/usr/bin/env Rscript
# Plot benchmark and test results
library(tidyverse)
library(ggpubr)
library(ggh4x)
dir.create("plots", showWarnings = FALSE)

theme_set(theme_pubclean() + theme(legend.position = 'right',
                                   plot.title = element_text(hjust = 0.5),
                                   plot.subtitle = element_text(hjust = 0.5),
                                   strip.background = element_blank(),
                                   legend.key = element_blank()))

# Correctness tests
test_pairs <- tibble(
  observed = str_remove(dir("data/tests/", pattern = "*\\.counts.tsv"), ".counts.tsv")
) %>%
  separate_wider_delim(observed, delim = ":", names = c("test", "mode", "distance", "expected", "end", "threads"),
                       cols_remove = FALSE, too_few = "align_start")

test_pair <- function(observed, expected) {
  true_counts <- read_tsv(str_c("data/tests/", expected, ".true_counts.tsv"))
  obs_counts <- read_tsv(str_c("data/tests/", observed, ".counts.tsv"))
  
  reg_names <- select(true_counts, -group, -ends_with("_nearest"), -combination_status, -combinations_in_library, -combination_indexes, -count) %>%
    colnames()
  
  true_exact <- apply(select(true_counts, one_of(reg_names)) == select(true_counts, ends_with("_nearest")), 1, all)
  
  # Test region extraction
  all_obs <- full_join(
    obs_counts,
    select(true_counts, -ends_with("_nearest"), -combination_status, -combinations_in_library, -combination_indexes) %>% rename(true_count = count)
  ) %>%
    mutate(observed = observed, expected = expected, id = str_c("i", 1:n()), split = "all_regions") %>%
    select(observed, expected, split, id, combination_status, count, true_count)
  
  lib_file <- str_c("data/tests/", observed, ".library_counts.tsv")
  if (file.exists(lib_file)) {
    all_obs <- full_join(
      read_tsv(lib_file),
      select(true_counts, group, ends_with("_nearest"), true_count = count) %>%
        rename_with(~str_remove(., "_nearest"), ends_with("_nearest")) %>%
        count(across(c(-true_count)), wt = true_count, name = "true_count")
    ) %>%
      mutate(observed = observed, expected = expected, id = str_c("i", 1:n()), split = "library") %>%
      select(observed, expected, split, id, combination_status, count, true_count) %>%
      bind_rows(all_obs, .)
  }
  
  summary_file <- str_c("data/tests/", observed, ".summary.tsv")
  if (file.exists(summary_file)) {
    all_obs <- full_join(
      read_tsv(summary_file),
      select(true_counts, combination_status, true_count = count) %>%
        mutate(combination_status = if_else(
          combination_status %in% c("match", "recombination"),
          str_c(if_else(true_exact, "exact", "nearest"), "_", combination_status),
          combination_status
        )) %>%
        count(combination_status, wt = true_count, name = "true_count"),
      by = join_by(metric == combination_status)
    ) %>%
      filter(!metric == "total") %>%
      mutate(observed = observed, expected = expected, id = metric, split = "summary") %>%
      select(observed, expected, split, id, count, true_count) %>%
      bind_rows(all_obs, .)
  }
  
  all_obs
}
quiet_test <- purrr::quietly(test_pair)

test_counts <- filter(test_pairs, threads == 1) %>%
  {map2(.$observed, .$expected, ~quiet_test(.x, .y)$result, .progress = TRUE)} %>%
  bind_rows() %>%
  replace_na(replace = list(count = 0, true_count = 0)) %>%
  select(-expected) %>%
  separate_wider_delim(observed, delim = ":", names = c("test", "mode", "distance", "library", "end", "threads"), cols_remove = TRUE) %>%
  mutate(threads = as.integer(threads))

category_colours <- c(
  "match" = "green", "exact_match" = "green", "nearest_match" = "darkgreen",
  "mismatch" = "orange", "nonmatch" = "red", 
  "recombination" = "blue", "exact_recombination" = "blue", "nearest_recombination" = "darkblue",
  "low_mean_quality" = "brown", "bad_alignment" = "grey", "multimatch" = "purple"
)

p_all_scatter <- filter(test_counts, split == "all_regions", threads == 1) %>%
  ggplot(., aes(x = true_count, y = count, colour = combination_status)) +
  facet_nested(rows = vars(mode, distance), cols = vars(library, end), render_empty = FALSE, solo_line = FALSE, nest_line = element_line(colour = "grey")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(shape = 20) +
  coord_fixed() +
  scale_x_continuous(transform = "pseudo_log") +
  scale_y_continuous(transform = "pseudo_log") +
  scale_colour_manual(values = category_colours) +
  labs(x = "Simulated Count", y = "Observed Count") +
  theme(legend.position = "bottom")
ggsave("plots/test/all_counts_scatter.png", p_all_scatter, units = "cm", height = 50, width = 50)

p_library_scatter <- filter(test_counts, split == "library", threads == 1) %>%
  ggplot(., aes(x = true_count, y = count, colour = combination_status)) +
  facet_nested(rows = vars(mode, distance), cols = vars(library, end), render_empty = FALSE, solo_line = FALSE, nest_line = element_line(colour = "grey")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(shape = 20) +
  coord_fixed() +
  scale_x_continuous(transform = "pseudo_log") +
  scale_y_continuous(transform = "pseudo_log") +
  scale_colour_manual(values = category_colours) +
  labs(x = "Simulated Count", y = "Observed Count") +
  theme(legend.position = "bottom")
ggsave("plots/test/library_scatter.png", p_library_scatter, units = "cm", height = 50, width = 50)

p_summary_bars <- filter(test_counts, split == "summary", threads == 1) %>%
  filter(!id == "uncompared") %>%
  pivot_longer(c(count, true_count), names_to = "group", values_to = "count") %>%
  mutate(group = c(count = "Observed", true_count = "Expected")[group],
         id = factor(id, levels = names(category_colours))) %>%
  ggplot(aes(x = group, y = count, fill = id)) +
  facet_nested(rows = vars(library), cols = vars(mode, distance, end), solo_line = FALSE, nest_line = element_line(colour = "grey")) +
  geom_col(position = position_stack(), width = 0.7) +
  scale_fill_manual(name = "", values = category_colours) +
  labs(x = "", y = "Count") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("plots/test/summary_bars.png", p_summary_bars, units = "cm", height = 40, width = 40)

count_cors <- filter(test_counts, split == "all_regions", threads == 1) %>%
  group_by(mode, distance, library, end) %>%
  group_modify(~broom::tidy(cor.test(.$count, .$true_count))) %>%
  ungroup()

p_test_cors <- ggplot(count_cors, aes(y = distance, x = estimate, xmin = conf.low, xmax = conf.high)) +
  facet_nested(rows = vars(library, mode), cols = vars(end), switch = "y", solo_line = FALSE, nest_line = element_line(colour = "grey")) +
  geom_col(fill = "#377eb8", width = 0.6) +
  geom_errorbarh(height = 0.3) +
  labs(x = "Pearson's r", y = "") +
  theme(panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.placement = "outside")
ggsave("plots/test/observed_expected_correlation.png", p_test_cors, units = "cm", height = 30, width = 16)

# Look at count correlation across thread counts
threaded_counts <- filter(test_pairs, mode == "align", distance == "bounded-levenshtein", expected == "mutant_pegrna", end == "single") %>%
  mutate(path = str_c("data/tests/", observed, ".counts.tsv")) %>%
  select(threads, path) %>%
  mutate(counts = map(path, read_tsv)) %>%
  unnest(counts) %>%
  select(-path) %>%
  pivot_wider(names_from = threads, names_prefix = "threads", values_from = count, values_fill = 0) %>%
  mutate(counts_equal = apply(across(starts_with("threads")), 1, n_distinct) == 1)

p_thread_cor <- select(threaded_counts, starts_with("threads")) %>%
  as.matrix() %>%
  cor() %>%
  as_tibble(rownames = "x") %>%
  pivot_longer(-x, names_to = "y", values_to = "cor") %>%
  mutate(x = str_remove(x, "threads"),
         y = str_remove(y, "threads")) %>%
  ggplot(aes(x = x, y = y, fill = cor, label = signif(cor))) +
  geom_tile() +
  geom_text(colour = "white") +
  labs(x = "Threads", y = "Threads") +
  scale_fill_distiller(name = "Count\nR", palette = "RdBu", direction = 1, limits = c(-1, 1))
ggsave("plots/test/thread_count_correlation.png", p_thread_cor, units = "cm", height = 10, width = 10)

# Benchmarks
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
  "name", "fwd", "rev", "lib_spec", "mode", "metric",
  "no_cache", "sort", "group", "library_counts",
  "library_size", "read_length", "additional_args",
  "total_time", "reads", "region_time", "region_rate",
  "unique_regions", "library_time", "library_rate",
  "summary_size", "summary_time", "summary_rate"
)
benchmark <- dir("data/benchmark", pattern = "bench_[0-9]*_[0-9]*.tsv", full.names = TRUE) %>%
  set_names() %>%
  map(read_tsv, col_names = bench_cols, skip = 1) %>%
  bind_rows(.id = "rep") %>%
  extract(rep, c("rep", "threads"), "data/benchmark/bench_([0-9]*)_([0-9]*)", convert = TRUE)

p_modes_time <- filter(benchmark, threads == 1) %>%
  mutate(mode = if_else(mode == "align", str_c(mode, if_else(no_cache, " (uncached)", " (cached)")), mode),
         mode = factor(mode, levels = c("full-read", "inframe", "pattern", "align (cached)", "align (uncached)")),
         target = str_c(str_match(name, "lib:([a-z]*)")[,2], " (", read_length, "bp)"),
         ends = if_else(rev == "None", "Single end", "Paired end")) %>%
  {
    ggplot(., aes(x = reads, y = region_time, colour = as.character(library_size), linetype = ends, shape = as.character(rep))) +
      facet_grid2(rows = vars(target), cols = vars(mode), render_empty = FALSE) +
      geom_point() +
      geom_line() +
      scale_colour_brewer(name = "Library size", palette = "Set1") +
      scale_shape_discrete(name = "Rep") +
      scale_linetype_discrete(name = "") +
      scale_x_log10(breaks = c(1e5, 1e6, 1e7), labels = c("100k", "1M", "10M")) +
      scale_y_continuous(breaks = c(0, 1, 30, 60, 300, 600, 1800, 3600, 7200, 14400), labels = time_label, transform = "pseudo_log") +
      labs(x = "Number of Reads", y = "Time")
  }
ggsave("plots/bench/modes_time.png", p_modes_time, units = "cm", height = 25, width = 25)

p_modes_rate <- filter(benchmark, threads == 1) %>%
  mutate(mode = if_else(mode == "align", str_c(mode, if_else(no_cache, " (uncached)", " (cached)")), mode),
         mode = factor(mode, levels = c("full-read", "inframe", "pattern", "align (cached)", "align (uncached)")),
         target = str_c(str_match(name, "lib:([a-z]*)")[,2], " (", read_length, "bp)"),
         ends = if_else(rev == "None", "Single end", "Paired end")) %>%
  {
    ggplot(., aes(x = reads, y = region_rate, colour = as.character(library_size), linetype = ends, shape = as.character(rep))) +
      facet_grid2(rows = vars(target), cols = vars(mode), render_empty = FALSE) +
      geom_point() +
      geom_line() +
      scale_colour_brewer(name = "Library size", palette = "Set1") +
      scale_shape_discrete(name = "Rep") +
      scale_linetype_discrete(name = "") +
      scale_x_log10(breaks = c(1e5, 1e6, 1e7), labels = c("100k", "1M", "10M")) +
      scale_y_continuous(breaks = c(1, 1e3, 1e4, 1e5, 1e6), transform = "pseudo_log", limits = c(0, 2e6)) +
      labs(x = "Number of Reads", y = "Processing rate (reads/s)")
  }
ggsave("plots/bench/modes_rate.png", p_modes_rate, units = "cm", height = 25, width = 25)

p_metrics_time <- filter(benchmark, threads == 1, mode == "align", !no_cache, reads < 1e7) %>%
  mutate(target = str_c(str_match(name, "lib:([a-z]*)")[,2], " (", read_length, "bp)"),
         metric = factor(metric, levels = c("exact", "hamming", "bounded-levenshtein", "levenshtein")),
         ends = if_else(rev == "None", "Single end", "Paired end")) %>%
  {
    ggplot(., aes(x = reads, y = total_time, colour = as.character(library_size), linetype = ends, shape = as.character(rep))) +
      facet_grid2(rows = vars(target), cols = vars(metric), render_empty = FALSE) +
      geom_point() +
      geom_line() +
      scale_colour_brewer(name = "Library size", palette = "Set1") +
      scale_shape_discrete(name = "Rep") +
      scale_linetype_discrete(name = "") +
      scale_x_log10(breaks = c(1e5, 1e6), labels = c("100k", "1M")) +
      scale_y_continuous(breaks = c(0, 1, 30, 60, 300, 600, 1800, 3600, 7200, 14400), labels = time_label, transform = "pseudo_log",
                         limits = c(0, 1.5e4)) +
      labs(x = "Number of Reads", y = "Time")
  }
ggsave("plots/bench/metrics_time.png", p_metrics_time, units = "cm", height = 25, width = 25)

p_metrics_rate <- filter(benchmark, threads == 1, mode == "align", !no_cache, reads < 1e7) %>%
  mutate(target = str_c(str_match(name, "lib:([a-z]*)")[,2], " (", read_length, "bp)"),
         metric = factor(metric, levels = c("exact", "hamming", "bounded-levenshtein", "levenshtein")),
         ends = if_else(rev == "None", "Single end", "Paired end")) %>%
  {
    ggplot(., aes(x = reads, y = library_rate, colour = as.character(library_size), linetype = ends, shape = as.character(rep))) +
      facet_grid2(rows = vars(target), cols = vars(metric), render_empty = FALSE) +
      geom_point() +
      geom_line() +
      scale_colour_brewer(name = "Library size", palette = "Set1") +
      scale_shape_discrete(name = "Rep") +
      scale_linetype_discrete(name = "") +
      scale_x_log10(breaks = c(1e5, 1e6), labels = c("100k", "1M")) +
      scale_y_continuous(breaks = c(1e6, 2e6, 3e6, 4e6)) +
      labs(x = "Number of Reads", y = "Regions matched/s")
  }
ggsave("plots/bench/metrics_rate.png", p_metrics_rate, units = "cm", height = 25, width = 25)

p_threads <- mutate(benchmark,
       mode = if_else(mode == "align", str_c(mode, if_else(no_cache, " (uncached)", " (cached)")), mode),
       mode = factor(mode, levels = c("full-read", "inframe", "pattern", "align (cached)", "align (uncached)")),
       target = str_c(str_match(name, "lib:([a-z]*)")[,2], " (", read_length, "bp)"),
       metric = factor(metric, levels = c("exact", "hamming", "bounded-levenshtein", "levenshtein")),
       ends = if_else(rev == "None", "Single end", "Paired end"),
       reads = factor(case_match(reads, 100000 ~ "100k", 1000000 ~ "1M", 10000000 ~ "10M"), levels = c("100k", "1M", "10M"))) %>%
  {
    ggplot(., aes(x = threads, y = total_time, colour = as.character(library_size), linetype = ends, shape = as.character(rep))) +
      facet_nested(rows = vars(mode, metric), cols = vars(target, reads), render_empty = FALSE, scales = "free_y") +
      geom_point() +
      geom_line() +
      scale_colour_brewer(name = "Library size", palette = "Set1") +
      scale_shape_discrete(name = "Rep") +
      scale_linetype_discrete(name = "") +
      scale_x_continuous(breaks = c(1, 2, 4, 6)) +
      scale_y_continuous(breaks = c(1, 30, 60, 300, 600, 1800, 3600, 7200, 14400, 28800), labels = time_label,
                         transform = "pseudo_log") +
      labs(x = "Number of Reads", y = "Regions matched/s")
  }
ggsave("plots/bench/threads.png", p_threads, units = "cm", height = 30, width = 25)
