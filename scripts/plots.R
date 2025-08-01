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
  separate_wider_delim(observed, delim = ":", names = c("mode", "distance", "expected", "end", "threads"), cols_remove = FALSE)

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
  separate_wider_delim(observed, delim = ":", names = c("mode", "distance", "library", "end", "threads"), cols_remove = TRUE) %>%
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
ggsave("plots/test_all_counts_scatter.png", p_all_scatter, units = "cm", height = 50, width = 50)

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
ggsave("plots/test_library_scatter.png", p_library_scatter, units = "cm", height = 50, width = 50)

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
ggsave("plots/test_summary_bars.png", p_summary_bars, units = "cm", height = 40, width = 40)

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
ggsave("plots/test_observed_expected_correlation.png", p_test_cors, units = "cm", height = 30, width = 16)

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
ggsave("plots/thread_count_correlation.png", p_thread_cor, units = "cm", height = 10, width = 10)

# Benchmarks
na_or_zero <- function(x) {
  if_else(is.na(x), 0, x)
}

time_label <- function(x) {
  out <- str_c(x, "s")
  out[x >= 60] <- str_c(round(x[x >= 60]/60, 1), "min")
  out[x >= 3600] <- str_c(round(x[x >= 3600]/3600, 1), "h")
  return(out)
}

# Captures any number of benchmark TSVs called bench1, bench2, ...
bench_cols <- c(
  "name", "fwd", "rev", "lib_spec", "mode", "metric",
  "no_cache", "sort", "group", "library_counts",
  "library_size", "read_length", "additional_args",
  "total_time", "reads", "region_time", "region_rate",
  "unique_regions", "library_time", "library_rate",
  "summary_size", "summary_time", "summary_rate"
)
benchmark <- dir("data/benchmark/", pattern = "bench[0-9]*.tsv", full.names = TRUE) %>%
  set_names(1:length(.)) %>%
  map(read_tsv, col_names = bench_cols, skip = 1) %>%
  bind_rows(.id = "rep")

p_format <- filter(benchmark, str_detect(name, "Format")) %>%
  select(rep, name, total_time, region_time) %>%
  mutate(format = c("Format (fq)"="Fastq", "Format (fa)"="Fasta")[name]) %>%
  pivot_longer(ends_with("time"), names_to = "type", values_to = "time") %>%
  ggplot(aes(x = format, y = time, fill = rep)) +
  facet_wrap(~type, nrow = 1, labeller = labeller(type = c(total_time="Total", region_time="Region Extraction"))) +
  geom_col(position = position_dodge(), width = 0.6) +
  scale_fill_brewer(name = "Rep", palette = "Dark2") +
  labs(x = "", y = "Time (s)", caption = "Using inframe matching with simulated 1M read gRNA files")
ggsave("plots/bench_file_format.png", p_format, units = "cm", height = 12, width = 16)

p_mode <- filter(benchmark, str_detect(name, "Extraction")) %>%
  mutate(paired = if_else(rev == "None", "Single-end", "Paired-end"),
         library = str_match(name, "\\((.*)/.*/.*\\)")[,2],
         library = str_c(str_to_lower(library), "\n(", read_length, "nt)"),
         mode = factor(mode, levels = c("full-read", "inframe", "pattern", "align"))) %>%
  select(rep, library, read_length, mode, paired, time = region_time) %>%
  group_by(library, read_length, mode, paired) %>%
  summarise(mean = mean(time), min = min(time), max = max(time), .groups = "drop") %>%
  ggplot(aes(x = library, y = mean, ymin = min, ymax = max, colour = mode, shape = paired, group = mode)) +
  geom_errorbar(position = position_dodge(width = 0.5), width = 0.5, show.legend = FALSE) +
  geom_point(position = position_dodge(width = 0.5)) +
  scale_colour_brewer(name = "", palette = "Set1") +
  scale_shape_discrete(name = "") +
  labs(x = "", y = "Region Extraction Time (s)",
       caption = "Based on simulated 1M read gRNA fastq files with no mutations")
ggsave("plots/bench_mode.png", p_mode, units = "cm", height = 12, width = 16)

p_metric <- filter(benchmark, str_detect(name, "Library comparison")) %>%
  select(rep, metric, library_time, library_rate) %>%
  pivot_longer(c(library_time, library_rate), names_to = "type", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("exact", "hamming", "bounded-levenshtein", "levenshtein"))) %>%
  {
    ggplot(., aes(x = metric, y = value, fill = metric)) +
      geom_boxplot(show.legend = FALSE) +
      facet_wrap(~type, scales = "free_y", strip.position = "left",
                 labeller = labeller(type = c(library_rate = "Items/s", library_time = "Time (s)"))) +
      scale_fill_brewer(palette = "Set2") +
      labs(x = "Distance Metric", y = "",
           caption = "Based on simulated 1M read gRNA fastq files from a 10k library") +
      theme(strip.placement = "outside")
  }
ggsave("plots/bench_metric.png", p_metric, units = "cm", height = 12, width = 30)

p_lib_size <- filter(benchmark, str_detect(name, "Library size")) %>%
  select(rep, name, library_size, library_time, library_rate) %>%
  pivot_longer(c(library_time, library_rate), names_to = "type", values_to = "time") %>%
  group_by(library_size, type) %>%
  summarise(mean = mean(time), min = min(time), max = max(time), .groups = "drop") %>%
  ggplot(aes(x = library_size, y = mean, ymin = min, ymax = max)) +
  facet_wrap(~type, scales = "free_y", strip.position = "left",
             labeller = labeller(type = c(library_rate = "Items/s", library_time = "Time (s)"))) +
  geom_errorbar(width = 2000, show.legend = FALSE) +
  geom_point() +
  labs(x = "Library Size", y = "",
       caption = "Based on simulated 1M read gRNA files with moderate mutation") +
  theme(strip.placement = "outside")
ggsave("plots/bench_lib_size.png", p_lib_size, units = "cm", height = 12, width = 16)

p_read_count <- filter(benchmark, str_detect(name, "Sequence file size")) %>%
  select(rep, name, reads, total_time, region_time, library_time, summary_time) %>%
  pivot_longer(ends_with("time"), names_to = "type", values_to = "time") %>%
  group_by(reads, type) %>%
  summarise(mean = mean(time), min = min(time), max = max(time), .groups = "drop") %>%
  filter(type != "summary_time") %>%
  {
    bar_width <- 0.05 * log10(max(.$reads))
    ggplot(., aes(x = reads, y = mean, ymin = min, ymax = max, colour = type)) +
      geom_line() +
      geom_point() +
      geom_errorbar(width = bar_width) +
      scale_y_log10() +
      scale_x_log10() +
      scale_colour_brewer(palette = "Set1", labels = c(total_time="Total", region_time="Region Extraction",
                                                       summary_time="Summarisation", library_time = "Library")) +
      labs(x = "Number of Reads", y = "Time (s)",
           caption = "Based on simulated gRNA files using a 10k library with moderate mutation")
  }
ggsave("plots/bench_read_count.png", p_read_count, units = "cm", height = 12, width = 16)
