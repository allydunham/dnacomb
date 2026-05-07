#!/usr/bin/env Rscript
# Plot test results
library(tidyverse)
library(ggpubr)
library(ggh4x)
dir.create("plots", showWarnings = FALSE)

theme_set(theme_pubclean() + theme(legend.position = 'right',
                                   plot.title = element_text(hjust = 0.5),
                                   plot.subtitle = element_text(hjust = 0.5),
                                   strip.background = element_blank(),
                                   legend.key = element_blank()))

## Import data
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
  
  mutate(all_obs, combination_status = replace_na(combination_status, "missed"))
}
quiet_test <- purrr::quietly(test_pair)

test_counts <- filter(test_pairs, threads == 1) %>%
  {map2(.$observed, .$expected, ~quiet_test(.x, .y)$result, .progress = TRUE)} %>%
  bind_rows() %>%
  replace_na(replace = list(count = 0, true_count = 0)) %>%
  select(-expected) %>%
  separate_wider_delim(observed, delim = ":", names = c("test", "mode", "distance", "library", "end", "threads"), cols_remove = TRUE) %>%
  mutate(threads = as.integer(threads))

## Observed vs expected
category_colours <- c(
  "match" = "green", "exact_match" = "green", "nearest_match" = "darkgreen",
  "mismatch" = "orange", "nonmatch" = "red", 
  "recombination" = "blue", "exact_recombination" = "blue", "nearest_recombination" = "darkblue",
  "low_mean_quality" = "brown", "bad_alignment" = "grey", "multimatch" = "purple",
  "missed" = "black"
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

## Corelations
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
