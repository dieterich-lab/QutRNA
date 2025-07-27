#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(scales)
library(optparse)
library(patchwork)

option_list <- list(
  make_option(c("-t", "--type"),
              type = "character",
              help = "BAM type(s)"),
  make_option(c("-c", "--cutoff"),
              type = "character",
              help = "Cutoff(s)"),
  make_option(c("-o", "--output"),
              type = "character",
              help = "Output")
)


opts <- parse_args(
  OptionParser(option_list = option_list),
  #c("--type", "mapped,mapped-rev",
  #  "--output", "~/tmp/results/plots/alignment/threshold_summary.pdf",
  #  "/beegfs/prj/tRNA_Francesca_Tuorto/qutrna_paper/test/adapter_length/test1/results/stats/alignment_score.txt"),
  positional_arguments = TRUE
)

stopifnot(length(opts$args) == 1)
stopifnot(!is.null(opts$options$cutoffs))
read_types <- strsplit(opts$options$type, ",")[[1]]
stopifnot(read_types != "")
stopifnot(length(read_types) > 1)

cutoffs <- read.table(opts$options$cutoffs, header = TRUE, sep = "\t")
scores <- read.table(opts$args,
                 sep = "\t",
                 header = TRUE) |>
  filter(read_type %in% read_types) |>
  mutate(base_calling = factor(base_calling, levels = c("pass", "fail", "merged", "unknown"), ordered = TRUE),
         category = case_match(
           read_type,
           "mapped" ~ "Real",
           "mapped-rev" ~ "Random"),
         category = as.factor(category),
         label = paste0("sample~", sample, ", subsample~", subsample, ", ", "base calling~", base_calling))

alignment_score_lim <- range(scores$alignment_score)
count_lim <- range(scores$count)
plots <- NULL
for (cond in unique(scores$condition)) {
  filtered_cutoffs <- cutoffs |>
    filter(condition == cond)
  p <- scores |>
    filter(condition == cond) |>
    ggplot(aes(x = alignment_score, y = count, fill = category)) +
    ylab("Frequency") +
    scale_fill_manual(
      name = "Alignment",
      values = c("Real" = "orange", "Random" = "blue"),
      breaks = c("Real", "Random")) +
    xlab("Alignment score") +
    geom_col() +
    geom_vline(data = filtered_cutoffs, aes(xintercept = cutoff), colour = "red") +
    geom_text(data = filtered_cutoffs, aes(x = cutoff, y = I(0.5), label = label), colour = "red", hjust = -.1) +
    theme_bw() +
    facet_grid(label ~ condition) +
    xlim(alignment_score_lim) + ylim(count_lim)
  if (is.null(plots)) {
    plots <- p
  } else {
    plots <- plots + p
  }
}
plots <- plots +
  plot_layout(
    guides = "collect",
    axes = "collect",
    ncol = length(unique(df$condition)))

ggsave(opts$options$output, p)
