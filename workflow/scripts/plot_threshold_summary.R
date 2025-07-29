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
  # TODO
  # c("--type", "mapped,mapped-rev",
  #   "--output", "~/tmp/results/plots/alignment/threshold_summary.pdf",
  #   "--cutoff", "/beegfs/prj/tRNA_Francesca_Tuorto/qutrna_paper/test/adapter_length/test_bam/results/stats/cutoff.txt",
  #   "/beegfs/prj/tRNA_Francesca_Tuorto/qutrna_paper/test/adapter_length/test_bam/results/stats/alignment_score.txt"),
  positional_arguments = TRUE
)

stopifnot(length(opts$args) == 1)
stopifnot(!is.null(opts$options$cutoff))
read_types <- strsplit(opts$options$type, ",")[[1]]
stopifnot(read_types != "")
stopifnot(length(read_types) > 1)

cutoffs <- read.table(opts$options$cutoff, header = TRUE, sep = "\t") |>
  mutate(label = paste0("Threshold: ", cutoff),
         category = case_match(
           read_type,
           "mapped" ~ "Real",
           "mapped-rev" ~ "Random"),
         category = as.factor(category))
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
         sample_desc = paste0("sample~", sample, "  subsample~", subsample, "  ", "base calling~", base_calling))

aln_score_lim <- range(scores$alignment_score)
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
    xlim(0, aln_score_lim[2] + 1) +
    ylim(0, count_lim[2] + 1) +
    geom_col(alpha = .6, position = "identity") +
    geom_vline(data = filtered_cutoffs, aes(xintercept = cutoff), colour = "red") +
    geom_text(data = filtered_cutoffs, aes(x = cutoff, y = I(0.5), label = label), colour = "red", hjust = -.1) +
    theme_bw() +
    facet_grid(sample_desc ~ condition)
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
    ncol = length(unique(scores$condition)))

ggsave(opts$options$output, plots)
