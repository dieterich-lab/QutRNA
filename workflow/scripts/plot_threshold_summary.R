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
read_types <- strsplit(opts$options$type, ",")[[1]]
stopifnot(read_types != "")
stopifnot(length(read_types) > 1)

df <- read.table(opts$args,
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

plots <- NULL
for (cond in unique(df$condition)) {
  p <- df |>
    filter(condition == cond) |>
    ggplot(aes(x = alignment_score, y = count, fill = category)) +
    ylab("Frequency") + 
    scale_fill_manual(
      name = "Alignment",
      values = c("Real" = "orange", "Random" = "blue"),
      breaks = c("Real", "Random")) +
    xlab("Alignment score") +
    geom_col() +
    # TODO add threshold
    # geom_vline(xintercept = cutoff, colour = "red") +
    # geom_text(data = annotation, aes(x = x, y = y, label = label), colour = "red", hjust = -.1) +
    theme_bw() +
    facet_grid(label ~ condition)
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
