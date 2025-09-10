#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(scales)
library(optparse)
library(patchwork)

source(file.path(Sys.getenv("QUTRNA2"), "scripts", "utils.R"))


option_list <- list(
  make_option(c("-t", "--type"),
              type = "character",
              help = "BAM type(s)"),
  make_option(c("-c", "--cutoff"),
              type = "character",
              help = "Cutoff(s)"),
  make_option(c("-o", "--output"),
              type = "character",
              help = "Output"),
  custom_ggsave_option()
)

if (exists("PLOT_ARGS")) {
  opts <- parse_args(
    OptionParser(option_list = option_list),
    args = PLOT_ARGS,
    positional_arguments = TRUE
  )
} else {
  opts <- parse_args(
    OptionParser(option_list = option_list),
    positional_arguments = TRUE
  )
}

stopifnot(length(opts$args) == 1)
stopifnot(!is.null(opts$options$cutoff))
read_types <- strsplit(opts$options$type, ",")[[1]]
stopifnot(read_types != "")
stopifnot(length(read_types) > 1)
ggsave_opts <- custom_ggsave_stopifnot(opts$options$ggsave_opts)

format_sample_desc <- function(df) {
  base_calling <- length(unique(df$base_calling)) > 1
  i <- df$sample != df$subsample
  sample_desc <- df$sample
  if (any(i)) {
    sample_desc[i] <- paste0(df$sample[i], "\n", df$subsample[i])
    if (base_calling) {
      sample_desc[i] <- paste0(sample_desc[i], "\n(", df$base_calling[i], ")")
    }
  } else {
    if (base_calling) {
      sample_desc <- paste0(sample_desc[i], "\n(", df$base_calling[i], ")")
    }
  }
  
  return(sample_desc)
}

cutoffs <- read.table(opts$options$cutoff, header = TRUE, sep = "\t") |>
  mutate(label = paste0("Threshold: ", cutoff),
         category = case_match(
           read_type,
           "mapped" ~ "Real",
           "mapped-random" ~ "Random"),
         category = as.factor(category))
cutoffs$sample_desc <- format_sample_desc(cutoffs)
scores <- read.table(opts$args,
                 sep = "\t",
                 header = TRUE) |>
  filter(is.element(read_type, read_types)) |>
  mutate(base_calling = factor(base_calling, levels = c("pass", "fail", "merged", "unknown"), ordered = TRUE),
         category = case_match(
           read_type,
           "mapped" ~ "Real",
           "mapped-random" ~ "Random"),
         category = as.factor(category))
scores$sample_desc <- format_sample_desc(scores)

aln_score_lim <- range(scores$alignment_score)
count_lim <- range(scores$count)
plots <- NULL
for (cond in unique(scores$condition)) {
  filtered_cutoffs <- cutoffs |>
    filter(condition == cond)
  p <- scores |>
    filter(condition == cond) |>
    ggplot(aes(x = alignment_score, y = count, fill = category)) +
    ylab("Reads") +
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
    theme_bw(base_size = 18) +
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


infer_width <- function(scores) {
  width <- length(unique(scores$condition))
  
  return(width)
}
infer_height <- function(scores) {
  return(
    scores |>
      summarise(n = length(unique(subsample)),
                .by = condition) |>
      pull(n) |>
      max() / 2)
}
infer_ggsave_opts <- function(df, ggsave_opts) {
  if (is.null(ggsave_opts)) {
    ggsave_opts <- list(
      "width" = infer_width(scores),
      "height" = infer_height(scores),
      scale = 4
    )
  } else {
    if (is.null(ggsave_opts$width)) {
      ggsave_opts$width = infer_width(scores)
    }
    if (is.null(ggsave_opts$height)) {
      ggsave_opts$height = infer_height(scores)
    }
  }
}

custom_ggsave(
  opts$options$output, 
  plots, 
  ggsave_opts = infer_ggsave_opts(scores, ggsave_opts))
