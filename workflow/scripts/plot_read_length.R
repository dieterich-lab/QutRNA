#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(scales)
library(optparse)

source(file.path(Sys.getenv("QUTRNA2"), "scripts", "utils.R"))


option_list <- list(
  make_option(c("-l", "--limits"),
              type = "character",
              help = "min,max"),
  make_option(c("-b", "--breaks"),
              type = "character",
              help = "Breaks separaterd by ','"),
  make_option(c("-t", "--type"),
            type = "character",
            default = "subsample",
            help = "sample, subsample, or condition"),
  make_option(c("-i", "--ignore_read_type"),
              type = "character",
              help = "Ignore read type"),
  make_option(c("-o", "--output"),
              type = "character",
              help = "Output"),
  rename_read_type_option(),
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
stopifnot(is.element(opts$options$type, c("condition", "sample", "subsample")))
ggsave_opts <- custom_ggsave_stopifnot(opts$options$ggsave_opts)
rename_read_type <- rename_read_type_stopifnot(opts$options$rename_read_type)

df <- read.table(opts$args,
                 sep = "\t",
                 header = TRUE)

df <- df |>
  dplyr::select(-c(unique_reads, multimapper_reads, records))

if (!is.null(opts$options$ignore_read_type)) {
  ignore_read_type <- strsplit(opts$options$ignore_read_type, ",")[[1]]

  df <- df |>
    filter(!is.element(read_type, ignore_read_type))
}

df <- df |>
  summarise(reads = sum(reads), .by = c(read_length, all_of(opts$options$type), base_calling, read_type)) |>
  mutate(base_calling = factor(base_calling, levels = c("pass", "fail", "merged", "unknown"), ordered = TRUE))

if (!is.null(rename_read_type)) {
  v <- tibble::deframe(rename_read_type)
  df <- df |>
    mutate(read_type =
             ifelse(read_type %in% names(v), v[read_type], read_type))
}
df <- df |>
  mutate(read_type = factor(read_type, levels = unique(read_type), ordered = TRUE))

plot <- df |>
  ggplot(aes(x = read_length, weight = reads, colour = !!sym(opts$options$type))) +
  geom_density() +
  labs(x = "read length [nt]") +
  guides(colour = guide_legend(opts$options$type)) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 18),
        strip.text.y.right = element_text(angle = 0))


if (is.null(opts$options$limits)) {
  limits <- NULL
} else {
  limits <- strsplit(opts$options$limits, ",")[[1]] |>
    as.integer()
}
if (is.null(opts$options$breaks)) {
  breaks <- waiver()
} else {
  breaks <- strsplit(opts$options$breaks, ",")[[1]] |>
    as.integer()
}

plot <- plot + 
  scale_x_continuous(limits = limits,
                     breaks = breaks)

if (length(unique(df$base_calling)) == 1) {
  plot <- plot + 
    facet_wrap(read_type ~ ., ncol = 1)
} else {
  plot <- plot + 
    facet_grid(read_type ~ base_calling, labeller = function(...) { return(label_both(sep = ":\n", ...)) })
}

infer_height <- function(df) {
  read_types <- df$read_type |>
    unique() |>
    length()
  
  height <- read_types * 1.1
  
  return(height)
}
infer_width <- function(df) {
  width <- 8
  
  return(width)
}
infer_ggsave_opts <- function(df, ggsave_opts) {
  if (is.null(ggsave_opts)) {
    ggsave_opts <- list(
      "width" = infer_width(df),
      "height" = infer_height(df),
      scale = 1
    )
  } else {
    if (is.null(ggsave_opts$width)) {
      ggsave_opts$width = infer_width(df)
    }
    if (is.null(ggsave_opts$height)) {
      ggsave_opts$height = infer_height(df)
    }
  }
}

ggsave_opts <- infer_ggsave_opts(df, ggsave_opts)
custom_ggsave(opts$options$output, plot, ggsave_opts)
