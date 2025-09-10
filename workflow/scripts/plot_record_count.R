#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(scales)
library(optparse)

source(file.path(Sys.getenv("QUTRNA2"), "scripts", "utils.R"))


option_list <- list(
  make_option(c("-t", "--type"),
              default = "subsample",
              type = "character",
              help = "Sample, subsample, or condition"),
  make_option(c("-i", "--ignore_read_type"),
              type = "character",
              help = "Ignore read type"),
  make_option(c("-p", "--show_percent"),
              action = "store_true",
              default = FALSE,
              help = "Show percentage"),
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
stopifnot(is.element(opts$options$type, c("sample", "subsample", "condition")))
ggsave_opts <- custom_ggsave_stopifnot(opts$options$ggsave_opts)
rename_read_type <- rename_read_type_stopifnot(opts$options$rename_read_type)

df <- read.table(opts$args,
                 sep = "\t",
                 header = TRUE)

if (!is.null(opts$options$ignore_read_type)) {
  ignore_read_type <- strsplit(opts$options$ignore_read_type, ",")[[1]]
  df <- df |>
    filter(!is.element(read_type, ignore_read_type))
}

df <- df |>
  mutate(base_calling = factor(base_calling, levels = c("pass", "fail", "merged", "unknown"), ordered = TRUE))

if (!is.null(rename_read_type)) {
  v <- tibble::deframe(rename_read_type)
  df <- df |>
    mutate(read_type =
             ifelse(read_type %in% names(v), v[read_type], read_type))
}
df <- df |>
  mutate(read_type = factor(read_type, levels = rev(unique(read_type)), ordered = TRUE))

df <- df |>
  summarise(reads = sum(reads), .by = c(all_of(opts$options$type), base_calling, read_type))

if (opts$options$show_percent) {
  n <- levels(df$read_type) |> length()
  ref_read_type <- levels(df$read_type)[n]
  df <- df |>
    mutate(percent = reads / reads[read_type == ref_read_type],
           .by = c(all_of(opts$options$type), base_calling)) |>
    mutate(percent_label =  label_percent(.1)(percent))
}

plot <- df |>
  ggplot(aes(x = read_type, y = reads, fill = !!sym(opts$options$type))) +
  xlab("Read type") +
  ylab("Reads") +
  guides(fill = guide_legend(opts$options$type)) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 18))

if (opts$options$show_percent) {
  plot <- plot +
    geom_text(aes(label = percent_label), hjust = "inward")
}

if (length(unique(df$base_calling)) == 1) {
  plot <- plot +
    facet_wrap(vars(!!sym(opts$options$type)),
               ncol = 1,
               labeller = function(...) { return(label_both(sep = ": ", ...)) })
} else {
  plot <- plot +
    facet_grid(rows = vars(!!sym(opts$options$type)),
               cols = vars(base_calling),
               labeller = function(...) { return(label_both(sep = ":\n", ...)) })
}

infer_height <- function(df) {
  types <- df[, opts$options$type] |>
    unique() |>
    length()
  
  read_types <- df$read_type |>
    unique() |>
    length()
  
  height <- types * read_types * 0.3
  
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
custom_ggsave(opts$options$output, plot, ggsave_opts = ggsave_opts)
