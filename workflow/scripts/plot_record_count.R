#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(scales)
library(optparse)


option_list <- list(
  make_option(c("-t", "--type"),
              type = "character",
              help = "Sample, subsample, or condition"),
  make_option(c("-i", "--ignore_read_type"),
              type = "character",
              help = "Ignore read type"),
  make_option(c("-w", "--width"),
              type = "character",
              default = NA,
              help = "Width of plot"),
  make_option(c("-e", "--height"),
              type = "character",
              default = NA,
              help = "Height of plot"),
  make_option(c("-p", "--add_percent"),
              type = "logical",
              default = FALSE,
              help = "Show percentage"),
  make_option(c("-o", "--output"),
                type = "character",
                help = "Output")
)

# FIXME remove c("--type", "sample", "--output", "~/tmp/test2.pdf", "/beegfs/prj/tRNA_Francesca_Tuorto/qutrna_paper/test/adapter_length/test1/results/stats/record_count.txt"),
#c("--type", "subsample", "--output", "~/tmp/test2.pdf", "~/scrap/qutrna/stats/record_count_pass_fail.txt"),
opts <- parse_args(
  OptionParser(option_list = option_list),
  positional_arguments = TRUE
)

stopifnot(length(opts$args) == 1)
stopifnot(is.element(opts$options$type, c("sample", "subsample", "condition")))

df <- read.table(opts$args,
                 sep = "\t",
                 header = TRUE)

if (!is.null(opts$options$ignore_read_type)) {
  ignore_read_type <- strsplit(opts$options$ignore_read_type, ",")[[1]]
  df <- df |>
    filter(!is.element(read_type, ignore_read_type))
}

df <- df |>
  mutate(base_calling = factor(base_calling, levels = c("pass", "fail", "merged", "unknown"), ordered = TRUE),
         read_type = factor(read_type, levels = rev(unique(read_type)), ordered = TRUE))

df <- df |>
  summarise(records = sum(records), .by = c(all_of(opts$options$type), base_calling, read_type))

p <- df |>
  ggplot(aes(x = read_type, y = records, fill = !!sym(opts$options$type))) +
  xlab("read type") +
  ylab("record count") +
  guides(fill = guide_legend(opts$options$type)) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 18))

if (opts$options$add_percentage) {
  ref_read_type <- levels(df$read_type)[1]
  df <- df |>
    mutate(percent = scales::percent(records / records[records == ref_read_type]),
           .by = c(all_of(opts$options$type), base_calling, read_type))
  p <- p +
    geom_text(aes(label = percent))
}

if (length(unique(df$base_calling)) == 1) {
  p <- p +
    facet_wrap(vars(!!sym(opts$options$type)),
               ncol = 1,
               labeller = function(...) { return(label_both(sep = ": ", ...)) })
} else {
  p <- p +
    facet_grid(rows = vars(!!sym(opts$options$type)),
               cols = vars(base_calling),
               labeller = function(...) { return(label_both(sep = ":\n", ...)) })
}

ggsave(opts$options$output, p)