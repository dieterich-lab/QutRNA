#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(scales)
library(optparse)

option_list <- list(
  make_option(c("-t", "--type"),
              type = "character",
              help = "sample, subsample, or condition"),
  make_option(c("-o", "--output"),
                type = "character",
                help = "Output")
)


opts <- parse_args(
  OptionParser(option_list = option_list),
  # FIXME remove c("--type", "sample", "--output", "~/tmp/test2.pdf", "/beegfs/prj/tRNA_Francesca_Tuorto/qutrna_paper/test/adapter_length/test1/results/stats/read_count.txt"),
  positional_arguments = TRUE
)

stopifnot(length(opts$args) == 1)
stopifnot(opts$options$type %in% c("sample", "subsample", "condition"))

df <- read.table(opts$args,
                 sep = "\t",
                 header = TRUE) |>
  mutate(base_calling = factor(base_calling, levels = c("pass", "fail", "merged", "unknown"), ordered = TRUE),
         read_type = factor(read_type, levels = rev(unique(read_type)), ordered = TRUE))

df <- df |>
  summarise(reads = sum(reads), .by = c(all_of(opts$options$type), base_calling, read_type))

p <- df |>
  ggplot(aes(x = read_type, y = reads, fill = !!sym(opts$options$type))) +
  xlab("read type") +
  ylab("reads") +
  guides(fill = guide_legend(opts$options$type)) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 18),
        strip.text.y.right = element_text(angle = 0))

if (length(unique(df$base_calling)) == 1) {
  p <- p +
    facet_grid(rows = vars(!!sym(opts$options$type)),
               labeller = function(...) { return(label_both(sep = ":\n", ...)) })
} else {
  p <- p +
    facet_grid(rows = vars(!!sym(opts$options$type)),
               cols = vars(base_calling),
               labeller = function(...) { return(label_both(sep = ":\n", ...)) })
}

# TODO width, height?
ggsave(opts$options$output, p, width = 8, height = 6)
