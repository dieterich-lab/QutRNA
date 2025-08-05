#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(scales)
library(optparse)

option_list <- list(
  make_option(c("-r", "--range"),
              type = "character",
              help = "min,max"),
  make_option(c("-b", "--breaks"),
              type = "character",
              help = "Breaks separaterd by ','"),
  make_option(c("-t", "--type"),
            type = "character",
            help = "sample, subsample, or condition"),
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
  make_option(c("-o", "--output"),
              type = "character",
              help = "Output")
)

opts <- parse_args(
  OptionParser(option_list = option_list),
  # FIXME remove 
  # c("--type", "sample", "--output", "~/tmp/test.pdf", "/beegfs/prj/tRNA_Francesca_Tuorto/qutrna_paper/test/adapter_length/test_bam/results/stats/samtools_RL.txt"),
  positional_arguments = TRUE
)

stopifnot(length(opts$args) == 1)
stopifnot(is.element(opts$options$type, c("condition", "sample", "subsample")))

df <- read.table(opts$args,
                 sep = "\t",
                 header = TRUE)

if (!is.null(opts$options$ignore_read_type)) {
  ignore_read_type <- strsplit(opts$options$ignore_read_type, ",")[[1]]
  df <- df |>
    filter(!is.element(read_type, ignore_read_type))
}

df <- df |>
  summarise(count = sum(count), .by = c(read_length, all_of(opts$options$type), base_calling, read_type)) |>
  mutate(read_type = factor(read_type, levels = unique(read_type), ordered = TRUE),
         base_calling = factor(base_calling, levels = c("pass", "fail", "merged", "unknown"), ordered = TRUE))


p <- df |>
  ggplot(aes(x = read_length, weight = count, colour = !!sym(opts$options$type))) +
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
  limits <- strsplit(opts$options$limits, ",")[[1]]
}
if (is.null(opts$options$breaks)) {
  breaks <- waiver()
} else {
  breaks <- strsplit(opts$options$breaks, ",")[[1]]
}

p <- p + 
  scale_x_continuous(limits = limits,
                     breaks = breaks)

if (length(unique(df$base_calling)) == 1) {
  p <- p + 
    facet_wrap(read_type ~ .,
               ncol = 1,
               labeller = function(...) { return(label_both(sep = ": ", ...)) })
} else {
  p <- p + 
    facet_grid(read_type ~ base_calling, labeller = function(...) { return(label_both(sep = ":\n", ...)) })
}

n <- df$read_type |> unique() |> length()
ggsave(opts$options$output, p)