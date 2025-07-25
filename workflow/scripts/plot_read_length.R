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
  # FIXME remove c("--type", "sample", "--output", "~/tmp/test.pdf", "/beegfs/prj/tRNA_Francesca_Tuorto/qutrna_paper/test/adapter_length/test1/results/stats/samtools_RL.txt"),
  positional_arguments = TRUE
)

stopifnot(length(opts$args) == 1)
stopifnot(opts$options$type %in% c("condition", "sample", "subsample"))

df <- read.table(opts$args,
                 sep = "\t",
                 header = TRUE)


df <- df |>
  summarise(count = sum(count), .by = c(read_length, all_of(opts$options$type), base_calling, read_type)) |>
  mutate(read_type = factor(read_type, levels = unique(read_type), ordered = TRUE),
         base_calling = factor(base_calling, levels = c("pass", "fail", "merged", "unknown"), ordered = TRUE))


p <- df |>
  ggplot(aes(x = read_length, weight = count, colour = !!sym(opts$options$type))) +
  geom_density() +
  scale_x_continuous(limits = c(20, 250),
                     breaks = c(20, 40, 60, 100, 150, 200, 250)) +
  labs(x = "read length [nt]") +
  guides(colour = guide_legend(opts$options$type)) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 18),
        strip.text.y.right = element_text(angle = 0))

if (length(unique(df$base_calling)) == 1) {
  p <- p + 
    facet_grid(read_type ~ ., labeller = function(...) { return(label_both(sep = ":\n", ...)) })
} else {
  p <- p + 
    facet_grid(read_type ~ base_calling, labeller = function(...) { return(label_both(sep = ":\n", ...)) })
}

n <- df$read_type |> unique() |> length()
ggsave(opts$options$output, p) # TODO, width = max(12, n * 1.5), height = 8)