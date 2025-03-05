#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(scales)

option_list <- list(
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output")
)

opts <- optparse::parse_args(
  optparse::OptionParser(option_list = option_list),
  # c("--output", "~/tmp/test2.pdf", "/beegfs/prj/tRNA_Francesca_Tuorto/data/20250115_FT_HCT116_tRNA_RNA004/qutrna/nuclear_trnas/output/pass_local_99-9_filter-reads-strict-no-marks/results/read_counts.tsv"),
  positional_arguments = TRUE
)

stopifnot(length(opts$args) == 1)

df <- read.table(opts$args,
                 sep = "\t",
                 header = TRUE) |>
  mutate(base_calling = factor(base_calling, levels = c("pass", "fail", "merged", "unknown"), ordered = TRUE))

p <- df |>
  ggplot(aes(x = subsample, y = read_count, fill = read_type)) +
  labs(x = "read type", y = "reads") +
  guides(fill = guide_legend("read type")) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_grid(base_calling ~ ., labeller = label_both)

ggsave(opts$options$output, p, width = 8, height = 6)
