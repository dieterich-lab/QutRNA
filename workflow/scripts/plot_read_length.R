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
  # c("--output", "~/tmp/test.pdf", "/beegfs/prj/tRNA_Francesca_Tuorto/data/20250115_FT_HCT116_tRNA_RNA004/qutrna/nuclear_trnas/output/pass_local_99-9_filter-reads-strict-no-marks/results/samtools/stats/RL.tsv"),
  positional_arguments = TRUE
)

stopifnot(length(opts$args) == 1)

df <- read.table(opts$args,
                 sep = "\t",
                 header = TRUE)

df_munged <- df |>
  select(-fname) |>
  mutate(proportion = count / sum(count),
         .by = c(read_type, subsample, base_calling, condition),
         base_calling = factor(base_calling, levels = c("pass", "fail", "merged", "unknown"), ordered = TRUE))

p <- df_munged |>
  ggplot(aes(x = read_length, y = proportion, group = subsample, colour = sample)) +
  geom_line() +
  scale_x_continuous(limits = c(20, 250),
                     breaks = c(20, 40, 60, 100, 150, 200, 250)) +
  labs(x = "read length [nt]", y = "proportion") +
  guides(colour = guide_legend("condition")) +
  theme_bw() +
  theme(legend.position = "bottom",
        text=element_text(size = 18)) +
  facet_grid(base_calling ~ read_type, labeller = label_both) 

ggsave(opts$options$output, p, width = 12, height = 8)