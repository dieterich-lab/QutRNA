#!/usr/bin/env Rscript

require(optparse)
require(yardstick)
require(ggplot2)
require(dplyr)

option_list <- list(
  make_option(c("-e", "--real"),
              type = "character",
              help = "Scores of real reads mapping"),
  make_option(c("-a", "--random"),
              type = "character",
              help = "Scores of reads mapping to random"),
  #make_option(c("-P", "--prc_plot"),
  #            type = "character",
  #            help = "PRC output"),
  make_option(c("-S", "--score_plot"),
              type = "character",
              help = "Score output"),
  make_option(c("-C", "--cutoff"),
              type = "character",
              help = "Cutoff  output"),
  make_option(c("-p", "--precision"),
              type = "double",
              default = 0.999,
              help = "Precision cutoff"),
  make_option(c("-t", "--title"),
              type = "character",
              default = "tRNA Alignment Score distributions",
              help = "Alignment score plot title")
)

opts <- parse_args(
  OptionParser(option_list = option_list),
  positional_arguments = TRUE
)

stopifnot(!is.null(opts$options$real))
stopifnot(!is.null(opts$options$random))
stopifnot(file.exists(opts$options$real))
stopifnot(file.exists(opts$options$random))
stopifnot(opts$options$score_plot != "")
#stopifnot(opts$options$prc_plot != "")
stopifnot(opts$options$cutoff != "")
stopifnot(opts$options$precision >= 0 & opts$options$precision <= 1)

read_score <- function(fname, category) {
 df <- read.table(fname, header = TRUE, sep = "\t")
 df <- data.frame(
   alignment_score = rep(df$alignment_score, times = df$count),
   category = category
 )

 return(df)
}

scores <- dplyr::bind_rows(
  read_score(opts$options$real, "POS"),
  read_score(opts$options$random, "NEG")) |>
  mutate(category = as.factor(category) |> relevel("POS"))

curve <- scores |>
  pr_curve(category, alignment_score)
#p_curve <- curve |>
#  autoplot()
#ggsave(opts$options$prc_plot, p_curve)
df_curve <- curve |>
  as.data.frame()

idx <- which(df_curve$precision > opts$options$precision)
cutoff <- df_curve[idx[length(idx)], 1]
if (as.character(cutoff) == "Inf") {
 stop("cutoff cannot be 'Inf'")
}
write.table(
  data.frame(cutoff = cutoff),
  opts$options$cutoff,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t")

annotation <- data.frame(
  x = cutoff,
  y = I(0.5),
  category = "",
  label = paste0("Threshold: ", cutoff)
)

p_score <- scores |>
  mutate(category = case_match(
    category,
    "POS" ~ "Real",
    "NEG" ~ "Random")) |>
  ggplot(aes(x = alignment_score, fill = category)) +
  ylab("Frequency") + 
  scale_fill_manual(
    name = "Alignment",
    values = c("Real" = "orange", "Random" = "blue"),
    breaks = c("Real", "Random")) +
  xlab("Alignment score") +
  geom_histogram(binwidth = 1, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = cutoff, colour = "red") +
  geom_text(data = annotation, aes(x = x, y = y, label = label), colour = "red", hjust = -.1) +
  theme_bw()
ggsave(opts$options$score_plot, p_score)