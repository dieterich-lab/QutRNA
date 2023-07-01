#!/usr/bin/env Rscript

require(optparse)
require(yardstick)
require(magrittr)
require(ggplot2)

option_list <- list(
  make_option(c("-f", "--forward"),
              type = "character",
              help = "Scores of reads mapping to forward"),
  make_option(c("-r", "--reverse"),
              type = "character",
              help = "Scores of reads mapping to reverse"),
  make_option(c("-o", "--output"),
              type = "character",
              help = "Output directory"),
  make_option(c("-p", "--precision"),
              type = "double",
              default = 0.95,
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

stopifnot(!is.null(opts$options$forward))
stopifnot(!is.null(opts$options$reverse))
stopifnot(file.exists(opts$options$forward))
stopifnot(file.exists(opts$options$reverse))
stopifnot(opts$options$output != "")
stopifnot(opts$options$precision >= 0 & opts$options$precision <= 1)

fwdScore <- scan(opts$options$forward, numeric())
revScore <- scan(opts$options$reverse, numeric())

dir.create(opts$options$output, showWarnings = FALSE)

df <- data.frame(
  Class = as.factor(
    c(rep("POS", length(fwdScore)),
      rep("NEG",length(revScore)))),
  Score = c(fwdScore, revScore))
df$Class <- relevel(df$Class, "POS")

p <- df %>%
  pr_curve(Class, Score) %>%
  autoplot()
ggsave(file.path(opts$options$output, "PRC_reads.pdf"), p)

score.df <- df %>%
  pr_curve(Class, Score) %>%
  as.data.frame()

idx <- which(score.df$precision > opts$options$precision)
CutOff <- score.df[idx[length(idx)], 1]
write(CutOff, file = file.path(opts$options$output, "CutOff.txt"))

pdf(file.path(opts$options$output, "Alignment_score_distributions.pdf"))
  hist(fwdScore, main = opts$options$title, sub = CutOff, col = "orange", xlab = "Score")
  hist(revScore, add = TRUE, col = "blue")
  abline(v = CutOff,col = "red", lwd = 2)
  legend("topright", c("Real", "Random"), fill = c("orange","blue"))
dev.off()
