#!/usr/bin/env Rscript

# Plor JACUSA2 score and existing modification info as a heatmap)

library(optparse)
library(pheatmap)

option_list <- list(
  make_option(c("-s", "--score"),
                type = "character",
                default = "Mis.Del.Ins",
                help = "Score column to use"),
  make_option(c("-o", "--output"),
                type = "character",
                help = "Output"),
  make_option(c("-l", "--left"),
                type = "character",
                default = 24,
                help = "Ignore positions left"),
  make_option(c("-t", "--length"),
                type = "character",
                default = 77,
                help = "Ignore positions right")
)

opts <- parse_args(
  OptionParser(option_list = option_list),
  positional_arguments = TRUE
)
stopifnot(!is.null(opts$options$output))
stopifnot(length(opts$args) == 1)

dir.create(opts$options$output, showWarnings = FALSE)
df <- read.csv(opts$args)
df$score <- df[, opts$options$score]
df$score[df$score < 0] <- 0

score <- tapply(df$score, list(df$Ref, paste(df$Pos, df$Mod, sep = ":")), max)
score[is.na(score)] <- 0

idxs <- order(as.numeric(sapply(sapply(colnames(score), strsplit, ":"), function(x){ x[1] } )))
Validxs <- sort(as.numeric(sapply(sapply(colnames(score), strsplit, ":"), function(x){ x[1] } )))
left <- as.numeric(opts$options$left)
len <- as.numeric(opts$options$length)
idxs <- idxs[which(Validxs > left & Validxs < left + len)]

pheatmap(score[, idxs],
         scale = "none",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.0f",
         width = 22,
         height = 14,
         file = file.path(opts$options$output, "main.pdf"))

score <- tapply(df$score, list(df$Ref, df$Pos3), max)
score[is.na(score)]<-0
mod <- sapply(tapply(df$Mod, df$Pos3, unique), paste, collapse = "_")
annotation_col = data.frame(Modification = ifelse(as.character(mod) == "", "Unmod", "Mod"))
rownames(annotation_col) <- names(mod)
pheatmap(score[, order(as.numeric(colnames(score)))],
         scale = "none",
         annotation_col = annotation_col,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.0f",
         width = 22,
         height = 14,
         file = file.path(opts$options$output, "small.pdf"))
