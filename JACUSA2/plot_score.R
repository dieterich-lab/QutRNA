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
  make_option(c("-c", "--coordinates"),
                type = "character",
                help = "Use column as coordinates, default: Pos3"),
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
df <- df[df$Pos3 > opts$options$left, ]
if (opts$options$len > 0) {
  df <- df[df$Pos3 < as.numeric(opts$options$left) + as.numeric(opts$options$len), ]
}
pos_col <- opts$options$coordinates
if (pos_col == "u_pos") {
  df <- df[df[[opts$options$coordinates]] != "-", ]
}

df$score <- df[, opts$options$score]
df$score[df$score < 0] <- 0

score <- tapply(df$score, list(df$Ref, paste(df[[pos_col]], df$Mod, sep = ":")), max)
score[is.na(score)] <- 0
idx <- order(as.numeric(sapply(sapply(colnames(score), strsplit, ":"), function(x){ x[1] } )))

pheatmap(score[, idx],
         scale = "none",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.0f",
         width = 22,
         height = 14,
         file = file.path(opts$options$output, "main.pdf"))

#score <- tapply(df$score, list(df$Ref, df[[pos_col]]), max)
#score[is.na(score)]<-0
#mod <- sapply(tapply(df$Mod, df[[pos_col]], unique), paste, collapse = "_")
#annotation_col = data.frame(Modification = ifelse(as.character(mod) == "", "Unmod", "Mod"))
#rownames(annotation_col) <- names(mod)
#pheatmap(score[, order(as.numeric(colnames(score)))],
#         scale = "none",
#         annotation_col = annotation_col,
#         cluster_rows = FALSE,
#         cluster_cols = FALSE,
#         display_numbers = TRUE,
#         number_format = "%.0f",
#         width = 22,
#         height = 14,
#         file = file.path(opts$options$output, "small.pdf"))
