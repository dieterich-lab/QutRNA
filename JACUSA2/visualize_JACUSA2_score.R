#!/usr/bin/env Rscript

library(optparse)
library(pheatmap)

option_list <- list(
  make_option(c("-s", "--score"),
                type = "character",
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
                help = "Ignore positions right"),
)

opts <- optparse::parse_args(
  optparse::OptionParser(option_list = option_list),
  positional_arguments = TRUE
)

df <- read.csv(opts$args)
df$lab <- paste(df$Pos, df$ModStatus, sep=":")
df$Mis.Del.Ins[df$Mis.Del.Ins < 0] <- 0

score <- tapply(df$Mis.Del.Ins, list(df$Ref, df$lab), max)
score[is.na(score)] <- 0

idxs <- order(as.numeric(sapply(sapply(colnames(scores), strsplit, ":"), function(x){ x[1] } )))
Validxs <- sort(as.numeric(sapply(sapply(colnames(scores), strsplit, ":"), function(x){ x[1] } )))
idxs <- idxs[which(Validxs > opts$options$left & Validxs < opts$options$left + opts$options$length)]

pheatmap(score[, idxs],
         scale="none",
         cluster_rows=F,
         cluster_cols=F,
         display_numbers = T,
         number_format = "%.0f",
         width=22,
         height=14,
         file=paste0(k,".pdf"))

score <- tapply(feat$Mis.Del.Ins, list(feat$Ref,feat$Pos3), max)
score[is.na(score)] <- 0
blu=sapply(tapply(feat$ModStatus,feat$Pos3,unique),paste,collapse="_")
annotation_col = data.frame(Modification = ifelse(as.character(blu) == "", "Unmod", "Mod"))
rownames(annotation_col) <- names(blu)
pheatmap(gaga[,order(as.numeric(colnames(gaga)))],
         scale = "none",
         annotation_col = annotation_col,
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = T,
         number_format = "%.0f",
         width = 22,
         height = 14,
         file = paste0(k,"_small.pdf"))
