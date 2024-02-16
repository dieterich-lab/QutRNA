#!/usr/bin/env Rscript

# Add additional data.
library(magrittr)

WIDTH <- 5

option_list <- list(
  optparse::make_option(c("-m", "--mods"),
                        type = "character",
                        help = "Known modifications"),
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output")
)

opts <- optparse::parse_args(
  optparse::OptionParser(option_list = option_list),
  positional_arguments = TRUE
)

stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$options$mods))
stopifnot(length(opts$args) == 1)

mods <- data.table::fread(opts$options$mods, header = TRUE) %>%
  as.data.frame() %>%
  dplyr::rename(seqnames = seq_id, start = u_pos) %>%
  dplyr::mutate(end = start) %>%
  GenomicRanges::GRanges()

result <- read.table(opts$args, header = TRUE)
suppressWarnings({result <- result %>% plyranges::join_overlap_left(mods)})

df %>% write.table(opts$options$output,
                   quote=FALSE, sep = ",",
                   col.names = TRUE, row.names = FALSE)
