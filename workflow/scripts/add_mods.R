#!/usr/bin/env Rscript

library(magrittr)

option_list <- list(
  optparse::make_option(c("-m", "--mods"),
                        type = "character",
                        help = "Known modifications"),
  optparse::make_option(c("-s", "--sprinzl"),
                        action="store_true",
                        default = FALSE,
                        help = "Sprinzl coordinates"),
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
  as.data.frame()

result <- read.table(opts$args, header = TRUE)

by <- dplyr::join_by(Ref == trna, Pos3 == pos)
if (opts$options$sprinzl) {
  by <- dplyr::join_by(Ref == trna, sprinzl == pos)
  mods$pos = as.character(mods$pos)
}
df <- dplyr::left_join(result, mods, by = by)

df %>% write.table(opts$options$output,
                   quote=FALSE, sep = "\t",
                   col.names = TRUE, row.names = FALSE)
