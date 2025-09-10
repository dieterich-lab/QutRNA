#!/usr/bin/env Rscript

library(dplyr)

option_list <- list(
  optparse::make_option(c("-m", "--mods"),
                        type = "character",
                        help = "Known modifications"),
  optparse::make_option(c("-a", "--abbrevs"),
                        type = "character",
                        help = "One character abbreviations for mods"),
  optparse::make_option(c("-s", "--sprinzl"),
                        action = "store_true",
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

get_contained_mods <- function(df, mods_fname, use_sprinzl) {
  mods <- read.table(mods_fname, header = TRUE, sep = "\t", quote = "", comment.char = "")
  
  if (use_sprinzl) {
    mods$pos <- as.character(mods$pos)
    contained_mods <- dplyr::left_join(df[, c("trna", "sprinzl")], mods, by = join_by(trna, sprinzl == pos)) |>
      rename(pos = sprinzl)
  } else {
    contained_mods <- dplyr::left_join(df[, c("trna", "sprinzl")], mods, by = join_by(trna, seq_position == pos)) |>
      rename(pos = seq_position)
  }

  contained_mods <- contained_mods |>
    dplyr::filter(!is.na(mod))
  
  return(contained_mods)
}

add_contained_mods <- function(df, contained_mods, use_sprinzl) {
  if (use_sprinzl) {
    df <- dplyr::left_join(df, contained_mods, by = join_by(trna, sprinzl == pos))
  } else {
    df <- dplyr::left_join(df, contained_mods, by = join_by(trna, seq_position == pos))
  }
  
  i <- is.na(df$mods)
  if (any(i)) {
    df$mods[i] <- ""
  }
  
  return(df)
}

################################################################################

df <- read.table(opts$args, header = TRUE, sep = "\t")
mods <- read.table(opts$options$mods, header = TRUE, sep = "\t", quote = "", comment.char = "")

contained_mods <- get_contained_mods(df, opts$options$mods, opts$options$sprinzl) 
if (nrow(contained_mods) > 0) {
  contained_mods <- contained_mods |>
    summarise(mods = paste(mod, collapse = " "),
              .by = c(trna, pos))
  df <- add_contained_mods(df, contained_mods, opts$options$sprinzl)
}

write.table(df,
            opts$options$output,
            quote=TRUE, sep = "\t",
            col.names = TRUE, row.names = FALSE)