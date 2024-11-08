#!/usr/bin/env Rscript

# Process JACUSA2 scores and add additional data.
library(GenomicRanges)
library(JACUSA2helper)
library(SummarizedExperiment)

WIDTH <- 5

option_list <- list(
  optparse::make_option(c("-f", "--fasta"),
                        type = "character",
                        help = "FASTA sequence"),
  optparse::make_option(c("-w", "--width"),
                        type = "numeric",
                        help = paste0("Width of sequence context. Default: ",
                                      WIDTH),
                        default = WIDTH),
  optparse::make_option(c("-s", "--stat"),
                        help = "Define statistics"),
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output")
)

opts <- optparse::parse_args(
  optparse::OptionParser(option_list = option_list),
  #args = c("--fasta=../data/S.pombe.fasta",
  #         "S.pombe_tRNAAsp_IVT-Q/JACUSA2.out"),
  positional_arguments = TRUE
)

stopifnot(!is.null(opts$options$fasta))
stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$options$stat))
stopifnot(length(opts$args) == 1)

summarise_ratio <- function(r, f) {
  by_conds <- split(colData(r), colData(r)$condition) |>
    lapply(rownames)
  cols1 <- by_conds[[names(by_conds)[1]]]
  cols2 <- by_conds[[names(by_conds)[2]]]

  ratios <- f(r)
  if (length(cols1) > 1) {
    ratios1 <- rowMeans(ratios[, cols1])
  } else {
    ratios1 <- ratios[, cols1]
  }
  if (length(cols2) > 1) {
    ratios2 <- rowMeans(ratios[, cols2])
  } else {
    ratios2 <- ratios[, cols2]
  }

  ratios1 - ratios2
}


norm_score <- function(r, score) {
  scores <- rowData(r)[[score]]
  scores_subsampled <- rowData(r)[[paste0("score_subsampled")]] |>
    lapply(median) |>
    unlist()

  i <- scores <= scores_subsampled
  if (any(i)) {
    scores[i] <- 0
  }

  scores
}

PARSE_COLS <- list("mismatch_score" = function(r) { return(rowData(r)$score) },
                   "insertion_score" = function(r) { return(rowData(r)$insertion_score) },
                   "deletion_score" = function(r) { return(rowData(r)$deletion_score) },
                   #
                   "norm_mismatch_score" = function(r) { return(norm_score(r, "score")) },
                   "norm_insertion_score" = function(r) { return(norm_score(r, "insertion_score")) },
                   "norm_deletion_score" = function(r) { return(norm_score(r, "deletion_score")) },
                   #
                   "non_ref_ratio" = function(r) { return(summarise_ratio(r, non_ref_ratio)) },
                   "insertion_ratio" = function(r) { return(summarise_ratio(r, insertion_ratio)) } ,
                   "deletion_ratio" = function(r) { return(summarise_ratio(r, deletion_ratio)) })

clean_score <- function(score) {
  tidyr::replace_na(score, 0)
}

parse_stats <- function(stat_opts) {
  l <- strsplit(stat_opts, ",") |>
    unlist() |>
    strsplit("::")
  df <- do.call(rbind, l) |>
    as.data.frame()
  colnames(df) <- c("new_col", "stat")
  df[, "stat"] <- gsub("norm:", "norm_", df[, "stat"])

  df
}

pick_cols <- function(stats) {
  stats <- strsplit(stats, split = "\\+") |>
    unlist()

  gsub("^norm:", "", stats) |>
    unique()
}

replace_prefix <- function(df, r, new_prefix) {
  new_cols <- colData(r) |>
    as.data.frame() |>
    dplyr::mutate(.by = condition,
                  new_prefix = paste0(new_prefix, "_", dplyr::cur_group_id(), "_", 1:dplyr::n())) |>
    dplyr::pull(new_prefix)

  df <- df[, rownames(colData(r))]
  colnames(df) <- new_cols

  df
}

parse_result <- function(r, stats) {
  cols <- c("ref", "ref_context")
  df <- rowRanges(r)[, c("ref", "ref_context")] |>
    as.data.frame() |>
    dplyr::rename(trna = seqnames,
                  pos = start) |>
    dplyr::select(trna, pos, strand, ref, ref_context)

  # add intermediate columns
  reads <- replace_prefix(assays(r)$reads, r, "reads")
  cols <- c(cols, colnames(reads))
  df <- cbind(df, reads)
  for (col in pick_cols(stats$stat)) {
    if (!col %in% colnames(df)) {
      f <- PARSE_COLS[[col]]
      df[, col] <- f(r)
    }
  }

  # add final columns
  for (i in rownames(stats)) {
    new_col <- stats[i, "new_col"]
    stat <- stats[i, "stat"]
    if (!new_col %in% colnames(df)) {
      df[, new_col] <- with(df, eval(parse(text = stat)))
    }
  }
  # keep what is needed
  cols <- c("trna", "pos", "strand", cols, stats$new_col)
  df <- df[, cols]

  df
}

# load JACUSA2 output and add... reference context
r <- JACUSA2helper::import_result(opts$args)
rowData(r) <- rowData(r) |>
  as.data.frame()
fasta <- Biostrings::readDNAStringSet(opts$options$fasta)
GenomeInfoDb::seqlevels(r) <- GenomeInfoDb::seqlevels(fasta)
seqinfo(r) <- seqinfo(fasta)
r <- JACUSA2helper::add_ref_context(r, fasta, width = opts$options$width)

# parse options and create output container
stats <- parse_stats(opts$options$stat)
# calculate stats and create output
df <- parse_result(r, stats)

write.table(df,
            opts$options$output,
            quote=FALSE,
            sep = "\t",
            col.names = TRUE, row.names = FALSE)
