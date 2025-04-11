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
  # args = c("-s", "MDI::mismatch_score+deletion_score+insertion_score,MDI_subsampled::norm_mismatch_score_subsampled+norm_deletion_score_subsampled+norm_insertion_score_subsampled",
  #         "-f", "data/ref.fasta",
  #         "-o", "~/tmp/test.tsv",
  #         "/beegfs/prj/tRNA_Berlin/isabel/20240820_AJ_Ecoli_ctrls_tRNA_RNA002_cust/qutrna/queF_2__vs__NC_plus_N3/output-qutrna2/results/jacusa2/cond1~queF_2/cond2~NC_plus_N3/JACUSA2.out"),
  positional_arguments = TRUE
)

stopifnot(!is.null(opts$options$fasta))
stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$options$stat))
stopifnot(length(opts$args) == 1)

summarise_ratio <- function(r, get_ratios, f) {
  by_conds <- split(colData(r), colData(r)$condition) |>
    lapply(rownames)
  cols1 <- by_conds[[names(by_conds)[1]]]
  cols2 <- by_conds[[names(by_conds)[2]]]

  ratios <- get_ratios(r)
  if (length(cols1) > 1) {
    ratios1 <- f(ratios[, cols1])
  } else {
    ratios1 <- ratios[, cols1]
  }
  if (length(cols2) > 1) {
    ratios2 <- f(ratios[, cols2])
  } else {
    ratios2 <- ratios[, cols2]
  }

  ratios1 - ratios2
}


norm_score <- function(r, score, score_sampled) {
  scores <- rowData(r)[[score]]
  scores_subsampled <- rowData(r)[[score_sampled]]
  runs <- lapply(scores_subsampled, length) |>
    unlist() |>
    unique()
  # stopifnot(length(runs) == 1 || runs == 0)

  new_scores <- mapply(function(observed, subsampled) { return(observed - max(subsampled, na.rm = TRUE)) },
               as.list(scores), scores_subsampled)

  i <- new_scores < 0
  if (any(i)) {
    new_scores[i] <- 0
  }

  new_scores
}

summarise_score <- function(r, score, f) {
  scores <- rowData(r)[[score]]
  runs <- lapply(scores, length) |>
    unlist() |>
    unique()
  stopifnot(length(runs) == 1 || runs == 0)

  scores <- lapply(scores, f) |>
    unlist()

  scores
}

PARSE_COLS <- list("mismatch_score" = function(r) { return(rowData(r)$score) },
                   "insertion_score" = function(r) { return(rowData(r)$insertion_score) },
                   "deletion_score" = function(r) { return(rowData(r)$deletion_score) },
                   #
                   "norm_mismatch_score_subsampled" = function(r) { return(norm_score(r, "score", "score_subsampled")) },
                   "norm_insertion_score_subsampled" = function(r) { return(norm_score(r, "insertion_score", "insertion_score_subsampled")) },
                   "norm_deletion_score_subsampled" = function(r) { return(norm_score(r, "deletion_score", "deletion_score_subsampled")) },
                   #
                   "norm_mismatch_score_downsampled" = function(r) { return(norm_score(r, "score", "score_downsampled")) },
                   "norm_insertion_score_downsampled" = function(r) { return(norm_score(r, "insertion_score", "insertion_score_downsampled")) },
                   "norm_deletion_score_subsampled" = function(r) { return(norm_score(r, "deletion_score", "deletion_score_downsampled")) },
                   #
                   "median_mismatch_score_downsampled" = function(r) { return(summarise_score(r, "score_downsampled", median)) },
                   "median_insertion_score_downsampled" = function(r) { return(summarise_score(r, "insertion_score_downsampled", median)) },
                   "median_deletion_score_downsampled" = function(r) { return(summarise_score(r, "deletion_score_downsampled", median)) },
                   #
                   "mean_mismatch_score_downsampled" = function(r) { return(summarise_score(r, "score_downsampled", mean)) },
                   "mean_insertion_score_downsampled" = function(r) { return(summarise_score(r, "insertion_score_downsampled", mean)) },
                   "mean_deletion_score_downsampled" = function(r) { return(summarise_score(r, "deletion_score_downsampled", mean)) },
                   #
                   "mean_non_ref_ratio" = function(r) { return(summarise_ratio(r, non_ref_ratio, rowMeans)) },
                   "mean_insertion_ratio" = function(r) { return(summarise_ratio(r, insertion_ratio, rowMeans)) } ,
                   "mean_deletion_ratio" = function(r) { return(summarise_ratio(r, deletion_ratio, rowMeans)) })

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
    dplyr::rename(trna = seqnames, pos = start) |>
    dplyr::select(trna, pos, strand, ref, ref_context)

  # add intermediate columns
  reads <- replace_prefix(assays(r)$reads, r, "coverage")
  cols <- c(cols, colnames(reads))
  df <- cbind(df, reads)
  for (col in pick_cols(stats$stat)) {
    if (!col %in% colnames(df)) {
      f <- PARSE_COLS[[col]]
      if (is.null(f)) {
        print(col)
      }
      df[, col] <- f(r)
    }
  }

  # score_cols <- c("score", "score_subsampled", "score_downsampled",
  #                 "deletion_score", "deletion_score_subsampled", "deletion_score_downsampled",
  #                 "insertion_score", "insertion_score_subsampled", "insertion_score_downsampled")
  #for (col in score_cols) {
  #  df[, col] <- rowRanges(r)[, col]
  #}

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

res_to_r <- function(res) {
  rowRanges <- res[, c("score", "filter", "ref",
                       "insertion_score", "deletion_score",
                       "score_subsampled", "deletion_score_subsampled", "insertion_score_subsampled",
                       "score_downsampled", "deletion_score_downsampled", "insertion_score_downsampled")]
  mcols(rowRanges)[, "ref_context"] <- "N"

  assays <- list()
  bases <- res$bases
  l <- list()
  condition <- c()
  for (cond in names(bases)) {
    s <- gsub("cond", "bam_", cond)
    r <- length(bases[[cond]])
    b <- bases[[cond]]
    s <- paste0(s, "_", seq_len(r))
    names(b) <- s
    l <- append(l, b)
    condition <- c(condition, rep(cond, r))
  }
  assays[["bases"]] <- l |>
    tibble::as_tibble()

  colData <- data.frame(condition = condition)
  row.names(colData) <- names(assays[["bases"]])

  scores <- c("score", "insertion_score", "deletion_score")
  scores <- c(paste0(scores, "_subsampled"),
              paste0(scores, "_downsampled"))
  for (score in scores) {
    mcols(rowRanges)[[score]] <- strsplit(mcols(rowRanges)[[score]], ",") |>
        lapply(as.numeric) |>
        lapply(unlist) |>
        clean_score()
  }

  l <- strsplit(res$reads, split = ",", fixed = TRUE) |>
    lapply(as.numeric)
  d <- do.call(rbind, l)
  colnames(d) <- rownames(colData)
  assays[["reads"]] <- as.data.frame(d)

  metadata <- list()
  metadata[["command"]] <- "UNKNOWN"

  se <- SummarizedExperiment(assays = assays,
                             rowRanges = rowRanges,
                             colData = colData,
                             metadata = metadata)

  se
}

# load JACUSA2 output and add... reference context
res <- JACUSA2helper::read_result(opts$args,
                                  unpack = c("score_subsampled", "score_downsampled",
                                             "reads",
                                             "insertion_score", "deletion_score",
                                             "insertion_score_subsampled", "deletion_score_subsampled",
                                             "deletion_score_downsampled", "insertion_score_downsampled"))

r <- res_to_r(res)



# parse options and create output container
stats <- parse_stats(opts$options$stat)
# calculate stats and create output
df <- parse_result(r, stats)
df$Ref <- df$trna
df$Pos3 <- df$pos

write.table(df,
            opts$options$output,
            quote = FALSE,
            sep = "\t",
            col.names = TRUE, row.names = FALSE)
