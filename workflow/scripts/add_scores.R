#!/usr/bin/env Rscript

# Process JACUSA2 scores and add additional data.
library(magrittr)
library(GenomicRanges)

# TODO move this part to JACUSA2helper
# process scores
# add kmer in JACUS2Ahelper

WIDTH <- 5

option_list <- list(
  optparse::make_option(c("-f", "--fasta"),
                        type = "character",
                        help = "FASTA sequence"),
  optparse::make_option(c("-n", "--normalize_score"),
                        action = "store_true",
                        help = "Normalize score with subsampled data",
                        default = FALSE),
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
stopifnot(length(opts$args) == 1)
fasta <- Biostrings::readDNAStringSet(opts$options$fasta)

# TODO normalize by sampled scores

common_cols <- c("Ref", "Pos3", "Mis", "Mis+Del+Ins", "Kmer")
if (opts$options$normalize_score) {
  result <- JACUSA2helper::read_result(opts$args, unpack = c("score_subsampled", "insertion_score_subsampled", "deletion_score_subsampled", "insertion_score", "deletion_score", "reads", "insertion_ratio", "deletion_ratio", "non_ref_ratio", "base_ratio")) |>
    IRanges::shift(-1) |>
    plyranges::select(bases, score,
                      score_subsampled,
                      insertion_score, insertion_score_subsampled, insertion_ratio,
                      deletion_score, deletion_score_subsampled, deletion_ratio,
                      reads,ref, non_ref_ratio, base_ratio)
  common_cols <- c(common_cols, "mean_delta_MisDelIns_subsampled")
} else {
  result <- JACUSA2helper::read_result(opts$args, unpack = c("insertion_score", "deletion_score", "reads")) %>%
    IRanges::shift(-1) %>%
    plyranges::select(bases, score,
                      ins, insertion_score,
                      del, deletion_score, reads, ref) %>%
    plyranges::mutate(Ref = seqnames,
                      Pos3 = start)
}


cond_repl <- c(names(result$bases$cond1) %>% gsub("rep", "", .) %>% paste0("1", .),
               names(result$bases$cond2) %>% gsub("rep", "", .) %>% paste0("2", .))

add_value <- function(result, key, label) {
  if (key %in% colnames(mcols(result))) {
    values <- lapply(data.table::tstrsplit(mcols(result)[, key], split = ",", fixed = TRUE), as.numeric) %>%
        as.data.frame()
    mcols(result)[, key] <- NULL
    colnames(values) <- paste0(label, cond_repl)
    common_cols <- c(common_cols, colnames(values))
    mcols(result) <- dplyr::bind_cols(mcols(result) |> as.data.frame(), values)
  }

  result
}

result <- add_value(result, "reads", "coverage_")
result <- add_value(result, "insertion_ratio", "insertion_ratio_")
browser()
result <- add_value(result, "deletion_ratio", "deletion_ratio_")
#result <- add_value(result, "non_ref_ratio", "non_ref_ratio_")

GenomeInfoDb::seqlevels(result) <- GenomeInfoDb::seqlevels(fasta)
seqinfo(result) <- seqinfo(fasta)

mcols(result, level = "within")[, "Mis"] <-
  mcols(result, level = "within")[, "score"]
mcols(result, level = "within")[, "Del"] <-
  mcols(result, level = "within")[, "deletion_score"]
mcols(result, level = "within")[, "Ins"] <-
  mcols(result, level = "within")[, "insertion_score"]

result <- result %>%
  plyranges::mutate(
    Mis = tidyr::replace_na(Mis, 0),
    Del = tidyr::replace_na(Del, 0),
    Ins = tidyr::replace_na(Ins, 0),
  )

# add kmer TODO add to JACUSA2helper
suppressWarnings({
  result <- result %>%
    plyranges::mutate(MisDelIns = Mis + Del + Ins) %>%
    IRanges::resize(width = WIDTH, fix = "center") %>%
    IRanges::shift(1) %>%
    plyranges::filter(start > 0 &
                      end < GenomeInfoDb::seqlengths(.)[as.character(GenomeInfoDb::seqnames(.))]) %>%
    plyranges::mutate(Kmer = BSgenome::getSeq(fasta, .) %>% as.character()) %>%
    IRanges::shift(-1) %>%
    IRanges::resize(width = IRanges::width(.) + 1)
})

if (opts$options$normalize_score) {
  unpack <- function(v) {
    strsplit(v, ",") |>
      lapply(function(x) { tidyr::replace_na(x, 0) } ) |>
      lapply(as.numeric)
  }

  score_subsampled <- unpack(result$score_subsampled)
  delta_score_subsampled <- mapply(`-`, result$score, score_subsampled, SIMPLIFY=FALSE)
  mean_delta_score_subsampled <- lapply(delta_score_subsampled, mean) |> unlist()

  insertion_score_subsampled <- unpack(result$insertion_score_subsampled)
  delta_insertion_score_subsampled <- mapply(`-`, result$insertion_score, insertion_score_subsampled, SIMPLIFY=FALSE)
  mean_delta_insertion_score_subsampled <- lapply(delta_insertion_score_subsampled, mean) |> unlist()

  deletion_score_subsampled <- unpack(result$deletion_score_subsampled)
  delta_deletion_score_subsampled <- mapply(`-`, result$deletion_score, deletion_score_subsampled, SIMPLIFY=FALSE)
  mean_delta_deletion_score_subsampled <- lapply(delta_deletion_score_subsampled, mean) |> unlist()

  common_cols <- c(common_cols,
                   "score", "score_subsampled",
                   "insertion_score", "insertion_score_subsampled",
                   "deletion_score", "deletion_score_subsampled")

  result <- result |>
    plyranges::mutate(Ref = seqnames,
                      Pos3 = start,
                      mean_delta_MisDelIns_subsampled = mean_delta_score_subsampled +
                      mean_delta_insertion_score_subsampled +
                      mean_delta_deletion_score_subsampled)
  i <- result$mean_delta_deletions_score_subsampled < 0 | is.na(result$mean_delta_deletions_score_subsampled)
  if (any(i)) {
    mcols(result)[i, "mean_delta_MisDelIns_subsampled"] <- 0
  }
}

df <- GenomicRanges::mcols(result) |>
  as.data.frame() |>
  dplyr::rename(`Mis+Del+Ins` = MisDelIns)

# good to go to files
df <- df %>% dplyr::select(dplyr::all_of(common_cols),
                           dplyr::starts_with("coverage_"))

write.table(df, opts$options$output, quote=FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
