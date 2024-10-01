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
  result <- JACUSA2helper::read_result(opts$args, unpack = TRUE) %>%
    IRanges::shift(-1) %>%
    plyranges::select(bases, score,
                      score_subsampled,
                      ins, insertion_score, insertion_score_subsampled,
                      del, deletion_score, deletions_score_subsampled,
                      ref) |>
    plyranges::mutate(Ref = seqnames,
                      Pos3 = start,
                      score_subsampled =
                        strsplit(score_subsampled,
                                 sep = ",") |> unlist() |> tidyr::replace_na(0),
                      insertion_score_subsampled =
                        strsplit(insertion_score_subsampled,
                                 sep = ",") |> unlist() |> tidyr::replace_na(0),
                      deletion_score_subsampled =
                        strsplit(deletions_score_subsampled,
                                 sep = ",") |> unlist() |> tidyr::replace_na(0))
  common_cols <- c(common_cols, MisDelIns_subsampled)
} else {
  result <- JACUSA2helper::read_result(opts$args, unpack = TRUE) %>%
    IRanges::shift(-1) %>%
    plyranges::select(bases, score,
                      ins, insertion_score,
                      del, deletion_score, ref) %>%
    plyranges::mutate(Ref = seqnames,
                      Pos3 = start)
}

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
    Mis = replace_na(Mis, 0),
    Del = replace_na(Del, 0),
    Ins = replace_na(Ins, 0),
  )

# add kmer TODO add to JACUSA2helper
suppressWarnings({
  result <- result |>
    dplyr::mutate(MisDelIns = Mis + Del + Ins) |>
    IRanges::resize(width = WIDTH, fix = "center") |>
    IRanges::shift(1) %>%
    plyranges::filter(start > 0 &
                      end < GenomeInfoDb::seqlengths(.)[as.character(GenomeInfoDb::seqnames(.))]) %>%
    plyranges::mutate(Kmer = BSgenome::getSeq(fasta, .) |> as.character()) |>
    IRanges::shift(-1) |>
    IRanges::resize(width = IRanges::width(.) + 1)
})

# TODO move to JACUSA2
non_ref <- function(bases, ref) {
  stopifnot(length(ref) ==  nrow(bases))
  n <- length(ref)

  cov <- rowSums(bases)
  non_zero_i <- cov > 0

  non_ref <- rep(0, n)
  # calculate ratio only for sites with cov > 0 -> x / cov

  # matrix indexing ref base in base call matrix
  ref_i <- as.matrix(
    cbind(
      seq_len(n),
      match(ref, JACUSA2helper:::.BASES)
    )
  )
  # coverage = ref bc + non-ref bc (base_calls[i] corresponds to ref bc)
  non_ref[non_zero_i] <- (cov[non_zero_i] - as.matrix(bases)[ref_i][non_zero_i])

  non_ref
}

bases_to_cols <- function(df, ref) {
  l <- list()
  for (cond in c(1:ncol(df))) {
    for (repl in c(1:ncol(df[[paste0("cond", cond)]]))) {
      k <- paste("mismatch_rate", cond, repl, sep = "_")
      bases <- df[[paste0("cond", cond)]][[paste0("rep", repl)]]
      l[[k]] <- non_ref(bases, ref) / rowSums(bases)
      k <- paste("coverage", cond, repl, sep = "_")
      l[[k]] <- rowSums(bases)
      for (b in JACUSA2helper:::.BASES) {
        k <- paste("bases", cond, repl, b, sep = "_")
        l[[k]] <- bases[[b]]
      }
    }
  }
  as.data.frame(l)
}

bases <- bases_to_cols(result$bases, result$ref)

indel_to_cols <- function(df, prefix) {
  l <- list()
  for (cond in seq_len(ncol(df))) {
    for (repl in seq_len(ncol(df[[paste0("cond", cond)]]))) {
      k <- paste(prefix, "rate", cond, repl, sep = "_")
      reads <-  df[[paste0("cond", cond)]][[paste0("rep", repl)]]$reads
      coverage <-  df[[paste0("cond", cond)]][[paste0("rep", repl)]]$coverage
      rate <- reads / coverage
      rate[is.na(rate)] <- 0
      l[[k]] <- rate
    }
  }

  as.data.frame(l)
}

ins <- indel_to_cols(result$ins, "insertion")
del <- indel_to_cols(result$del, "deletion")

df <- dplyr::bind_cols(
  GenomicRanges::mcols(result) %>%
    as.data.frame() %>%
    dplyr::select(Ref, Pos3,
                  Kmer,
                  Mis, MisDelIns),
  bases,
  ins,
  del
)

# TODO move to JACUSA2
# final polish
GenomicRanges::mcols(result) <- df
df <- GenomicRanges::mcols(result) %>%
  as.data.frame() %>%
  dplyr::select(
    Ref, Pos3,
    Mis, MisDelIns,
    Kmer,
    dplyr::starts_with("coverage_"),
    dplyr::starts_with("mismatch_rate_"),
    dplyr::starts_with("bases_"),
    dplyr::starts_with("insertion_"),
    dplyr::starts_with("deletion_")
  ) %>%
  dplyr::rename(`Mis+Del+Ins` = MisDelIns)

# good to go to files

df %>% dplyr::select(dplyr::all_of(common_cols),
                     dplyr::starts_with("coverage_"),
                     dplyr::starts_with("mismatch_rate"),
                     dplyr::starts_with("bases_"),
                     dplyr::starts_with("insertion_"),
                     dplyr::starts_with("deletion_")) %>%
  write.table(opts$options$output, quote=FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
