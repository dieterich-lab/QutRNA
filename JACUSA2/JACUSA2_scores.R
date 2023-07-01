#!/usr/bin/env Rscript

# TODO remove
options(error = traceback)
options(show.error.locations = TRUE)

#readr
#library(optparse)
#library(JACUSA2helper)
#library(GenomicRanges)
#library(IRanges)
#library(plyranges)
#library(dplyr)
#library(BSgenome)
#library(data.table)
library(magrittr)

WIDTH <- 5

option_list <- list(
  optparse::make_option(c("-f", "--fasta"),
                        type = "character",
                        help = "FASTA sequence"),
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

fasta <- Biostrings::readDNAStringSet(opts$options$fasta)
mods <- NULL
if (!is.null(opts$options$mods)) {
mods <- data.table::fread(opts$options$mods, header = TRUE) %>%
  as.data.frame() %>%
  dplyr::rename(seqnames = Chr, start = Position) %>%
  dplyr::mutate(end = start) %>%
  GenomicRanges::GRanges()
}

result <- JACUSA2helper::read_result(opts$args, unpack = TRUE) %>%
  IRanges::shift(-1) %>%
  plyranges::select(bases, score, ins, insertion_score, del, deletion_score, ref) %>%
  plyranges::mutate(
    Ref = seqnames,
    Pos3 = start,
    Ref_Pos = paste0(seqnames, "_", Pos3))

if (is.null(mods)) {
  result <- result %>%
    plyranges::mutate(Mod)
} else {
  result <- result %>%
    plyranges::join_overlap_left(mods)
}

GenomeInfoDb::seqlevels(result) <- GenomeInfoDb::seqlevels(fasta)
GenomicRanges::seqinfo(result) <- GenomicRanges::seqinfo(fasta)

GenomicRanges::mcols(result, level = "within")[, "Mis"] <- GenomicRanges::mcols(result, level = "within")[, "score"]
GenomicRanges::mcols(result, level = "within")[, "Del"] <- GenomicRanges::mcols(result, level = "within")[, "deletion_score"]
GenomicRanges::mcols(result, level = "within")[, "Ins"] <- GenomicRanges::mcols(result, level = "within")[, "insertion_score"]

result <- result %>%
  plyranges::mutate(
    Mis = replace(Mis, is.na(Mis), 0),
    Del = replace(Del, is.na(Del), 0),
    Ins = replace(Ins, is.na(Ins), 0),
    Mod = replace(Mod, is.na(Mod), ""),
  )

# add context scores
result <- result %>%
  dplyr::mutate(MisDelIns = Mis + Del + Ins) %>%
  IRanges::resize(width = WIDTH, fix = "center") %>%
  IRanges::shift(1) %>%
  plyranges::filter(start > 0 & end < GenomeInfoDb::seqlengths(.)[as.character(GenomeInfoDb::seqnames(.))]) %>%
  plyranges::mutate(
    Kmer = BSgenome::getSeq(fasta, .) %>% as.character(),
  ) %>%
  IRanges::shift(-1) %>%
  IRanges::resize(width = IRanges::width(.) + 1) %>%
  plyranges::mutate(
    Context = paste0(start , "_", end)
  )

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
  for (cond in c(1:ncol(df))) {
    for (repl in c(1:ncol(df[[paste0("cond", cond)]]))) {
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
    dplyr::select(Ref, Context, Ref_Pos, Pos3,
                  Mod,
                  Kmer,
                  Mis, MisDelIns),
  bases,
  ins,
  del
)


# final polish
GenomicRanges::mcols(result) <- df
df <- GenomicRanges::mcols(result) %>%
  as.data.frame() %>%
  dplyr::select(
    Ref, Context, Ref_Pos, Pos3,
    Mis, MisDelIns,
    Mod,
    Kmer,
    dplyr::starts_with("coverage_"),
    dplyr::starts_with("mismatch_rate_"),
    dplyr::starts_with("bases_"),
    dplyr::starts_with("insertion_"),
    dplyr::starts_with("deletion_")
  ) %>%
  dplyr::rename(`Mis+Del+Ins` = MisDelIns)

# good to go to files
common_cols <- c("Ref", "Context", "Ref_Pos", "Pos3")
common_cols <- c(common_cols, "Mis", "Mis+Del+Ins",
                 "Mod",
                 "Kmer")

df %>% dplyr::select(dplyr::all_of(common_cols),
                     dplyr::starts_with("coverage_"),
                     dplyr::starts_with("mismatch_rate"),
                     dplyr::starts_with("bases_"),
                     dplyr::starts_with("insertion_"),
                     dplyr::starts_with("deletion_")) %>%
  write.table(opts$options$output, quote=FALSE, sep = ",", col.names = TRUE, row.names = FALSE)
