#!/usr/bin/env Rscript

# Process JACUSA2 scores and add additional data.
library(GenomicRanges)
library(dplyr)
library(SummarizedExperiment)

##

.EMPTY <- "*"
.BASES <- c("A", "C", "G", "T")

.find_conditions <- function(cond_repl) {
  cond_count <- 1
  max_conds <- length(cond_repl)

  while (cond_count <= max_conds) {
    cond_count_nchar <- nchar(cond_count)
    cond <- substring(cond_repl, first = 1, last = cond_count_nchar)
    repl <- substring(cond_repl, first = cond_count_nchar + 1)
    if (all(nchar(cond) > 0 & nchar(repl) > 0) &&
        length(unique(cond)) == cond_count) {
      return(cond_count)
    }
    cond_count <- cond_count + 1
  }

  return(NULL)
}

guess_sample_info <- function(raw_header_names) {
  prefixes <- c("bases", "arrest_bases")

  for (prefix in prefixes) {
    prefix_regex <- paste0("^", prefix, "[0-9]+")
    i <- grepl(prefix_regex, raw_header_names)
    conditions <- .find_conditions(gsub(prefix, "", raw_header_names[i]))
    if (is.null(conditions)) {
      next
    }
    prefix_regex = paste0("^(", paste0(prefix, collapse = "|"), ")")
    condition_regex = paste0("([0-9]{", nchar(conditions), "})")
    replicate_regex = "([0-9]+)"
    m <- do.call(
      rbind,
      regmatches(raw_header_names[i],
                 regexec(paste0(prefix_regex,
                                condition_regex,
                                replicate_regex),
                         raw_header_names[i]))
    )
    df <- as.data.frame(m)
    colnames(df) <- c("matched_column", "matched_prefix", "condition", "replicate")
    df$condition <- as.integer(df$condition)
    df$replicate <- as.integer(df$replicate)
    df$sample <- paste0("bam_", df$condition, "_", df$replicate)

    return(df)
  }

  stop("Could not guess 'sample_info'!")
}

.fill_empty <- function(df, cols, new_cols) {
  new_value <- paste0(rep(0, length(new_cols)), collapse = ",")

  for (col in cols) {
    i <- df[, col] == .EMPTY | df[, col] == ""
    if (length(i)) {
      df[i, col] <- new_value
    }
  }

  return(df)
}

.unpack <- function(s, new_cols, sep = ",") {
  l <- lapply(data.table::tstrsplit(s, split = sep, fixed = TRUE), as.numeric)
  names(l) <- new_cols

  return(tibble::as_tibble(l))
}

.unpack_cols <- function(df, cols, samples, new_cols) {
  df <- .fill_empty(df, cols, new_cols)
  unpacked <- lapply(df[cols], .unpack, new_cols = new_cols)
  names(unpacked) <- samples

  return(tibble::as_tibble(unpacked))
}

# read extended JACUSA2 output and convert to data.frame
read_result <- function(file, tmpdir = tempdir(), showProgress = FALSE, ...) {
  # if url download first
  if (grepl("ftp://|http://", file)) {
    # partly adopted from data.table::fread
    tmpFile <- tempfile(tmpdir = tmpdir)
    if (!requireNamespace("curl", quietly = TRUE))
      stop("Input URL requires https:// connection for which fread() requires 'curl' package which cannot be found. Please install 'curl' using 'install.packages('curl')'.")

    tmpFile = tempfile(fileext = paste0(".", tools::file_ext(file)), tmpdir = tmpdir)
    curl::curl_download(file, tmpFile, mode = "wb", quiet = !showProgress)
    file = tmpFile
    on.exit(unlink(file), add = TRUE)
  }

  # pre process file
  # parse comments ^(#|##) to determine result/method type and number of conditions
  fun <- ifelse(grepl('\\.gz$', file), base::gzfile, base::file)
  con <- fun(file, "r")
  skip_lines <- 0
  header_names <- NULL
  jacusa_header <- c()
  while (TRUE) {
    line = readLines(con, n = 1)
    # quit reading: nothing to read or first no header line
    if (length(line) == 0 || length(grep("^#", line)) == 0) {
      break
    }
    # count header lines to ignore
    skip_lines <- skip_lines + 1

    if (length(grep("^#contig", line)) > 0) {
      # parse and store header
      # fix header: #contig -> contig
      header_names <- sub("^#", "", line);
      header_names <- unlist(base::strsplit(header_names, split = "\t", fixed = TRUE))
    } else if (length(grep("^##", line)) > 0) {
      jacusa_header <- c(gsub("^##", "", line), jacusa_header)
    }

  }
  # finished pre-processing
  close(con)

  # check that a header could be parsed
  if (is.null(header_names)) {
    stop("No header line for file: ", file)
  }

  # read data
  data <- data.table::fread(
    file,
    skip = skip_lines,
    sep = "\t",
    header = FALSE,
    ...
  )
  colnames(data) <- header_names
  # convert to numeric
  i <- data[, 5] == "*"
  if (any(i)) {
    data[i, 5] <- NA
    data[, 5] <- as.numeric(data[, 5])
  }

  return(as.data.frame(data))
}

.to_se <- function(df) {
  header_names <- colnames(df)
  sample_info <- guess_sample_info(colnames(df))
  assay_cols <-
  assays <- list()
  # add bases as an assay
  cols <- paste0("bases", sample_info$condition, sample_info$replicate)
  assays[["bases"]] <- .unpack_cols(df, cols, sample_info$sample, .BASES)
  assay_cols <- c(cols)

  for (prefix in c("reads", "insertions", "deletions")) {
    cols <- paste0(prefix, sample_info$condition, sample_info$replicate)
    assays[[prefix]] <- df[, cols] |>
      lapply(as.integer) |>
      as.data.frame()
    colnames(assays[[prefix]]) <- sample_info$sample
    assay_cols <- c(assay_cols, cols)
  }
  for (prefix in c("nonref_ratio", "insertion_ratio", "deletion_ratio")) { # FIXME add more keys
    cols <- paste0(prefix, sample_info$condition, sample_info$replicate)
    assays[[prefix]] <- df[, cols] |>
      lapply(as.numeric) |>
      as.data.frame()
    colnames(assays[[prefix]]) <- sample_info$sample
    assay_cols <- c(assay_cols, cols)
  }

  colData <- data.frame(condition = paste0("condition~", sample_info$condition))
  row.names(colData) <- names(sample_info$sample)

  # process scores
  scores <- c("score", "insertion_score", "deletion_score")
  scores <- c(paste0(scores, "_subsampled"))
  for (score in intersect(colnames(df), scores)) {
    df[[score]] <- strsplit(df[[score]], ",") |>
      lapply(as.numeric) |>
      lapply(unlist) |>
      clean_score()
  }

  metadata <- list()
  metadata[["command"]] <- "UNKNOWN"

  gr <- GenomicRanges::GRanges(
    seqnames = df$contig,
    ranges = IRanges::IRanges(start = df$start + 1, end = df$end),
    strand = GenomicRanges::strand(gsub("\\.", "*", df$strand))
  )
  df[c("contig", "start", "end", "strand", assay_cols)] <- NULL
  GenomicRanges::mcols(gr) <- df

  se <- SummarizedExperiment(assays = assays,
                             rowRanges = gr,
                             colData = colData,
                             metadata = metadata)

  se
}



# print missing keys

##

option_list <- list(
  optparse::make_option(c("-s", "--stat"),
                        help = "Define statistics"),
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output")
)

opts <- optparse::parse_args(
  optparse::OptionParser(option_list = option_list),
  #args = c("-s", "MDI::mismatch_score+deletion_score+insertion_score,MDI_subsampled::norm_mismatch_score_subsampled+norm_deletion_score_subsampled+norm_insertion_score_subsampled",
  #         "-o", "~/tmp/test.tsv",
  #         "/beegfs/prj/tRNA_Francesca_Tuorto/data/20250115_FT_HCT116_tRNA_RNA004/qutrna/mt_trnas/test_jacusa2_local_hac500/results/jacusa2/cond1~WT/cond2~DNMT2/JACUSA2.out"),
  positional_arguments = TRUE
)

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
  # runs <- lapply(scores_subsampled, length) |>
  #   unlist() |>
  #   unique()
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

# mapping to functions
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
                   "mean_deletion_score_downsampled" = function(r) { return(summarise_score(r, "deletion_score_downsampled", mean)) })
                   # TODO check
                   #"mean_non_ref_ratio" = function(r) { return(summarise_ratio(r, non_ref_ratio, rowMeans)) },
                   #"mean_insertion_ratio" = function(r) { return(summarise_ratio(r, insertion_ratio, rowMeans)) } ,
                   #"mean_deletion_ratio" = function(r) { return(summarise_ratio(r, deletion_ratio, rowMeans)) })

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

parse_result <- function(r, stats) {
  cols <- c("ref")
  df <- rowRanges(r)[, c("ref")] |>
    as.data.frame() |>
    rename(trna = seqnames, seq_position = start) |>
    select(all_of(c("trna", "seq_position", "strand", "ref")))

  # add intermediate columns
  for (a in c("reads", "insertion_ratio", "deletion_ratio", "nonref_ratio")) {
    tmp <- assays(r)[[a]]
    colnames(tmp) <- gsub("^bam", a, colnames(tmp))
    cols <- c(cols, colnames(tmp))
    df <- cbind(df, tmp)
  }
  for (col in pick_cols(stats$stat)) {
    if (!col %in% colnames(df)) {
      f <- PARSE_COLS[[col]]
      if (is.null(f)) {
        print(col)
      }
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
  cols <- c("trna", "seq_position", "strand", cols, stats$new_col)
  df <- df[, cols]

  df
}

df <- read_result(opts$args)
se <- .to_se(df)

# # parse options and create output container
stats <- parse_stats(opts$options$stat)
# calculate stats and create output
df <- parse_result(se, stats)

write.table(df,
            opts$options$output,
            quote = FALSE,
            sep = "\t",
            col.names = TRUE, row.names = FALSE)
