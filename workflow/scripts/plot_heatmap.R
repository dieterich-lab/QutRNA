#!/usr/bin/env Rscript

# Plot JACUSA2 score and existing modification info as a heatmap
# options(error = function() { traceback(3) ; stop("Error - could not create plots!")})

library(optparse)
library(ggplot2)
library(patchwork)
library(dplyr)

source(file.path(Sys.getenv("QUTRNA2"), "scripts", "utils.R"))
source(file.path(Sys.getenv("QUTRNA2"), "scripts", "plot_heatmap_utils.R"))


################################################################################
# debug / testing
################################################################################

args <- c(
  "--condition1=WT",
  "--condition2=QTRT1",
  "--bam_type=filtered-random_alignment",
  "--trna_annotation=results/data/trna_annotation.tsv",
  "--output_dir=results/plots/scores/cond1~WT/cond2~QTRT1/QTRT1_all-trnas_MDI-subsampled/bam~filtered-random_alignment",
  "--position_column=sprinzl",
  #"--group=amino_acid,anti_codon",
  "--five_adapter=24",
  "--remove_prefix=tRNA-",
  "--score_column=MDI_subsampled",
  "--sprinzl=data/sprinzl.txt",
  "--abbrevs=data/mod_abbrevs.tsv",
  "--max_scores=results/jacusa2/max_scores.tsv",
  "--harmonize_scaling=CONTRASTS",
  "--estimate_coverage",
  "--no_crop",
  "--trnas=Arg|Ala",
  #"--flag_positions='tRNA-Arg-CCG-2-1:37'",
  "--sort_by=trna",
  "--show_introns",
  "--mark_positions=:(34|38)$",
  "--title={condition1}, {condition2}, {score}, {bam_type}, {position_column}, {anti_codon}, {amino_acid}",
  "results/jacusa2/cond1~WT/cond2~QTRT1/bam~filtered-random_alignment/scores_sprinzl-mods.tsv"
)

################################################################################
# create command line options
################################################################################

option_list <- list(
  # adapter size
  make_option(c("--five_adapter"),
              type = "integer",
              help = "Ignore positions on the 5' side"),
  make_option(c("--three_adapter"),
              type = "integer",
              help = "Ignore positions on the 3' side"),
  # other files to read
  make_option(c("--abbrevs"),
              type = "character",
              help = "Filename mod abbrevs."),
  make_option(c("--trna_annotation"),
              type = "character",
              help = "Filename to tRNA annotation"),
  make_option(c("--max_scores"),
              type = "character",
              help = "File with max scores"),
  # which columns to use  
  make_option(c("--position_column"),
              type = "character",
              default = "seq_position",
              help = "Use column as coordinates, default: 'seq_position'"),
  make_option(c("--score_column"),
              type = "character",
              default = "MDI",
              help = "Score column to use"),
  make_option(c("--sprinzl"),
              type = "character",
              help = "Path to sprinzl file with labels"),
  # data set description
  make_option(c("--title"),
              type = "character",
              help = "Title of the plot. Patterns: {condition1}, {condition2}, {score}, {bam_type}, {position_column}, {anti_codon}, {amino_acid}"),
  make_option(c("--bam_type"),
              type = "character",
              help = "Set bam type",
              default = ""),
  make_option(c("--condition1"),
              type = "character",
              default = "condition 1",
              help = "Condition 1"),
  make_option(c("--condition2"),
              type = "character",
              default = "condition 2",
              help = "Condition 2"),
  # hide/show features  
  make_option(c("--estimate_coverage"),
              action = "store_true",
              default = FALSE,
              help = "Estimate coverage for tRNAs and show"),
  # make_option(c("--show_ref"),
  #             action = "store_true",
  #             default = FALSE,
  #             help = "Show reference base"),
  make_option(c("--hide_varm"),
              action = "store_true",
              default = FALSE,
              help = "Hide variable arm coordinates"),
  make_option(c("--hide_mods"),
              action = "store_true",
              default = FALSE,
              help = "Hide modifications"),
  make_option(c("--show_introns"),
              action = "store_true",
              default = FALSE,
              help = "Show intron positions"),
  # output    
  make_option(c("--group"),
              type = "character",
              help = "Group plots by: 'amino_acid', 'anti_codon', or 'amino_acid,anti_codon'. If empty, plot all tRNAs in one plot."),
  make_option(c("--remove_prefix"),
              type = "character",
              default = "",
              help = "Remove prefix from tRNA label for visualization"),
  make_option(c("--output_dir"),
              type = "character",
              help = "Output directory for plots"),
  make_option(c("--sort_by"),
              default = "trna",
              help = "Sort by 'coverage' or 'tRNA', default: trna"),
  make_option(c("--no_crop"),
              action = "store_true",
              default = FALSE,
              help = "Do NOT crop final pdf"),
  # filter positions or tRNAs  
  make_option(c("--flag_positions"),
              type = "character",
              help = "Flag positions, e.g.: tRNA:position, use ','"),
  make_option(c("--max_score"),
              type = "numeric",
              help = "Cap scores at max_score"),
  make_option(c("--trnas"),
              type = "character",
              help = "Regular expression to restrict tRNAs to visualize"),
  make_option(c("--ignore_trnas"),
              type = "character",
              help = "Regular expression to restrict tRNAs to visualize"),
  make_option(c("--min_score"),
              type = "numeric",
              help = "Min JACUSA2 score required to display a tRNA"),
  make_option(c("--min_reads"),
              type = "numeric",
              help = "Min reads per sample required to display a tRNA"),
  make_option(c("--score_quantile"),
              type = "numeric",
              help = "Min JACUSA2 score defined by quantile required to display a tRNA"),
  make_option(c("--mark_positions"),
              type = "character",
              help = "Regular expresison to mark specific positions in plot, e.g.: tRNA:position"),
  make_option(c("--positions"),
              type = "character",
              help = "Show only the following positions, use ',' to separate."),
  make_option(c("--harmonize_scaling"),
              type = "character",
              help = "Harmonize scaling, 'NONE', 'CONTRAST', 'FILTER', 'CONTRAST+FILTER' or a number"),
  # ggplot related
  make_option(c("--base_size"),
              type = "numeric",
              default = 9,
              help = "Base size for ggplot2."),
  custom_ggsave_option(),
  # debug related
  make_option(c("--debug"),
              action = "store_true",
              default = FALSE,
              help = "Show debug info")
)

################################################################################
# parse command line
################################################################################

if (exists("PLOT_ARGS")) {
  opts <- parse_args(
    OptionParser(option_list = option_list),
    args = PLOT_ARGS,
    positional_arguments = TRUE
  )
} else {
  opts <- parse_args(
    OptionParser(option_list = option_list),
    positional_arguments = TRUE
  )
}

BASE_SIZE <- opts$options$base_size

################################################################################
# check command line

# parse options and check if valid
if (opts$options$debug) {
  options(error = traceback) # add debug options
}

# check required values
stopifnot(!is.null(opts$options$output_dir))
stopifnot(length(opts$args) == 1)
stopifnot(file.exists(opts$args))
# add trna annotation
stopifnot(!is.null(opts$options$trna_annotation))
stopifnot(file.exists(opts$options$trna_annotation))

# check min_score and score_quantile
if (!is.null(opts$options$min_score) && !is.null(opts$options$score_quantile)) {
  stop("Set only 'min_score' OR 'quantile_score'")
}
if (!is.null(opts$options$min_score)) {
  stopifnot(opts$options$min_score >= 0)
}
if (!is.null(opts$options$score_quantile)) {
  stopifnot(opts$options$score_quantile >= 0 && opts$options$score_quantile <= 1)
}

# check 5' and 3' adapter settings
if (!is.null(opts$options$five_adapter)) {
  stopifnot(opts$options$five_adapter > 0)
}
if (is.null(opts$options$three_adapter)) {
  stopifnot(opts$options$three_adapter > 0)
  # ensure that trna_anntotation is provided to determine tRNA length
  stopifnot(!is.null(opts$options$trna_annotation))
}

# check if files provided that they exist
if (!is.null(opts$options$sprinzl)) {
  stopifnot(file.exists(opts$options$sprinzl))
}

# check that picked position are valid
stopifnot(is.element(opts$options$position_column, c("sprinzl", "seq_position")))

# how to group tRNAs in plots
if (!is.null(opts$options$group)) {
  stopifnot(is.element(opts$options$group,
                       c("amino_acid", "anti_codon", "amino_acid,anti_codon")))
}

# process optional ggsave opts
if (!is.null(opts$options$ggsave_opts)) {
  custom_ggsave_stopifnot(opts$options$ggsave_opts)
}

# read process JACUSA2 scores
df <- read.table(opts$args, sep = "\t", header = TRUE)

# keep optional subset of tRNAs
trna_annotation <- read.table(opts$options$trna_annotation, header = TRUE, sep = "\t")
if (!is.null(opts$options$trnas)) {
  pattern <- opts$options$trnas
  i <- grepl(pattern, trna_annotation$trna)
  if (!any(i)) {
    browser()
    stop("Pattern resulted in empty set of tRNAs!")
  }
  trna_annotation <- trna_annotation[i, ]
  df |>
    filter(is.element(trna, trna_annotation$trna))
}
if (!is.null(opts$options$ignore_trnas)) {
  pattern <- opts$options$ignore_trnas
  i <- grepl(pattern, trna_annotation$trna)
  if (!any(i)) {
    message(paste0("Pattern '--ignore_trnas=", pattern,"' resulted in empty set of tRNAs. No tRNA were filtered!"))
  }
  trna_annotation <- trna_annotation[!i, ]
  df |>
    filter(is.element(trna, trna_annotation$trna))
}
if (!is.null(opts$options$min_reads)) {
  trnas <- df |> 
    filter(if_all(starts_with("reads"), ~ . > !!opts$options$min_reads)) |>
    pull(trna) |>
    unique()
  df <- df |>
    filter(trna %in% trnas)
}
# add amino_acid and anti_codon
df <- inner_join(df, trna_annotation[, c("trna", "amino_acid", "anti_codon")],
                 by = join_by(trna))

# ensure chosen columns are present in data
stopifnot(is.element(opts$options$position_column, colnames(df)))
stopifnot(is.element(opts$options$score_column, colnames(df)))
stopifnot(is.element("seq_position", colnames(df)))

# remove 5' adapter, if any
if (!is.null(opts$options$five_adapter)) {
  df <- df[df$seq_position > opts$options$five_adapter, ]
}
# remove 3' adapter, if any
if (!is.null(opts$options$three_adapter)) {
  df <- remove_3adapter(df, trna_annotation,
                        opts$options$five_adapter, opts$options$three_adapter)
}

# process chosen coordinates system
if (opts$options$position_column == "sprinzl") {
  df <- process_sprinzl_coords(df, opts$options$sprinzl, opts$options$hide_varm, opts$options$show_introns)
} else if (opts$options$position_column == "seq_position") {
  df <- process_seq_coords(df)
} else {
  stop("Unknown position column: ", opts$options$position_column)
}

# add descriptive trna labels and coordinates
df$trna_coords <- trna_coords(df)
df$trna_label <- trna_label(df, opts$options$remove_prefix)

# hide/flag some positions
df$flag_position <- FALSE
if (!is.null(opts$options$flag_positions)) {
  positions <- strsplit(opts$options$flag_positions, ",") |>
    unlist()
  df <- flag_positions(df, positions)
}

# mark specific positions
df$mark_position <- FALSE
if (!is.null(opts$options$mark_position)) {
  df <- mark_positions(df, opts$options$mark_position)
}

# process score columns
df$score <- df[, opts$options$score_column]
df$score[df$score < 0] <- 0
df$exceed_max_score <- FALSE
if (!is.null(opts$options$max_score)) {
  i <- df$score > opts$options$max_score
  if (any(i, na.rm = TRUE)) {
    df$score[i] <- NA
    df$exceed_max_score[i] <- TRUE
  }
}

# create coverage summary: estimate from observed median read counts per tRNA
if (!is.null(opts$options$estimate_coverage)) {
  coverage_summary <- estimate_coverage_info(df, opts$options$condition1, opts$options$condition2)
  coverage_summary <- inner_join(coverage_summary, trna_annotation[, c("trna", "amino_acid", "anti_codon")],
                                 by = join_by(trna))
  coverage_summary$trna_label <- trna_label(coverage_summary, opts$options$remove_prefix)
} else {
  coverage_summary <- NULL
}


# add label from mods or mod. abbreviations
if (!is.null(df$mods) && any(df$mods != "") && !opts$options$hide_mods) {
  if (!is.null(opts$options$abbrevs)) {
    mod_abbrevs <- read_mod_abbrevs(opts$options$abbrevs)
    abbrevs <- strsplit(df$mods, " ") |>
      lapply(function(mods) {
        mod_abbrevs <- mod_abbrevs |>
          dplyr::filter(is.element(short_name, mods)) |>
          pull(abbrev) |>
          paste0(collapse = " ")
        
        return(mod_abbrevs)
      })
    df$mod_abbrevs <- unlist(abbrevs)
    i <- is.na(df$mod_abbrevs)
    if (any(i)) {
      df$mod_abbrevs[i] <- ""
    }
    df$mod_label <- df$mod_abbrevs
  } else {
    df$mod_label <- df$mods
  }
}

# add missing values
df <- add_missing_values(df)

# keep specific positions
if (!is.null(opts$options$positions)) {
  positions <- strsplit(opts$options$positions, ",") |>
    unlist()
  
  src <- as.character(df[[opts$options$position_column]])
  dst <- as.character(positions)
  
  i <- is.element(src, dst)
  if (any(i)) {
    df <- df[i, ]
  }
  
  return(df)
}

# filter by min_score or quantile_score
if (!is.null(opts$options$min_score)) {
  trnas <- df |>
    filter(score >= opts$options$min_score) |>
    pull(trna) |>
    unique()
  df <- df |>
    filter(trna %in% trnas)
}
if (!is.null(opts$options$score_quantile)) {
  scores <- df |>
    filter(!is.na(score)) |>
    pull(score)
  
  min_score <- quantile(scores, prob = opts$options$score_quantile, names = FALSE)
  trnas <- df |>
    filter(score >= min_score) |>
    pull(trna) |>
    unique()
  df <- df |>
    filter(trna %in% trnas)
}

################################################################################

# group data and create plots  
plots <- group_trna(opts$options$group,
                    df, coverage_summary,
                    opts$options$output_dir,
                    opts)


# print additional debug info
if (opts$options$debug) {
  warnings()
}
