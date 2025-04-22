#!/usr/bin/env Rscript

# Plot JACUSA2 score and existing modification info as a heatmap)


library(optparse)
library(ggplot2)
library(patchwork)
library(magrittr)
library(dplyr)

# FIXES - https://github.com/krassowski/complex-upset/issues/158
guide_axis_label_trans <- function(label_trans = identity, ...) {
  axis_guide <- guide_axis(...)
  axis_guide$label_trans <- rlang::as_function(label_trans)
  class(axis_guide) <- c("guide_axis_trans", class(axis_guide))
  axis_guide
}

guide_train.guide_axis_trans <- function(x, ...) {
  trained <- NextMethod()
  trained$key$.label <- x$label_trans(trained$key$.label)
  trained
}


get_ref_fasta <- function(fname) {
  ref_fasta <- Biostrings::readBStringSet(fname)

  return(ref_fasta)
}

harmonize_scaling <- function(df) {
  return(max(df$score))
}

trna_label <- function(df, remove_prefix = "") {
  new_label <- df$trna
  if (remove_prefix != "") {
    new_label <- gsub(remove_prefix, "", new_label)
  }
  
  return(new_label)
}

trna_coords <- function(df) {
  return(paste0(df$trna, ":", df$position))
}

remove_3adapter <- function(df, ref_fasta, three_adapter) {
  # convert biostring to data frame: trna, seq_length
  df_ref_fasta <- ref_fasta |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "trna") |>
    rename(seq = "x") |>
    mutate(seq_length = nchar(seq)) |>
    select(all_of(c("trna", "seq_length")))

  df <- left_join(df, df_ref_fasta, by = join_by(trna)) |>
    filter(seq_position <= (seq_length - three_adapter))
  stopifnot(!any(is.na(df$seq_length)))

  return(df)
}

process_sprinzl_coords <- function(df, sprinzl_fname, hide_varm, show_introns, intron_start = "37") {
  # remove positions with no sprinzl mapping
  df <- df[!df$sprinzl %in% c("-", "."), ] # - -> gap, . unmatched
  sprinzl_coords <- read.table(sprinzl_fname, header = FALSE)$V1
  
  if (hide_varm) {
    # remove from sequenced trnas
    df <- df[!grepl("e[0-9]+", df$sprinzl), ]
    # remove from sprinzl model
    sprinzl_coords <- sprinzl_coords[!grepl("e[0-9]+", sprinzl_coords)]
  }
  if (!show_introns) {
    # remove from sequenced trnas
    df <- df[!grepl("i[0-9]+", df$sprinzl), ]
    # remove from sprinzl model
    sprinzl_coords <- sprinzl_coords[!grepl("i[0-9]+", sprinzl_coords)]
  } else {
    # add intron coordinates to sprinzl coordinates
    i <- grepl("i[0-9]+", df$sprinzl)
    if (any(i)) {
      intron_start_i <- which(sprinzl_coords == intron_start)
      sprinzl_coords <- c(sprinzl_coords[1:intron_start_i],
                          df[i, "sprinzl"] |>
                            unique() |>
                            stringr::str_sort(numeric = TRUE),
                          sprinzl_coords[(intron_start_i + 1):length(sprinzl_coords)])
    }
  }

  df$position <- factor(df$sprinzl, levels = sprinzl_coords, ordered = TRUE)

  return(df)
}

process_seq_coords <- function(df) {
  df$position <- factor(df$position,
                        levels = stringr::str_sort(
                          unique(df$position),
                          numeric = TRUE), ordered = TRUE)

  return(df)
}

flag_positions <- function(df, flas_positions) {
  i <- df$trna_coords %in% flag_positions
  if (any(i)) {
    df[i, "score"] <- NA
    df[i, "flag_position"] <- TRUE
  }

  return(df)
}
mark_positions <- function(df, mark_positions) {
  i <- df$trna_coords %in% mark_positions
  if (any(i)) {
    df[i, "mark_position"] <- TRUE
  }

  return(df)
}

# extract amino_acid anti_codon
add_trna_details <- function(df) {
  df$amino_acid <- ""
  df$anti_codon <- ""
  i <- grepl("tRNA", df$trna)
  if (any(i)) {
    df[i, "amino_acid"] <- stringr::str_extract(df[i, "trna"], ".*tRNA-([A-Za-z]+)-([A-Za-z]{3}).*", group = 1)
    df[i, "anti_codon"] <- stringr::str_extract(df[i, "trna"], ".*tRNA-([A-Za-z]+)-([A-Za-z]{3}).*", group = 2)
  }

  return(df)
}

estimate_coverage_info <- function(df, condition1, condition2) {
  coverage_summary <- df |>
    select("trna", starts_with("reads_")) |>
    group_by(trna) |>
    summarise_all(median) |>
    ungroup() |>
    tidyr::pivot_longer(starts_with("reads_"),
                        names_to = "condition",
                        values_to = "reads") |>
    mutate(replicate = gsub("^reads_\\d+_", "", condition),
           condition = gsub("^reads_|_\\d+$", "", condition),
           condition_i = condition,
           condition = case_match(
             condition,
             "1" ~ condition1,
             "2" ~ condition2)) |>
    as.data.frame()

  # reorder condition
  df_condition <- coverage_summary |>
    select(all_of(c("condition_i", "condition"))) |>
    distinct() |>
    arrange(condition_i)
  coverage_summary$condition <- factor(coverage_summary$condition, levels = df_condition$condition, ordered = TRUE)
  
  return(coverage_summary |> select(-all_of("condition_i")))
}

# FIXME - change in Snakemake
# expect format trna, condition, replicate, reads
get_coverage_info <- function(coverage_fname) {
  df_coverage <- read.table(coverage_fname, header = TRUE, sep = "\t") |>
    #dplyr::rename(Ref = "rname", condition = "sample", coverage = "numreads") |>
    mutate(replicate = gsub(".+_(\\d+)$", "\\1", condition),
           condition = gsub("_\\d+$", "", condition)) |>
    filter(condition %in% c(opts$options$cond1, opts$options$cond2),
           Ref %in% df$Ref)

  # cov_cond1 <- cov |>
  #   filter(condition == opts$options$cond1) |>
  #   select(Ref, condition, coverage) |>
  #   group_by(Ref, condition) |>
  #   summarise(coverage1 = mean(coverage))
  # 
  # cov_summary <- cov %>%
  #   select(Ref, coverage) %>%
  #   summarise(total_coverage = sum(coverage),
  #             expression = median(coverage),
  #             .by = Ref)
  browser()
  return(df_coverage)
}

order_by_helper <- function(df, coverage_summary, trna_labels) {
  df$trna_label <- factor(df$trna_label, levels = trna_labels, ordered = TRUE)
  if (!is.null(coverage_summary)) {
    coverage_summary$trna_label <- factor(coverage_summary$trna_label, levels = trna_labels, ordered = TRUE)
  }

  return(list("df" = df, coverage_summary = coverage_summary))
}
order_by_coverage <- function(df, coverage_summary) {
  total_reads_info <- coverage_summary |> 
    summarise(.by = c("trna_label"), total_reads = sum(reads)) |>
    select(all_of(c("trna_label", "total_reads"))) |>
    distinct()
  trna_labels <- total_reads_info$trna_label[order(total_reads_info$total_reads)]

  return(order_by_helper(df, coverage_summary, trna_labels))
}
order_by_trna <- function(df, coverage_summary) {
  trna_labels <- unique(df$trna_label) |>
    rev()

  df$trna <- factor(df$trna_label, levels = trna_labels, ordered = TRUE)
  if (!is.null(coverage_summary)) {
    coverage_summary$trna <- factor(coverage_summary$trna, levels = trna_labels, ordered = TRUE)
  }

  return(order_by_helper(df, coverage_summary, trna_labels))
}

add_mods <- function(df, mods_fname) {
  mods <- read.table(mods_fname, header = TRUE, sep = "\t", quote = "", comment.char = "")
  mods$trna_coords <- paste0(df$trna, ":", mods$pos)

  df <- left_join(df, mods,
                  by = join_by(trna_coords))
  i <- is.na(df$mod)
  if (any(i)) {
    df$mod <- ""
  }

  return(df)
}

add_mod_abbrevs <- function(df, mod_abbrevs_fname) {
  mod_abbrevs <- read.table(mod_abbrevs_fname, header = TRUE, sep = "\t", quote = "", comment.char = "")
  colnames(mod_abbrevs) <- paste0("mod_", colnames(mod_abbrevs))
  df <- left_join(df, mod_abbrevs,
                  by = join_by(mod == mod_abbrev))

  return(mod_abbrevs)
}

add_mod_label <- function(df) {
  df$mod_label <- ""
  if ("mod" %in% colnames(df)) {
    df$mod_label <- df$mod
  }
  if ("mod_abbrev" %in% colnames(df)) {
    i <- is.na(df$mod_abbrev) || df$mod_abbrev != ""
    df$mod_label[i] <- df$mod_abbrev[i]
  }

  return(df)
}


add_missing_values <- function(df, opts, cols = c("trna_label", "position")) {
  fill <- list("exceed_max_score" = FALSE,
               "mark_position" = FALSE,
               "flag_position" = FALSE)
  if (!is.null(df$mod_label)) {
    i <- is.na(df$mod_label)
    if (any(i)) {
      df[i, "mod_label"] <- ""
    }
  }

  df <- tidyr::complete(df, !!! rlang::syms(cols), fill = fill)

  df
}

## plots

# debug code
# plot.background = element_rect(fill = "yellow")

plot_coverage <- function(coverage_summary, xlab = "reads") {
  p <- coverage_summary |>
    ggplot(aes(x = reads, y = trna_label, fill = condition)) +
    scale_colour_manual(values = c("#bf568b", "#bcd357")) +
    scale_fill_manual(values = c("#bf568b", "#bcd357")) +
    scale_x_continuous(breaks = scales::breaks_pretty(3),
                       labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    ylab("") +
    xlab(xlab) +
    theme_bw() +
    theme(plot.margin = margin(0, 0, 0, 0),
          #axis.text.y = element_blank(),
          legend.position = "top",
          text = element_text(family = "mono"),
          legend.title = element_blank(),
          legend.direction = "vertical",
          axis.title.y = element_blank())

  if (unique(coverage_summary$replicate) %>% length() == 1) {
    p <- p + geom_col(position = position_dodge2())
  } else {
    p <- p + geom_point(position = position_dodge2(width = 0.5))
  }

  p
}

plot_mod_table <- function(df) {
  df <- select(df, all_of(c("mod_short_name", "mod_abbrev"))) |>
    filter(!is.na(mod_abbrev)) |>
    filter(mod_short_name != mod_abbrev) |>
    distinct() |>
    select(all_of(c("mod_abbrev", "mod_short_name"))) |>
    dplyr::rename(`RNAMods\ncode` = "mod_abbrev", `Mod.` = "mod_short_name")

  return(gridExtra::tableGrob(df, rows = NULL))
}

plot_heatmap <- function(df, xlab = "position", harmonize_scaling = NULL) {
  p <- df |>
    ggplot(aes(x = position, y = trna_label, fill = score)) +
    geom_tile(colour = "white") +
    ylab("") +
    xlab(xlab) +
    coord_equal() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white", color = "white"),
          legend.position = "top",
          panel.border = element_blank(),
          axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(family = "mono"),
          axis.title.y = element_blank(),
          plot.margin = margin(0, .5, 0, 0))

  # ensure same scaling for all plots
  if (!is.null(harmonize_scaling)) {
    p <- p +
      scale_fill_gradient(low = "yellow", high = "blue",
                          na.value = "grey",
                          limits = c(0, harmonize_scaling), oob = scales::squish)
  } else {
    p <- p +
      scale_fill_gradient(low = "yellow", high = "blue",
                          na.value = "grey")
  }

  # mark flagged positions
  if (any(df$flag_position)) {
    p <- p + geom_point(data = df |>
                          filter(flag == TRUE),
                        shape = 4, color = "darkred")
  }
  # mark score exceeding max_score
  if (any(df$exceed_max_score)) {
    p <- p + geom_point(data = df |>
                          filter(exceed_max_score == TRUE),
                        shape = 4, color = "red")
  }
  
  # FIXME add rows -> show label
  # change size of tile or add more rows
  # if (opts$options$show_ref) {
  #     p <- p + geom_text(aes(label = ref_base))
  #     if (show_mods) {
  #       p <- p + geom_text(aes(label = mod), position = position_nudge(x = 0.5, y = 0.5))
  #     }
  # } else {
  #     if (!opts$options$hide_mods) {
  #       p <- p + geom_text(aes(label = mod))
  #     }
  # }

  # mark special positions
  if (any(df$mark_position)) {
    p <- p + geom_point(data = df |>
                        filter(mark_position == TRUE),
                        shape = 4, color = "green")
  }

  p
}

plot_combined <- function(df, coverage_summary, plot_args) {
  ncol <- 1
  widths <- c(8)
  p <- plot_heatmap(df,
                    xlab = plot_args[["position_xlab"]],
                    harmonize_scaling = plot_args[["harmonize_scaling"]])
  if (!is.null(plot_args$title)) {
    title <- plot_args$title
    title <- gsub("\\{anti_codon\\}", paste0(unique(na.omit(df$anti_codon)), collapse = ", "), title)
    title <- gsub("\\{amino_acid\\}", paste0(unique(na.omit(df$amino_acid)), collapse = ", "), title)
    p <- p + ggtitle(title)
  }

  if (!is.null(coverage_summary)) {
    ncol <- ncol + 1
    widths <- c(widths, 2)
    p <- p + guides(y.sec=guide_axis_label_trans(identity)) +
      plot_coverage(coverage_summary, plot_args[["coverage_xlab"]]) +
      theme(axis.text.y = element_blank())
  }

  if ("mod_label" %in% colnames(df) && any(df$mod_label != "")) {
    ncol <- ncol + 1
    widths <- c(widths, 2)
    p <- p | plot_mod_table(df)
  }

  p <- p + plot_layout(ncol = ncol, widths = unit(c(-1, 4), c("null", "cm")))

  p
}

# FIXME width and height
save_plot <- function(df, p, output) {
  # use opts
  # infer width and height from df
  # knitr::plot_crop()
  ggsave(output, p, width = 2 * 8.56, height = 5.49)
}

## split

split_isoacceptor <- function(df, coverage_summary, output_dir, plot_args, process_args) {
  # check if no amino_acid or anticodon
  i <- is.na(df$amino_acid) || df$amino_acid == ""
  if (any(i)) {
    warning("Missing 'amino_acid': ", output_dir)
  }

  l_df <- split(df, df$amino_acid)
  amino_acids <- names(l_df)

  l_coverage_summary <- list
  if (!is.null(coverage_summary)) {
    stopifnot(length(setdiff(unique(coverage_summary$amino_acid), unique(df$amino_acid))))
    l_coverage_summary <- split(coverage_summary, coverage_summary$amino_acid)
  }

  for (amino_acid in amino_acids) {
    tmp_df <- l_df[[amino_acid]]
    tmp_coverage_summary <- l_coverage_summary[[amino_acid]]
    if (!is.null(process_args$sort_by_callback)) {
      res <- process_args$sort_by_callback(tmp_df, tmp_coverage_summary)
      tmp_df <- res$df
      tmp_coverage_summary <- res$coverage_summary
    }

    p <- plot_combined(tmp_df, tmp_coverage_summary, plot_args)
    save_plot(df, p, file.path(output_dir, paste0(amino_acid, ".pdf")))
  }
}

split_isodecoder <- function(df, coverage_summary, output_dir, plot_args, process_args) {
  # check if no amino_acid or anticodon
  i <- is.na(df$amino_acid) | df$amino_acid == "" |
           is.na(df$anti_codon) | df$anti_codon == ""
  if (any(i)) {
    warning("Missing 'amino_acid' or 'anti_codon' for: ", output_dir)
  }

  l_df <- split(df, paste0(df$amino_acid, "-", df$anti_codon))
  isodecoders <- names(l_df)

  l_coverage_summary <- list()
  if (!is.null(coverage_summary)) {
    l_coverage_summary <- split(
      coverage_summary,
      paste0(coverage_summary$amino_acid, "-", coverage_summary$anti_codon))
  }

  for (isodecoder in isodecoders) {
    tmp_df <- l_df[[isodecoder]]
    tmp_coverage_summary <- l_coverage_summary[[isodecoder]]
    if (!is.null(process_args$sort_by_callback)) {
      res <- process_args$sort_by_callback(tmp_df, tmp_coverage_summary)
      tmp_df <- res$df
      tmp_coverage_summary <- res$coverage_summary
    }
    p <- plot_combined(tmp_df, tmp_coverage_summary, plot_args)
    save_plot(df, p, file.path(output_dir, paste0(isodecoder, ".pdf")))
  }
}
split_all <- function(df, coverage_summary, output_dir, plot_args, process_args) {
  if (!is.null(process_args$sort_by_callback)) {
    res <- process_args$sort_by_callback(df, coverage_summary)
    df <- res$df
    coverage_summary <- res$coverage_summary
  }
  p <- plot_combined(df, coverage_summary, plot_args)
  output <- file.path(output_dir, "all.pdf")
  save_plot(df, p, output)
}

split_extended <- function(df, coverage_summary, output_dir, plot_args, process_args) {
  l_df <- split(df, df$trna)
  trnas <- names(l_df)

  l_coverage_summary <- list()
  if (!is.null(coverage_summary)) {
    # all trnas in df must have a match in coverage summary

    l_coverage_summary <- split(
      coverage_summary,
      coverage_summary$trna)
  }

  for (trna in trnas) {
    # TODO
  }
}

## CLI

# built arguments for processing score and coverage data frame
get_process_args <- function(opts) {
  sort_by_callback <- NULL
  if (!is.null(opts$options$sort_by)) {
    if (opts$options$sort_by == "coverage") {
      sort_by_callback <- order_by_coverage
    } else if (opts$options$sort_by == "trna") {
      sort_by_callback <- order_by_trna
    } else {
      stop("Unknown sort_by: ", opts$options$sort_by)
    }
  }
  
  return(
    list(
      "sort_by_callback" = sort_by_callback
    ))
}

get_plot_args <- function(df, opts) {
  title <- NULL
  if (!is.null(opts$options$title)) {  
    title <- opts$options$title
    title <- gsub("\\{condition1\\}", opts$options$condition1, title)
    title <- gsub("\\{condition2\\}", opts$options$condition2, title)
    title <- gsub("\\{score\\}", opts$options$score_column, title)
    title <- gsub("\\{position_column\\}", opts$options$position_column, title)
  }

  position_xlab <- "position"
  if (opts$options$position_column == "sprinzl") {
    position_xlab <- "Sprinzl pos."
  } else if (opts$options$position_column == "seq_position") {
    position_xlab <- "pos. in tRNA seq."
  } else {
    stop("Unknown position column: ", opts$options$position_column)
  }
  
  coverage_xlab <- NULL
  if (!is.null(opts$options$estimate_coverage)) {
    coverage_xlab <- "median coverage"
  } else if (opts$options$coverage_info) {
    coverage_xlab <- "reads"
  }
  
  harmonize_scaling <- opts$options$harmonize_scaling
  if (!is.null(harmonize_scaling) && harmonize_scaling == -1) {
    harmonize_scaling <- max(df$score)
  }
  return(list(
    "title" = title,
    "position_xlab" = position_xlab,
    "coverage_xlab" = coverage_xlab,
    "harmonize_scaling" = harmonize_scaling
  ))
}

option_list <- list(
  make_option(c("--score_column"),
              type = "character",
              default = "MDI",
              help = "Score column to use"),
  make_option(c("--condition1"),
              type = "character",
              default = "condition 1",
              help = "Condition 1"),
  make_option(c("--condition2"),
              type = "character",
              default = "condition 2",
              help = "Condition 2"),
  make_option(c("--output_dir"),
              type = "character",
              help = "Output"),
  make_option(c("--five_adapter"),
              type = "numeric",
              help = "Ignore positions on the 5' side"),
  make_option(c("--position_column"),
              type = "character",
              default = "seq_position",
              help = "Use column as coordinates, default: seq_position"),
  make_option(c("--three_adapter"),
              type = "numeric",
              help = "Ignore positions on the 3' side"),
  make_option(c("--ref_fasta"),
              type = "character",
              help = "Filename to sequence header"),
  make_option(c("--title"),
              type = "character",
              help = "Title of the plot"),
  make_option(c("--coverage_info"),
              type = "character",
              help = "Path to file with coverages"),
  make_option(c("--sprinzl"), # _mapping _info
              type = "character",
              help = "Path to sprinzl file with labels"),
  make_option(c("--hide_varm"),
              action = "store_true",
              default = FALSE,
              help = "Hide variable arm coordinates"),
  make_option(c("--hide_mods"),
              action = "store_true",
              default = FALSE,
              help = "Hide mods"),
  make_option(c("--sort_by"),
              default = "trna",
              help = "Sort by 'coverage' or 'tRNA', default: trna"),
  make_option(c("--show_introns"),
              action = "store_true",
              default = FALSE,
              help = "Show intron positions"),
  make_option(c("--estimate_coverage"),
              action = "store_true",
              default = FALSE,
              help = "Show coverage"),
  make_option(c("--debug"),
              action = "store_true",
              default = FALSE,
              help = "Show debug info"),
  make_option(c("--show_ref"),
              action = "store_true",
              default = FALSE,
              help = "Show reference base TODO"),
  # FIXME - remove this
  make_option(c("--crop"),
              action = "store_true",
              default = FALSE,
              help = "Do  crop pdf"),
  make_option(c("--modmap"),
              type = "character",
              help = "File mapping of modomics short name to abbrevations"),
  make_option(c("--split"),
              type = "character",
              help = "Split plots by: 'isodecoder', 'isoacceptor', 'all'"), # TODO 'extended'
  make_option(c("--remove_prefix"),
              type = "character",
              default = "",
              help = "Remove prefix from tRNA label"),
  make_option(c("--flag_positions"),
              type = "character",
              default = NULL,
              help = "Flag positions, e.g.: tRNA:position, use ','"),
  make_option(c("--amino_acids"),
              type = "character",
              default = NULL,
              help = "tRNAs to show, e.g.: tRNA, use ','"),
  make_option(c("--mark_positions"),
              type = "character",
              default = NULL,
              help = "Mark specific positions in plot, e.g.: tRNA:position, use ','"),
  make_option(c("--max_score"),
              type = "numeric",
              help = "Hide all scores exceeding MAX_SCORE"),
  make_option(c("--harmonize_scaling"),
              type = "numeric",
              help = "Set to -1 to harmonize scaling between all plots; or set to some value."),
  make_option(c("--positions"),
              type = "character",
              default = NULL,
              help = "Positions to show, use ','")
)

## parse CLI

opts <- parse_args(
  OptionParser(option_list = option_list),
   # args = c("--condition1=WT",
   #          "--condition2=DNMT2",
   #          "--score_column=MDI_subsampled",
   #          "--split=isodecoder",
   #          "--position_column=sprinzl",
   #          "--ref_fasta=/beegfs/prj/tRNA_Francesca_Tuorto/index/all_human_tRNAs_v2_linkerNovoa_final.fa",
   #          "--show_introns",
   #          "--title={condition1} vs. {condition2}; score: {score}",
   #          "--modmap=/beegfs/homes/mpiechotta/git/QutRNA/data/Hsapi38/human_modomics_sprinzl.tsv",
   #          "--sort_by=coverage",
   #          "--estimate_coverage",
   #          "--five_adapter=23",
   #          "--three_adapter=33",
   #          "--harmonize_scaling=-1",
   #          "--sprinzl=/beegfs/homes/mpiechotta/git/QutRNA/data/euk_sprinzl.txt",
   #          "--output_dir=~/tmp/plot_score",
   #          "--remove_prefix=Homo_sapiens_",
   # "~/tmp/new_test_data.tsv"),
  positional_arguments = TRUE
)

# parse options and check if valid
if (opts$options$debug) {
  options(error = traceback) # add debug options
}
stopifnot(!is.null(opts$options$output))
stopifnot(length(opts$args) == 1)
stopifnot(file.exists(opts$args))

# check 5' and 3' adapter settings
if (!is.null(opts$options$fife_adapter)) {
  stopifnot(opts$options$five_adapter > 0)
}
if (is.null(opts$options$three_adapter)) {
  stopifnot(opts$options$three_adapter > 0)
  # ensure that ref_fasta is provided to determine tRNA length
  stopifnot(!is.null(opts$options$ref_fasta))
}
# ensure that ref_fasta is provided to determine tRNA length
if (!is.null(opts$options$ref_fasta)) {
  stopifnot(file.exists(opts$options$ref_fasta))
}
if (!is.null(opts$options$coverage_info)) {
  stopifnot(file.exists(opts$options$coverage_info))
}
if (!is.null(opts$options$sprinzl)) {
  stopifnot(file.exists(opts$options$sprinz))
}
if (!is.null(opts$options$modmap)) {
  stopifnot(file.exists(opts$options$modmap))
}

stopifnot(opts$options$position_column %in% c("sprinzl", "seq_position"))

if (!is.null(opts$options$mod_abbrevs)) {
  stopifnot(!is.null(opts$options$mods))
}

## 

# read process JACUSA2 scores
df <- read.table(opts$args, sep = "\t", header = TRUE)
# ensure chosen columns are present in data
stopifnot(opts$options$position_column %in% colnames(df))
stopifnot(opts$options$score_column %in% colnames(df))
stopifnot("seq_position" %in% colnames(df))

# remove 5' adapter, if any
if (!is.null(opts$options$five_adapter)) {
  df <- df[df$seq_position > opts$options$five_adapter, ]  
}
# remove 3' adapter, if any
if (!is.null(opts$options$three_adapter)) {
  ref_fasta <- get_ref_fasta(opts$options$ref_fasta)
  df <- remove_3adapter(df, ref_fasta, opts$options$three_adapter)
}

# process chosen coordinates system
if (opts$options$position_column == "sprinzl") {
  df <- process_sprinzl_coords(df, opts$options$sprinzl, opts$options$hide_varm, opts$options$show_introns)
} else if (opts$options$position_column == "seq_position") {
  df <- process_seq_coords(df)
} else {
  stop("Unknown position column: ", opts$options$position_column)
}

# keep specific positions
if (!is.null(opts$options$positions)) {
  positions = strsplit(opts$options$position, ",") |>
    unlist()
  i <- df[[opts$options$position_column]] %in% positions
  if (any(i)) {
    df <- df[i, ]
  }
  
  return(df)
}

# add descriptive trna labels and coordinates
df$trna_coords <- trna_coords(df)
df$trna_label <- trna_label(df, opts$options$remove_prefix)

# split tRNA-XXX-YYY to XXX=amino acid and YYY=anti codon
df <- add_trna_details(df)

# filter by amino acid(s)
if (!is.null(opts$options$amino_acids)) {
  amino_acids <- strsplit(opts$options$amino_acids, ",") |>
    unlist()
  df <- df[df$amino_acid %in% amino_acids, ]
}

# hide/flag some positions
df$flag_position <- FALSE
if (!is.null(opts$options$flag_position)) {
  flag_positions <- strsplit(opts$options$flag_positions, ",") |>
    unlist()
  df <- flag_positions(df, flag_positions)
}

# mark specific positions
df$mark_position <- FALSE
if (!is.null(opts$options$mark_position)) {
  mark_positions <- strsplit(opts$options$mark_positions, ",") |>
    unlist()
  df <- mark_positions(df, mark_positions)
}

# process score columns
df$score <- df[, opts$options$score_column]
df$score[df$score < 0] <- 0
df$exceed_max_score <- FALSE
if (!is.null(opts$options$max_score)) {
  i <- df$score > opts$max_score
  if (any(i)) {
    df$score[i] <- NA
    df$exceed_max_score[i] <- TRUE
  }
}

# create coverage summary: estimate from observed median read counts per tRNA or
# use external quantification
if (!is.null(opts$options$coverage_info)) {
  coverage_summary <- get_coverage_info(opts$options$coverage_info)
} else if (!is.null(opts$options$estimate_coverage)) {
  coverage_summary <- estimate_coverage_info(df, opts$options$condition1, opts$options$condition2)
} else {
  coverage_summary <- NULL
}
# keep only observed trnas
coverage_summary <- coverage_summary |>
  filter(trna %in% unique(df$trna))
coverage_summary <- add_trna_details(coverage_summary)
coverage_summary$trna_label <- trna_label(coverage_summary, opts$options$remove_prefix)

# add mods and abbreviations
if (!is.null(opts$options$mods)) {
  df <- add_mods(df, opts$options$mods)

  if (!is.null(opts$options$mod_abbrevs)) {
    df <- add_mods(df, opts$options$mods)
  }
}

# add missing values
df <- add_missing_values(df, opts)

# split data, plot, and save
if (opts$options$split == "isoacceptor") {
  split_isoacceptor(df, coverage_summary, opts$options$output_dir, get_plot_args(df, opts), get_process_args(opts))
} else if (opts$options$split == "isodecoder") {
  split_isodecoder(df, coverage_summary, opts$options$output_dir, get_plot_args(df, opts), get_process_args(opts))
} else if (opts$options$split == "all") {
  split_all(df, coverage_summary, opts$options$output_dir, get_plot_args(df, opts), get_process_args(opts))
} else if (opts$options$split == "extended") {
  split_extended(df, coverage_summary, opts$options$output_dir, get_plot_args(df, opts), get_process_args(opts))
} else {
  stop("Unknown split: ", opts$options$split)
}

# print additional debug info
if (opts$options$debug) {
  warnings()
}
