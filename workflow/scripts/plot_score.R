# TODO remove

#!/usr/bin/env Rscript

# Plot JACUSA2 score and existing modification info as a heatmap

library(optparse)
library(ggplot2)
library(patchwork)
library(dplyr)

source(file.path(Sys.getenv("QUTRNA2"), "scripts", "utils.R"))
options(error = function() traceback(3))

args <- c("--condition1=WT",
          "--condition2=DNMT2",
          "--score_column=MDI_subsampled",
          "--split=all",
          "--position_column=sprinzl",
          "--ref_fasta=/beegfs/prj/tRNA_Francesca_Tuorto/index/all_human_tRNAs_v2_linkerNovoa_final.fa",
          "--hide_varm",
          #"--show_introns",
          "--title={condition1} vs. {condition2}; score: {score}",
          "--modmap=/beegfs/homes/mpiechotta/git/QutRNA/data/Hsapi38/human_modomics_sprinzl.tsv",
          "--sort_by=trna",
          "--estimate_coverage",
          #"--mark_positions=34",
          #"--amino_acids=Asp,Asn,Tyr,Gln,Glu,His,Leu,Thr,Gly,Val",
          #"--amino_acids=Asn", # Asp,Asn,Tyr,Gln,Glu,His,Leu,Thr,Gly,Val",
          "--five_adapter=23",
          "--three_adapter=33",
          "--show_ref",
          #"--harmonize_scaling=-1",
          "--sprinzl=/beegfs/homes/mpiechotta/git/QutRNA/data/euk_sprinzl.txt",
          "--output_dir=~/tmp/plot_score/all",
          "--remove_prefix=Homo_sapiens_",
          "/beegfs/prj/tRNA_Francesca_Tuorto/data/HCT116_tRNA_RNA004/2_biolrepl/output-strict-basq1//results/jacusa2/cond1~WT/cond2~DNMT2/scores_sprinzl.tsv")


BASE_SIZE <- 9 # TODO make this an argument

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


# TODO
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

remove_missing_positions <- function(df) {
  tbl <- table(df$position, is.na(df$ref))
  if (is.element("FALSE", colnames(tbl))) {
    tbl[, "FALSE"] == 0
    i <- tbl[, "FALSE"] == 0
    if (any(i)) {
      missing_positions <- rownames(tbl[i, ])
      df <- df |>
        filter(!is.element(position, missing_positions))
    }
  }

  return(df)
}


add_missing_values <- function(df, cols = c("trna_label", "position")) {
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

  # FIX reference
  i <- is.na(df$ref)
  if (any(i)) {
    df[i, "ref"] <- ""
  }
  
  # trna
  # amino_acid
  # anti_codon

  fixes <- df |>
    select(trna_label, trna, strand, amino_acid, anti_codon) |>
    filter(!if_any(everything(), is.na)) |>
    distinct()
  df <- df |>
    select(-c(trna, strand, amino_acid, anti_codon)) |>
    inner_join(
      fixes,
      by = join_by(trna_label))

  df
}

###
# Plots
###

# debug code
# plot.background = element_rect(fill = "yellow")

plot_coverage <- function(coverage_summary, xlab = "reads") {
  p <- coverage_summary |>
    ggplot(aes(x = reads, y = trna_label, colour = condition, fill = condition)) +
    scale_colour_manual(values = c("#bf568b", "#bcd357")) +
    scale_fill_manual(values = c("#bf568b", "#bcd357")) +
    scale_x_continuous(breaks = scales::breaks_pretty(3),
                       labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    ylab("") +
    xlab(xlab) +
    theme_bw(base_size = BASE_SIZE) +
    theme(plot.margin = margin(0, 0, 0, 0),
          legend.position = "top",
          legend.title = element_blank(),
          legend.direction = "vertical",
          axis.title.y = element_blank())

  if (unique(coverage_summary$replicate) |> length() == 1) {
    p <- p + geom_col(position = position_dodge2())
  } else {
    p <- p + geom_point(position = position_dodge2(width = 0.5), size = rel(.75))
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

plot_heatmap <- function(df, xlab = "position", harmonize_scaling = NULL, title = NULL) {
  p <- df |>
    ggplot(aes(x = position, y = trna_label, fill = score)) +
    geom_tile(colour = "white", width = rel(.95), height = rel(.95)) +
    ylab("") +
    xlab(xlab) +
    theme_bw(base_size = BASE_SIZE) +
    theme(panel.background = element_rect(fill = "white", color = "white"),
          legend.position = "top",
          panel.border = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          #axis.title.x = element_text(size = 10),
          #axis.text.y = element_text(size = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
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

  # mark special positions
  if (any(df$mark_position)) {
    p <- p + geom_tile(data = df |>
                       filter(mark_position == TRUE),
                       color = "red")
  }
    
  rel_size <- 1.7
  show_mods <- is.element("mod_label", colnames(df))
  if (opts$options$show_ref && show_mods) {
    df <- df |>
      mutate(combined_label = paste0(ref, mod_label))
    p <- p + geom_text(aes(label = combined_label), size = rel(rel_size)) +
      coord_equal(ratio = 1 / 2)
  } else if (opts$options$show_ref) {
    p <- p + geom_text(aes(label = ref), size = rel(rel_size)) +
      coord_equal()
  } else if (show_mods) {
    p <- p + geom_text(aes(label = mod_label), size = rel(rel_size)) +
      coord_equal()
  } else {
    p <- p + coord_equal()
  }

  if (!is.null(title)) {
    title <- gsub("\\{anti_codon\\}", paste0(unique(na.omit(df$anti_codon)), collapse = ", "), title)
    title <- gsub("\\{amino_acid\\}", paste0(unique(na.omit(df$amino_acid)), collapse = ", "), title)
    p <- p + ggtitle(title)
  }
  
  p
}

plot_combined <- function(df, coverage_summary, plot_args) {
  ncol <- 1
  widths <- c(NA)
  p <- plot_heatmap(df,
                    xlab = plot_args[["position_xlab"]],
                    harmonize_scaling = plot_args$harmonize_scaling,
                    title = plot_args$title)
  
  if (!is.null(coverage_summary)) {
    ncol <- ncol + 1
    widths <- c(widths, .1)
    p <- p + guides(y.sec = guide_axis_label_trans(identity))
    p <- p | plot_coverage(coverage_summary, plot_args[["coverage_xlab"]]) +
      theme(axis.text.y = element_blank())
  }
  
  if (is.element("mod_label", colnames(df)) && any(df$mod_label != "")) {
    ncol <- ncol + 1
    p <- p | plot_mod_table(df)
  }
  widths <- unit(c(-1, 100), c("null", "pt"))
  p <- p + plot_layout(ncol = ncol, widths = widths)

  p
}

###################
# Save plots and utility functions
###################

deduce_height <- function(df) {
  height <- df$trna |>
    unique() |>
    length()
  
  height <- height / 2
  
  return(height)
}


deduce_width <- function(df) {
  width <- df$position |>
    unique() |>
    length()
  width <- width + df$trna |>
    unique() |>
    as.character() |>
    nchar() |>
    max()
  
  width <- width / 10
  
  return(width)
}



save_plot <- function(p, output, width = NA, height = NA) {
  ggsave(output, p, width = width, height = height, limitsize = FALSE)
  if (!opts$options$no_crop) {
    knitr::plot_crop(output, quiet = FALSE)
  }
}

####
# Split and output tRNAs isoacceptor, isodecoder, or all
####

split_isoacceptor <- function(df, coverage_summary, output_dir, plot_args, process_args) {
  # check if no amino_acid or anticodon
  i <- is.na(df$amino_acid) | df$amino_acid == ""
  if (any(i)) {
    warning("Missing 'amino_acid': ", output_dir)
  }
  
  l_df <- split(df, df$amino_acid)
  amino_acids <- names(l_df)
  
  l_coverage_summary <- list
  if (!is.null(coverage_summary)) {
    l_coverage_summary <- split(coverage_summary, coverage_summary$amino_acid)
  }
  
  for (amino_acid in amino_acids) {
    tmp_df <- l_df[[amino_acid]]
    tmp_df <- remove_missing_positions(tmp_df)
    tmp_coverage_summary <- l_coverage_summary[[amino_acid]]
    if (!is.null(process_args$sort_by_callback)) {
      res <- process_args$sort_by_callback(tmp_df, tmp_coverage_summary)
      tmp_df <- res$df
      tmp_coverage_summary <- res$coverage_summary
    }
    
    p <- plot_heatmap(tmp_df,
                      xlab = plot_args[["position_xlab"]],
                      harmonize_scaling = plot_args$harmonize_scaling,
                      title = plot_args$title)
    save_plot(p, file.path(output_dir, paste0(amino_acid, "_heatmap.pdf")), width = deduce_width(df), height = deduce_height(df))
    save(p, file = file.path(output_dir, paste0(amino_acid, "_heatmap.rdata")))
    if (!is.null(tmp_coverage_summary)) {
      p <- plot_combined(tmp_df, NULL, plot_args)
      save_plot(p, file.path(output_dir, paste0(amino_acid, "_combined.pdf")), width = deduce_width(df) + 5, height = deduce_height(df))
    }
  }
}

split_isodecoder <- function(df, coverage_summary, output_dir, plot_args, process_args) {
  # check if no amino_acid or anticodon
  i <- is.na(df$anti_codon) | df$anti_codon == ""
  if (any(i)) {
    warning("Missing 'anti_codon' for: ", output_dir)
  }
  
  l_df <- split(df, df$anti_codon)
  isodecoders <- names(l_df)
  
  l_coverage_summary <- list()
  if (!is.null(coverage_summary)) {
    l_coverage_summary <- split(
      coverage_summary,
      coverage_summary$anti_codon)
  }
  
  for (isodecoder in isodecoders) {
    tmp_df <- l_df[[isodecoder]]
    tmp_df <- remove_missing_positions(tmp_df)
    tmp_coverage_summary <- l_coverage_summary[[isodecoder]]
    if (!is.null(process_args$sort_by_callback)) {
      res <- process_args$sort_by_callback(tmp_df, tmp_coverage_summary)
      tmp_df <- res$df
      tmp_coverage_summary <- res$coverage_summary
    }
    
    p <- plot_heatmap(tmp_df,
                      xlab = plot_args$position_xlab,
                      harmonize_scaling = plot_args$harmonize_scaling,
                      title = plot_args$title)
    save_plot(p, file.path(output_dir, paste0(isodecoder, "_heatmap.pdf")), width = deduce_width(df), height = deduce_height(df))
    save(p, file = file.path(output_dir, paste0(isodecoder, "_heatmap.rdata")))
    if (!is.null(coverage_summary)) {
      p <- plot_combined(tmp_df, tmp_coverage_summary, plot_args)
      save_plot(p, file.path(output_dir, paste0(isodecoder, "_combined.pdf")), width = deduce_width(df) + 5, height = deduce_height(df))
    }
  }
}

split_all <- function(df, coverage_summary, output_dir, plot_args, process_args) {
  if (!is.null(process_args$sort_by_callback)) {
    res <- process_args$sort_by_callback(df, coverage_summary)
    df <- res$df
    coverage_summary <- res$coverage_summary
  }
  df <- remove_missing_positions(df)
  
  p <- plot_heatmap(df,
                    xlab = plot_args$position_xlab,
                    harmonize_scaling = plot_args$harmonize_scaling,
                    title = plot_args$title)
  output <- file.path(output_dir, "all_heatmap.pdf")
  save_plot(p, output, height = deduce_height(df), width = deduce_width(df))
  save(p, file = file.path(output_dir, "all_heatmap.rdata"))
  if (!is.null(coverage_summary)) {
    p <- plot_combined(df, coverage_summary, plot_args)
    output <- file.path(output_dir, "all_combined.pdf")
    save_plot(p, output, height = deduce_height(df), width = deduce_width(df) + 5)
  }
}

##
## CLI
##

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
    title <- gsub("\\{bam_type\\}", opts$options$bam_type, title)
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
  if (!is.null(harmonize_scaling)) {
    max_scores <- read.table(opts$options$max_scores, header = TRUE, sep = "\t")
    
    if (harmonize_scaling == "-1") {
      harmonize_scaling <- max(df$score, na.rm = TRUE)
    } else if (harmonize_scaling == "-2") {
      harmonize_scaling <- max_scores |>
        filter(cond1 == opts$options$condition1, cond2 == opts$options$condition2) |>
        pull(opts$options$score) |>
        max()
    } else if (harmonize_scaling == "-3") {
      harmonize_scaling <- max(max_scores[[opts$options$score]])
    }
    harmonize_scaling <- as.numeric(harmonize_scaling)
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
              help = "Show reference base"),
  make_option(c("--no_crop"),
              action = "store_false",
              default = TRUE,
              help = "Do  crop pdf"),
  make_option(c("--modmap"),
              type = "character",
              help = "File mapping of modomics short name to abbrevations"),
  make_option(c("--split"),
              type = "character",
              help = "Split plots by: 'isodecoder', 'isoacceptor', 'extended', 'all'"),
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
  make_option(c("--min_scores"),
              type = "numeric",
              default = NULL,
              help = "Show positions with minimal score"),
  make_option(c("--max_scores"),
              type = "character",
              help = "File to max scores"),
  make_option(c("--harmonize_scaling"),
              type = "numeric",
              help = "Set to -1 to harmonize scaling between all plots; or set to some value."),
  make_option(c("--bam_type"),
              type = "character",
              help = "Set to the bam type used.",
              default = ""),
  make_option(c("--positions"),
              type = "character",
              default = NULL,
              help = "Positions to show, use ','")
)

## parse CLI

opts <- parse_args(
  OptionParser(option_list = option_list),
  # args = args,
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
  stopifnot(file.exists(opts$options$sprinzl))
}
if (!is.null(opts$options$modmap)) {
  stopifnot(file.exists(opts$options$modmap))
}

stopifnot(is.element(opts$options$position_column, c("sprinzl", "seq_position")))

if (!is.null(opts$options$mod_abbrevs)) {
  stopifnot(!is.null(opts$options$modmap))
}

##

# read process JACUSA2 scores
df <- read.table(opts$args, sep = "\t", header = TRUE)
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
  positions <- strsplit(opts$options$position, ",") |>
    unlist()
  i <- is.element(df[[opts$options$position_column]], positions)
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
  df <- df[is.element(df$amino_acid, amino_acids), ]
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
  positions <- strsplit(opts$options$mark_positions, ",") |>
    unlist()
  
  i <- grepl(":", positions)
  if (any(!i)) {
    j <- df$position[!i] == positions
    positions <- paste0(unique(df$trna[j]), ":", positions)
  }
  
  df <- mark_positions(df, positions)
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

# create coverage summary: estimate from observed median read counts per tRNA
if (!is.null(opts$options$estimate_coverage)) {
  coverage_summary <- estimate_coverage_info(df, opts$options$condition1, opts$options$condition2)
} else {
  coverage_summary <- NULL
}
# keep only observed trnas
coverage_summary <- coverage_summary |>
  filter(is.element(trna, unique(df$trna)))
coverage_summary <- add_trna_details(coverage_summary)
coverage_summary$trna_label <- trna_label(coverage_summary, opts$options$remove_prefix)

# add mods and abbreviations
if (!is.null(opts$options$modmap)) {
  df <- add_mods(df, opts$options$modmap)
  
  if (!is.null(opts$options$mod_abbrevs)) {
    df <- add_mods(df, opts$options$modmap)
  }
}

# add missing values
df <- add_missing_values(df)

# QUICK FIX
#df <- add_trna_details(df, col = "trna")
#coverage_summary <- add_trna_details(coverage_summary, col = "trna")

# split data, plot, and save
if (opts$options$split == "isoacceptor") {
  split_isoacceptor(df, coverage_summary, opts$options$output_dir, get_plot_args(df, opts), get_process_args(opts))
} else if (opts$options$split == "isodecoder") {
  split_isodecoder(df, coverage_summary, opts$options$output_dir, get_plot_args(df, opts), get_process_args(opts))
} else if (opts$options$split == "all") {
  split_all(df, coverage_summary, opts$options$output_dir, get_plot_args(df, opts), get_process_args(opts))
} else {
  stop("Unknown split: ", opts$options$split)
}

# print additional debug info
if (opts$options$debug) {
  warnings()
}
