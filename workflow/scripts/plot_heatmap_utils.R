library(dplyr)
library(ggplot2)


################################################################################
# command line tools
################################################################################

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
    if (harmonize_scaling == "NONE") {
      harmonize_scaling <- NULL
    } else if (harmonize_scaling == "MAX_SCORE") {
      harmonize_scaling <- max(df$score, na.rm = TRUE)
    } else if (harmonize_scaling == "CONTRASTS") {
      max_scores <- read.table(opts$options$max_scores, header = TRUE, sep = "\t")
      harmonize_scaling <- max_scores |>
        filter(cond1 == opts$options$condition1,
               cond2 == opts$options$condition2,
               bam_type == opts$options$bam_type) |>
        pull(opts$options$score_column) |>
        max()
    } else if (harmonize_scaling == "OVERALL") {
      max_scores <- read.table(opts$options$max_scores, header = TRUE, sep = "\t")
      harmonize_scaling <- max(max_scores[[opts$options$score]])
    } else {
      harmonize_scaling <- as.numeric(harmonize_scaling)
    } 
  }

  # quick fix
  if (is.element("mods", colnames(df))) {
    i <- is.na(df$mods)
    df$mods[i] <- ""
  }
  show_mods <- is.element("mods", colnames(df)) &&
    any(df$mods != "") &&
    !opts$options$hide_mods

  return(list(
    "title" = title,
    "position_xlab" = position_xlab,
    "coverage_xlab" = coverage_xlab,
    "harmonize_scaling" = harmonize_scaling,
    "show_mods" = show_mods,
    "show_ref" = FALSE #opts$options$show_ref
  ))
}

################################################################################
# estimate coverage by median from observed base calls per tRNA
################################################################################

estimate_coverage_info <- function(df, condition1, condition2) {
  coverage_summary <- df |>
    dplyr::select("trna", starts_with("reads_")) |>
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
    dplyr::select(all_of(c("condition_i", "condition"))) |>
    distinct() |>
    arrange(condition_i)
  coverage_summary$condition <- factor(coverage_summary$condition,
                                       levels = df_condition$condition,
                                       ordered = TRUE)
  
  return(coverage_summary |> dplyr::select(-all_of("condition_i")))
}

################################################################################
# plots
################################################################################

plot_coverage <- function(coverage_summary, opts, xlab = "reads") {
  p <- coverage_summary |>
    ggplot(aes(x = reads, y = trna_label,
               colour = condition, fill = condition)) +
    scale_colour_manual(values = c("#bf568b", "#bcd357")) +
    scale_fill_manual(values = c("#bf568b", "#bcd357")) +
    scale_x_continuous(breaks = scales::breaks_pretty(3),
                       labels = scales::label_number(
                         scale_cut = scales::cut_short_scale())) +
    ylab("") +
    xlab(xlab) +
    theme_bw(base_size = opts$options$base_size) +
    theme(plot.margin = margin(0, 0, 0, 0),
          legend.position = "top",
          legend.title = element_blank(),
          legend.direction = "vertical",
          axis.title.y = element_blank())

  if (opts$options$debug) {
    p <- p + theme(plot.background = element_rect(fill = "yellow"))
  }

  if (unique(coverage_summary$replicate) |> length() == 1) {
    p <- p + geom_col(position = position_dodge2())
  } else {
    p <- p + geom_point(position = position_dodge2(width = 0.5), size = rel(.75))
  }

  p
}


read_mod_abbrevs <- function(mod_abbrevs_fname) {
  mod_abbrevs <- read.table(mod_abbrevs_fname,
                            header = TRUE,
                            sep = "\t",
                            quote = "",
                            comment.char = "")

  return(mod_abbrevs)
}


plot_mod_abbrevs_table <- function(mod_abbrevs) {
  mod_abbrevs <- mod_abbrevs |>
    dplyr::select(all_of(c("short_name", "abbrev"))) |>
    filter(!is.na(abbrev)) |>
    distinct() |>
    dplyr::select(all_of(c("abbrev", "short_name"))) |>
    dplyr::rename(`RNAMods\ncode` = "abbrev", `Mod.` = "short_name")

  return(gridExtra::tableGrob(mod_abbrevs, rows = NULL))
}


plot_heatmap <- function(df, base_size = 9, xlab = "position",
                         harmonize_scaling = NULL,
                         title = NULL,
                         show_ref = FALSE, show_mods = FALSE) {

  show_ref <- FALSE
  if (show_ref && show_mods) {
    df <- df |>
      mutate(combined_label = paste0(ref, mod_label))
  }

  p <- df |>
    ggplot(aes(x = position, y = trna_label, fill = score)) +
    geom_tile(colour = "white", width = rel(.95), height = rel(.95)) +
    ylab("") +
    xlab(xlab) +
    theme_bw(base_size = base_size) +
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
                        filter(flag_position == TRUE),
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
  if (show_ref && show_mods) {
    p <- p + geom_text(aes(label = combined_label), size = rel(rel_size)) +
      coord_equal(ratio = 1 / 2)
  } else if (show_ref) {
    p <- p + geom_text(aes(label = ref), size = rel(rel_size)) +
      coord_equal()
  } else if (show_mods) {
    p <- p + geom_text(aes(label = mod_label), size = rel(rel_size)) +
      coord_equal()
  } else {
    p <- p + coord_equal()
  }

  if (!is.null(title)) {
    title <- gsub("\\{anti_codon\\}",
                  paste0(unique(na.omit(df$anti_codon)), collapse = ", "), title)
    title <- gsub("\\{amino_acid\\}",
                  paste0(unique(na.omit(df$amino_acid)), collapse = ", "), title)
    p <- p + ggtitle(title)
  }

  p
}


################################################################################
## infer dimensions


infer_height <- function(df) {
  trna_height <- df$trna |>
    unique() |>
    length()

  height <- trna_height + 15

  return(height / 10)
}


infer_width <- function(df, coverage_summary = NULL, abbrevs = NULL) {
  position_width <- df$position |>
    unique() |>
    length()

  trna_width <- df$trna |>
    unique() |>
    as.character() |>
    nchar() |>
    max()

  width <- position_width + trna_width
  if (!is.null(coverage_summary)) {
    width <- width + trna_width
  }
  if (!is.null(abbrevs)) {
    width <- width +
      nchar(abbrevs$short_name) |> max()
      nchar(abbrevs$abbrev) |> max()
  }

  return(width / 10)
}


infer_ggsave_opts <- function(df, coverage_summary = NULL,
                              mod_abbrevs = NULL,
                              ggsave_opts) {

  if (is.null(ggsave_opts)) {
    ggsave_opts <- list(
      "width" = infer_width(df, coverage_summary, mod_abbrevs),
      "height" = infer_height(df)
    )
  } else {
    if (is.null(ggsave_opts$width)) {
      ggsave_opts$width <- infer_width(df, coverage_summary)
    }
    if (is.null(ggsave_opts$height)) {
      ggsave_opts$height <- infer_height(df)
    }
  }

  return(ggsave_opts)
}

#################################################################################
## group tRNA and output
#################################################################################

add_groups <- function(group_opts, df, coverage_summary = NULL) {
  if (!is.null(group_opts)) {
    group_opts <- strsplit(group_opts, ",")[[1]]
    df <- df |>
      tidyr::unite(group, all_of(group_opts), sep = "-", remove = FALSE)
    if (!is.null(coverage_summary)) {
      coverage_summary <- coverage_summary |>
        tidyr::unite(group, all_of(group_opts), sep = "-", remove = FALSE)
    }
  } else {
    df$group <- "_"
    if (!is.null(coverage_summary)) {
      coverage_summary$group <- "_"
    }
  }

  return(
    list(
      df = df,
      coverage_summary = coverage_summary))
}

split_groups <- function(df, coverage_summary = NULL) {
  # split based on group
  df <- split(df, df$group)

  # container for split coverage_summary
  if (!is.null(coverage_summary)) {
    coverage_summary <- split(coverage_summary, coverage_summary$group)
  } else {
    coverage_summary <- list()
  }

  return(
    list(
      df = df,
      coverage_summary = coverage_summary))
}

process_group <- function(df, coverage_summary, process_args) {
  df <- remove_missing_positions(df)

  # optional sort df and coverage
  if (!is.null(process_args$sort_by_callback)) {
    sorted <- process_args$sort_by_callback(df, coverage_summary)
    df <- sorted$df
    coverage_summary <- sorted$coverage_summary
  }

  return(
    list(
      df = df,
      coverage_summary = coverage_summary))
}

group_trna <- function(group_opts,
                       df, coverage_summary = NULL,
                       output_dir,
                       opts) {
  plot_args <- get_plot_args(df, opts)
  process_args <- get_process_args(opts)
  
  plots <- list()
  
  if (!is.null(group_opts)) {
    i <- is.na(df$amino_acid)
    if (any(i)) {
      warning("Missing 'amino_acid'")
    }
    i <- is.na(df$anti_codon)
    if (any(i)) {
      warning("Missing 'anti_codon'")
    }
  }

  l <- add_groups(group_opts, df, coverage_summary)
  groups <- split_groups(l$df,
                         l$coverage_summary)

  for (group in names(groups$df)) {
    l <- process_group(groups$df[[group]],
                       groups$coverage_summary[[group]],
                       process_args)

    group_plots <- list(
      "df" = l$df,
      "coverage_summary" = l$coverage_summary
    )
        
    sep <- ifelse(group == "_", "", "-")
    fname_group <- ifelse(group == "_", "", group)

    ncol <- 1
    heatmap <- plot_heatmap(l$df,
                            base_size = opts$options$base_size,
                            xlab = plot_args[["position_xlab"]],
                            harmonize_scaling = plot_args$harmonize_scaling,
                            title = plot_args$title,
                            show_ref = plot_args$show_ref,
                            show_mods = plot_args$show_mods)
    custom_ggsave(file.path(output_dir, paste0(fname_group, sep, "heatmap.pdf")),
                  heatmap,
                  ggsave_opts = infer_ggsave_opts(df, ggsave_opts = opts$options$ggsave_opts),
                  opts$options$no_crop)
    group_plots[["heatmap"]] <- heatmap

    # add coverage summary
    if (!is.null(l$coverage_summary)) {
      ncol <- ncol + 1
      widths <- unit(c(-1, .10), c("null", "npc"))
      heatmap <- heatmap +
        guides(y.sec = guide_axis_label_trans(identity)) +
        plot_coverage(l$coverage_summary, opts, plot_args[["coverage_xlab"]]) +
        theme(axis.text.y = element_blank())
      
      p <- heatmap +
        plot_layout(ncol = ncol, widths = widths)

      ggsave_opts <- infer_ggsave_opts(l$df, l$coverage_summary,
                                       ggsave_opts = opts$options$ggsave_opts)
      ggsave_opts$width <- ggsave_opts$width + 3
      
      custom_ggsave(file.path(output_dir, paste0(fname_group, sep, "heatmap_cov.pdf")),
                    p,
                    ggsave_opts = ggsave_opts,
                    opts$options$no_crop)
      group_plots[["heatmap_cov"]] <- p
      group_plots[["coverage_summary"]] <- l$coverage_summary
    }

    # quick fix
    if (is.element("mod_abbrevs", colnames(l$df))) {
      i <- is.na(l$df$mod_label)
      l$df$mod_label <- ""
    }
    # add mod abbrevs table
    if (is.element("mod_abbrevs", colnames(l$df)) && any(l$df$mod_label != "")) {
      ncol <- ncol + 1
      widths <- unit(c(-1, .10, .10), c("null", "npc", "npc"))
      
      mod_abbrevs <- read_mod_abbrevs(opts$options$abbrevs)
      mods <- paste0(l$df$mods, " ") |>
        unique() |>
        strsplit(" ")

      mod_abbrevs <- mod_abbrevs |>
        dplyr::filter(is.element(short_name, mods))

      heatmap <- heatmap + plot_mod_abbrevs_table(mod_abbrevs)
      p <- heatmap +
        plot_layout(ncol = ncol, widths = widths)

      ggsave_opts <- infer_ggsave_opts(l$df, l$coverage_summary, mod_abbrevs,
                                      ggsave_opts = opts$options$ggsave_opts)
      ggsave_opts$width <- ggsave_opts$width + 3
      
      custom_ggsave(file.path(output_dir, paste0(fname_group, sep, "heatmap_mods.pdf")),
                    p,
                    ggsave_opts = ggsave_opts,
                    opts$options$no_crop)
      group_plots[["heatmap_mods"]] <- p
    }
    plots[[group]] <- group_plots
  }

  return(plots)
}

flag_positions <- function(df, flag_positions) {
  i <- is.element(df$trna_coords, flag_positions)
  if (any(i)) {
    df[i, "score"] <- NA
    df[i, "flag_position"] <- TRUE
  }

  return(df)
}


mark_positions <- function(df, pattern) {
  i <- grepl(pattern, df$trna_coords)
  if (!any(i)) {
    stop("Pattern did not mark any position!")
  }
  df[i, "mark_position"] <- TRUE

  return(df)
}


remove_3adapter <- function(df, trna_annotation, five_adapter, three_adapter) {
  five_adapter <- ifelse(is.null(five_adapter), 0, five_adapter)
  df_trna_annotation <- trna_annotation |>
    mutate(trna_length = nchar(seq)) |>
    dplyr::select(trna, trna_length)

  df <- left_join(df, df_trna_annotation, by = join_by(trna)) |>
    filter(seq_position <= (trna_length + five_adapter)) |>
    dplyr::select(-trna_length)
  stopifnot(!any(is.na(df$trna_length)))

  return(df)
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

  df <- tidyr::complete(df, !!! rlang::syms(cols), fill = fill)

  # fixing trna and strand
  fixes <- df |>
    dplyr::select(trna_label, trna, strand) |>
    filter(!if_any(everything(), is.na)) |>
    distinct()
  df <- df |>
    dplyr::select(-c(trna, strand)) |>
    inner_join(
      fixes,
      by = join_by(trna_label))

  if (any(!is.na(df$amino_acid))) {
    fixes <- df |>
      dplyr::select(trna_label, amino_acid) |>
      filter(!if_any(everything(), is.na)) |>
      distinct()
    if (nrow(fixes) > 0) {
      df <- df |>
        dplyr::select(-amino_acid) |>
        left_join(
          fixes,
          by = join_by(trna_label))
    }
  }
  if (any(!is.na(df$anti_codon))) {
    fixes <- df |>
      dplyr::select(trna_label, anti_codon) |>
      filter(!if_any(everything(), is.na)) |>
      distinct()
    if (nrow(fixes) > 0) {
      df <- df |>
        dplyr::select(-anti_codon) |>
        left_join(
          fixes,
          by = join_by(trna_label))
    }
  }

  # fixing reference
  i <- is.na(df$ref)
  if (any(i)) {
    df[i, "ref"] <- "*"
  }
  if (!is.null(df$mods)) {
    i <- is.na(df$mods)
    if (any(i)) {
      df[i, "mods"] <- ""
    }
  }
  if (!is.null(df$mod_label)) {
    i <- is.na(df$mod_label)
    if (any(i)) {
      df[i, "mod_label"] <- ""
    }
  }
  
  # fixing coords
  df$coords <- trna_coords(df)

  df
}


process_sprinzl_coords <- function(df, sprinzl_fname, hide_varm, show_introns, intron_start = "37") {
  # remove positions with no sprinzl mapping
  df <- df[!is.element(df$sprinzl, c("-", ".")), ] # - -> gap, . unmatched
  sprinzl_coords <- read.table(sprinzl_fname, header = FALSE)$V1
  i <- sprinzl_coords == "-"
  sprinzl_coords <- sprinzl_coords[!i]

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

save_plot <- function(p, filename, opts) {
  output <- file.path(opts$options$output_dir, filename)
  
  save(p, file = output)
  
  return(output)
}

##
# Order data frame by different keys
##

order_by_helper <- function(df, coverage_summary, trna_labels) {
  df$trna_label <- factor(df$trna_label, levels = trna_labels, ordered = TRUE)
  if (!is.null(coverage_summary)) {
    coverage_summary$trna_label <- factor(coverage_summary$trna_label,
                                          levels = trna_labels, ordered = TRUE)
  }

  return(list("df" = df, coverage_summary = coverage_summary))
}
order_by_coverage <- function(df, coverage_summary) {
  total_reads_info <- coverage_summary |>
    summarise(.by = c("trna_label"), total_reads = sum(reads)) |>
    dplyr::select(all_of(c("trna_label", "total_reads"))) |>
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
