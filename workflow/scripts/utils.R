flag_positions <- function(df, flag_positions) {
  i <- is.element(df$trna_coords, flag_positions)
  if (any(i)) {
    df[i, "score"] <- NA
    df[i, "flag_position"] <- TRUE
  }

  return(df)
}


mark_positions <- function(df, mark_positions) {
  i <- is.element(df$trna_coords, mark_positions)
  if (any(i)) {
    df[i, "mark_position"] <- TRUE
  }

  return(df)
}
get_ref_fasta <- function(fname) {
  ref_fasta <- Biostrings::readBStringSet(fname)

  return(ref_fasta)
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
    select(-seq)

  df <- left_join(df, df_ref_fasta, by = join_by(trna)) |>
    filter(seq_position <= (seq_length - three_adapter))
  stopifnot(!any(is.na(df$seq_length)))

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
    dplyr::select(trna_label, trna, strand, amino_acid, anti_codon) |>
    filter(!if_any(everything(), is.na)) |>
    distinct()
  df <- df |>
    dplyr::select(-c(trna, strand, amino_acid, anti_codon)) |>
    inner_join(
      fixes,
      by = join_by(trna_label))
  
  df
}

process_sprinzl_coords <- function(df, sprinzl_fname, hide_varm, show_introns, intron_start = "37") {
  # remove positions with no sprinzl mapping
  df <- df[!is.element(df$sprinzl, c("-", ".")), ] # - -> gap, . unmatched
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

process_sprinzl_coords <- function(df, sprinzl_fname, hide_varm, show_introns, intron_start = "37") {
  # remove positions with no sprinzl mapping
  df <- df[!is.element(df$sprinzl, c("-", ".")), ] # - -> gap, . unmatched
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


# extract amino_acid anti_codon
add_trna_details <- function(df, col = "trna") {
  df$amino_acid <- ""
  df$anti_codon <- ""
  i <- grepl("tRNA", df[[col]])
  if (any(i)) {
    df[i, "amino_acid"] <- stringr::str_extract(df[i, col], ".*tRNA-([A-Za-z]+)-([A-Za-z]{3}).*", group = 1)
    df[i, "anti_codon"] <- stringr::str_extract(df[i, col], ".*tRNA-([A-Za-z]+)-([A-Za-z]{3}).*", group = 2)
  }

  return(df)
}
add_mods <- function(df, mods_fname) {
  mods <- read.table(mods_fname, header = TRUE, sep = "\t", quote = "", comment.char = "")
  mods$trna_coords <- paste0(mods$trna, ":", mods$pos)
  mods <- mods |>
    dplyr::select(-trna)

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

  return(df)
}

add_mod_label <- function(df) {
  df$mod_label <- ""
  if (is.element("mod", colnames(df))) {
    df$mod_label <- df$mod
  }
  if (is.element("mod_abbrev", colnames(df))) {
    i <- is.na(df$mod_abbrev) || df$mod_abbrev != ""
    df$mod_label[i] <- df$mod_abbrev[i]
  }

  return(df)
}

##
# Order data frame by different keys
##

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


save_plot <- function(p, output, width = NA, height = NA) {
  ggsave(output, p, width = width, height = height, limitsize = FALSE)
  if (!opts$options$no_crop) {
    knitr::plot_crop(output, quiet = FALSE)
  }
}


########################################################################################################################

# CLI
# stopifnot
# serialize, des

custom_ggsave_test <- function(opts) {
  tryCatch(
    {
      l <- eval(parse(text=opts))
      stopifnot(!is.null(l))
      return(l)
    }, error = function(msg){
      
    })
}

custom_ggsave_options <- function() {
  return(
    make_option(c("--ggsave_opts"),
                type = "character",
                help = "Text representation of list."))

custom_ggsave <- function(filename, plot, opts) {
  ggsave_opts <- list(
    "filename" = filename,
    "plot" = plot
  )
  if (!is.null(opts$options$ggsave_opts)) {
    # use opts$options$ggsave
  }
  do.call(ggsave, ggsave_opts)
}