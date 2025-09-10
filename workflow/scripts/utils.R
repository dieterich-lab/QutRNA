library(dplyr)
library(optparse)

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

########################################################################################################################

rename_read_type_stopifnot <- function(opts) {
  if (is.null(opts)) {
    return(NULL)
  }

  rename_read_types <- do.call(rbind, strsplit(strsplit(opts, ";")[[1]], "=")) |>
    as.data.frame()
  colnames(rename_read_types) <- c("old", "new")

  return(rename_read_types)
}


rename_read_type_option <- function() {
  return(
    make_option(c("--rename_read_type"),
                type = "character",
                help = "Separate by ';'"))
}

########################################################################################################################


custom_ggsave_stopifnot <- function(opts) {
  if (is.null(opts)) {
    return(NULL)
  }
  
  tryCatch(
    {
      l <- eval(parse(text = paste0("list(", opts, ")")))
      stopifnot(!is.null(l))
      return(l)
    }, error = function(msg){
      stop(msg)
    })
}


custom_ggsave_option <- function() {
  return(
    make_option(c("--ggsave_opts"),
                type = "character",
                help = "Text representation of list."))
}

custom_ggsave <- function(filename, plot, ggsave_opts = NULL, no_crop = TRUE) {
  if (is.null(ggsave_opts)) {
    ggsave_opts <- list()
  }
  
  ggsave_opts[["filename"]] <- filename
  ggsave_opts[["plot"]] <- plot
  ggsave_opts[["limitsize"]] <- FALSE

  do.call(ggsave, ggsave_opts)
  if (!no_crop) {
    knitr::plot_crop(filename, quiet = FALSE)
  }

}

################################################################################

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
