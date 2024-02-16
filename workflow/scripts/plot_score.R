#!/usr/bin/env Rscript

# Plot JACUSA2 score and existing modification info as a heatmap)

library(optparse)
library(ggplot2)
library(cowplot)
library(ggpubr)


option_list <- list(
  make_option(c( "--score"),
                type = "character",
                default = "Mis.Del.Ins",
                help = "Score column to use"),
  make_option(c( "--output"),
                type = "character",
                help = "Output"),
  make_option(c("--left"),
                type = "numeric",
                default = 24,
                help = "Ignore positions left"),
  make_option(c("--column"),
                type = "character",
                help = "Use column as coordinates, default: Pos3"),
  make_option(c("--length"),
                type = "numeric",
                default = 77,
                help = "Ignore positions right"),
  make_option(c("--title"),
                type = "character",
                help = "Title of the plot"),
  make_option(c("--hide_varm"),
                action = "store_true",
                default = FALSE,
                help = "Hide variable arm coordinates"),
  make_option(c("--show_introns"),
                action = "store_true",
                default = FALSE,
                help = "Show intron positions"),
  make_option(c("--show_ref"),
                action = "store_true",
                default = FALSE,
                help = "Show reference base"),
  make_option(c("--modmap"),
                type = "character",
                help = "File mapping of modomics short name to abbrevations"),
  make_option(c("--split"),
                type = "character",
                help = "Split by tRNA"),
  make_option(c("--positions"),
                type = "character",
                default = NULL,
                help = "Positions to show, use ','")
)

opts <- parse_args(
  OptionParser(option_list = option_list),
  positional_arguments = TRUE
)
stopifnot(!is.null(opts$options$output))
stopifnot(length(opts$args) == 1)

dir.create(opts$options$output, showWarnings = FALSE)
dir.create(paste0(opts$options$output, "/final"), showWarnings = FALSE)
df <- read.csv(opts$args)
df <- df[df$Pos3 >= opts$options$left, ]

pos_col <- opts$options$column
if (pos_col == "u_pos") {
  df <- df[df[[pos_col]] != "-", ]
  df <- df[df[[pos_col]] != ".", ]
  if (opts$options$hide_varm) {
    i <- !grepl("e[0-9]+", df[[pos_col]])
    df <- df[i, ]
  }
  if (!opts$options$show_introns) {
    i <- !grepl("i[0-9]+", df[[pos_col]])
    df <- df[i, ]
  }
} else {
  df[[pos_col]] <- factor(df[[pos_col]],
                        levels = stringr::str_sort(unique(df[[pos_col]]), numeric = TRUE))
}

df$trna <- df$Ref
i <- grepl("tRNA", df$Ref)
df[i, "trna"] <- do.call(rbind, stringr::str_match_all(df[i, "Ref"], ".*tRNA-([A-Za-z]{3}).*"))[, 2]
df[["Ref"]] <- factor(df[["Ref"]],
                      levels = stringr::str_sort(unique(df[["Ref"]]), numeric = TRUE, decreasing = TRUE))

df[["ref_base"]] <- substr(df[["Kmer"]], 3, 3)

df$score <- df[, opts$options$score]
df$score[df$score < 0] <- 0
if (!is.null(opts$options$modmap)) {
  modmap <- read.table(opts$options$modmap, header = TRUE, sep = "\t", quote = "", comment.char = "")
  stopifnot(length(modmap$shortname) == length(unique(modmap$shortname)))
  df <- merge(df, modmap, by.x = "mod", by.y = "short_name", all.x = TRUE, sort = FALSE)
  i <- !is.na(df$abbrev)
  df$short_name <- df$mod
  df[i, "mod"] <- df[i, "abbrev"]
  df <- dplyr::select(df, "Ref", tidyr::all_of(pos_col), "trna", "score", "mod", "short_name", "abbrev", "Pos3", "ref_base")
} else {
  df <- dplyr::select(df, "Ref", tidyr::all_of(pos_col), "trna", "score", "mod", "Pos3", "ref_base")
}
if (is.null(opts$options$positions)) {
  df <- tidyr::complete(df, Ref, .data[[pos_col]])
  df[is.na(df$mod), "mod"] <- ""
  df[is.na(df$score), "score"] <- NA
} else {
  positions = strsplit(opts$options$position, ",")[[1]]
  i <- df[[pos_col]] %in% positions
  df <- df[i, ]
}

plot_table <- function(df) {
  df <- dplyr::select(df, short_name, abbrev) |>
      dplyr::filter(!is.na(abbrev)) |>
      dplyr::filter(short_name != abbrev) |>
      dplyr::distinct() |>
      dplyr::rename(`Mod.` = short_name, `RNAMods code` = abbrev)

  if (nrow(df) > 0) {
    return(ggtexttable(df, rows = NULL))
  }

  return(NULL)
}

plot_main <- function(df) {
  p <- df |>
    dplyr::mutate(Ref = gsub("Homo_sapiens_", "", Ref)) |>
    ggplot(aes(x = .data[[pos_col]], y = Ref)) +
      ylab("") +
      scale_fill_gradient(low = "yellow", high = "blue", na.value = "grey") +
      coord_equal() +
      theme_bw() +
      theme(panel.background = element_rect(fill = "white", color = "white"),
            panel.border = element_blank(),
            axis.text.x = element_text(angle = 45),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = unit(c(0,0,0,0), "cm"))


   if (opts$options$column == "u_pos") {
     if (opts$options$show_introns) {
      p <- p + xlab("clover leaf pos. (+i_ntron)")
     } else {
      p <- p + xlab("clover leaf pos.")
     }
   } else {
      p <- p + xlab("pos. in tRNA seq.")
   }

   if (opts$options$show_ref) {
      # TODO position
      p <- p +
        geom_tile(aes(fill = score), width = 0.9, height = 0.9) +
        geom_text(aes(label = mod, nudge_x = 0.5, y_nudge=0.5)) +
        geom_text(aes(label = ref_base))
   } else {
      p <- p +
        geom_tile(aes(fill = score), width = 0.9, height = 0.9) +
        geom_text(aes(label = mod))
   }

  if (!is.null(opts$options$title)) {
    p <- p + ggtitle(opts$options$title)
  }

  if (!is.null(opts$options$modmap)) {
    tbl <- plot_table(df)
    if (!is.null(tbl)) {
      p <- plot_grid(plotlist = list(p + theme(legend.position = "none"),
                                     get_legend(p),
                                     plot_table(df)),
                      ncol = 3,
                      nrow = 1,
                      rel_widths = c(80, 4, 10),
                      greedy = FALSE)
    } else {
      p <- plot_grid(plotlist = list(p + theme(legend.position = "none"),
                                     get_legend(p)),
                      ncol = 2,
                      nrow = 1,
                      rel_widths = c(80, 4),
                      greedy = FALSE)
    }
  } else {
    p <- p + theme(legend.position = "right")
  }

  p
}

if (opts$options$hide_varm) {
  df[[pos_col]] <- factor(df[[pos_col]],
                          levels = c(1:17, "17a", 18:20, "20a", "20b", 21:37, paste0("i", 1:24), 38:76 ),
                          ordered = TRUE)
} else {
  df[[pos_col]] <- factor(df[[pos_col]],
                          levels = c(1:17, "17a", 18:20, "20a", "20b", 21:37, paste0("i", 1:24), 38:45, paste0("e", 1:19), 46:76 ),
                          ordered = TRUE)
}

if (is.null(opts$options$split)) {
  df <- df |>
    dplyr::filter(any(!is.na(score)), .by = u_pos)

  df[[pos_col]] <- droplevels(df[[pos_col]])
  p <- plot_main(df)
  ggsave(file.path(opts$options$output, "all.pdf"), p, width = 22, height = 33)
} else if (opts$options$split == "trna") {
  l <- split(df, df$trna)
  for (trna in names(l)) {
    df <- l[[trna]]

    df$Ref <- droplevels(df$Ref)
    df[[pos_col]] <- droplevels(df[[pos_col]])
    df <- tidyr::complete(df, Ref, .data[[pos_col]])
    df[is.na(df$mod), "mod"] <- ""
    df[is.na(df$score), "score"] <- NA

    p <- plot_main(df)
    ggsave(file.path(opts$options$output, paste0(trna, ".pdf")), p, width = 22, height = 11)
  }
}

warnings()
