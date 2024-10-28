#!/usr/bin/env Rscript

# Plot JACUSA2 score and existing modification info as a heatmap)

library(optparse)
library(ggplot2)
library(patchwork)
library(magrittr)
library(dplyr)


option_list <- list(
  make_option(c("--score"),
              type = "character",
              default = "Mis.Del.Ins",
              help = "Score column to use"),
  make_option(c("--cond1"),
              type = "character",
              default = "condition 1",
              help = "Condition 1"),
  make_option(c("--cond2"),
              type = "character",
              default = "condition 2",
              help = "Condition 2"),
  make_option(c("--output"),
              type = "character",
              help = "Output"),
  make_option(c("--left"),
              type = "numeric",
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
  make_option(c("--hide_mods"),
              action = "store_true",
              default = FALSE,
              help = "Hide mods"),
  make_option(c("--sort"),
              action = "store_true",
              default = FALSE,
              help = "Sort by coverage,tRNA"),
  make_option(c("--show_introns"),
              action = "store_true",
              default = FALSE,
              help = "Show intron positions"),
  make_option(c("--show_coverage"),
              action = "store_true",
              default = FALSE,
              help = "Show coverage"),
  make_option(c("--show_ref"),
              action = "store_true",
              default = FALSE,
              help = "Show reference base"),
  make_option(c("--crop"),
              action = "store_true",
              default = FALSE,
              help = "Do  crop pdf"),
  make_option(c("--modmap"),
              type = "character",
              help = "File mapping of modomics short name to abbrevations"),
  make_option(c("--split"),
              type = "character",
              help = "Split by: isodecoder, isoacceptor"),
  make_option(c("--remove_prefix"),
              type = "character",
              default = "",
              help = "Remove prefix from tRNA"),
  make_option(c("--positions"),
              type = "character",
              default = NULL,
              help = "Positions to show, use ','")
)

opts <- parse_args(
  OptionParser(option_list = option_list),
  # args = c("--cond1=wt",
  #          "--cond2=test",
  #          "--split=isodecoder",
  #          "--column=sprinzl",
  # #          "--show_introns",
  # #          "--title=isodecoder: {amino_acid}-{anti_codon}",
  # #          "--sort",
  # #          "--show_coverage",
  #          "--left=23",
  #            "--modmap=/beegfs/homes/mpiechotta/git/QutRNA/data/Hsapi38/human_mods_map.tsv",
  # #          #"--hide_mods",
  # #         "--crop",
  #          "--output=~/tmp/plot_score",
  # #          "--remove_prefix=Homo_sapiens_",
  # #          #"--remove_prefix=Mus_musculus_",
  #          "scores_sprinzl-mods.tsv"),
  positional_arguments = TRUE
)
stopifnot(!is.null(opts$options$output))
stopifnot(length(opts$args) == 1)

df <- read.table(opts$args, sep = "\t", header = TRUE)
df <- df[df$Pos3 > opts$options$left, ]

if (!"mod" %in% colnames(df)) {
  df$mod <- ""
}

pos_col <- opts$options$column
if (pos_col == "sprinzl") {
  df <- df[!df[[pos_col]] %in% c("-", "."), ]

  coords <- c(1:17, "17a", 18:20, "20a", "20b", 21:37,
              paste0("i", 1:24),
              38:45,
              paste0("e", c(11:17, 1:5, 27:21)), 46:76)

  if (opts$options$hide_varm) {
    i <- !grepl("e[0-9]+", df[[pos_col]])
    df <- df[i, ]
    i <- !grepl("e[0-9]+", coords)
    coords <- coords[i]
  }
  if (!opts$options$show_introns) {
    i <- !grepl("i[0-9]+", df[[pos_col]])
    df <- df[i, ]
    i <- !grepl("i[0-9]+", coords)
    coords <- coords[i]
  }
  df[[pos_col]] <- factor(df[[pos_col]], levels = coords, ordered = TRUE)
} else {
  df[[pos_col]] <- factor(df[[pos_col]],
                          levels = stringr::str_sort(unique(df[[pos_col]]), numeric = TRUE), ordered = TRUE)
}

df$trna <- df$Ref
i <- grepl("tRNA", df$Ref)
df$amino_acid <- ""
df[i, "amino_acid"] <- stringr::str_extract(df[i, "Ref"], ".*tRNA-([A-Za-z]+)-([A-Za-z]{3}).*", group = 1)
df$anti_codon <- ""
df[i, "anti_codon"] <- stringr::str_extract(df[i, "Ref"], ".*tRNA-([A-Za-z]+)-([A-Za-z]{3}).*", group = 2)

cov <- df %>%
  select(Ref, amino_acid, anti_codon, starts_with("coverage_")) %>%
  group_by(Ref, amino_acid, anti_codon) %>%
  summarise_all(median) %>%
  ungroup() %>%
  tidyr::pivot_longer(starts_with("coverage_"),
                      names_to = "condition",
                      values_to = "coverage") %>%
  mutate(replicate = gsub("^coverage_\\d+_", "", condition),
                condition = gsub("^coverage_|_\\d+$", "", condition),
                condition = case_match(
                  condition,
                  "1" ~ opts$options$cond1,
                  "2" ~ opts$options$cond2,
                ))

cov_summary <- cov %>%
  select(Ref, coverage) %>%
  summarise(median_coverage = median(coverage),
            total_coverage = sum(coverage),
            expression = median_coverage,
            .by = Ref)

df <- df %>%
  full_join(cov_summary, by = join_by(Ref))
cov <- cov %>%
  full_join(cov_summary, by = join_by(Ref))

if (!is.null(opts$options$remove_prefix)) {
  df <- df %>% mutate(Ref = gsub(opts$options$remove_prefix, "", Ref))
  cov <- cov %>% mutate(Ref = gsub(opts$options$remove_prefix, "", Ref))
}
df[["Ref"]] <- factor(df[["Ref"]],
                      levels = stringr::str_sort(unique(df[["Ref"]]), numeric = TRUE, decreasing = TRUE))
df[["ref_base"]] <- substr(df[["Kmer"]], 3, 3)

df$score <- df[, opts$options$score]
df$score[df$score < 0] <- 0
cols <- c("Ref", pos_col, "anti_codon", "amino_acid", "score", "mod", "Pos3", "ref_base", "total_coverage", "median_coverage", "expression")
if (!is.null(opts$options$modmap)) {
  modmap <- read.table(opts$options$modmap, header = TRUE, sep = "\t", quote = "", comment.char = "")
  stopifnot(length(modmap$shortname) == length(unique(modmap$shortname)))
  df <- merge(df, modmap, by.x = "mod", by.y = "short_name", all.x = TRUE, sort = FALSE)
  i <- !is.na(df$abbrev)
  df$short_name <- df$mod
  df[i, "mod"] <- df[i, "abbrev"]
  cols <- c(cols, "short_name", "abbrev")
}
df <- select(df, tidyr::all_of(cols))

add_missing <- function(df) {
  helper <- function(e) {
    unique(e[!is.na(e)])
  }

  tidyr::complete(df, Ref, .data[[pos_col]]) %>%
    mutate(expression = helper(expression), .by = Ref) %>%
    tidyr::replace_na(list(mod = ""))
}

if (!is.null(opts$options$positions)) {
  positions = strsplit(opts$options$position, ",")[[1]]
  i <- df[[pos_col]] %in% positions
  df <- df[i, ]
}

plot_table <- function(df) {
  df <- select(df, short_name, abbrev) |>
      filter(!is.na(abbrev)) |>
      filter(short_name != abbrev) |>
      distinct() |>
      select(abbrev, short_name) |>
      rename(`RNAMods\ncode` = abbrev, `Mod.` = short_name)

  if (nrow(df) > 0) {
    return(gridExtra::tableGrob(df, rows = NULL))
  }

  return(NULL)
}

plot_main <- function(df) {
  sprinzl <- levels(df$sprinzl)
  vjust <- rep(0, length(sprinzl))
  vjust[grep("^(i|e)", sprinzl)] <- 0.25

  p <- df |>
    ggplot(aes(x = .data[[pos_col]], y = Ref)) +
      ylab("") +
      scale_fill_gradient(low = "yellow", high = "blue", na.value = "grey") +
      theme_bw() +
      theme(panel.background = element_rect(fill = "white", color = "white"),
            legend.position = "top",
            panel.border = element_blank(),
            axis.text.x = element_text(angle = 45),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(family = "mono"),
            plot.margin = margin(0, .5, 0, 0))

  if (opts$options$column == "sprinzl") {
    p <- p + xlab("Sprinzl pos.")
  } else {
    p <- p + xlab("pos. in tRNA seq.")
  }

  p <- p + geom_tile(aes(fill = score), colour = "white", width = 1, height = 1)
  if (opts$options$show_ref) {
      p <- p + geom_text(aes(label = ref_base))
      if (!opts$options$hide_mods) {
        p <- p + geom_text(aes(label = mod), position = position_nudge(x = 0.5, y = 0.5))
      }
  } else {
      if (!opts$options$hide_mods) {
        p <- p + geom_text(aes(label = mod))
      }
  }
  p <- p + coord_equal()

  if (!is.null(opts$options$title)) {
    title <- gsub("\\{anti_codon\\}", paste0(unique(na.omit(df$anti_codon)), collapse = ", "), opts$options$title)
    title <- gsub("\\{amino_acid\\}", paste0(unique(na.omit(df$amino_acid)), collapse = ", "), title)
    p <- p + ggtitle(title)
  }

  p
}

plot_coverage <- function(cov, max_cov = NULL) {
  limits = NULL
  if (!is.null(max_cov)) {
    limits <- c(0, max_cov)
  }

  p <- cov |>
    ggplot(aes(x = coverage, y = Ref, fill = condition, colour = condition)) +
    scale_colour_manual(values = c("#bf568b", "#bcd357")) +
    scale_fill_manual(values = c("#bf568b", "#bcd357")) +
    scale_y_discrete(labels = NULL) +
    scale_x_continuous(breaks = scales::breaks_pretty(3),
                       limits = limits,
                       labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    ylab("") +
    xlab("median cov.") +
    theme_bw() +
    theme(plot.margin = margin(0, 0, 0, 0),
          legend.position = "top",
          legend.title = element_blank(),
          legend.direction = "vertical")

  if (unique(cov$replicate) %>% length() == 1) {
    p <- p + geom_col(width = 0.8, position = "dodge")
  } else {
    p <- p + geom_point(position = position_dodge2(width = 0.5))
  }

  p
}

sort_ref <- function(df, cov) {
  tmp_cov <- cov[, c("Ref", "total_coverage")] %>%
    distinct()
  o <- order(tmp_cov$total_coverage)
  ref <- tmp_cov$Ref[o]

  df$Ref <- factor(df$Ref, levels = ref, ordered = TRUE)
  cov$Ref <- factor(cov$Ref, levels = ref, ordered = TRUE)

  list(df = df, cov = cov)
}
#
plot_complex <- function(df, cov) {
  padding_helper <- function(s) {
    stringr::str_pad(s, max(nchar(s)), side = "right", pad = "-")
  }

  if (opts$options$sort) {
    tmp <- sort_ref(df, cov)
    df <- tmp$df
    cov <- tmp$cov
  }

  p1 <- plot_main(df) +
    scale_y_discrete(labels = padding_helper, position = "right")
  p2 <- plot_coverage(cov)

  if (!is.null(opts$options$modmap)) {
    tbl <- plot_table(df)
    if (!is.null(tbl)) {
      p <- (p1 | (p2 + theme(plot.margin = margin(0, 20, 0, 0))) | tbl) +
        plot_layout(ncol = 3,
                    widths = unit(c(-1, 4 , 4), c("null", "cm", "cm"))) +
        theme(plot.margin = margin(0, 0, 0, 0))
    } else {
      p <- (p1 | p2) +
        plot_layout(ncol = 2,
                    widths = unit(c(-1, 4), c("null", "cm"))) +
        theme(plot.margin = margin(0, 0, 0, 0))
    }
  } else {
    p <- p1 | p2 +
      plot_layout(ncol = 2,
                  widths = unit(c(-1, 4), c("null", "cm"))) +
      theme(plot.margin = margin(0, 0, 0, 0))
  }

  p
}


save_plot <- function(df, cov, e) {
  df$Ref <- droplevels(df$Ref)
  df[[pos_col]] <- droplevels(df[[pos_col]])
  df <- add_missing(df)

  if (opts$options$show_coverage) {
    p <- plot_complex(df, cov)
  } else {
    p <- plot_main(df)
    if (!is.null(opts$options$modmap)) {
      tbl <- plot_table(df)
      if (!is.null(tbl)) {
        p <- ((p + theme(plot.margin = margin(0, 20, 0, 0))) | tbl) +
          plot_layout(ncol = 3,
                      widths = unit(c(-1, 4 , 4), c("null", "cm", "cm"))) +
          theme(plot.margin = margin(0, 0, 0, 0))
      }
    }
  }

  if (opts$options$crop) {
    ggsave(tmp, p, width = 33, height = 11)
    tmp <- file.path(opts$options$output, paste0(e, "_tmp.pdf"))
    output <- file.path(opts$options$output, paste0(e, ".pdf"))
    system2("pdfcrop.pl", c("--margin=5", tmp, output))
    unlink(tmp)
  } else {
    output <- file.path(opts$options$output, paste0(e, ".pdf"))
    ggsave(output, p, width = 33, height = 11)
  }
}

# TODO what abot no isoacceptor or isodecoder
if (opts$options$split == "isoacceptor") {
  l1 <- split(df, df$amino_acid)
  l2 <- split(cov, cov$amino_acid)
  for (aa in names(l1)) {
    save_plot(l1[[aa]], l2[[aa]], aa)
  }
} else if (opts$options$split == "isodecoder") {
  l1 <- split(df, paste0(df$amino_acid, "-", df$anti_codon))
  l2 <- split(cov, paste0(cov$amino_acid, "-", cov$anti_codon))
  for (isodec in names(l1)) {
    save_plot(l1[[isodec]], l2[[isodec]], isodec)
  }
} else {
  stop("Unknown option: --split")
}

warnings()
