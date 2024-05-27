#!/usr/bin/env Rscript

options(error = traceback)
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
  make_option(c("--scale_by_cov"),
              action = "store_true",
              default = FALSE,
              help = "Scale tiles by coverage"),
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
  make_option(c("--flag"),
              type = "character",
              default = NULL,
              help = "Flag positions, e.g.: tRNA:pos, use ','"),
  make_option(c("--positions"),
              type = "character",
              default = NULL,
              help = "Positions to show, use ','")
)

opts <- parse_args(
  OptionParser(option_list = option_list),
  # args = c("--cond1=wt",
  #         "--cond2=NSUN2",
  #         "--split=isoacceptor",
  #         "--column=sprinzl",
  #         "--show_introns",
  # #          "--title=isodecoder: {amino_acid}-{anti_codon}",
  #         "--sort",
  # #         "--show_coverage",
  #         "--left=23",
  #         "--scale_by_cov",
  #         "--modmap=/beegfs/homes/mpiechotta/git/QutRNA/data/Hsapi38/human_mods_map.tsv",
  # #          #"--hide_mods",
  # #         "--crop",
  #         "--output=~/tmp/plot_score",
  #        "--remove_prefix=Homo_sapiens_",
  # #          #"--remove_prefix=Mus_musculus_",
  #         "scores_sprinzl-mods.tsv"),
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

df$flag <- FALSE
if (!is.null(opts$options$flag)) {
  df$tmp_coords <- paste0(df$Ref, ":", df[[pos_col]])
  flag <- strsplit(opts$options$flag, ",")[[1]]
  i <- df$tmp_coords %in% flag
  df[i, opts$options$score] <- NA
  df[i, "flag"] <- TRUE
  df$tmp_coords <- NULL
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
cols <- c("Ref", pos_col, "anti_codon", "amino_acid", "score", "mod", "Pos3", "ref_base", "total_coverage", "median_coverage", "expression", "flag")
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

  if (opts$options$scale_by_cov) {
    s <- "median cov. Q1\n(1st. quartile)"
    
    heights <- df |>
      distinct(Ref, median_coverage) |>
      filter(!is.na(median_coverage)) |>
      mutate(quartile = paste0("Q", ntile(median_coverage, 4))) |>
      group_by(quartile) |>
      arrange(Ref) |>
      mutate(height = median_coverage / sum(median_coverage, na.rm = TRUE),
             y_end = cumsum(height),
             y_pos = y_end - height / 2) |>
      ungroup() |>
      select(-median_coverage) |>
      mutate(
        quartile = case_when(
          quartile == "Q1" ~ s,
          .default = quartile),
        quartile = factor(quartile, levels = c(paste0("Q", 4:2), s)))

    df <- df |>
      left_join(heights, by = join_by(Ref))

    helper_y_pos <- function(q) {
      heights |>
        filter(quartile == q) |>
        pull(y_pos)
    }
    helper_Ref <- function(q) {
      heights |>
        filter(quartile == q) |>
        pull(Ref)
    }

    p <- df |> 
      ggplot(aes(x = .data[[pos_col]], y = y_pos, fill = score)) +
      facet_grid(quartile ~ ., scales = "free_y", space = "free_y") +
      geom_tile(aes(height = height), width = 1, colour = "white") +
      ylab("median cov.")
      
      q = unique(heights$quartile)
      n = length(q)
      
      l <- list()
      l <- append(l, quartile == s ~ scale_y_continuous(breaks = helper_y_pos(s), labels = helper_Ref(s)))
      if (n > 1) {
        l <- append(l, quartile == "Q2" ~ scale_y_continuous(breaks = helper_y_pos("Q2"), labels = helper_Ref("Q2")))
      }
      if (n > 2) {
        l <- append(l, quartile == "Q3" ~ scale_y_continuous(breaks = helper_y_pos("Q3"), labels = helper_Ref("Q3")))
      }
      if (n > 3) {
        l <- append(l, quartile == "Q4" ~ scale_y_continuous(breaks = helper_y_pos("Q4"), labels = helper_Ref("Q4")))
      }
      # FIXME less than 4 scales
      p <- p + ggh4x::facetted_pos_scales(y = l)
  } else {
    p <- df |> 
      ggplot(aes(x = .data[[pos_col]], y = Ref, fill = score)) +
      geom_tile(colour = "white", width = 1, height = 1) +
      ylab("") +
      coord_equal()
  }

  p <- p +
      scale_fill_gradient(low = "yellow", high = "blue", na.value = "grey") +
      theme_bw() +
      theme(panel.background = element_rect(fill = "white", color = "white"),
            legend.position = "top",
            panel.border = element_blank(),
            axis.text.x = element_text(angle = 45, size = 8),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(family = "mono"),
            plot.margin = margin(0, .5, 0, 0))

  if (opts$options$column == "sprinzl") {
    p <- p + xlab("Sprinzl pos.")
  } else {
    p <- p + xlab("pos. in tRNA seq.")
  }

  if (!is.null(opts$options$flag)) {
    if ("flag" %in% names(df)) {
      data <- df |>
        filter(flag == TRUE)
      p <- p + geom_point(data = data,
                          shape = 4, color = "darkred")
    }
  }

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

  if (opts$options$sort) {
    tmp <- sort_ref(df, cov)
    df <- tmp$df
    cov <- tmp$cov
  }
  
  if (opts$options$show_coverage) {
    p <- plot_complex(df, cov)
  } else {
    p <- plot_main(df)
    if (!is.null(opts$options$modmap)) {
      tbl <- plot_table(df)
      if (!is.null(tbl)) {
        p <- ((p + theme(plot.margin = margin(0, 20, 0, 0))) | tbl) +
          plot_layout(ncol = 3,
                      widths = unit(c(-1, 4 , 2), c("null", "cm", "cm"))) +
          theme(plot.margin = margin(0, 0, 0, 0))
      }
    }
  }

  if (opts$options$crop) {
    tmp <- file.path(opts$options$output, paste0(e, "_tmp.pdf"))
    ggsave(tmp, p, width = 33, height = 11)
    output <- file.path(opts$options$output, paste0(e, ".pdf"))
    system2("pdfcrop.pl", c("--margin=5", tmp, output))
    unlink(tmp)
  } else {
    output <- file.path(opts$options$output, paste0(e, ".pdf"))
    ggsave(output, p, width = 33, height = 11)
  }
}

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
  save_plot(df, cov, "all")
}

warnings()
