options(error = function() traceback(3))

`%nin%` = Negate(`%in%`)
library(optparse)
library(ggplot2)
library(dplyr)

option_list <- list(
  make_option(c("--score_columns"),
              type = "character",
              help = "Score column to use"),
  make_option(c("--output"),
              type = "character",
              help = "Output"),
  make_option(c("--object"),
              type = "character",
              help = "R"),
  make_option(c("--names"),
              type = "character",
              help = "Names"),
  make_option(c("--position_column"),
              type = "character",
              default = "seq_position",
              help = "Use column as coordinates, default: seq_position"),
  make_option(c("--title"),
              type = "character",
              help = "Title of the plot"),
  make_option(c("--remove_prefix"),
              type = "character",
              default = "",
              help = "Remove prefix from tRNA label"),
  make_option(c("--amino_acids"),
              type = "character",
              help = "Amino acids"),
  make_option(c("--merge"),
              type = "character",
              help = "Merge")
)

args_mesc <- c("--score_columns=MDI,MDI_subsampled,MDI_subsampled,MDI_subsampled,MDI_subsampled",
              "--output=~/slides/DFG/plots/mouse-QTRT1/mESC.pdf",
              "--names=raw,subsampled,trim-cigar,adapter-overlap,remove-multimapper",
              "--position_column=sprinzl",
              "--remove_prefix=Mus_musculus_tRNA-",
              "--amino_acids=Asn,Asp,His,Tyr",
              "--merge=other",
              "/beegfs/prj/tRNA_Francesca_Tuorto/data/20240312_FT_mESC_mHC_tRNA_RNA002_custom/qutrna/output-final-intermediate-plots/results/jacusa2/cond1~mESC-wt/cond2~mESC-Qtrt1/bams~samtools/scores_sprinzl.tsv",
              "/beegfs/prj/tRNA_Francesca_Tuorto/data/20240312_FT_mESC_mHC_tRNA_RNA002_custom/qutrna/output-final-intermediate-plots/results/jacusa2/cond1~mESC-wt/cond2~mESC-Qtrt1/bams~samtools/scores_sprinzl.tsv",
              "/beegfs/prj/tRNA_Francesca_Tuorto/data/20240312_FT_mESC_mHC_tRNA_RNA002_custom/qutrna/output-final-intermediate-plots/results/jacusa2/cond1~mESC-wt/cond2~mESC-Qtrt1/bams~trim_cigar/scores_sprinzl.tsv",
              "/beegfs/prj/tRNA_Francesca_Tuorto/data/20240312_FT_mESC_mHC_tRNA_RNA002_custom/qutrna/output-final-intermediate-plots/results/jacusa2/cond1~mESC-wt/cond2~mESC-Qtrt1/bams~overlap//scores_sprinzl.tsv",
              "/beegfs/prj/tRNA_Francesca_Tuorto/data/20240312_FT_mESC_mHC_tRNA_RNA002_custom/qutrna/output-final-intermediate-plots/results/jacusa2/cond1~mESC-wt/cond2~mESC-Qtrt1/bams~remove_multimappers/scores_sprinzl.tsv")

args_mhc <- c("--score_columns=MDI,MDI_subsampled,MDI_subsampled,MDI_subsampled,MDI_subsampled",
              "--output=~/slides/DFG/plots/mouse-QTRT1/mHC.pdf",
              "--names=raw,subsampled,trim-cigar,adapter-overlap,remove-multimapper",
              "--position_column=sprinzl",
              "--remove_prefix=Mus_musculus_tRNA-",
              "--amino_acids=Asn,Asp,His,Tyr",
              "--merge=other",
              "/beegfs/prj/tRNA_Francesca_Tuorto/data/20240312_FT_mESC_mHC_tRNA_RNA002_custom/qutrna/output-final-intermediate-plots/results/jacusa2/cond1~mHC-wt/cond2~mHC-Qtrt1/bams~samtools/scores_sprinzl.tsv",
              "/beegfs/prj/tRNA_Francesca_Tuorto/data/20240312_FT_mESC_mHC_tRNA_RNA002_custom/qutrna/output-final-intermediate-plots/results/jacusa2/cond1~mHC-wt/cond2~mHC-Qtrt1/bams~samtools/scores_sprinzl.tsv",
              "/beegfs/prj/tRNA_Francesca_Tuorto/data/20240312_FT_mESC_mHC_tRNA_RNA002_custom/qutrna/output-final-intermediate-plots/results/jacusa2/cond1~mHC-wt/cond2~mHC-Qtrt1/bams~trim_cigar/scores_sprinzl.tsv",
              "/beegfs/prj/tRNA_Francesca_Tuorto/data/20240312_FT_mESC_mHC_tRNA_RNA002_custom/qutrna/output-final-intermediate-plots/results/jacusa2/cond1~mHC-wt/cond2~mHC-Qtrt1/bams~overlap//scores_sprinzl.tsv",
              "/beegfs/prj/tRNA_Francesca_Tuorto/data/20240312_FT_mESC_mHC_tRNA_RNA002_custom/qutrna/output-final-intermediate-plots/results/jacusa2/cond1~mHC-wt/cond2~mHC-Qtrt1/bams~remove_multimappers/scores_sprinzl.tsv")
opts <- parse_args(
  OptionParser(option_list = option_list),
  args = args_mhc,
  positional_arguments = TRUE
)

stopifnot(!is.null(opts$options$output))
stopifnot(length(opts$args) >= 1)
names <- strsplit(opts$options$names, ",") |> unlist()
stopifnot(length(opts$args) == length(names))
score_cols <- strsplit(opts$options$score_columns, ",") |> unlist()
stopifnot(length(opts$args) == length(score_cols))

amino_acids <- strsplit(opts$options$amino_acids, ",") |> unlist()

# title
# merge

trna_label <- function(df, remove_prefix = "") {
  new_label <- df$trna
  if (remove_prefix != "") {
    new_label <- gsub(remove_prefix, "", new_label)
  }
  
  return(new_label)
}

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

l <- mapply(function(f, n, s) {
  df <- read.table(f, header = TRUE)
  df$score <- df[[s]]
  df$position <- df[[opts$options$position_column]]
  df$name <- n
  df <- df |>
    filter(position != ".")
  
  df <- add_trna_details(df)
  if (!is.null(opts$options$remove_prefix)) {
    df$trna_label <- trna_label(df, opts$options$remove_prefix)
  } else {
    df$trna_label <- df$trna
  }
  
  df <- df |>
    select(name, trna_label, position, score, amino_acid, anti_codon)
  
  return(df)
}, f = opts$args, n = names, s = score_cols, SIMPLIFY = FALSE)
df <- bind_rows(l)
df$name <- factor(df$name, levels = unique(df$name), ordered = TRUE)

# plot
#geom_point(position = position_jitter(seed = 1, width = 0.2)) +
#filter(amino_acid %in% amino_acids) |>

i <- df$amino_acid %in% amino_acids

df <- df |>
  mutate(
    trna_label2 = case_when(
      amino_acid %in% amino_acids ~ trna_label,
      amino_acid %nin% amino_acids ~ "other"),
    amino_acid_label = case_when(
      amino_acid %in% amino_acids ~ amino_acid,
      amino_acid %nin% amino_acids ~ "other"),
    amino_acid_label = factor(amino_acid_label, levels = c(amino_acids, "other"), ordered = TRUE),
    trna_label2 = factor(trna_label2, levels = c("other", rev(unique(trna_label))), ordered = TRUE),
    name_label = gsub("-", "\n", name),
    name_label = factor(df$name, levels = unique(df$name), ordered = TRUE))

max_scores <- df |> filter(trna_label2 == "other") |>
  summarise(.by = name_label,
            median_score = max(score),
            mean_score = max(score),
            max_score = max(score))

p <- df |>
  ggplot(aes(x = trna_label2, y = score, fill = amino_acid_label)) +
  labs(fill = "Amino acid") +
  xlab("") +
  ylab("score") +
  geom_violin(scale = "width") +
  geom_hline(
    mapping = aes(yintercept = max_score),
    linetype = "dashed",
    colour = "gray",
    data = max_scores) +
  facet_wrap(name_label ~ ., ncol = 1, strip.position = "right") +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, size = 8, hjust = 1))
