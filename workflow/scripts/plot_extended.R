# TODO remove

plot_extended <- function(df, coverage_summary, plot_args, condition1, condition2) {
  title <- df$trna_label |>
    unique() |>
    as.character()
  
  vlines <- data.frame(x = c(1, 10, 20, 30, 40, 50, 60, 70))
  p_score <- df |>
    mutate(trna_label = "Score") |>
    plot_heatmap(xlab = plot_args[["position_xlab"]],
                 harmonize_scaling = plot_args[["harmonize_scaling"]]) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 10, colour = "black" ),
          axis.title.x = element_blank())
  
  # p_read_cov <- df |> 
  #   select(position, starts_with("reads_")) |>
  #   tidyr::pivot_longer(cols = starts_with("reads_"), values_to = "reads") |>
  #   mutate(replicate = gsub("^reads_\\d+_", "", name),
  #          condition = gsub("^reads_|_\\d+$", "", name),
  #          condition_i = condition,
  #          condition = case_match(condition,
  #                                 "1" ~ condition1,
  #                                 "2" ~ condition2),
  #          condition = factor(condition, levels = c(condition1, condition2))) |>
  #   ggplot(aes(x = position, y = reads, group = name, colour = condition, fill = condition)) +
  #   ylab("Read coverage") +
  #   geom_line() +
  #   geom_vline(data = vlines, mapping = aes(xintercept = x), linewidth = 0.25, colour = "gray", linetype = 2) +
  #   scale_y_continuous(breaks = scales::breaks_pretty(3),
  #                      labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  #   theme_bw() +
  #   theme(panel.background = element_rect(fill = "white", color = "white"),
  #         panel.border = element_blank(),
  #         axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5, hjust = 1),
  #         axis.title.x = element_blank(),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         plot.margin = margin(0, .5, 0, 0)) +
  #   ggtitle(title)
  
  p_cov <- df |>
    select(position, starts_with("ref_base_calls_"), starts_with("nonref_base_calls_"), starts_with("insertions_"), starts_with("deletions_")) |>
    tidyr::pivot_longer(cols = c(starts_with("ref_base_calls_"), starts_with("nonref_base_calls_"), starts_with("insertions_"), starts_with("deletions_")), values_to = "count") |>
    mutate(replicate = gsub("^(ref_base_calls|nonref_base_calls|insertions|deletions)_\\d+_", "", name),
           condition = gsub("^(ref_base_calls|nonref_base_calls|insertions|deletions)_|_\\d+$", "", name),
           type = gsub("(_\\d+_\\d+)$", "", name),
           condition_i = condition,
           condition = case_match(condition,
                                  "1" ~ condition1,
                                  "2" ~ condition2),
           condition = factor(condition, levels = c(condition1, condition2)),
           sample = paste0(condition, "_", replicate)) |>
    ggplot(aes(x = position, y = count, group = sample, colour = type, fill = type)) +
    ylab("Base coverage") +
    geom_col() +
    geom_vline(data = vlines, mapping = aes(xintercept = x), linewidth = 0.25, colour = "gray", linetype = 2) +
    scale_y_continuous(breaks = scales::breaks_pretty(3),
                       labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white", color = "white"),
          panel.border = element_blank(),
          axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(0, .5, 0, 0)) +
    ggtitle(title) +
    facet_wrap(sample ~ ., ncol = 1, strip.position = "right")
  
  ratios_long <- df |> 
    select(position, starts_with("nonref_ratio_"), starts_with("insertion_ratio_"), starts_with("deletion_ratio_")) |>
    tidyr::pivot_longer(cols = c(starts_with("nonref_ratio_"), starts_with("insertion_ratio_"), starts_with("deletion_ratio_")), values_to = "ratio")
  ratio_details <- as.data.frame(stringr::str_match(ratios_long$name, "(nonref|insertion|deletion)_ratio_(\\d+)_(\\d+)"))[, -1]
  colnames(ratio_details) <- c("error", "condition", "replicate")
  ratios <- cbind(ratios_long, ratio_details)
  
  p_ratios <- ratios |>
    summarise(.by = c(position, error, condition),
              mean_ratio = mean(ratio)) |>
    summarise(.by = c(position, error),
              diff_mean_ratio = mean_ratio[condition == 1] - mean_ratio[condition == 2]) |>
    ggplot(aes(x = position, y = diff_mean_ratio, group = error, colour = error, fill = error)) +
    ylab(expression(Delta ~ Error)) +
    geom_line() +
    geom_abline(linetype = "dashed", colour = "gray", slope = 0, intercept = 0) +
    geom_vline(data = vlines, mapping = aes(xintercept = x), linewidth = 0.25, colour = "gray", linetype = 2) +
    scale_y_continuous(labels = scales::label_percent()) +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white", color = "white"),
          panel.border = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(0, .5, 0, 0))
  
  p_ref <- df |> 
    mutate(trna_label = "Reference") |>
    ggplot(aes(x = position, y = trna_label, label = ref)) +
    geom_text(size = 3) +
    coord_equal() +
    theme_bw() +
    xlab(plot_args[["position_xlab"]]) + 
    theme(panel.background = element_rect(fill = "white", color = "white"),
          legend.position = "none",
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 10, colour = "black" ),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5, hjust = 1),
          plot.margin = margin(0, .5, 0, 0))
  
  # TODO
  # add indels-sequence logo
  
  p <- 
    p_cov / 
    p_ratios /
    p_score / 
    p_ref + 
    plot_layout(nrow = 4, heights = unit(c(6, 3, -1, -1), c("cm", "cm", "null", "null"))) +
    plot_layout(guides = "collect")
  
  return(p)
}


split_extended <- function(df, coverage_summary, output_dir, plot_args, process_args, condition1, condition2) {
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
    tmp_df <- l_df[[trna]]
    tmp_coverage_summary <- l_coverage_summary[[trna]]
    if (!is.null(process_args$sort_by_callback)) {
      res <- process_args$sort_by_callback(tmp_df, tmp_coverage_summary)
      tmp_df <- res$df
      tmp_coverage_summary <- res$coverage_summary
    }
    p <- plot_extended(tmp_df, tmp_coverage_summary, plot_args, condition1, condition2)
    trna_label <- unique(tmp_df$trna_label)
    output <- file.path(output_dir, paste0(trna_label, ".pdf"))
    save_plot(p, output)
  }
}

split_extended(df, coverage_summary,
               opts$options$output_dir,
               get_plot_args(df, opts),
               get_process_args(opts), opts$options$condition1, opts$options$condition2)

split_extended <- function(df, coverage_summary, output_dir, plot_args, process_args, condition1, condition2) {
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
    tmp_df <- l_df[[trna]]
    tmp_coverage_summary <- l_coverage_summary[[trna]]
    if (!is.null(process_args$sort_by_callback)) {
      res <- process_args$sort_by_callback(tmp_df, tmp_coverage_summary)
      tmp_df <- res$df
      tmp_coverage_summary <- res$coverage_summary
    }
    p <- plot_extended(tmp_df, tmp_coverage_summary, plot_args, condition1, condition2)
    trna_label <- unique(tmp_df$trna_label)
    output <- file.path(output_dir, paste0(trna_label, ".pdf"))
    save_plot(p, output)
  }
}
