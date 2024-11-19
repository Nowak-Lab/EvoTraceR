library(ggplot2)
library(gridExtra)
library(ggrepel)

# Ensure directory exists
prepare_directory <- function(figure_dir) {
  if (!dir.exists(figure_dir)) {dir.create(figure_dir)}
  return(figure_dir)
}


# 1. Unique ASV Count
plot_unique_asv <- function(
  EvoTraceR_object,
  figure_dir
  )
{
  figure_dir <- prepare_directory(figure_dir)
  seqtab_history <- EvoTraceR_object$seqtab_history
  
  # Calculate unique ASV counts and percentage change
  unique_asv_animal <- data.frame(
    name = factor(names(seqtab_history), levels = c("Starting Seqs", "Hamming Merging", "Flanking Seq. Filter", "Substitutions Merging", "Frequency Filter", "Final ASVs")),
    num = sapply(seqtab_history, nrow)
  )
  unique_asv_animal <- unique_asv_animal %>%
    mutate(diff_perc = ifelse(row_number() == 1, "", 
                             paste0(ceiling((num / lag(num) - 1) * 100), " %"))) # Calculate percentage change
  
  max_value <- max(unique_asv_animal$num) * 1.3
  
  bar_plot <- ggplot(unique_asv_animal, aes(x = name, y = num, fill = name)) +
    geom_bar(stat = "identity", width = 0.8) +
    geom_text(aes(label = format(num, big.mark = ",", scientific = FALSE)), vjust = -1.8, size = 3) + # Main number above
    geom_text(aes(label = diff_perc), vjust = -0.3, size = 3, color = "blue") + # Percentage below main number
    labs(title = "Unique ASV Count", x = "Filtering Steps", y = "Unique ASV Count") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_value), labels = function(x) format(x, big.mark = ",", scientific = FALSE)) + # Format y-axis with commas
    scale_fill_manual(values = c("#444c5c", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#78a5a3")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
    ggsave(filename = file.path(figure_dir, "unique_asv_count.pdf"), plot = bar_plot, width = 15, height = 10, units = "cm")
}

# 2. Total ASV Count per Animal
plot_total_asv <- function(
  EvoTraceR_object,
  figure_dir
) {
  figure_dir <- prepare_directory(figure_dir)
  seqtab_history <- EvoTraceR_object$seqtab_history
  sample_columns <- setdiff(colnames(seqtab_history[[1]]), c("seq_names", "seq", "totalCounts"))

  total_asv_animal <- data.frame(
    name = factor(names(seqtab_history), levels = c("Starting Seqs", "Hamming Merging", "Flanking Seq. Filter", "Substitutions Merging", "Frequency Filter", "Final ASVs")),
    num = sapply(seqtab_history, function(df) sum(rowSums(df[, sample_columns])))
  )
  
  # Calculate percentage change
  total_asv_animal <- total_asv_animal %>%
    mutate(diff_perc = ifelse(row_number() == 1, "", paste0(ceiling((num / lag(num) - 1) * 100), " %")))
  
  max_value <- max(total_asv_animal$num) * 1.3
  
  bar_plot <- ggplot(total_asv_animal, aes(x = name, y = num, fill = name)) +
    geom_bar(stat = "identity", width = 0.8) +
    geom_text(aes(label = format(num, big.mark = ",", scientific = FALSE)), vjust = -1.8, size = 3) + # Main number above
    geom_text(aes(label = diff_perc), vjust = -0.3, size = 3, color = "blue") + # Percentage below main number
    labs(title = "Total ASV Count", x = "Filtering Steps", y = "Total ASV Count") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_value), labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    scale_fill_manual(values = c("#444c5c", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#78a5a3")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  ggsave(filename = file.path(figure_dir, "total_asv_count.pdf"), plot = bar_plot, width = 15, height = 10, units = "cm")
}

# 3. Unique ASV Count per Organ with optional donut plots for contamination
plot_unique_asv_per_organ <- function(
  EvoTraceR_object,
  figure_dir
) {
  figure_dir <- prepare_directory(figure_dir)
  seqtab_history <- EvoTraceR_object$seqtab_history
  sample_columns <- setdiff(colnames(seqtab_history[[1]]), c("seq_names", "seq", "totalCounts"))

  unique_asv_organ <- do.call(rbind, lapply(names(seqtab_history), function(step_name) {
    df <- seqtab_history[[step_name]]
    data.frame(
      name = step_name,
      organ = sample_columns,
      unique_count = colSums(df[, sample_columns] > 0)
    )
  }))
  
  unique_asv_organ$name <- factor(unique_asv_organ$name, levels = c("Starting Seqs", "Hamming Merging", "Flanking Seq. Filter", "Substitutions Merging", "Frequency Filter", "Final ASVs"))
  
  unique_asv_organ <- unique_asv_organ %>%
    group_by(organ) %>%
    mutate(diff_perc = ifelse(row_number() == 1, "", paste0(ceiling((unique_count / lag(unique_count) - 1) * 100), " %")))
  
  max_value <- max(unique_asv_organ$unique_count) * 1.3
  
  bar_plot <- ggplot(unique_asv_organ, aes(x = name, y = unique_count, fill = name)) +
    geom_bar(stat = "identity", width = 0.8) +
    geom_text(aes(label = format(unique_count, big.mark = ",", scientific = FALSE)), vjust = -1.8, size = 3) +
    geom_text(aes(label = diff_perc), vjust = -0.3, size = 3, color = "blue") +
    facet_wrap(~organ, scales = "free_y") +
    labs(title = "Unique ASV Count per Organ", x = "Filtering Steps", y = "Unique ASV Count") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_value)) +
    scale_fill_manual(values = c("#444c5c", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#78a5a3")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
    ggsave(filename = file.path(figure_dir, "unique_asv_count_per_organ.pdf"), plot = bar_plot, width = 15, height = 10, units = "cm")
}
plot_total_asv_per_organ <- function(
  EvoTraceR_object, 
  figure_dir
  )
  {
  figure_dir <- prepare_directory(figure_dir)
  seqtab_history <- EvoTraceR_object$seqtab_history
  sample_columns <- setdiff(colnames(seqtab_history[[1]]), c("seq_names", "seq", "totalCounts"))

  total_asv_organ <- do.call(rbind, lapply(names(seqtab_history), function(step_name) {
    df <- seqtab_history[[step_name]]
    organ_totals <- colSums(df[, sample_columns])
    data.frame(
      name = step_name,
      organ = names(organ_totals),
      total_count = organ_totals
    )
  }))
  
  total_asv_organ$name <- factor(total_asv_organ$name, levels = c("Starting Seqs", "Hamming Merging", "Flanking Seq. Filter", "Substitutions Merging", "Frequency Filter", "Final ASVs"))
  
  # Calculate percentage change
  total_asv_organ <- total_asv_organ %>%
    group_by(organ) %>%
    mutate(diff_perc = ifelse(row_number() == 1, "", paste0(ceiling((total_count / lag(total_count) - 1) * 100), " %")))
  
  max_value <- max(total_asv_organ$total_count) * 1.3
  
  bar_plot <- ggplot(total_asv_organ, aes(x = name, y = total_count, fill = name)) +
    geom_bar(stat = "identity", width = 0.8) +
    geom_text(aes(label = format(total_count, big.mark = ",", scientific = FALSE)), vjust = -1.8, size = 3) + # Main count above
    geom_text(aes(label = diff_perc), vjust = -0.3, size = 3, color = "blue") + # Percentage below main count
    facet_wrap(~organ, scales = "free_y") +
    labs(title = "Total ASV Count per Organ", x = "Filtering Steps", y = "Total ASV Count") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_value), labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    scale_fill_manual(values = c("#444c5c", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#78a5a3")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
    ggsave(filename = file.path(figure_dir, "total_asv_count_per_organ.pdf"), plot = bar_plot, width = 15, height = 10, units = "cm")
}

# Wrapper function to generate all ASV plots at once
generate_all_asv_plots <- function(EvoTraceR_object, figure_dir = EvoTraceR_object$output_directory) {
  print("Generating ASV count plots")
  plot_unique_asv(EvoTraceR_object, figure_dir)
  plot_total_asv(EvoTraceR_object, figure_dir)
  plot_unique_asv_per_organ(EvoTraceR_object, figure_dir)
  plot_total_asv_per_organ(EvoTraceR_object, figure_dir)
}