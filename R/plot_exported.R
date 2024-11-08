# Ensure directory exists
prepare_directory <- function(figure_dir) {
  figure_dir = paste0(figure_dir, '/asv_analysis/')
  if (!dir.exists(figure_dir)) {dir.create(figure_dir)}
  return(figure_dir)
}

# 1. Unique ASV Count per Animal (Total unique ASVs across all organs)
plot_unique_asv_per_animal <- function(EvoTraceR_object, figure_dir) {
  figure_dir <- prepare_directory(figure_dir)
  seqtab_history = EvoTraceR_object$seqtab_history
  
  unique_asv_animal = data.frame(
    name = factor(names(seqtab_history), levels = c("Starting ASVs", "Hamming Merging", "Flanking Seq. Filter", "Substitutions Merging", "Frequency Filter", "Final ASVs")),
    num = sapply(seqtab_history, nrow)
  )
  
  max_value <- max(unique_asv_animal$num) * 1.3  # Add 10% buffer to max value
  
  plot <- ggplot(unique_asv_animal, aes(x=name, y=num, fill=name)) +
    geom_bar(stat="identity", width=0.8) +
    geom_text(aes(label = num), vjust = -0.5, size = 3) +
    labs(title = "Unique ASV Count", x = "Filtering Steps", y = "Unique ASV Count") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_value)) +
    scale_fill_manual(values=c("#444c5c", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#78a5a3")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  ggsave(filename = file.path(figure_dir, "unique_asv_count.pdf"), plot = plot, width = 15, height = 10, units = "cm")
}

# 2. Total ASV Count per Animal (Sum across all organs)
plot_total_asv_per_animal <- function(EvoTraceR_object, figure_dir) {
  figure_dir <- prepare_directory(figure_dir)
  seqtab_history = EvoTraceR_object$seqtab_history
  sample_columns = setdiff(colnames(seqtab_history[[1]]), c("seq_names", "seq", "totalCounts"))

  total_asv_animal = data.frame(
    name = factor(names(seqtab_history), levels = c("Starting ASVs", "Hamming Merging", "Flanking Seq. Filter", "Substitutions Merging", "Frequency Filter", "Final ASVs")),
    num = sapply(seqtab_history, function(df) sum(rowSums(df[, sample_columns])))
  )
  
  max_value <- max(total_asv_animal$num) * 1.3  # Add 10% buffer to max value
  
  plot <- ggplot(total_asv_animal, aes(x=name, y=num, fill=name)) +
    geom_bar(stat="identity", width=0.8) +
    geom_text(aes(label = num), vjust = -0.5, size = 3) +
    labs(title = "Total ASV Count", x = "Filtering Steps", y = "Total ASV Count") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_value)) +
    scale_fill_manual(values=c("#444c5c", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#78a5a3")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  ggsave(filename = file.path(figure_dir, "total_asv_count.pdf"), plot = plot, width = 15, height = 10, units = "cm")
}

# 3. Unique ASV Count per Organ (Count of unique ASVs in each organ at each step)
plot_unique_asv_per_organ <- function(EvoTraceR_object, figure_dir) {
  figure_dir <- prepare_directory(figure_dir)
  seqtab_history = EvoTraceR_object$seqtab_history
  sample_columns = setdiff(colnames(seqtab_history[[1]]), c("seq_names", "seq", "totalCounts"))

  unique_asv_organ = do.call(rbind, lapply(names(seqtab_history), function(step_name) {
    df = seqtab_history[[step_name]]
    data.frame(
      name = step_name,
      organ = sample_columns,
      unique_count = colSums(df[, sample_columns] > 0)  # Count non-zero entries as unique ASVs
    )
  }))
  
  unique_asv_organ$name <- factor(unique_asv_organ$name, levels = c("Starting ASVs", "Hamming Merging", "Flanking Seq. Filter", "Substitutions Merging", "Frequency Filter", "Final ASVs"))
  
  max_value <- max(unique_asv_organ$unique_count) * 1.3  # Add 10% buffer to max value
  
  plot <- ggplot(unique_asv_organ, aes(x=name, y=unique_count, fill=name)) +
    geom_bar(stat="identity", width=0.8) +
    geom_text(aes(label = unique_count), vjust = -0.5, size = 3) +
    facet_wrap(~organ, scales = "free_y") +
    labs(title = "Unique ASV Count per Organ", x = "Filtering Steps", y = "Unique ASV Count") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_value)) +
    scale_fill_manual(values=c("#444c5c", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#78a5a3")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  ggsave(filename = file.path(figure_dir, "unique_asv_count_per_organ.pdf"), plot = plot, width = 15, height = 10, units = "cm")
}

# 4. Total ASV Count per Organ (Sum of ASV counts in each organ at each step)
plot_total_asv_per_organ <- function(EvoTraceR_object, figure_dir) {
  figure_dir <- prepare_directory(figure_dir)
  seqtab_history = EvoTraceR_object$seqtab_history
  sample_columns = setdiff(colnames(seqtab_history[[1]]), c("seq_names", "seq", "totalCounts"))

  total_asv_organ = do.call(rbind, lapply(names(seqtab_history), function(step_name) {
    df = seqtab_history[[step_name]]
    organ_totals = colSums(df[, sample_columns])
    data.frame(
      name = step_name,
      organ = names(organ_totals),
      total_count = organ_totals
    )
  }))
  
  total_asv_organ$name <- factor(total_asv_organ$name, levels = c("Starting ASVs", "Hamming Merging", "Flanking Seq. Filter", "Substitutions Merging", "Frequency Filter", "Final ASVs"))
  
  max_value <- max(total_asv_organ$total_count) * 1.3  # Add 10% buffer to max value
  
  plot <- ggplot(total_asv_organ, aes(x=name, y=total_count, fill=name)) +
    geom_bar(stat="identity", width=0.8) +
    geom_text(aes(label = total_count), vjust = -0.5, size = 3) +
    facet_wrap(~organ, scales = "free_y") +
    labs(title = "Total ASV Count per Organ", x = "Filtering Steps", y = "Total ASV Count") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_value)) +
    scale_fill_manual(values=c("#444c5c", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#78a5a3")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  ggsave(filename = file.path(figure_dir, "total_asv_count_per_organ.pdf"), plot = plot, width = 15, height = 10, units = "cm")
}

# Wrapper function to generate all plots at once
generate_all_asv_plots <- function(EvoTraceR_object, figure_dir = EvoTraceR_object$output_directory) {
  plot_unique_asv_per_animal(EvoTraceR_object, figure_dir)
  plot_total_asv_per_animal(EvoTraceR_object, figure_dir)
  plot_unique_asv_per_organ(EvoTraceR_object, figure_dir)
  plot_total_asv_per_organ(EvoTraceR_object, figure_dir)
}