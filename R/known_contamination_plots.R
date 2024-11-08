# =============================================================================
# These plots are used in our known contamination analysis, where a lab 
# may be aware of specific contamination present in their sequencing data 
# and is working to filter it out.
# =============================================================================


# Ensure directory exists
prepare_directory <- function(figure_dir) {
  figure_dir = paste0(figure_dir, '/asv_analysis/')
  if (!dir.exists(figure_dir)) {dir.create(figure_dir)}
  return(figure_dir)
}

# Helper function to calculate contamination percentage based on substring matching
calculate_contamination <- function(df, contaminant_seqs) {
  # Check if any contaminant sequence is present as a substring within each ASV sequence
  is_contaminant <- sapply(df$seq, function(asv_seq) any(sapply(contaminant_seqs, function(contam_seq) grepl(contam_seq, asv_seq))))
  
  # Select only numeric columns to avoid errors with rowSums
  numeric_df <- df[is_contaminant, sapply(df, is.numeric), drop = FALSE]
  
  # Calculate the total contamination count and percentage
  list(
    total_contamination = sum(rowSums(numeric_df) > 0),
    percentage = sum(is_contaminant) / nrow(df) * 100
  )
}

# 1. Unique ASV Count per Animal with Contamination
plot_unique_asv_per_animal_contamination <- function(EvoTraceR_object, figure_dir, contaminant_seqs) {
  figure_dir <- prepare_directory(figure_dir)
  seqtab_history = EvoTraceR_object$seqtab_history
  
  unique_asv_animal = data.frame(
    name = factor(names(seqtab_history), levels = c("Starting ASVs", "Hamming Merging", "Flanking Seq. Filter", "Substitutions Merging", "Frequency Filter", "Final ASVs")),
    num = sapply(seqtab_history, nrow),
    contamination_percentage = sapply(seqtab_history, function(df) calculate_contamination(df, contaminant_seqs)$percentage)
  )
  
  max_value <- max(unique_asv_animal$num) * 1.1
  
  plot <- ggplot(unique_asv_animal, aes(x=name, y=num, fill=name)) +
    geom_bar(stat="identity", width=0.8) +
    geom_text(aes(label = num), vjust = -0.5, size = 3) +
    geom_line(aes(x=name, y=contamination_percentage * max_value / 100, group = 1, color = "Contamination"), size = 1) +
    geom_point(aes(x=name, y=contamination_percentage * max_value / 100, color = "Contamination"), size = 3) +
    labs(title = "Unique ASV Count with Contamination", x = "Filtering Steps", y = "Unique ASV Count") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_value), sec.axis = sec_axis(~./max_value*100, name = "Contamination %")) +
    scale_fill_manual(values=c("#444c5c", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#78a5a3")) +
    scale_color_manual(values = c("Contamination" = "red")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  ggsave(filename = file.path(figure_dir, "unique_asv_count_contamination.pdf"), plot = plot, width = 15, height = 10, units = "cm")
}

# 2. Total ASV Count per Animal with Contamination
plot_total_asv_per_animal_contamination <- function(EvoTraceR_object, figure_dir, contaminant_seqs) {
  figure_dir <- prepare_directory(figure_dir)
  seqtab_history = EvoTraceR_object$seqtab_history
  sample_columns = setdiff(colnames(seqtab_history[[1]]), c("seq_names", "seq", "totalCounts"))

  total_asv_animal = data.frame(
    name = factor(names(seqtab_history), levels = c("Starting ASVs", "Hamming Merging", "Flanking Seq. Filter", "Substitutions Merging", "Frequency Filter", "Final ASVs")),
    num = sapply(seqtab_history, function(df) sum(rowSums(df[, sample_columns]))),
    contamination_percentage = sapply(seqtab_history, function(df) calculate_contamination(df, contaminant_seqs)$percentage)
  )
  
  max_value <- max(total_asv_animal$num) * 1.1  # Add 10% buffer to max value
  
  plot <- ggplot(total_asv_animal, aes(x=name, y=num, fill=name)) +
    geom_bar(stat="identity", width=0.8) +
    geom_text(aes(label = num), vjust = -0.5, size = 3) +
    geom_line(aes(x=name, y=contamination_percentage * max_value / 100, group = 1, color = "Contamination"), size = 1) +
    geom_point(aes(x=name, y=contamination_percentage * max_value / 100, color = "Contamination"), size = 3) +
    labs(title = "Total ASV Count with Contamination", x = "Filtering Steps", y = "Total ASV Count") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_value), sec.axis = sec_axis(~./max_value*100, name = "Contamination %")) +
    scale_fill_manual(values=c("#444c5c", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#78a5a3")) +
    scale_color_manual(values = c("Contamination" = "red")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  ggsave(filename = file.path(figure_dir, "total_asv_count_contamination.pdf"), plot = plot, width = 15, height = 10, units = "cm")
}

# 3. Unique ASV Count per Organ with Contamination
plot_unique_asv_per_organ_contamination <- function(EvoTraceR_object, figure_dir, contaminant_seqs) {
  figure_dir <- prepare_directory(figure_dir)
  seqtab_history = EvoTraceR_object$seqtab_history
  sample_columns = setdiff(colnames(seqtab_history[[1]]), c("seq_names", "seq", "totalCounts"))

  unique_asv_organ = do.call(rbind, lapply(names(seqtab_history), function(step_name) {
    df = seqtab_history[[step_name]]
    contamination = calculate_contamination(df, contaminant_seqs)
    data.frame(
      name = step_name,
      organ = sample_columns,
      unique_count = colSums(df[, sample_columns] > 0),  # Count non-zero entries as unique ASVs
      contamination_percentage = contamination$percentage
    )
  }))
  
  unique_asv_organ$name <- factor(unique_asv_organ$name, levels = c("Starting ASVs", "Hamming Merging", "Flanking Seq. Filter", "Substitutions Merging", "Frequency Filter", "Final ASVs"))
  max_value <- max(unique_asv_organ$unique_count) * 1.1  # Add 10% buffer to max value
  
  plot <- ggplot(unique_asv_organ, aes(x=name, y=unique_count, fill=name)) +
    geom_bar(stat="identity", width=0.8) +
    geom_text(aes(label = unique_count), vjust = -0.5, size = 3) +
    geom_line(aes(x=name, y=contamination_percentage * max_value / 100, group = 1, color = "Contamination"), size = 1) +
    geom_point(aes(x=name, y=contamination_percentage * max_value / 100, color = "Contamination"), size = 3) +
    facet_wrap(~organ, scales = "free_y") +
    labs(title = "Unique ASV Count per Organ with Contamination", x = "Filtering Steps", y = "Unique ASV Count") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_value), sec.axis = sec_axis(~./max_value*100, name = "Contamination %")) +
    scale_fill_manual(values=c("#444c5c", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#78a5a3")) +
    scale_color_manual(values = c("Contamination" = "red")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  ggsave(filename = file.path(figure_dir, "unique_asv_count_per_organ_contamination.pdf"), plot = plot, width = 15, height = 10, units = "cm")
}

# 4. Total ASV Count per Organ with Contamination
plot_total_asv_per_organ_contamination <- function(EvoTraceR_object, figure_dir, contaminant_seqs) {
  figure_dir <- prepare_directory(figure_dir)
  seqtab_history = EvoTraceR_object$seqtab_history
  sample_columns = setdiff(colnames(seqtab_history[[1]]), c("seq_names", "seq", "totalCounts"))

  total_asv_organ = do.call(rbind, lapply(names(seqtab_history), function(step_name) {
    df = seqtab_history[[step_name]]
    contamination = calculate_contamination(df, contaminant_seqs)
    organ_totals = colSums(df[, sample_columns])
    data.frame(
      name = step_name,
      organ = names(organ_totals),
      total_count = organ_totals,
      contamination_percentage = contamination$percentage
    )
  }))
  
  total_asv_organ$name <- factor(total_asv_organ$name, levels = c("Starting ASVs", "Hamming Merging", "Flanking Seq. Filter", "Substitutions Merging", "Frequency Filter", "Final ASVs"))
  max_value <- max(total_asv_organ$total_count) * 1.1  # Add 10% buffer to max value
  
  plot <- ggplot(total_asv_organ, aes(x=name, y=total_count, fill=name)) +
    geom_bar(stat="identity", width=0.8) +
    geom_text(aes(label = total_count), vjust = -0.5, size = 3) +
    geom_line(aes(x=name, y=contamination_percentage * max_value / 100, group = 1, color = "Contamination"), size = 1) +
    geom_point(aes(x=name, y=contamination_percentage * max_value / 100, color = "Contamination"), size = 3) +
    facet_wrap(~organ, scales = "free_y") +
    labs(title = "Total ASV Count per Organ with Contamination", x = "Filtering Steps", y = "Total ASV Count") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_value), sec.axis = sec_axis(~./max_value*100, name = "Contamination %")) +
    scale_fill_manual(values=c("#444c5c", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#78a5a3")) +
    scale_color_manual(values = c("Contamination" = "red")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  ggsave(filename = file.path(figure_dir, "total_asv_count_per_organ_contamination.pdf"), plot = plot, width = 15, height = 10, units = "cm")
}

# Wrapper function to generate all plots with contamination at once
generate_all_asv_plots_with_contamination <- function(EvoTraceR_object, figure_dir = EvoTraceR_object$output_directory, contaminant_seqs) {
  plot_unique_asv_per_animal_contamination(EvoTraceR_object, figure_dir, contaminant_seqs)
  plot_total_asv_per_animal_contamination(EvoTraceR_object, figure_dir, contaminant_seqs)
  plot_unique_asv_per_organ_contamination(EvoTraceR_object, figure_dir, contaminant_seqs)
  plot_total_asv_per_organ_contamination(EvoTraceR_object, figure_dir, contaminant_seqs)
}