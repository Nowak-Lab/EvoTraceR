library(ggplot2)
library(gridExtra)

# Ensure directory exists
prepare_directory <- function(figure_dir) {
  if (!dir.exists(figure_dir)) {dir.create(figure_dir)}
  return(figure_dir)
}

# Helper function to create donut plot for sequence classification for a specified filtering step
create_donut_plot <- function(seqtab, step_name, right_flank, known_contaminations) {
  # Classify each ASV sequence
  classification = sapply(seqtab$seq, function(asv_seq) {
    if (!is.null(right_flank) && grepl(paste0(right_flank, "$"), asv_seq)) {
      return("Identified ASVs")
    } else if (!is.null(known_contaminations) && any(sapply(known_contaminations, grepl, asv_seq))) {
      return("Known Contaminations")
    } else {
      return("Unknown Contaminations")
    }
  })
  
  # Summarize counts for each classification
  classification_counts = table(classification)
  classification_df = data.frame(
    category = names(classification_counts),
    count = as.numeric(classification_counts)
  )
  
  # Define colors consistent with other plots
  color_mapping <- c("Identified ASVs" = "#78a5a3", 
                     "Known Contaminations" = "#444c5c", 
                     "Unknown Contaminations" = "#aaaaaa")
  
  # Create donut plot
  plot <- ggplot(classification_df, aes(x = 2, y = count, fill = category)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    geom_text(aes(label = paste0(round(count / sum(count) * 100, 1), "%")), 
              position = position_stack(vjust = 0.5)) +
    labs(title = paste("Sequence Classification -", step_name), fill = "Category") +
    scale_fill_manual(values = color_mapping) +
    theme_void() +
    theme(legend.position = "right")
  return(invisible(plot))  # Prevent auto-printing of the plot
}

# Helper function to create combined bar and donut plots
create_combined_plot <- function(bar_plot, donut_starting, donut_hamming) {
  library(grid)
  library(ggplot2)
  library(gridExtra)  
  # Arrange the plots without drawing them
  combined_grob <- arrangeGrob(
    bar_plot,
    arrangeGrob(donut_starting, donut_hamming, ncol = 2),
    nrow = 2,
    heights = c(5, 1)
  )
  grid.draw(combined_grob)
  invisible(NULL)  # Ensure the function doesn't return a value that could be auto-printed
}

# 1. Unique ASV Count per Animal with optional donut plots for contamination
plot_unique_asv_per_animal <- function(EvoTraceR_object, figure_dir, right_flank = NULL, known_contaminations = NULL) {
  figure_dir <- prepare_directory(figure_dir)
  seqtab_history = EvoTraceR_object$seqtab_history
  
  unique_asv_animal = data.frame(
    name = factor(names(seqtab_history), levels = c("Starting ASVs", "Hamming Merging", "Flanking Seq. Filter", "Substitutions Merging", "Frequency Filter", "Final ASVs")),
    num = sapply(seqtab_history, nrow)
  )
  max_value <- max(unique_asv_animal$num) * 1.3
  bar_plot <- ggplot(unique_asv_animal, aes(x=name, y=num, fill=name)) +
    geom_bar(stat="identity", width=0.8) +
    geom_text(aes(label = num), vjust = -0.5, size = 3) +
    labs(title = "Unique ASV Count", x = "Filtering Steps", y = "Unique ASV Count") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_value)) +
    scale_fill_manual(values=c("#444c5c", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#78a5a3")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  # Check if contamination plots are needed
  if (!is.null(right_flank) || !is.null(known_contaminations)) {
    donut_starting = create_donut_plot(seqtab_history[["Starting ASVs"]], "Starting ASVs", right_flank, known_contaminations)
    donut_hamming = create_donut_plot(seqtab_history[["Hamming Merging"]], "Hamming Merging", right_flank, known_contaminations)
    
    pdf(file = file.path(figure_dir, "unique_asv_count_with_donuts.pdf"), width = 10, height = 12)
    create_combined_plot(bar_plot, donut_starting, donut_hamming)
    dev.off()
  } else {
    ggsave(filename = file.path(figure_dir, "unique_asv_count.pdf"), plot = bar_plot, width = 15, height = 10, units = "cm")
  }
}

# 2. Total ASV Count per Animal
plot_total_asv_per_animal <- function(EvoTraceR_object, figure_dir, right_flank = NULL, known_contaminations = NULL) {
  figure_dir <- prepare_directory(figure_dir)
  seqtab_history = EvoTraceR_object$seqtab_history
  sample_columns = setdiff(colnames(seqtab_history[[1]]), c("seq_names", "seq", "totalCounts"))

  total_asv_animal = data.frame(
    name = factor(names(seqtab_history), levels = c("Starting ASVs", "Hamming Merging", "Flanking Seq. Filter", "Substitutions Merging", "Frequency Filter", "Final ASVs")),
    num = sapply(seqtab_history, function(df) sum(rowSums(df[, sample_columns])))
  )
  
  max_value <- max(total_asv_animal$num) * 1.3
  
  bar_plot <- ggplot(total_asv_animal, aes(x=name, y=num, fill=name)) +
    geom_bar(stat="identity", width=0.8) +
    geom_text(aes(label = num), vjust = -0.5, size = 3) +
    labs(title = "Total ASV Count", x = "Filtering Steps", y = "Total ASV Count") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_value)) +
    scale_fill_manual(values=c("#444c5c", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#78a5a3")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  if (!is.null(right_flank) || !is.null(known_contaminations)) {
    donut_starting = create_donut_plot(seqtab_history[["Starting ASVs"]], "Starting ASVs", right_flank, known_contaminations)
    donut_hamming = create_donut_plot(seqtab_history[["Hamming Merging"]], "Hamming Merging", right_flank, known_contaminations)
    
    pdf(file = file.path(figure_dir, "total_asv_count_with_donuts.pdf"), width = 10, height = 12)
    create_combined_plot(bar_plot, donut_starting, donut_hamming)
    dev.off()
  } else {
    ggsave(filename = file.path(figure_dir, "total_asv_count.pdf"), plot = bar_plot, width = 15, height = 10, units = "cm")
  }
}

# 3. Unique ASV Count per Organ with optional donut plots for contamination
plot_unique_asv_per_organ <- function(EvoTraceR_object, figure_dir, right_flank = NULL, known_contaminations = NULL) {
  figure_dir <- prepare_directory(figure_dir)
  seqtab_history = EvoTraceR_object$seqtab_history
  sample_columns = setdiff(colnames(seqtab_history[[1]]), c("seq_names", "seq", "totalCounts"))

  unique_asv_organ = do.call(rbind, lapply(names(seqtab_history), function(step_name) {
    df = seqtab_history[[step_name]]
    data.frame(
      name = step_name,
      organ = sample_columns,
      unique_count = colSums(df[, sample_columns] > 0)
    )
  }))
  
  unique_asv_organ$name <- factor(unique_asv_organ$name, levels = c("Starting ASVs", "Hamming Merging", "Flanking Seq. Filter", "Substitutions Merging", "Frequency Filter", "Final ASVs"))
  
  max_value <- max(unique_asv_organ$unique_count) * 1.3
  
  bar_plot <- ggplot(unique_asv_organ, aes(x=name, y=unique_count, fill=name)) +
    geom_bar(stat="identity", width=0.8) +
    geom_text(aes(label = unique_count), vjust = -0.5, size = 3) +
    facet_wrap(~organ, scales = "free_y") +
    labs(title = "Unique ASV Count per Organ", x = "Filtering Steps", y = "Unique ASV Count") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_value)) +
    scale_fill_manual(values=c("#444c5c", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#78a5a3")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  if (!is.null(right_flank) || !is.null(known_contaminations)) {
    donut_starting = create_donut_plot(seqtab_history[["Starting ASVs"]], "Starting ASVs", right_flank, known_contaminations)
    donut_hamming = create_donut_plot(seqtab_history[["Hamming Merging"]], "Hamming Merging", right_flank, known_contaminations)
    
    pdf(file = file.path(figure_dir, "unique_asv_count_per_organ_with_donuts.pdf"), width = 10, height = 12)
    create_combined_plot(bar_plot, donut_starting, donut_hamming)
    dev.off()
  } else {
    ggsave(filename = file.path(figure_dir, "unique_asv_count_per_organ.pdf"), plot = bar_plot, width = 15, height = 10, units = "cm")
  }
}

# 4. Total ASV Count per Organ with optional donut plots for contamination
plot_total_asv_per_organ <- function(EvoTraceR_object, figure_dir, right_flank = NULL, known_contaminations = NULL) {
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
  
  max_value <- max(total_asv_organ$total_count) * 1.3  # Add 30% buffer to max value
  
  bar_plot <- ggplot(total_asv_organ, aes(x=name, y=total_count, fill=name)) +
    geom_bar(stat="identity", width=0.8) +
    geom_text(aes(label = total_count), vjust = -0.5, size = 3) +
    facet_wrap(~organ, scales = "free_y") +
    labs(title = "Total ASV Count per Organ", x = "Filtering Steps", y = "Total ASV Count") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_value)) +
    scale_fill_manual(values=c("#444c5c", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#78a5a3")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  if (!is.null(right_flank) || !is.null(known_contaminations)) {
    donut_starting = create_donut_plot(seqtab_history[["Starting ASVs"]], "Starting ASVs", right_flank, known_contaminations)
    donut_hamming = create_donut_plot(seqtab_history[["Hamming Merging"]], "Hamming Merging", right_flank, known_contaminations)
    
    pdf(file = file.path(figure_dir, "total_asv_count_per_organ_with_donuts.pdf"), width = 10, height = 12)
    create_combined_plot(bar_plot, donut_starting, donut_hamming)
    dev.off()
  } else {
    ggsave(filename = file.path(figure_dir, "total_asv_count_per_organ.pdf"), plot = bar_plot, width = 15, height = 10, units = "cm")
  }
}

# Wrapper function to generate all ASV plots at once, including donut plots if specified
generate_all_asv_plots <- function(EvoTraceR_object, figure_dir = EvoTraceR_object$output_directory, right_flank = NULL, known_contaminations = NULL) {
  print("Generating ASV count plots")
  plot_unique_asv_per_animal(EvoTraceR_object, figure_dir, right_flank, known_contaminations)
  plot_total_asv_per_animal(EvoTraceR_object, figure_dir, right_flank, known_contaminations)
  plot_unique_asv_per_organ(EvoTraceR_object, figure_dir, right_flank, known_contaminations)
  plot_total_asv_per_organ(EvoTraceR_object, figure_dir, right_flank, known_contaminations)
}