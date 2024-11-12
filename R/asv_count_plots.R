library(ggplot2)
library(gridExtra)

# Ensure directory exists
prepare_directory <- function(figure_dir) {
  if (!dir.exists(figure_dir)) {dir.create(figure_dir)}
  return(figure_dir)
}

# Enhanced function to create donut plot for ASV sequence classification
create_unique_count_donut <- function(seqtab, step_name, right_flank = NULL, left_flank = NULL, known_contaminations = NULL, flanking_filtering = "right") {
  # Classify each ASV sequence based on flanking and contamination criteria
  classification = sapply(seqtab$seq, function(asv_seq) {
    is_right_match <- !is.null(right_flank) && grepl(right_flank, asv_seq)
    is_left_match <- !is.null(left_flank) && grepl(left_flank, asv_seq)
    
    # Classification logic based on flanking_filtering
    if (flanking_filtering == "right" && is_right_match) {
      return("Identified ASVs")
    } else if (flanking_filtering == "left" && is_left_match) {
      return("Identified ASVs")
    } else if (flanking_filtering == "either" && (is_right_match || is_left_match)) {
      return("Identified ASVs")
    } else if (flanking_filtering == "both" && is_right_match && is_left_match) {
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
  
  # Define consistent colors
  color_mapping <- c("Identified ASVs" = "#4C9F70", 
                     "Known Contaminations" = "#FF6F61", 
                     "Unknown Contaminations" = "#A0C4FF")
  
  # Create the donut plot
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
  
  return(invisible(plot))  # Return the plot without printing automatically
}

# Enhanced function to create a donut plot for total ASV counts
create_total_count_donut_plot <- function(seqtab, step_name, right_flank = NULL, left_flank = NULL, known_contaminations = NULL, flanking_filtering = "right") {
  sample_columns <- setdiff(names(seqtab), c("seq_names", "seq", "totalCounts"))
  
  # Classify each ASV sequence based on flanking and contamination criteria
  classification = sapply(seqtab$seq, function(asv_seq) {
    is_right_match <- !is.null(right_flank) && grepl(right_flank, asv_seq)
    is_left_match <- !is.null(left_flank) && grepl(left_flank, asv_seq)
    
    # Classification logic based on flanking_filtering
    if (flanking_filtering == "right" && is_right_match) {
      return("Identified ASVs")
    } else if (flanking_filtering == "left" && is_left_match) {
      return("Identified ASVs")
    } else if (flanking_filtering == "either" && (is_right_match || is_left_match)) {
      return("Identified ASVs")
    } else if (flanking_filtering == "both" && is_right_match && is_left_match) {
      return("Identified ASVs")
    } else if (!is.null(known_contaminations) && any(sapply(known_contaminations, grepl, asv_seq))) {
      return("Known Contaminations")
    } else {
      return("Unknown Contaminations")
    }
  })
  
  # Calculate total counts per classification
  seqtab$class <- classification
  total_counts <- aggregate(rowSums(seqtab[, sample_columns]), by = list(seqtab$class), FUN = sum)
  colnames(total_counts) <- c("category", "count")
  
  # Define consistent colors
  color_mapping <- c("Identified ASVs" = "#4C9F70", 
                     "Known Contaminations" = "#FF6F61", 
                     "Unknown Contaminations" = "#A0C4FF")
  
  # Create the donut plot for total ASV counts
  plot <- ggplot(total_counts, aes(x = 2, y = count, fill = category)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    geom_text(aes(label = paste0(round(count / sum(count) * 100, 1), "%")), 
              position = position_stack(vjust = 0.5)) +
    labs(title = paste("Total ASV Count Classification -", step_name), fill = "Category") +
    scale_fill_manual(values = color_mapping) +
    theme_void() +
    theme(legend.position = "right")
  
  return(invisible(plot))  # Return the plot without printing automatically
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
plot_unique_asv <- function(
  EvoTraceR_object,
  figure_dir,
  right_flank = NULL,
  left_flank=NULL,
  known_contaminations = NULL,
  flanking_filtering) {

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
    donut_starting = create_unique_count_donut(seqtab_history[["Starting ASVs"]], "Starting ASVs", right_flank, left_flank = left_flank, known_contaminations=known_contaminations, flanking_filtering=flanking_filtering)
    donut_hamming = create_unique_count_donut(seqtab_history[["Hamming Merging"]], "Hamming Merging", right_flank, left_flank = left_flank, known_contaminations=known_contaminations, flanking_filtering=flanking_filtering)
    
    pdf(file = file.path(figure_dir, "unique_asv_count_with_donuts.pdf"), width = 10, height = 12)
    create_combined_plot(bar_plot, donut_starting, donut_hamming)
    dev.off()
  } else {
    ggsave(filename = file.path(figure_dir, "unique_asv_count.pdf"), plot = bar_plot, width = 15, height = 10, units = "cm")
  }
}

# 2. Total ASV Count per Animal
plot_total_asv_per <- function(
  EvoTraceR_object,
  figure_dir,
  right_flank = NULL,
  left_flank=NULL,
  known_contaminations = NULL,
  flanking_filtering) {

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
    donut_starting = create_total_count_donut_plot(seqtab_history[["Starting ASVs"]], "Starting ASVs", right_flank, left_flank = left_flank, known_contaminations=known_contaminations, flanking_filtering=flanking_filtering)
    donut_hamming = create_total_count_donut_plot(seqtab_history[["Hamming Merging"]], "Hamming Merging", right_flank, left_flank = left_flank, known_contaminations=known_contaminations, flanking_filtering=flanking_filtering)
    
    pdf(file = file.path(figure_dir, "total_asv_count_with_donuts.pdf"), width = 10, height = 12)
    create_combined_plot(bar_plot, donut_starting, donut_hamming)
    dev.off()
  } else {
    ggsave(filename = file.path(figure_dir, "total_asv_count.pdf"), plot = bar_plot, width = 15, height = 10, units = "cm")
  }
}

# 3. Unique ASV Count per Organ with optional donut plots for contamination
plot_unique_asv_per_organ <- function(
  EvoTraceR_object,
  figure_dir,
  right_flank = NULL,
  left_flank = NULL,
  known_contaminations = NULL,
  flanking_filtering) {

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
    donut_starting = create_unique_count_donut(seqtab_history[["Starting ASVs"]], "Starting ASVs", right_flank, left_flank = left_flank, known_contaminations=known_contaminations, flanking_filtering=flanking_filtering)
    donut_hamming = create_unique_count_donut(seqtab_history[["Hamming Merging"]], "Hamming Merging", right_flank, left_flank = left_flank, known_contaminations=known_contaminations, flanking_filtering=flanking_filtering)
    
    pdf(file = file.path(figure_dir, "unique_asv_count_per_organ_with_donuts.pdf"), width = 10, height = 12)
    create_combined_plot(bar_plot, donut_starting, donut_hamming)
    dev.off()
  } else {
    ggsave(filename = file.path(figure_dir, "unique_asv_count_per_organ.pdf"), plot = bar_plot, width = 15, height = 10, units = "cm")
  }
}

# 4. Total ASV Count per Organ with optional donut plots for contamination
plot_total_asv_per_organ <- function(
  EvoTraceR_object, 
  figure_dir, 
  right_flank = NULL, 
  left_flank = NULL,
  known_contaminations = NULL, 
  flanking_filtering) {
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
    donut_starting = create_total_count_donut_plot(seqtab_history[["Starting ASVs"]], "Starting ASVs", right_flank, left_flank = left_flank, known_contaminations=known_contaminations, flanking_filtering=flanking_filtering)
    donut_hamming = create_total_count_donut_plot(seqtab_history[["Hamming Merging"]], "Hamming Merging", right_flank, left_flank = left_flank, known_contaminations=known_contaminations, flanking_filtering=flanking_filtering)
    
    pdf(file = file.path(figure_dir, "total_asv_count_per_organ_with_donuts.pdf"), width = 10, height = 12)
    create_combined_plot(bar_plot, donut_starting, donut_hamming)
    dev.off()
  } else {
    ggsave(filename = file.path(figure_dir, "total_asv_count_per_organ.pdf"), plot = bar_plot, width = 15, height = 10, units = "cm")
  }
}

# Wrapper function to generate all ASV plots at once, including donut plots if specified
generate_all_asv_plots <- function(EvoTraceR_object, figure_dir = EvoTraceR_object$output_directory, right_flank = NULL, left_flank=NULL, known_contaminations = NULL, flanking_filtering="right") {
  print("Generating ASV count plots")
  print(left_flank)
  plot_unique_asv(EvoTraceR_object, figure_dir, right_flank,left_flank, known_contaminations, flanking_filtering)
  plot_total_asv_per(EvoTraceR_object, figure_dir, right_flank,left_flank, known_contaminations, flanking_filtering)
  plot_unique_asv_per_organ(EvoTraceR_object, figure_dir, right_flank, left_flank, known_contaminations, flanking_filtering)
  plot_total_asv_per_organ(EvoTraceR_object, figure_dir, right_flank, left_flank, known_contaminations, flanking_filtering)
}