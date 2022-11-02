
seq_filtering_plot = function(EvoTraceR_object, figure_dir = EvoTraceR_object$output_directory) {
  figure_dir = paste0(figure_dir, '/asv_analysis/')
  if (!dir.exists(figure_dir)) {dir.create(figure_dir)}
  track_data = EvoTraceR_object$preprocessing$seq_filters

  track_data <- dplyr::mutate(track_data, name = fct_relevel(name, c("Starting ASVs", "Hamming Merging", "Substitutions Merging", "Frequency Filter", "Flanking Seq. Filter",  "Final ASVs")))
  track_data$num_names <- paste0(track_data$num, " x ASVs") # numbers of ASV included on the top of bar  
  
  # calculate percent change from previous filter
  track_data <-
    track_data %>%
    mutate(track_data, diff_perc = ceiling(x=(num/lag(num)-1)*100)) %>%
    mutate(diff_perc = paste0(diff_perc, " %")) %>%
    mutate(diff_perc = if_else(name == "Starting ASVs" | name == "Final ASVs", "", diff_perc)) %>%
    mutate(num_names_sf = if_else(name == "Starting ASVs" | name == "Final ASVs", num_names, "")) %>%
    mutate(num_names_ins = if_else(name == "Starting ASVs" | name == "Final ASVs", "", num_names))
  
  # start graph plotting 
  seqtab_df_clean_track <-
    ggplot(data=track_data) +
    geom_bar(aes(x=name, y=num, fill=name), position = "dodge", stat = "identity", width=0.8, size=0.2, show.legend = FALSE) +
    geom_text(aes(x=name, y=num, label=num_names_sf), vjust=-0.25, size=3) + # change order to have up whatever you choose, opposite to order
    geom_text(aes(x=name, y=num, label=num_names_ins), vjust=-1.75, size=3) + # change order to have up whatever you choose, opposite to order
    geom_text(aes(x=name, y=num, label=diff_perc), vjust=-0.25, size=3, col="blue") + # change order to have up whatever you choose, opposite to order
    scale_y_continuous(expand = c(0, 0), 
                       limits= c(0, plyr::round_any(max(track_data$num), 1000, f = ceiling)+1000)) +
    scale_fill_manual(values=c("#444c5c", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#78a5a3")) +
    labs(x = "ASVs Filtering Steps", y = element_blank()) + 
    barplot_nowaklab_theme() + # add theme 
    theme(plot.margin = unit(c(0, 0, 0, 0), "mm"), # update theme specifically 
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # , hjust = 1, vjust = 1
          axis.line.x = element_blank(), # disable x axis lines
          axis.line.y=element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) # disable x axis ticks lines
  # save pdf
  ggsave(filename=file.path(figure_dir, "asv_filtering_freq.pdf"), plot=seqtab_df_clean_track, width=15, height=15, units = "cm")
}












