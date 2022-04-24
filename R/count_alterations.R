count_alterations <- function(EvoTraceR_object, output_dir_files, output_dir_figures) {
  
  alignment_tidy_ref_alt = EvoTraceR_object$alignment$asv_barcode_alignment
  
  percentages = EvoTraceR_object$statistics$asv_df_percentages %>%
    ungroup() %>% add_row(data.frame(asv_names = EvoTraceR_object$reference$ref_name,
                                     stringsAsFactors = F))
  
  # Join with the sequences df to have the frequency of each ASV in each sample
  alignment_tidy_ref_alt_mrg_final <- inner_join(alignment_tidy_ref_alt, 
                                                 percentages, 
                                                 by="asv_names") %>%
    dplyr::select(c(asv_names, sample, perc_in_sample, position_bc260, ref_asv, read_asv, alt))
  
  
  # For each ASV, compute the total number of alterations for each type
  alignment_tidy_ref_alt_mrg_final_width_summ <-
    alignment_tidy_ref_alt %>%
    ungroup() %>% group_by(asv_names, alt) %>%
    tally(name="width") %>%
    tidyr::pivot_wider(names_from = alt, values_from = width, values_fill = 0) %>%
    rename_at(vars(!matches("asv_names")), ~ paste0("width_total_", .))
  
  utils::write.csv(alignment_tidy_ref_alt_mrg_final_width_summ, 
                   file.path(output_dir_files, "asv_alterations_width.csv"), 
                   row.names = FALSE)
  
  EvoTraceR_object$alignment$ASV_alterations_width = alignment_tidy_ref_alt_mrg_final_width_summ
  
  # Select only deletions and substitutions
  del_sub_df <- 
    alignment_tidy_ref_alt_mrg_final %>%
    filter(!alt %in% c("ins", "ins_smwr")) %>% 
    dplyr::select(asv_names, sample, position_bc260, alt, perc_in_sample)
  
  # Leave only insertion with ASV
  ins_df <- 
    alignment_tidy_ref_alt_mrg_final %>%
    filter(alt == "ins") %>%
    group_by(asv_names, sample, perc_in_sample) %>% 
    mutate(cons_bin = c(0, abs(diff(position_bc260)) == 1)) %>% # find if number is consecutive = 0, if not = 1
    #filter(cons_bin == 0) %>%
    dplyr::select(asv_names, sample, position_bc260, alt, perc_in_sample)
  
  # Plotting All Alterations with Insertions as one coordinate (i.e. 17 means insertion is here but no length info)
  del_sub_ins_df <- 
    rbind(del_sub_df, ins_df) #%>%
  
  sample_columns = setdiff(colnames(EvoTraceR_object$asv_prefilter), c("seq_names", "seq"))
  
  sample_order = EvoTraceR_object$sample_order
  
  # prepare levels and orders
  del_sub_ins_df <- 
    del_sub_ins_df %>%
    dplyr::mutate(sample = forcats::fct_relevel(sample, sample_order)) %>%
    arrange(match(sample, sample_columns))
  
  
  # Plot CNA Frequency based on "del_sub_ins_df"
  # summarize stat for "del" and "sub" -> position is not stacked but added as one; i.e. pos 10 & freq: 12%, 10%, will be pos: 10 freq: 22%
  del_sub_df_data_to_plot_sum_perc <-
    del_sub_ins_df %>% ungroup() %>%
    group_by(sample, alt, position_bc260) %>% #
    dplyr::summarise(sum_perc = sum(perc_in_sample), .groups = 'drop') %>%
    dplyr::filter(alt != "wt" & alt != "ins") # don't plot "wt" and "ins
  
  # "ins" only -> position is not stacked but added as one; i.e. pos 10 & freq: 12%, 10%, will be pos: 10 freq: 22%
  ins_df_data_to_plot_sum_perc <-
    del_sub_ins_df %>%
    dplyr::select(asv_names, sample, position_bc260, alt, perc_in_sample) %>%
    dplyr::filter(alt == "ins") %>% # get only "ins"
    unique() %>%
    group_by(sample, alt, position_bc260) %>%
    dplyr::summarise(sum_perc = sum(perc_in_sample), .groups = 'drop')
  # bind "ins" after recalculation with "del_sub"
  del_sub_ins_df_data_to_plot_sum_perc <- rbind(del_sub_df_data_to_plot_sum_perc, ins_df_data_to_plot_sum_perc)
  
  # add full names for labeling of plots
  del_sub_ins_df_data_to_plot_sum_perc <-
    del_sub_ins_df_data_to_plot_sum_perc %>%
    mutate(alt_long_names = ifelse(alt == "i", "Insertions", 
                                   ifelse(alt == "d", "Deletions", 
                                          ifelse(alt == "s", "Substitutions", "No Edits"))))
  # save file
  utils::write.csv(del_sub_ins_df_data_to_plot_sum_perc, file.path(output_dir_files, "mutations_frequency.csv"), row.names = FALSE)
  
  # Histogram Graph ------------------------------------------------------
  # auxiliary data for plotting histogram graph   
  # position of PAM in guides
  pam_pos <- EvoTraceR_object$reference$ref_cut_sites
  # length of barcode
  bc_len <- nchar(EvoTraceR_object$reference$ref_seq)
  # annotating rectangles for target sites = 26 bp -> 20x bp (target site) + 3x bp (PAM) + 3x bp (spacer)
  
  backup = del_sub_ins_df_data_to_plot_sum_perc
  del_sub_ins_df_data_to_plot_sum_perc = del_sub_ins_df_data_to_plot_sum_perc %>% filter(alt_long_names != 'No Edits')
  alt_count_annot_rect <-
    ggplot(data = del_sub_ins_df_data_to_plot_sum_perc, aes(x=position_bc260, y=sum_perc, fill=alt_long_names, group=sample))
  # add rectangles
  annot_rect <- EvoTraceR_object$reference$ref_border_sites
  for (i in seq(1, length(annot_rect), by = 2)) {
    alt_count_annot_rect = alt_count_annot_rect + annotate("rect", xmin=annot_rect[i], xmax=annot_rect[i+1], ymin=-Inf, max=Inf, fill="black", alpha=0.1) 
  }
  
  # plot rest of histogram
  alt_count_bc <-
    alt_count_annot_rect + # add earlier prepared rectangles to graph
    # geom bar
    geom_bar(stat="identity", width=1.5) + #, size=1
    scale_fill_manual(values=c("Substitutions"="#329932", "Insertions" = "#FF0033", "Deletions" = "#3366FF"), 
                      breaks=c("Deletions", "Insertions", "Substitutions")) +
    scale_x_continuous(labels=scales::comma, 
                       breaks=c(1, seq(ceiling(bc_len/10), bc_len, ceiling(bc_len/10))), 
                       limits=c(-4, bc_len + 5), 
                       expand = c(0.001, 0.001)) +
    scale_y_continuous(labels=function(x) paste0(x, "%"),
                       limits=c(0, plyr::round_any(max(del_sub_ins_df_data_to_plot_sum_perc$sum_perc), 10, f = ceiling)), # automatic
                       breaks=seq(from=0, to=plyr::round_any(max(del_sub_ins_df_data_to_plot_sum_perc$sum_perc), 10, f = ceiling), by=10),
                       #breaks=c(0, 25, 50, 75, 100), # manual 
                       expand = c(0, 0)) +
    geom_vline(xintercept=pam_pos, linetype="dashed", size=0.3, col="orange") + # Cas9 Cleavage
    lemon::coord_capped_cart(left="both", bottom="both") +
    lemon::facet_rep_grid(rows = vars(sample), cols=vars(alt_long_names), repeat.tick.labels = TRUE) +
    labs(y = "Editing Frequency", x = "Position of Nucleotides", fill = "Type of Editing")
  
  # add theme
  alt_count_bc <- 
    alt_count_bc + 
    theme(plot.margin = unit(c(1, 1, 1, 1), "mm"),
          axis.ticks = element_blank(), # disable ticks lines
          axis.line.y = element_line(colour="black", size=0.3), # axis y line only
          axis.line.x = element_line(colour="black", size=0.3), # axis x line only
          panel.border = element_blank(), # disable panel border
          panel.grid.major = element_blank(), # disable lines in grid on X-axis
          panel.grid.minor = element_blank(), # disable lines in grid on X-axis
          axis.text.y = element_text(size=8, angle=0, hjust=1, vjust=0.5),
          axis.text.x = element_text(size=8, angle=0, hjust=0.5, vjust=0.5),
          axis.ticks.x = element_line(colour="black", size=0.3),
          axis.ticks.y = element_line(colour="black", size=0.3),
          legend.position="bottom", legend.box = "horizontal",
          strip.background=element_blank(), 
          strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-7.5,-5,-5,-5),
          panel.background = element_rect(fill="white")) 
  
  # Save PDF
  #ggsave(filename=file.path(output_dir_figures, "hist_del_sub_ins_perc.pdf"), plot=alt_count_bc, width=25, height=5*length(sample_columns), units = "cm")
  # temp
  ggsave(filename=file.path(output_dir_figures, "/hist_del_sub_ins_data.pdf"), 
         plot=alt_count_bc, width=25, height=5*length(sample_columns), units = "cm") 
  # write.csv(del_sub_ins_df_data_to_plot_sum_perc, file.path(output_dir_figures, "/hist_del_sub_ins_data.csv"), row.names = FALSE, quote = FALSE)
  
  EvoTraceR_object$alignment$mutations_df = del_sub_ins_df
  return(EvoTraceR_object)
  
}
