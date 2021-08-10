# Input: "dnastringset_msa.fasta" -> "alignment_tidy_ref_alt_final.csv"   ------------------------------------------------------
# Output 1: "alignment_tidy_ref_alt_final.csv" -> "del_sub_ins_df.csv" ------------------------------------------------------
# Output 2: "del_sub_ins_df.csv" -> "alt_count_bc.pdf   ------------------------------------------------------
count_alterations <- function(REvoBC_object, alignment_tidy, output_dir) {
  
  # Create tidy alignment df where, for each position and each ASV you store the 
  # reference and the observed nucleotide
  alignment_tidy_ref_alt <- 
    merge(alignment_tidy, filter(alignment_tidy, stringr::str_detect(name, ".ORG")), by="position") %>%
    dplyr::select(1, 2, 5, 3)
  
  # Adjust Column Names
  colnames(alignment_tidy_ref_alt) <- c("position", "asv_names", "ref_asv", "read_asv")
  
  percentages = REvoBC_object$statistics$asv_df_percentages %>%
    ungroup() %>%
    add_row(dplyr::select(REvoBC_object$barcode, !c(seq_start, seq_end, seq))) %>%
    filter(!stringr::str_detect(asv_names, c("NMBC")))
  
  # Join with the sequences df to have the frequency of each ASV in each sample
  alignment_tidy_ref_alt <- inner_join(alignment_tidy_ref_alt, 
                                       percentages, 
                                       by="asv_names") %>%
    dplyr::select(c(asv_names, day_organ, perc_in_sample, position, ref_asv, read_asv))
  
  
  # Assign Alterations types; wt - wild type, del - deletions, ins - insertion, sub - substitution
  # In cases where ref_asv == "-" & read_asv == "-", then the alteration type is set to ins_smwr 
  # (insertion sites, it causes separation of the barcode sequence)
  alignment_tidy_ref_alt_mrg <-
    alignment_tidy_ref_alt %>%
    mutate(alt = ifelse(ref_asv == read_asv, "wt",
                        ifelse(read_asv == "-", "del", 
                               ifelse(ref_asv == "-", "ins", "sub")))) %>%
    mutate(alt = ifelse(ref_asv == "-" & read_asv == "-", "ins_smwr", alt)) %>%
    mutate(alt_bin = ifelse(alt == "wt" , 0, 1)) # present/absent mutation
  
  # Now translate barcode to barcode 260 nts long.
  # In case of sites of insertion, the last position not affected by the insertion is propagated
  # for all the nucleotides affected by that insertion.
  alignment_tidy_ref_alt_mrg_final <-
    alignment_tidy_ref_alt_mrg %>%
    mutate(position_bc260 = ifelse(ref_asv == "-", NA, position)) %>% # add "bc260" scale
    group_by(asv_names, day_organ) %>%
    mutate(position_bc260 = data.table::nafill(position_bc260, "locf")) %>% # fill NA with consecutive numbers
    dplyr::select(asv_names, day_organ, position, position_bc260, ref_asv, read_asv, alt, perc_in_sample, alt_bin)  
  
  # recount for consecutive
  alignment_tidy_ref_alt_mrg_final$position_bc260 <- as.numeric(as.character(factor(alignment_tidy_ref_alt_mrg_final$position_bc260, 
                                                                                    levels = unique(alignment_tidy_ref_alt_mrg_final$position_bc260), 
                                                                                    labels = seq_along(unique(alignment_tidy_ref_alt_mrg_final$position_bc260)))))
  # For each ASV, compute the total number of alterations for each type
  alignment_tidy_ref_alt_mrg_final_width_summ <-
    alignment_tidy_ref_alt_mrg_final %>%
    ungroup() %>% 
    dplyr::select(asv_names, position, alt) %>%
    distinct(asv_names, position, .keep_all = TRUE) %>%
    group_by(asv_names, alt) %>%
    tally(name="width") %>%
    tidyr::pivot_wider(names_from = alt, values_from = width, values_fill = 0) %>%
    rename_at(vars(!matches("asv_names")), ~ paste0("width_total_", .))
  
  utils::write.csv(alignment_tidy_ref_alt_mrg_final_width_summ, 
            file.path(output_dir, "ASV_alterationType_frequency.csv"), 
            row.names = FALSE)
  
  REvoBC_object$alignment$asv_alterationType_frequency = alignment_tidy_ref_alt_mrg_final_width_summ
  
  # Select only deletions and substitutions
  del_sub_df <- 
    alignment_tidy_ref_alt_mrg_final %>%
    filter(!alt %in% c("ins", "ins_smwr")) %>% 
    group_by(asv_names, day_organ) %>% 
    dplyr::select(asv_names, day_organ, position, position_bc260, alt, perc_in_sample)
  
  # Leave only insertion with ASV
  ins_df <- 
    alignment_tidy_ref_alt_mrg_final %>%
    filter(alt == "ins") %>%
    group_by(asv_names, day_organ, perc_in_sample) %>% 
    mutate(cons_bin = c(0, abs(diff(position_bc260)) == 1)) %>% # find if number is consecutive = 0, if not = 1
    filter(cons_bin == 0) %>%
    dplyr::select(asv_names, day_organ, position, position_bc260, alt, perc_in_sample)
  # Add Alteration Type
  ins_df$alt <- "ins"
  # Plotting All Alterations with Insertions as one coordinate (i.e. 17 means insertion is here but no length info)
  del_sub_ins_df <- 
    rbind(del_sub_df, ins_df) %>%
    dplyr::select(asv_names, day_organ, position, position_bc260, alt, perc_in_sample)
  
  sample = setdiff(colnames(REvoBC_object$dada2_asv_prefilter), c("seq_names", "seq"))
  # prepare levels and orders
  del_sub_ins_df <- 
    del_sub_ins_df %>%
    dplyr::mutate(day_organ = forcats::fct_relevel(day_organ, sample)) %>%
    arrange(match(day_organ, sample))
  
  # Save File
  utils::write.csv(del_sub_ins_df, 
            file.path(output_dir, "mutations_df.csv"), 
            row.names = FALSE)
  
  ### Plot mutations Frequency based on "del_sub_ins_df"
  # sumarise stat for Del, Sub and Ins -> postion is not stacked but added as one; i.e. pos 10 & freq: 12%, 10%, will be pos: 10 freq: 22% 
  del_sub_ins_df_data_to_plot_sum_perc <-
    del_sub_ins_df %>% 
    group_by(day_organ, alt, position_bc260) %>%
    dplyr::summarise(sum_perc = sum(perc_in_sample)) %>% # all deletions that happened in the barcode
    dplyr::filter(!alt == "wt") # don't plot wt 
  # Save File
  utils::write.csv(del_sub_ins_df_data_to_plot_sum_perc, 
            file.path(output_dir, "mutations_frequency.csv"), 
            row.names = FALSE)
  
  
  ### "msa Plot" - Auxiliary Plot Decorations ------------------------------------------------------
  # Position of PAM in guides
  pam_pos <- c(17.5, 42.5, 68.5, 94.5, 120.5, 146.5, 171.5, 198.5, 224.5, 251.5)
  # plot
  alt_count_bc <-
    ggplot(data= del_sub_ins_df_data_to_plot_sum_perc, aes(x=position_bc260, y=sum_perc, fill=alt, group=day_organ)) +
    annotate("rect", xmin=1, xmax=26, ymin=-Inf, max=Inf, fill="black", alpha=.1) +
    annotate("rect", xmin=52, xmax=78, ymin=-Inf, max=Inf, fill="black", alpha=.1) +
    annotate("rect", xmin=104, xmax=130, ymin=-Inf, max=Inf, fill="black", alpha=.1) +
    annotate("rect", xmin=156, xmax=182, ymin=-Inf, max=Inf, fill="black", alpha=.1) +
    annotate("rect", xmin=208, xmax=234, ymin=-Inf, max=Inf, fill="black", alpha=.1) +
    # geom bar
    geom_bar(stat="identity", width=1.5, size=1) +
    scale_fill_manual(values=c("sub"="#329932", "ins" = "#FF0033", "ins_smwr" = "pink", "del" = "#3366FF", "wt" = "#f2f2f2"), breaks=c("wt", "del", "sub", "ins", "ins_smwr")) +
    scale_x_continuous(labels=scales::comma, breaks=c(1, seq(26, 260, 26)), limits=c(-4, 265), expand = c(0.001, 0.001)) +
    scale_y_continuous(labels=function(x) paste0(x, "%"),
                       limits=c(0, 3+max(del_sub_ins_df_data_to_plot_sum_perc$sum_perc)), 
                       expand = c(0, 0)) +
    geom_vline(xintercept=pam_pos, linetype="dashed", size=0.3, col="orange") + # Cas9 Cleavage
    lemon::coord_capped_cart(left="both", bottom="both") +
    lemon::facet_rep_grid(rows = vars(day_organ), cols=vars(alt), repeat.tick.labels = TRUE) 
  
  # Add Theme
  alt_count_bc <- 
    alt_count_bc + 
    theme_bw() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
          axis.ticks = element_blank(), # disable ticks lines
          axis.line.y = element_line(colour="black", size=0.3), # axis y line only
          axis.line.x = element_line(colour="black", size=0.3), # axis x line only
          axis.title = element_blank(), # disable panel border,
          panel.border = element_blank(), # disable panel border
          #panel.grid.major.y = element_blank(), # y grid line 
          panel.grid.major = element_blank(), # disable lines in grid on X-axis
          #panel.grid.minor.x = element_blank(), # disable lines in grid on X-axis
          panel.grid.minor = element_blank(), # disable lines in grid on X-axis
          axis.text.y = element_text(size=6, angle=0, hjust=1, vjust=0.5),
          axis.text.x = element_text(size=6, angle=0, hjust=0.5, vjust=0.5),
          axis.ticks.x = element_line(colour="black", size=0.3),
          axis.ticks.y = element_line(colour="black", size=0.3),
          legend.position="bottom", legend.box = "horizontal",
          strip.background=element_blank(),
          panel.background = element_rect(fill="white"))
  
  # Save PDF
  ggsave(filename=file.path(output_dir, "gghist_del_sub_ins_perc.pdf"), 
         plot=alt_count_bc, 
         #device=grDevices::cairo_pdf, 
         width=25, 
         height=5*length(sample), 
         units = "cm") #17.5 for 4x
  
  REvoBC_object$alignment$mutations_df = dplyr::select(del_sub_ins_df, -c('position'))
  return(REvoBC_object)
  
}
