# Input: "dnastringset_msa.fasta" -> "alignment_tidy_ref_alt_final.csv"   ------------------------------------------------------
# Output 1: "alignment_tidy_ref_alt_final.csv" -> "del_sub_ins_df.csv" ------------------------------------------------------
# Output 2: "del_sub_ins_df.csv" -> "alt_count_bc.pdf   ------------------------------------------------------
count_alterations <- function(REvoBC_object, alignment_tidy, output_dir_files, output_dir_figures) {
  
  # Create tidy alignment df where, for each position and each ASV you store the 
  # reference and the observed nucleotide
  alignment_tidy_ref_alt <- 
    merge(alignment_tidy, filter(alignment_tidy, !stringr::str_detect(name, "ASV")), by="position") %>%
    dplyr::select(1, 2, 5, 3)
  
  # Adjust Column Names
  colnames(alignment_tidy_ref_alt) <- c("position", "asv_names", "ref_asv", "read_asv")
  
  percentages = REvoBC_object$statistics$asv_df_percentages %>%
    ungroup() %>%
    add_row(dplyr::select(REvoBC_object$barcode, !c(seq_start, seq_end, seq)))# %>%
    #filter(!stringr::str_detect(asv_names, c("NMBC")))
  
  # Join with the sequences df to have the frequency of each ASV in each sample
  alignment_tidy_ref_alt <- inner_join(alignment_tidy_ref_alt, 
                                       percentages, 
                                       by="asv_names") %>%
    dplyr::select(c(asv_names, sample, perc_in_sample, position, ref_asv, read_asv))
  
  
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
    group_by(asv_names, sample) %>%
    mutate(position_bc260 = data.table::nafill(position_bc260, "locf")) %>% # fill NA with consecutive numbers
    dplyr::select(asv_names, sample, position, position_bc260, ref_asv, read_asv, alt, perc_in_sample, alt_bin)  
  
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
            file.path(output_dir_files, "ASV_alterations_width.csv"), 
            row.names = FALSE)
  
  REvoBC_object$alignment$ASV_alterations_width = alignment_tidy_ref_alt_mrg_final_width_summ
  
  # Select only deletions and substitutions
  del_sub_df <- 
    alignment_tidy_ref_alt_mrg_final %>%
    filter(!alt %in% c("ins", "ins_smwr")) %>% 
    group_by(asv_names, sample) %>% 
    dplyr::select(asv_names, sample, position, position_bc260, alt, perc_in_sample)
  
  # Leave only insertion with ASV
  ins_df <- 
    alignment_tidy_ref_alt_mrg_final %>%
    filter(alt == "ins") %>%
    group_by(asv_names, sample, perc_in_sample) %>% 
    mutate(cons_bin = c(0, abs(diff(position_bc260)) == 1)) %>% # find if number is consecutive = 0, if not = 1
    filter(cons_bin == 0) %>%
    dplyr::select(asv_names, sample, position, position_bc260, alt, perc_in_sample)
  # Add Alteration Type
  ins_df$alt <- "ins"
  # Plotting All Alterations with Insertions as one coordinate (i.e. 17 means insertion is here but no length info)
  del_sub_ins_df <- 
    rbind(del_sub_df, ins_df) %>%
    dplyr::select(asv_names, sample, position, position_bc260, alt, perc_in_sample)
  
  sample_columns = setdiff(colnames(REvoBC_object$dada2_asv_prefilter), c("seq_names", "seq"))
  # prepare levels and orders
  del_sub_ins_df <- 
    del_sub_ins_df %>%
    dplyr::mutate(sample = forcats::fct_relevel(sample, sample_columns)) %>%
    arrange(match(sample, sample_columns))
  
  # Save File
  utils::write.csv(del_sub_ins_df, 
            file.path(output_dir_files, "mutations_df.csv"), 
            row.names = FALSE)
  
  
  ################
  ### Plot CNA Frequency based on "del_sub_ins_df" ------------------------------------------------------
  # summarize stat for "del" and "sub" -> position is not stacked but added as one; i.e. pos 10 & freq: 12%, 10%, will be pos: 10 freq: 22%
  del_sub_df_data_to_plot_sum_perc <-
    del_sub_ins_df %>%
    group_by(sample, alt, position_bc260) %>% #
    dplyr::summarise(sum_perc = sum(perc_in_sample)) %>%
    dplyr::filter(!alt == "wt") %>% # don't plot "wt"
    dplyr::filter(!alt == "ins") # don't plot "ins"
  # "ins" only -> position is not stacked but added as one; i.e. pos 10 & freq: 12%, 10%, will be pos: 10 freq: 22%
  ins_df_data_to_plot_sum_perc <-
    del_sub_ins_df %>%
    dplyr::select(asv_names, sample, position_bc260, alt, perc_in_sample) %>%
    dplyr::filter(alt == "ins") %>% # get only "ins"
    unique() %>%
    group_by(sample, alt, position_bc260) %>%
    dplyr::summarise(sum_perc = sum(perc_in_sample))
  # bind "ins" after recalculation with "del_sub"
  del_sub_ins_df_data_to_plot_sum_perc <- rbind(del_sub_df_data_to_plot_sum_perc, ins_df_data_to_plot_sum_perc)

  ### Change 10-12-21 ###
  
  # ### Plot mutations Frequency based on "del_sub_ins_df"
  # # sumarise stat for Del, Sub and Ins -> postion is not stacked but added as one; i.e. pos 10 & freq: 12%, 10%, will be pos: 10 freq: 22% 
  # del_sub_ins_df_data_to_plot_sum_perc <-
  #   del_sub_ins_df %>% 
  #   group_by(sample, alt, position_bc260) %>%
  #   dplyr::summarise(sum_perc = sum(perc_in_sample)) %>% # all deletions that happened in the barcode
  #   dplyr::filter(!alt == "wt") # don't plot wt 
  
  # Save File
  utils::write.csv(del_sub_ins_df_data_to_plot_sum_perc, 
            file.path(output_dir_files, "mutations_frequency.csv"), 
            row.names = FALSE)
  
  
  ### "msa Plot" - Auxiliary Plot Decorations ------------------------------------------------------
  # Position of PAM in guides
  pam_pos <- c(17.5, 42.5, 68.5, 94.5, 120.5, 146.5, 171.5, 198.5, 224.5, 251.5)
  # plot
  alt_count_bc <-
    ggplot(data= del_sub_ins_df_data_to_plot_sum_perc, aes(x=position_bc260, y=sum_perc, fill=alt, group=sample)) +
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
    lemon::facet_rep_grid(rows = vars(sample), cols=vars(alt), repeat.tick.labels = TRUE) 
  
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
  ggsave(filename=file.path(output_dir_figures, "gghist_del_sub_ins_perc.pdf"), 
         plot=alt_count_bc, 
         #device=grDevices::cairo_pdf, 
         width=25, 
         height=5*length(sample_columns), 
         units = "cm") 
  write.csv(del_sub_ins_df_data_to_plot_sum_perc,  
            file.path(output_dir_figures, "/gghist_del_sub_ins_data.csv"),
            row.names = FALSE, quote = FALSE)
  
  
  REvoBC_object$alignment$mutations_df = dplyr::select(del_sub_ins_df, -c('position'))
  return(REvoBC_object)
  
}

# Create the binary mutations matrix
binary_mutation_matrix = function(REvoBC_object, output_dir, output_dir_figures) {
  output_dir_figures = file.path(REvoBC_object$output_directory, "msa_figures")
  if (!dir.exists(output_dir_figures)) {dir.create(output_dir_figures)}
  
  output_dir = file.path(REvoBC_object$output_directory, "msa")
  if (!dir.exists(output_dir)) {dir.create(output_dir)}
  
  # Adjust Alignment
  dnastringset_msa = REvoBC_object$alignment$msa_stringset 
  dnastringset_msa <- Biostrings::unmasked(dnastringset_msa)
  
  # Set references to "perfect_match_single" 
  ref <- dnastringset_msa[REvoBC_object$barcode$asv_names]
  # Reference Sequence (Barcode Plus)
  ref <- as.vector(as.matrix(ref))
  
  ### Convert to Binary Alignment File (for Camin-Sokal MP)
  
  # Call Subs, Insertions and Deletions, Based on msa
  arle.list <- list()
  asv_sequences = names(dnastringset_msa)
  for(ai in 1:length(dnastringset_msa)) {
    # another sequence
    seq <- as.vector(as.matrix(dnastringset_msa[ai]))
    mut.df <- identifyMutations(seq, ref)
    if (nrow(mut.df)>0) {
      mut.df$seq.id <- names(dnastringset_msa)[ai]
      arle.list[[names(dnastringset_msa)[ai]]] <- mut.df 
    }
  }
  
  # bind mutations from all sequences
  mut.df.samples <- do.call('rbind', arle.list)
  mut.df.samples$seq.id <- factor(mut.df.samples$seq.id, levels=names(dnastringset_msa))
  # define mutation ID
  mut.df.samples$mut.id <- paste0(mut.df.samples$type, ".", mut.df.samples$start.barcode, "_",mut.df.samples$end.barcode)
  
  
  # for data frame and remove sub in columns ---------------------------------------------------
  # choose for data tree building: c("del", "ins", "sub", "complex") # what will be included for analysis
  mut.df.samples <-  dplyr::filter(mut.df.samples, type %in% c("del", "ins", "sub", "complex"))
  
  # convert to wide format. make sure only 0 and 1 are allowed]
  arle.all.df.wide= (table(mut.df.samples$seq.id, mut.df.samples$mut.id)>0)*1
  arle.all.df.wide[is.na(arle.all.df.wide)] <- 0
  
  # binary variants 
  asv_bin_var <- arle.all.df.wide
  
  # Show heatmap
  #Heteatmap(asv_bin_var, cluster_columns = T)
  p = pheatmap::pheatmap(asv_bin_var, color = c("#042e61", "#e29892"), fontsize = 6)
  ggsave(filename=file.path(output_dir_figures, "mutations_heatmap.pdf"), 
         plot=p, 
         #device=grDevices::cairo_pdf, 
         width=4*ncol(asv_bin_var), 
         height=5*nrow(asv_bin_var), 
         units = 'mm')
  # Save csv
  write.csv(asv_bin_var, file.path(output_dir, "/binary_mutation_matrix.csv"))
  
  
  # assign and organise BC10 variant (for plotting purposes)
  asv_var <- tibble::tibble(mut.df.samples) %>%
    dplyr::add_row(seq.id  = REvoBC_object$barcode$asv_names, type = 'wt', lengths = 0)
  
  # Save csv
  write.csv(asv_var, file.path(output_dir, "/mutations_coordinates.csv"))
  
  REvoBC_object$alignment$mutations_coordinates = asv_var
  REvoBC_object$alignment$binary_mutation_matrix = asv_bin_var
  
  return(REvoBC_object)
  
}

