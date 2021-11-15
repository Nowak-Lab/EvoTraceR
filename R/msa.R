# Input: "dnastringset_msa.fasta" -> "alignment_tidy_ref_alt_final.csv"   ------------------------------------------------------
# Output 1: "alignment_tidy_ref_alt_final.csv" -> "del_sub_ins_df.csv" ------------------------------------------------------
# Output 2: "del_sub_ins_df.csv" -> "alt_count_bc.pdf   ------------------------------------------------------
count_alterations <- function(REvoBC_object, output_dir_files, output_dir_figures) {
  
  alignment_tidy_ref_alt = REvoBC_object$alignment$asv_barcode_alignment
  
  percentages = REvoBC_object$statistics$asv_df_percentages %>%
    ungroup() %>% add_row(dplyr::select(REvoBC_object$barcode, !c(seq_start, seq_end, seq)))
  
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
            file.path(output_dir_files, "ASV_alterations_width.csv"), 
            row.names = FALSE)
  
  REvoBC_object$alignment$ASV_alterations_width = alignment_tidy_ref_alt_mrg_final_width_summ
  
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
    
  
  sample_columns = setdiff(colnames(REvoBC_object$dada2_asv_prefilter), c("seq_names", "seq"))
  # prepare levels and orders
  del_sub_ins_df <- 
    del_sub_ins_df %>%
    dplyr::mutate(sample = forcats::fct_relevel(sample, sample_columns)) %>%
    arrange(match(sample, sample_columns))
  
  
  ################
  ### Plot CNA Frequency based on "del_sub_ins_df" ------------------------------------------------------
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

  # Save File
  utils::write.csv(del_sub_ins_df_data_to_plot_sum_perc, 
            file.path(output_dir_files, "mutations_frequency.csv"), 
            row.names = FALSE)
  
  
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
  
  
  REvoBC_object$alignment$mutations_df = del_sub_ins_df
  return(REvoBC_object)
  
}

align_asv = function(REvoBC_object, output_dir_files, output_dir_figures, ...) {
  dots = list(...)
  df_to_plot_org_tree <- 
    dplyr::select(REvoBC_object$statistics$asv_toBarcode_similarity, asv_names, seq) %>%
    filter(!str_detect(asv_names, "NMBC")) %>%
    distinct() %>%
    arrange(asv_names)
  
  ### Create Biostring and FASTA Files Data
  
  dnastringset <- Biostrings::DNAStringSet(df_to_plot_org_tree$seq) 
  names(dnastringset) <- df_to_plot_org_tree$asv_names
  
  # Output as FASTA files
  outputFileFASTA <- file.path(output_dir_files, "dnastringset.fa")
  Biostrings::writeXStringSet(dnastringset, outputFileFASTA)
  
  utils::write.csv(dnastringset, file.path(output_dir_files, "dnastringset.csv"), quote=FALSE)
  
  # Perform Multiple Sequence Alignment with muscle
  dnastringset_msa <- do.call(muscle::muscle, c(list(stringset = dnastringset),
                                                get_args_from_dots(dots, muscle::muscle))) # suggestion: gapopen = -400
  Biostrings::writeXStringSet(methods::as(dnastringset_msa, "DNAStringSet"), 
                              file.path(output_dir_files, "dnastringset_muscle-muscle_msa.fasta"))
  
  # Store MSA result
  msa = Biostrings::DNAMultipleAlignment(dnastringset_msa)
  REvoBC_object$alignment$msa_stringset = msa
  
  # Transform Alignment to "Tidy Data"
  alignment_tidy <- ggmsa::tidy_msa(msa=msa, 
                                    start = 1, 
                                    end = ncol(msa))
  
  # Create tidy alignment df where, for each position and each ASV you store the 
  # reference and the observed nucleotide
  alignment_tidy_ref_alt <- 
    merge(alignment_tidy, filter(alignment_tidy, !stringr::str_detect(name, "ASV")), by="position") %>%
    dplyr::select(1, 2, 5, 3)
  
  # Adjust Column Names
  colnames(alignment_tidy_ref_alt) <- c("position", "asv_names", "ref_asv", "read_asv")
  
  # Assign Alterations types; wt - wild type, del - deletions, ins - insertion, sub - substitution
  # In cases where ref_asv == "-" & read_asv == "-", then the alteration type is set to ins_smwr 
  alignment_tidy_ref_alt <-
    alignment_tidy_ref_alt %>%
    mutate(alt = ifelse(ref_asv == read_asv, "wt",
                        ifelse(read_asv == "-", "del", 
                               ifelse(ref_asv == "-", "ins", "sub")))) %>%
    mutate(alt = ifelse(ref_asv == "-" & read_asv == "-", "ins_smwr", alt)) %>%
    mutate(alt_bin = ifelse(alt == "wt" , 0, 1)) # present/absent mutation
  
  # Now translate barcode to barcode 260 nts long.
  # In case of sites of insertion, the last position not affected by the insertion is propagated
  # for all the nucleotides affected by that insertion.
  alignment_tidy_ref_alt <-
    alignment_tidy_ref_alt %>%
    mutate(position_bc260 = ifelse(ref_asv == "-", NA, position)) %>% # add "bc260" scale
    group_by(asv_names) %>%
    mutate(position_bc260 = data.table::nafill(position_bc260, "locf")) %>% # fill NA with consecutive numbers
    dplyr::select(asv_names,  position, position_bc260, ref_asv, read_asv, alt, alt_bin)#,  sample, perc_in_sample,)  
  
  # recount for consecutive
  alignment_tidy_ref_alt$position_bc260 <- as.numeric(as.character(factor(alignment_tidy_ref_alt$position_bc260, 
                                                                          levels = unique(alignment_tidy_ref_alt$position_bc260), 
                                                                          labels = seq_along(unique(alignment_tidy_ref_alt$position_bc260)))))
  
  alignment_tidy_ref_alt = dplyr::filter(alignment_tidy_ref_alt, alt != 'ins_smwr')
  
  REvoBC_object$alignment$asv_barcode_alignment = alignment_tidy_ref_alt %>% 
    dplyr::select(-c('alt_bin', 'position'))
  
  utils::write.csv(REvoBC_object$alignment$asv_barcode_alignment, 
                   file.path(output_dir_files, "asv_barcode_alignment.csv"), 
                   quote=F,
                   row.names = FALSE)
  return(REvoBC_object)
  
}
# Create the binary mutations matrix
binary_mutation_matrix = function(REvoBC_object, output_dir_files, output_dir_figures) {
  
  
  # Start from the dataframe that contains a line for each position in each ASV.
  # This dataframe contains also cases where we have the so-called 'insertion somewhere' events
  # which are cases where in the reference and mutated sequences we find a dash '-'.
  # There is a dash on the reference for those position where in some of the ASV an insertion happened.
  # When we count for mutations, we don't care about them
  # Before removing any mutation we create the column contig alt
  # This code creates (with head and tail) two vectors where the one (head) contains all 
  # alterations except for the last one, and the second (tail) refers to all mutations
  # except for the first one. In this way, by comparing head and tail you know which 
  # positions are different from the one right below them. 
  # By computing cumsum of the different position you get a vector of consecutive numbers
  # where the quantity is incremented only in those positions that contain a different value
  # from the one next to them. For example, if the column "alt" is the following:
  # ins, ins, ins, sub, del, del,del, then the new column contig_alt would be 1,1,1,2,3,3,3.
  # We need this so we can group consecutive deletions and inertions and consider them
  # as one whole event.
  
  mut_df = REvoBC_object$alignment$asv_barcode_alignment %>% ungroup() %>% 
    arrange(asv_names, position_bc260) %>%
    group_by(asv_names) %>% 
    mutate(contig_alt = cumsum(c(1, head(alt, -1) != tail(alt, -1)))) %>%
    filter(alt!= 'wt' & alt != 'ins_smwr') %>%
    dplyr::rename(mutation_type = alt, ref_seq = ref_asv, alt_seq = read_asv) 
  
  # Temporary remove substitutions, prepare them in the same format as we will have the other mutations
  sub_df = mut_df %>% filter(mutation_type == 'sub') %>% 
    mutate(start = position_bc260, end = position_bc260, n_nucleotides = 1) %>% 
    dplyr::select(asv_names, mutation_type, start, end, n_nucleotides, ref_seq, alt_seq)
  
  
  # For insertions, the position refers to the nucleotide on the left on the insertion
  # For example, if the insertion is assigned to position 20, this means that after 20 nts 
  # of the original barcode, a new sequence is inserted (so, it starts from the position 21).

  mut_df = dplyr::filter(mut_df, mutation_type!= 'sub') 

  # Now we group alterations first based on the ASV, and then based on the new column contig_alt
  # in this way we can count the rows in the group to obtain the length of the insertion
  # We finally insert back the substitutions.
  mut_df = mut_df %>% group_by(asv_names, mutation_type, contig_alt) %>%
    dplyr::summarise(start = min(position_bc260), 
                     end = max(position_bc260), 
                     n_nucleotides = n(), 
                     ref_seq = paste(ref_seq, collapse = ''),
                     alt_seq = paste(alt_seq, collapse = ''),
                     .groups='drop') %>%
    dplyr::select(asv_names, mutation_type, start, end, n_nucleotides, ref_seq, alt_seq) %>% bind_rows(sub_df)
  
  mut_df$mut_id = paste0(mut_df$mutation_type, "_", mut_df$start, "_", mut_df$n_nucleotides, "nts")
  
  
  # Now smooth deletions and insertions (assign start and end to the closest cutting site)
  smoothed_del_ins = smooth_deletions(mut_df %>% filter(mutation_type != 'sub'))
  
  # After smoothing, we insert back substitutions, but we need to remove those that
  # fall within a smoothed deletion
  sub_df = mut_df %>% ungroup %>% filter(mutation_type == 'sub')

  for (asv in unique(as.character(sub_df$asv_names))) {
    tmp_sub = sub_df %>% filter(asv_names == asv)
    sub_df = sub_df %>% filter(asv_names != asv)
    tmp_smoothed = smoothed_del_ins  %>% filter(asv_names == asv)
    tmp_sub = tmp_sub[!sapply(tmp_sub$start, function(x) 
      return(any(tmp_smoothed$start_smoothed <= x & tmp_smoothed$end_smoothed >=x))),]
    
    sub_df = dplyr::bind_rows(sub_df, tmp_sub)
  }

  smoothed_mutations = sub_df %>%
    mutate(start_smoothed = start, end_smoothed = end, smoothed_id = mut_id) %>%
    bind_rows(smoothed_del_ins)
  
  # Now that we have mutations smoothed, compute the tidy alignment tibble.
  alignment_smoothed_tidy = tidy_alignment_smoothed(REvoBC_object, smoothed_mutations)
  
  REvoBC_object$smoothed_deletions_insertions$asv_barcode_alignment = alignment_smoothed_tidy
  # convert to binary mutation matrix. Compute the count of each mutation of each ASV.
  # and then convert the tibble from long (one row is ASV mutation count) to wide, in order 
  # to have a column for each mutation
  mut_df_wide = plyr::count(mut_df, vars = c("asv_names", "mut_id")) %>% 
    pivot_wider(names_from = mut_id, values_from = freq) %>% column_to_rownames("asv_names")
  
  mut_df_wide[is.na(mut_df_wide)] <- 0
  mut_df_wide[mut_df_wide >= 1] = 1
  
  # Binarize smoothed deletions
  smoothed_df_wide = plyr::count(smoothed_mutations, vars = c("asv_names", "smoothed_id")) %>% 
    pivot_wider(names_from = smoothed_id, values_from = freq) %>% column_to_rownames("asv_names")
  
  smoothed_df_wide[is.na(smoothed_df_wide)] <- 0
  smoothed_df_wide[smoothed_df_wide >= 1] = 1
  
  # Show heatmap
  
  p = pheatmap::pheatmap(mut_df_wide, color = c("#042e61", "#e29892"), fontsize = 6)
  ggsave(filename=file.path(output_dir_figures, "mutations_heatmap.pdf"), 
         plot=p, 
         #device=grDevices::cairo_pdf, 
         width=4*ncol(mut_df_wide), 
         height=5*nrow(mut_df_wide), 
         units = 'mm')
  
  p = pheatmap::pheatmap(smoothed_df_wide, color = c("#042e61", "#e29892"), fontsize = 6)
  ggsave(filename=file.path(output_dir_figures, "smoothed_deletions_insertions_heatmap.pdf"), 
         plot=p, 
         #device=grDevices::cairo_pdf, 
         width=4*ncol(smoothed_df_wide), 
         height=5*nrow(smoothed_df_wide), 
         units = 'mm')

  # assign BC10 variant (for plotting purposes)
  
  mut_df <- tibble::tibble(mut_df) %>%
    dplyr::add_row(asv_names  = REvoBC_object$barcode$asv_names, mutation_type = 'wt', n_nucleotides = 0)
  
  # In case some ASVs are not affected by any mutation, add a row to the binary mutation matrix 
  # which contains all zeros
  #wt_asv = c(setdiff(mut_df$asv_names, smoothed_del$asv_names))
  smoothed_mutations <- tibble::tibble(smoothed_mutations) %>%
    dplyr::add_row(asv_names  = REvoBC_object$barcode$asv_names, mutation_type = 'wt', n_nucleotides = 0)
  
  # Add row to binary mutation matrix corresponding to the original barcode (i.e. all mutations = 0)
  bc_mut = as.list(rep(0, ncol(mut_df_wide)))
  names(bc_mut) = colnames(mut_df_wide)
  mut_df_wide = dplyr::bind_rows(data.frame(bc_mut, row.names = REvoBC_object$barcode$asv_names),
                                 mut_df_wide)
  
  # Add row to binary mutation matrix corresponding to the original barcode (i.e. all mutations = 0)
  bc_mut = as.list(rep(0, ncol(smoothed_df_wide)))
  names(bc_mut) = colnames(smoothed_df_wide)

  smoothed_df_wide = dplyr::bind_rows(data.frame(bc_mut, row.names = REvoBC_object$barcode$asv_names),
                                      data.frame(smoothed_df_wide))
  
  write.csv(mut_df_wide, file.path(output_dir_files, "/binary_mutation_matrix.csv"))
  utils::write.csv(mut_df, 
                   file.path(output_dir_files, "mutations_coordinates.csv"), 
                   row.names = FALSE)
  
  REvoBC_object$alignment$binary_mutation_matrix = mut_df_wide
  REvoBC_object$alignment$mutations_coordinates = mut_df
  
  REvoBC_object$smoothed_deletions_insertions$binary_matrix = smoothed_df_wide
  REvoBC_object$smoothed_deletions_insertions$coordinate_matrix = smoothed_mutations
  

  
  return(REvoBC_object)
}

# mut_df is a dataframe with start and end of insertions and deletions in each ASV (no substitutions)
smooth_deletions = function(mut_df) {
  
  deletions_insertions = mut_df #%>% filter(mutation_type != 'ins')
  orange_lines = data.frame(site = c(17, 42, 68, 94, 120, 146, 171, 198, 224, 251),
                            index_cut = c(1:10))
  data.table::setDT(orange_lines)
  data.table::setDT(deletions_insertions)    
  
  # Create time column by which to do a rolling join
  orange_lines[, start_site := site]
  deletions_insertions[, start_site := start]
  
  deletions_insertions = orange_lines[deletions_insertions, on = "start_site", roll = "nearest"] %>%
    mutate(index_cut = ifelse(start - site > 5, index_cut + 1, index_cut)) %>%
    mutate(start_smoothed = orange_lines$site[index_cut]) %>%
    #dplyr::mutate(start_smoothed = ifelse(start - site > 5, orange_lines$site[index_cut +1], site)) %>%
    dplyr::select(-c(site, index_cut, start_site))
  
  orange_lines = dplyr::rename(orange_lines, end_site = start_site)
  deletions_insertions[, end_site := end]
  deletions_insertions = orange_lines[deletions_insertions, on = "end_site", roll = "nearest"] %>%
    mutate(index_cut = ifelse(site - end > 5, index_cut - 1, index_cut)) %>%
    mutate(end_smoothed = orange_lines$site[index_cut]) %>%
    #dplyr::mutate(end_smoothed = if_else(site - end > 5, orange_lines$site[index_cut - 1], site))
    dplyr::select(-c(site, index_cut, end_site))
  
  deletions_insertions = deletions_insertions %>% dplyr::mutate(end_smoothed = 
                                                                  ifelse(mutation_type == 'ins', start_smoothed, end_smoothed))
  
  deletions_insertions = mutate(deletions_insertions, 
                                smoothed_id = ifelse(mutation_type == 'del', 
                                                     paste0("del_", start_smoothed, "-", end_smoothed),
                                                     paste0("ins_", start_smoothed, "-", n_nucleotides, "nts")))
  
  return(deletions_insertions)
}

tidy_alignment_smoothed = function(REvoBC_object, smoothed_df) {
  del_sub_ins_df = REvoBC_object$alignment$asv_barcode_alignment
  # In case of smoothed deletions, we need to re-create a tibble where each row corresponds
  # to a position in each ASV. This tibble already exists for non-smoothed mutations
  # and now we want to replace its column "alt" with the status of each nucleotide 
  # in case we consider smoothed deletions and/or insertions.
  # Insertions and substitutions are characterized by only one position, so in order to create the tidy 
  # alignment tibble we start by removing them from the smoothed dataframe, and then we re-insert them. 
 
  smoothed_sub = filter(smoothed_df, mutation_type =='sub')
  smoothed_ins = filter(smoothed_df, mutation_type =='ins')
  smoothed_del = filter(smoothed_df, mutation_type == 'del')
  # Non dovrebbe servire, io assegno la mutazione smooothed prima in base a se trovo una
  # delezione smoothed. Poi, a tutti gli altri nucleotidi, se compaiono tra inserzioni o
  # sotituzioni smoothed allora gli assegno il tipo
  #del_sub_ins_df = filter(del_sub_ins_df, alt == 'del') 
  
  to_plot_df = list()
  for (asv in unique(as.character(smoothed_df$asv_names))) {
    tmp_smoothed_del = dplyr::filter(smoothed_del, asv_names == asv)
    # Remove insertinons from tidy alignment, as they are a duplicated position
    # of a wildtype. We remove them and then add to the tidy dataframe the smoothed insertions positions.
    tmp_tidy = dplyr::filter(del_sub_ins_df, asv_names == asv & alt != 'ins')
    tmp_smoothed_ins = dplyr::filter(smoothed_ins, asv_names == asv) 
    tmp_smoothed_sub = dplyr::filter(smoothed_sub, asv_names == asv) 
    
    tmp_tidy$alt <- ifelse(sapply(tmp_tidy$position_bc260, function(p) 
      any(tmp_smoothed_del$start_smoothed <= p & tmp_smoothed_del$end_smoothed >= p)),
      "del", 'smoothed_wt')
    
    tmp_tidy = ungroup(tmp_tidy) %>% 
      mutate(alt = ifelse(position_bc260 %in% tmp_smoothed_sub$start, 'sub', alt))
    
    tmp_tidy = tmp_smoothed_ins %>% 
      rename(position_bc260 = start_smoothed,
             read_asv = alt_seq, alt = mutation_type) %>%
      mutate(ref_asv = '-') %>% select(colnames(tmp_tidy)) %>%
      bind_rows(tmp_tidy) %>%
      dplyr::arrange(position_bc260)
    
    # tmp_tidy = ungroup(tmp_tidy) %>% 
    #   mutate(alt = ifelse(position_bc260 %in% tmp_smoothed_ins$start, 'ins', 
    #                                ifelse(position_bc260 %in% tmp_smoothed_sub$start, 'sub', alt)))

    to_plot_df[[asv]] = tmp_tidy
  }
  for (asv in setdiff(del_sub_ins_df$asv_names, names(to_plot_df))) {
    tmp_tidy = dplyr::filter(del_sub_ins_df, asv_names == asv)
    tmp_tidy$alt = 'non_smoothed_del'
    to_plot_df[[asv]] = tmp_tidy
  }
  to_plot_df = dplyr::bind_rows(to_plot_df)
  return(to_plot_df)
}




