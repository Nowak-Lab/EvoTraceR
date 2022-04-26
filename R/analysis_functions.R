perform_flanking_filtering = function(barcodes_info, seqtab_df, flanking_filtering) {
  # Removal of contamination: keep only those ASV that start (5 nucleotides) and end (10 nuceotides) like the original barcode.
  # The first 5 nts are the same for all barcodes, while the end is barcode-specific
  RD1_10 <- barcodes_info$ref_flank_left
  RD2_10 <- barcodes_info$ref_flank_right
  
  nmbc <- paste0( barcodes_info$ref_name, ".NMBC")
  # filter based on 5' and 3' 10x nts of
  if (flanking_filtering == 'both') {
    seqtab_df <- dplyr::filter(seqtab_df,
                               stringr::str_detect(string = seq, pattern = !!RD1_10) & stringr::str_detect(string = seq, pattern = !!RD2_10)) # the same for different barcodes: 1.0 - site less affected
    
  } else if (flanking_filtering == 'right') {
    seqtab_df <- dplyr::filter(seqtab_df,
                               stringr::str_detect(string = seq, pattern = !!RD2_10)) # the same for different barcodes: 1.0 - site less affected
    
  } else if(flanking_filtering == 'left') {
    seqtab_df <- dplyr::filter(seqtab_df,
                               stringr::str_detect(string = seq, pattern = !!RD1_10)) # the same for different barcodes: 1.0 - site less affected
    
  } else if (flanking_filtering == 'either'){
    seqtab_df <- dplyr::filter(seqtab_df,
                               stringr::str_detect(string = seq, pattern = !!RD1_10) | stringr::str_detect(string = seq, pattern = !!RD2_10)) # the same for different barcodes: 1.0 - site less a
  } else {
    stop('Flanking filtering must be one of "left", "right", "both" or "either". Exiting.')
  }
  
  return(seqtab_df)
  
}



# This function takes the ASVs and cleans them, by collapsing the ones that differ from one another only by substitutions.
asv_collapsing = function(seqtab, 
                          barcode,
                          pwa_match,
                          pwa_mismatch,
                          pwa_gapOpening,
                          pwa_gapExtension, sample.names,
                          barcode_name,
                          cut_sites,
                          cleaning_window = c(3,3)) {
  
  dnastringset <- Biostrings::DNAStringSet(seqtab$seq) 
  names(dnastringset) <- seqtab$seq_names
  
  
  # seqtab = seqtab %>% mutate(pairwise_sim =  list(Biostrings::pid(Biostrings::pairwiseAlignment(subject=seq, pattern=dnastringset))))
  # 
  # pid_matrix = seqtab$pairwise_sim
  # names(pid_matrix) = seqtab$seq_names
  # pid_matrix = bind_rows(pid_matrix)
  # Output as FASTA files
  
  mx_crispr <- Biostrings::nucleotideSubstitutionMatrix(match = pwa_match, mismatch = pwa_mismatch, baseOnly = TRUE)
  
  mpwa <- Biostrings::pairwiseAlignment(subject = barcode, 
                                        pattern = dnastringset, 
                                        substitutionMatrix = mx_crispr,
                                        gapOpening = pwa_gapOpening,
                                        gapExtension = pwa_gapExtension,
                                        type = 'global')
  cli::cli_alert_info('Computing pattern alignment')
  aligned_sequences = as.data.frame(Biostrings::alignedPattern(mpwa))
  cli::cli_alert_info('Computing subject alignment')
  aligned_reference = as.data.frame(Biostrings::alignedSubject(mpwa))
  
  
  cli::cli_alert_info('Computing tidy alignment for cleaning')
  alignment_tidy_ref_alt <- foreach::foreach(i = seq(1, length(mpwa)), .combine=rbind) %do% {
    
    # postaligned_seqs <- Biostrings:::.makePostalignedSeqs(mpwa[i])[[1]]
    # aligned_sequence = as.character(postaligned_seqs[1])
    # aligned_reference = as.character(postaligned_seqs[2])
    
    data.frame("seq_names" = rownames(aligned_sequences)[i], #names(postaligned_seqs)[1], #names(dnastringset)[i],#
               "read_asv" = strsplit(aligned_sequences$x[i], split='')[[1]], #strsplit(aligned_sequence, split='')[[1]],#
               "ref_asv" = strsplit(aligned_reference$x[i], split='')[[1]],#strsplit(aligned_reference, split='')[[1]], #
               "position" = seq(1, nchar(aligned_reference$x[i])))#seq(1, nchar(aligned_reference))) #
  }
  alignment_tidy_ref_alt = alignment_tidy_ref_alt %>% arrange(position, seq_names)
  
  
  
  # Assign Alterations types; wt - wild type, del - deletions, ins - insertion, sub - substitution
  # In cases where ref_asv == "-" & read_asv == "-", then the alteration type is set to ins_smwr 
  alignment_tidy_ref_alt <-
    alignment_tidy_ref_alt %>%
    mutate(alt = ifelse(ref_asv == read_asv, "w",
                        ifelse(read_asv == "-", "d", 
                               ifelse(ref_asv == "-", "i", "s")))) %>%
    mutate(alt = ifelse(ref_asv == "-" & read_asv == "-", "ins_smwr", alt)) 
  
  # algn = perform_collapsing(alignment_tidy_ref_alt, seqtab, sample.names)
  # 
  # alignment_tidy_ref_alt = alignment_tidy_ref_alt %>% filter(seq_names %in% algn$seq_names)
  
  # Now translate barcode to barcode 260 nts long.
  # In case of sites of insertion, the last position not affected by the insertion is propagated
  # for all the nucleotides affected by that insertion.
  alignment_tidy_ref_alt <-
    alignment_tidy_ref_alt %>%
    mutate(position_bc260 = ifelse(ref_asv == "-", NA, position)) %>% # add "bc260" scale
    group_by(seq_names) %>%
    mutate(position_bc260 = data.table::nafill(position_bc260, "locf")) %>% # fill NA with consecutive numbers
    mutate(position_bc260 = tidyr::replace_na(position_bc260, 0)) %>% #locf doesn't fill NAs that are at the beginning of the sequence
    dplyr::select(seq_names,  position, position_bc260, ref_asv, read_asv, alt)#,  sample, perc_in_sample,)  
  
  alignment_tidy_ref_alt = alignment_tidy_ref_alt %>% group_by(seq_names) %>%
    mutate(position_bc260 = ifelse(position_bc260 == 0, position_bc260, position_bc260 - min(position_bc260[position_bc260 != 0]))) %>%
    mutate(position_bc260 = as.numeric(as.character(factor(x = position_bc260, levels = unique(position_bc260),
                                                           labels = seq_along(unique(position_bc260))))))
  
  # algn_to_long = algn %>% tidyr::pivot_longer(all_of(position_columns), names_to = 'position', values_to = 'alt') %>%
  #   filter(!is.na(alt))
  # 
  # algn_to_long = algn_to_long %>% select(c(seq_names, position, alt)) %>%
  #   mutate(position = as.integer(position)) %>%
  #   arrange(seq_names, position, alt)
  # subset_tidy = subset_tidy %>% select(c(seq_names, position, alt)) %>%
  #   arrange(seq_names, position, alt)
  
  coord = mutation_coordinate_matrix(alignment_tidy_ref_alt, barcode_name)
  # Select sequences that have no indels, they are collapsed to the original barcode
  no_indels =  setdiff(seqtab$seq_names, coord$seq_names)
  no_indels = seqtab %>% filter(seq_names %in% no_indels)
  no_indels = no_indels %>% summarise(seq = seq[which.max(totalCounts)], #compute_consensus_sequence(seq, sum_counts), 
                                      seq_names = seq_names[which.max(totalCounts)],
                                      across(sample.names, sum))
  
  
  # I take the coordinate matrix, I clean it and I binarize it
  # and then I collapse the sequences that are the same
  
  # Clean the mutations coordinate matrix
  coord_cleaned = clean_mutations(coord, 
                                  orange_lines = cut_sites, 
                                  left_right_window = cleaning_window)
  coord_cleaned = coord_cleaned %>% mutate(asv_names = seq_names)
  binary_matrix = coordinate_to_binary(coord_cleaned, barcode_name)
  
  # Join the binary matrix with the sequences dataframe, in order to have counts
  mutations = unique(coord_cleaned$mut_id)
  binary_matrix$seq_names = rownames(binary_matrix) 
  seqtab_collapsed = binary_matrix %>% dplyr::left_join(seqtab) %>% filter(!is.na(totalCounts))
  # Collapse sequences that have the same cleaned mutations
  seqtab_collapsed = seqtab_collapsed %>% group_by(across(all_of(mutations))) %>% 
    summarize(consensus_seq = seq[which.max(totalCounts)], #compute_consensus_sequence(seq, sum_counts), 
              seq_names = seq_names[which.max(totalCounts)],
              across(sample.names, sum)) %>% ungroup() %>%
    rename(seq = consensus_seq) %>%
    select(-all_of(mutations))
  
  seqtab_collapsed = dplyr::bind_rows(seqtab_collapsed, no_indels)
  
  alignment_tidy_ref_alt = alignment_tidy_ref_alt %>% filter(seq_names %in% seqtab_collapsed$seq_names)
  # Ora devo anche pulire le mutazioni nel tidy
  clean_tidy_alignment = tidy_alignment_cleaned(alignment_tidy_ref_alt, 
                                                coord_cleaned, 
                                                no_indels %>% pull(seq_names))
  
  return(list(seqtab_df = seqtab_collapsed, tidy_alignment = clean_tidy_alignment,
              mutations_coordinates = coord_cleaned %>% select(-c(asv_names)), binary_matrix = binary_matrix))
}


# perform_collapsing <- function(alignment_tidy_ref_alt, seqtab, sample.names) {
#   algn = alignment_tidy_ref_alt %>% ungroup() %>% select(c(seq_names, position, alt)) %>%
#     tidyr::pivot_wider(names_from = 'position', values_from = 'alt') 
#   
#   position_columns = setdiff(colnames(algn), "seq_names")
#   
#   
#   algn = algn %>% 
#     inner_join(seqtab %>% mutate(seq_names = gsub(pattern='.NMBC', replacement = '', x = seq_names)), by = "seq_names") %>%
#     #tibble::column_to_rownames("seq_names")  %>% 
#     mutate(sum_counts = seqtab %>% select(all_of(sample.names)) %>% rowSums(na.rm = TRUE))
#   
#   
#   algn = algn %>% group_by(across(all_of(position_columns))) %>% 
#     summarise(consensus_seq = seq[which.max(sum_counts)], #compute_consensus_sequence(seq, sum_counts), 
#               seq_names = seq_names[which.max(sum_counts)],
#               across(sample.names, sum)) %>%
#     ungroup() %>% 
#     select(-all_of(position_columns))
#   
#   # algn$seq_names = paste0("SEQ", formatC(c(1:(nrow(algn))), 
#   #                                        width = nchar(trunc(nrow(algn))), 
#   #                                        format = "d", flag = "0")) # -1 to start 00 with no changes sequence
#   
#   algn = algn %>% rename(seq = consensus_seq) %>% tibble::column_to_rownames("seq")
#   algn$seq = rownames(algn)
#   
#   return(algn)
# }



# Compute the frequency of the different ASV in each organ/day and the frequency of counts in each organ for every ASV.
# Compute also, for each ASV the sample in which the frequency is maximum (store this information in a column \code{perc_fold_to_max} of tibble seqtab_df_clean_asv_long in the EvoTraceR object)
# Then, compute for each sample, indices of diversity of ASVs in each samples (Shannon Entropy, Simpson Index).
# It stores all results in the field \code{statistics} of the EvoTraceR object.
# 
# EvoTraceR_object where the statistics on ASV will be computed. 
# sample. List of the column names containing the organs/days.
# asv_count_cutoff. Cutoff on the minimum number of counts for an ASV to be 
# pwa. Object resulted from the pairwise alignment performed on the ASVs.
asv_statistics <- function(EvoTraceR_object, sample_columns, asv_count_cutoff, figure_dir, nmbc, output_dir) {
  seqtab_df_clean_asv = EvoTraceR_object$clean_asv_dataframe
  
  seqtab_df_clean_asv_long <-
    tibble(seqtab_df_clean_asv) %>%
    dplyr::select(-seq) %>% 
    tidyr::gather(sample, count, sample_columns, factor_key=TRUE) %>% # from wide to long -> - 1 for seq
    filter(count > asv_count_cutoff) %>%
    group_by(sample) %>%
    dplyr::mutate(count = as.numeric(count)) %>%
    dplyr::mutate(perc_in_sample = prop.table(count)*100) %>% # calculate percentage during i.e. Day or Organ
    dplyr::ungroup() %>% dplyr::mutate_if(is.numeric, round, digits=2) # round percent to two digits
  
  seqtab_df_clean_asv_long <-
    tibble(seqtab_df_clean_asv_long) %>%
    group_by(asv_names) %>%
    dplyr::mutate(count = as.numeric(count)) %>%
    dplyr::mutate(perc_asv = prop.table(count)*100) %>% 
    ungroup() %>%
    dplyr::mutate_if(is.numeric, round, digits=2)
  
  ### Order of Rows and Columns
  # order for columns summary bar based on total counts
  seq_names_ord <- 
    seqtab_df_clean_asv_long %>%
    group_by(asv_names) %>%
    summarise(sum = sum(count)) %>% # sum of all of ASVs
    arrange(-sum) %>%
    pull(asv_names) %>%
    as.character()
  
  # prepare levels and orders for days or organs
  sample_order = EvoTraceR_object$sample_order
  
  seqtab_df_clean_asv_long <- 
    seqtab_df_clean_asv_long %>%
    dplyr::mutate(sample = forcats::fct_relevel(sample, sample_order)) %>%
    arrange(match(sample, sample_columns))
  
  # prepare levels and orders of asv_names (ASVs)
  seqtab_df_clean_asv_long <- 
    seqtab_df_clean_asv_long %>%
    dplyr::mutate(asv_names = forcats::fct_relevel(asv_names, seq_names_ord)) %>%
    arrange(match(asv_names, seq_names_ord))
  
  # Normalization of Counts in ASVs (%) 
  # Rule 1.: max number perc calculation based on count - the highest number of reads will be 100% or 1.0
  # Rule 2.: log2 calculation based on count - first colony appearance value numbers of reads different than 0 in Day type experiment
  # pick the highest count in asv_names group (ASVs)
  
  seqtab_df_clean_asv_long <-
    seqtab_df_clean_asv_long %>%
    group_by(asv_names) %>%
    mutate(perc_fold_to_max = 100 * (count/max(count))) %>%
    arrange(-count)
  EvoTraceR_object$statistics$asv_df_percentages = seqtab_df_clean_asv_long
  
  # save as Data csv
  write.csv(seqtab_df_clean_asv_long, 
            paste0(output_dir, "/asv_df_percentages.csv"),
            row.names = F)
  
  ### Abundance and Richness Calculations
  # asv_names summary
  asv_names_stat <-
    seqtab_df_clean_asv_long %>%
    group_by(asv_names) %>%
    dplyr::summarize(
      abundance_asv_total = sum(count), 
      richness_asv_total = sum(count > 0, na.rm = TRUE)) 
  # day or organ summary
  sample_stat <-
    seqtab_df_clean_asv_long %>%  
    group_by(sample) %>%
    dplyr::summarize(
      abundance_asv_persample = sum(count), 
      richness_asv_persample = sum(count > 0, na.rm = TRUE))
  # Store in object and in csv file
  EvoTraceR_object$statistics$asv_totalCounts = asv_names_stat
  EvoTraceR_object$statistics$sample_totalcounts = sample_stat
  
  write.csv(asv_names_stat, 
            file.path(output_dir, "asv_totalCounts.csv"),
            row.names = F)
  
  write.csv(sample_stat, 
            file.path(output_dir, "sample_totalcounts.csv"),
            row.names = F)
  
  # merge data
  seqtab_df_clean_asv_long <- tibble(merge(seqtab_df_clean_asv_long, asv_names_stat, "asv_names")) 
  seqtab_df_clean_asv_long <- tibble(merge(seqtab_df_clean_asv_long, sample_stat, "sample")) 
  
  ### Diversity Calculations (Package "benthos") 
  # index calcuations
  diversity <-
    seqtab_df_clean_asv_long %>%
    mutate(asv_names = as.character(asv_names)) %>%
    dplyr::group_by(sample, .drop=FALSE) %>%  
    dplyr::summarise(
      # Measures of clonal richness
      abundance_asv_persample = benthos::total_abundance(count = count), 
      richness_asv_persample = benthos::species_richness(taxon = asv_names, count = count), 
      shannons_index_persample = benthos::shannon(taxon = asv_names, count = count, base=exp(1)), 
      richness_asv_persample = unique(richness_asv_persample)
    ) %>%
    mutate(pielous_evenness_persample = shannons_index_persample / log(richness_asv_persample) )
  
  EvoTraceR_object$statistics$asv_diversity_persample = diversity
  write.csv(diversity, 
            file.path(output_dir, "asv_diversity_persample.csv"))
  
  # matrix rows: organ_day and coumns ASV names
  mx_freq <- 
    seqtab_df_clean_asv_long %>%
    dplyr::select(sample, asv_names, count) %>%
    distinct() %>%
    tidyr::pivot_wider(names_from = asv_names, values_from = count, values_fill = 0) %>%
    tibble::column_to_rownames(var="sample") %>% 
    as.matrix()
  # Save as Data csv
  mx_freq= mx_freq[,c(nmbc, sort(colnames(mx_freq)[colnames(mx_freq)!=nmbc]))]
  write.csv(mx_freq, file.path(output_dir, "asv_persample_frequency.csv"))
  
  mx_freq_bin <- as.matrix(ifelse(mx_freq == 0, 0, 1))
  # Save as Data csv
  write.csv(mx_freq_bin, file.path(output_dir, "asv_persample_detection.csv"))
  
  EvoTraceR_object$statistics$asv_persample_frequency = mx_freq
  EvoTraceR_object$statistics$asv_persample_detection = mx_freq_bin
  
  df_to_plot_perf_match <- dplyr::inner_join(x=seqtab_df_clean_asv_long, 
                                             y=dplyr::select(EvoTraceR_object$clean_asv_dataframe, seq, asv_names), 
                                             by="asv_names") # %>%
  # add_row(data.frame(asv_names = EvoTraceR_object$reference$ref_name, 
  #                    seq = EvoTraceR_object$reference$ref_seq, 
  #                    stringsAsFactors = F))
  
  
  # Count length of barcode seq
  df_to_plot_perf_match$seq_n <- nchar(df_to_plot_perf_match$seq)
  
  pwa <- Biostrings::pairwiseAlignment(subject = EvoTraceR_object$reference$ref_seq, #df_to_plot_perf_match[stringr::str_detect(df_to_plot_perf_match$asv_names,'ORG|NMBC'),'seq']),
                                       pattern = df_to_plot_perf_match$seq,
                                       type="global", gapOpening = 20, gapExtension = 1)
  # TOLTO PERCHE FORSE NON SERVE
  # Perform pairwise Alignment Stat
  # pid = Computes the percent sequence identity
  df_to_plot_perf_match$pid <- Biostrings::pid(pwa)
  # nedit = Computes the Levenshtein edit distance of the alignments
  df_to_plot_perf_match$nedit <- Biostrings::nedit(pwa)
  # score = Extracts the pairwise sequence alignment scores
  df_to_plot_perf_match$alignment_score <- Biostrings::score(pwa)
  #
  # # Save as Data csv
  write.csv(df_to_plot_perf_match[c('asv_names','seq', 'pid', 'nedit', 'alignment_score')],
            file.path(output_dir, "asv_toBarcode_similarity.csv"))
  
  EvoTraceR_object$statistics$asv_toBarcode_similarity = df_to_plot_perf_match[c('asv_names','seq', 'pid', 'nedit', 'alignment_score')]
  
  EvoTraceR_object$statistics$all_asv_statistics = df_to_plot_perf_match
  # Histogram of Length    
  # skip nmbc
  data_for_hist = df_to_plot_perf_match %>% filter(stringr::str_detect(asv_names, "ASV"))
  if (!is.null(figure_dir)) {
    hist_seq_count <- 
      ggplot(data=data_for_hist, 
             aes(y=perc_in_sample, x=seq_n)) + # remove NMBC (usually too big and masking smaller ASVs) and ORG (because of NAs)
      geom_bar(stat="identity", position = "stack", fill="#B484A9", width=1) +
      scale_x_continuous(labels=scales::comma, breaks=c(1, seq(floor(nchar(EvoTraceR_object$reference$ref_seq)*2/10), 
                                                               nchar(EvoTraceR_object$reference$ref_seq)*2, 
                                                               floor(nchar(EvoTraceR_object$reference$ref_seq)*2/10))), 
                         limits=c(0, nchar(EvoTraceR_object$reference$ref_seq)*2), 
                         expand = c(0.01, 0.01)) +
      scale_y_continuous(labels=function(x) paste0(x, "%"),
                         limits=c(0, 1.25* max(data_for_hist$perc_in_sample)),  
                         expand = c(0.01, 0.01)) +
      xlab("ASV Length") +
      ylab("ASV Frequency") +
      geom_vline(xintercept=nchar(EvoTraceR_object$reference$ref_seq), linetype="dotted", size=0.25, col="#84B48F") + # expected size
      lemon::coord_capped_cart(left="both", bottom="left") + # axis with lemon
      lemon::facet_rep_grid(rows = vars(sample), repeat.tick.labels = TRUE) + 
      # add theme
      barplot_nowaklab_theme() +
      theme(aspect.ratio = 1/4,
            plot.margin = unit(c(0, 2, -2, -2), "mm"),
            legend.position = "None")
    
    # Save PDF
    ggsave(filename=file.path(figure_dir, "histogram_sequenceLength.pdf"), 
           plot=hist_seq_count, 
           #device=cairo_pdf, 
           width=25, 
           height=5*length(sample_columns), 
           units = "cm") #17.5 for 4x
    #pg <- ggplot_build(hist_seq_count)
    write.csv(data_for_hist, #pg$data[[1]],  
              file.path(figure_dir, "/histogram_sequenceLength_data.csv"),
              row.names = FALSE, quote = FALSE)
  }
  
  return(EvoTraceR_object)
}
