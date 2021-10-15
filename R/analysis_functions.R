# Compute the frequency of the different ASV in each organ/day and the frequency of counts in each organ for every ASV.
# Compute also, for each ASV the sample in which the frequency is maximum (store this information in a column \code{perc_fold_to_max} of tibble seqtab_df_clean_asv_long in the REvoBC object)
# Then, compute for each sample, indices of diversity of ASVs in each samples (Shannon Entropy, Simpson Index).
# It stores all results in the field \code{statistics} of the REvoBC object.
# 
# REvoBC_object where the statistics on ASV will be computed. 
# sample. List of the column names containing the organs/days.
# asv_count_cutoff. Cutoff on the minimum number of counts for an ASV to be 
# pwa. Object resulted from the pairwise alignment performed on the ASVs.
asv_statistics <- function(REvoBC_object, sample_columns, asv_count_cutoff, figure_dir, nmbc, output_dir) {
  seqtab_df_clean_asv = REvoBC_object$clean_asv_dataframe
  
  seqtab_df_clean_asv_long <-
    tibble(seqtab_df_clean_asv) %>%
    dplyr::select(-seq) %>% # remove sequence => not useful for census analysis
    tidyr::gather(sample, count, sample_columns, factor_key=TRUE) %>% # from wide to long -> - 1 for seq
    filter(count > asv_count_cutoff) %>% # remove data that percentage less than i.e. 0.5% -> visualization
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
    dplyr::select(asv_names) %>%
    pull() %>%
    as.character()
  
  # prepare levels and orders for days or organs
  seqtab_df_clean_asv_long <- 
    seqtab_df_clean_asv_long %>%
    dplyr::mutate(sample = forcats::fct_relevel(sample, sample_columns)) %>%
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
  REvoBC_object$statistics$asv_df_percentages = seqtab_df_clean_asv_long
  
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
  REvoBC_object$statistics$asv_totalCounts = asv_names_stat
  REvoBC_object$statistics$sample_totalcounts = sample_stat
  
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
      # D_Margalef = benthos::margalef(taxon = asv_names, count = count), 
      # Ryggs = benthos::rygg(taxon = asv_names, count = count), 
      # Ryggs_adjusted = benthos::rygg(taxon = asv_names, count = count, adjusted = TRUE), 
      # # Measures of heterogeneity/evenness
      # Simpsons_Index_D = benthos::simpson(taxon = asv_names, count = count), 
      # HurlbertsProbability_PIE = benthos::hpie(taxon = asv_names, count = count), 
      # Hills_Diversity_Num1 = benthos::hill1(taxon = asv_names, count = count), 
      # Hills_Diversity_Num2 = benthos::hill2(taxon = asv_names, count = count), 
      shannons_index_persample = benthos::shannon(taxon = asv_names, count = count, base=exp(1)), 
      # Simpsons_Reciprocal_Index = 1/benthos::simpson(taxon = asv_names, count = count), 
      richness_asv_persample = unique(richness_asv_persample)
    ) %>%
    mutate(pielous_evenness_persample = shannons_index_persample / log(richness_asv_persample) )
  
  REvoBC_object$statistics$asv_diversity_persample = diversity
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
  
  REvoBC_object$statistics$asv_persample_frequency = mx_freq
  REvoBC_object$statistics$asv_persample_detection = mx_freq_bin
  
  df_to_plot_perf_match <- dplyr::inner_join(x=seqtab_df_clean_asv_long, 
                                             y=dplyr::select(REvoBC_object$clean_asv_dataframe, seq, asv_names), 
                                             by="asv_names") %>%
    add_row(dplyr::select(REvoBC_object$barcode, -c("seq_start", "seq_end")))
  

  # Count length of barcode seq
  df_to_plot_perf_match$seq_n <- nchar(df_to_plot_perf_match$seq)
  
  pwa <- Biostrings::pairwiseAlignment(subject = REvoBC_object$barcode$seq, #df_to_plot_perf_match[stringr::str_detect(df_to_plot_perf_match$asv_names,'ORG|NMBC'),'seq']), 
                                       pattern = df_to_plot_perf_match$seq, 
                                       type="global", gapOpening = 20, gapExtension = 1)
  
  # Perform pairwise Alignment Stat
  # pid = Computes the percent sequence identity
  df_to_plot_perf_match$pid <- Biostrings::pid(pwa)
  # nedit = Computes the Levenshtein edit distance of the alignments
  df_to_plot_perf_match$nedit <- Biostrings::nedit(pwa)
  # score = Extracts the pairwise sequence alignment scores
  df_to_plot_perf_match$alignment_score <- Biostrings::score(pwa)
  
  # Save as Data csv
  write.csv(df_to_plot_perf_match[c('asv_names','seq', 'pid', 'nedit', 'alignment_score')], 
            file.path(output_dir, "asv_toBarcode_similarity.csv"))
  
  REvoBC_object$statistics$asv_toBarcode_similarity = df_to_plot_perf_match[c('asv_names','seq', 'pid', 'nedit', 'alignment_score')]
  
  REvoBC_object$statistics$all_asv_statistics = df_to_plot_perf_match
  # Histogram of Length    
  # skip nmbc
  data_for_hist = df_to_plot_perf_match %>% filter(!stringr::str_detect(asv_names, "NMBC|ORG"))
  if (!is.null(figure_dir)) {
    hist_seq_count <- 
      ggplot(data=data_for_hist, 
             aes(y=perc_in_sample, x=seq_n)) + # remove NMBC (usually too big and masking smaller ASVs) and ORG (because of NAs)
      geom_bar(stat="identity", position = "stack", fill="#B484A9", width=1) +
      scale_x_continuous(labels=scales::comma, breaks=c(1, seq(52, 520, 52)), 
                         limits=c(0, 520), 
                         expand = c(0.01, 0.01)) +
      scale_y_continuous(labels=function(x) paste0(x, "%"),
                         limits=c(0, 1.25* max(data_for_hist$perc_in_sample)),  
                         expand = c(0.01, 0.01)) +
      xlab("ASV Length") +
      ylab("ASV Frequency") +
      geom_vline(xintercept=260, linetype="dotted", size=0.25, col="#84B48F") + # expected size
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
  
  return(REvoBC_object)
}
