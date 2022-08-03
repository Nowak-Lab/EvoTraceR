count_alterations <- function(EvoTraceR_object) {
  
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
  
  EvoTraceR_object$alignment$ASV_alterations_width = alignment_tidy_ref_alt_mrg_final_width_summ
  
  # Select only deletions and substitutions
  del_sub_df <- 
    alignment_tidy_ref_alt_mrg_final %>%
    filter(alt != "i") %>% 
    dplyr::select(asv_names, sample, position_bc260, alt, perc_in_sample)
  
  # Leave only insertion with ASV
  ins_df <- 
    alignment_tidy_ref_alt_mrg_final %>%
    filter(alt == "i") %>%
    group_by(asv_names, sample, perc_in_sample) %>% 
    mutate(cons_bin = c(0, abs(diff(position_bc260)) == 1)) %>% # find if number is consecutive = 0, if not = 1
    #filter(cons_bin == 0) %>%
    dplyr::select(asv_names, sample, position_bc260, alt, perc_in_sample)
  
  del_sub_ins_df <- 
    rbind(del_sub_df, ins_df) #%>%
  
  sample_columns = setdiff(colnames(EvoTraceR_object$asv_prefilter), c("seq_names", "seq"))
  
  sample_order = EvoTraceR_object$sample_order
  
  # prepare levels and orders
  del_sub_ins_df <- 
    del_sub_ins_df %>%
    dplyr::mutate(sample = forcats::fct_relevel(sample, sample_order)) %>%
    arrange(match(sample, sample_columns))
  
  EvoTraceR_object$alignment$mutations_df = del_sub_ins_df
  return(EvoTraceR_object)
  
}
