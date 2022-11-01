
# Create the binary mutations matrix
mutation_coordinate_matrix = function(tidy_alignment, barcode_name) {
  
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
  
  mut_df = tidy_alignment %>% ungroup() %>% 
    arrange(seq_names, position_bc260) %>%
    group_by(seq_names) %>% 
    mutate(contig_alt = cumsum(c(1, head(alt, -1) != tail(alt, -1)))) %>%
    filter(alt!= 'w' & alt != 'ins_smwr' & alt != 's') %>%
    dplyr::rename(mutation_type = alt, ref_seq = ref_asv, alt_seq = read_asv) 
  
  # Now we group alterations first based on the ASV, and then based on the new column contig_alt
  # in this way we can count the rows in the group to obtain the length of the insertion
  # We finally insert back the substitutions.
  mut_df = mut_df %>% group_by(seq_names, mutation_type, contig_alt) %>%
    dplyr::summarise(start = min(position_bc260), 
                     end = max(position_bc260), 
                     n_nucleotides = n(), 
                     ref_seq = paste(ref_seq, collapse = ''),
                     alt_seq = paste(alt_seq, collapse = ''),
                     .groups='drop') %>%
    dplyr::select(seq_names, mutation_type, start, end, n_nucleotides, ref_seq, alt_seq) 
  
  mut_df$mut_id = paste0(mut_df$mutation_type, "_", mut_df$start, "_", mut_df$n_nucleotides, "nts")
  
  return(mut_df)
}

# mut_df is a dataframe with start and end of insertions and deletions in each ASV (no substitutions)
# This function removes the small events that happen too far from any cut site and counts the number of 
# target sites spanned by the deletion
# For each indel, we expand it to the left and to the right by left_right_window[1] and left_right_window[2] positions
# respectively, and then we check if the extended event spans any of the cut sites.
clean_mutations = function(mut_df, orange_lines, left_right_window = c(3,3)) {
  
  deletions_insertions = mut_df %>% 
    mutate(extended_start = start - left_right_window[1], 
           extended_end = end + left_right_window[2]) %>%
    rowwise() %>%
    mutate(spanned_cutSites = list(intersect(seq(extended_start, extended_end), orange_lines))) %>%
    mutate(n_sites = length(spanned_cutSites)) %>%
    filter(n_sites > 0)
  
  return(deletions_insertions)
}

coordinate_to_binary = function(mut_df, barcode) {
  mut_df_wide = plyr::count(mut_df, vars = c("asv_names", "mut_id")) %>% 
    pivot_wider(names_from = mut_id, values_from = freq) %>% tibble::column_to_rownames("asv_names")
  
  mut_df_wide[is.na(mut_df_wide)] <- 0
  mut_df_wide[mut_df_wide >= 1] = 1
  
  # Add row to binary mutation matrix corresponding to the original barcode (i.e. all mutations = 0)
  bc_mut = as.list(rep(0, ncol(mut_df_wide)))
  names(bc_mut) = colnames(mut_df_wide)
  mut_df_wide = dplyr::bind_rows(data.frame(bc_mut, row.names = barcode),
                                 mut_df_wide)
  
  return(mut_df_wide)
  
}

tidy_alignment_cleaned = function(tidy_alignment_full, cleaned_df, barcode) {
  cleaned_df = cleaned_df %>% 
    filter(!is.na(start)) %>% 
    rowwise() %>% mutate(positions = list(seq(start, end)))
  
  paste_fun = function(name, pos, alt_type) paste(name, pos, alt_type)
  valid_positions = unlist(mapply(paste_fun, cleaned_df$seq_names, cleaned_df$positions, cleaned_df$mutation_type))
  
  tidy_alignment_full = tidy_alignment_full %>% mutate(tmp_col = paste(seq_names, position_bc260, alt))
  tidy_alignment_full_new = tidy_alignment_full %>% mutate(alt = ifelse(tmp_col %in% valid_positions, alt, 'w') )
  tidy_alignment_full_new = tidy_alignment_full_new %>% 
    filter(seq_names %in% c(barcode, cleaned_df$seq_names))
  tidy_alignment_full_new = tidy_alignment_full_new %>% filter(!(ref_asv == '-' & alt == 'w')) %>%
    select(-c(tmp_col))
  
  return(tidy_alignment_full_new)
}




