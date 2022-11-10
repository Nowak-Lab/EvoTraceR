# This function calls Trimmomatic and Flash.
alignment_pipeline = function(fnFs,
                              fnRs,
                              output_dir,
                              map_file_sample,
                              sample.names,
                              trimmomatic_path,
                              flash_path
                              ) {
  fnFs <- gsub(" ", "\ ", fnFs)
  fnRs <- gsub(" ", "\ ", fnRs)
  # Make directories for the filtered fastqs and (optionally) for figures
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = T)

  fastq_dir = file.path(output_dir, "fastq_analysis/")
  trimmed_dir = file.path(fastq_dir, "fastq_trimmed/")
  flash_input_dir = file.path(fastq_dir, "fastq_flash_merged/")
  figure_dir = file.path(output_dir, "asv_analysis/")
  
  if (!dir.exists(fastq_dir)) dir.create(fastq_dir, recursive = T)
  if (!dir.exists(trimmed_dir)) dir.create(trimmed_dir, recursive = T)
  if (!dir.exists(flash_input_dir)) dir.create(flash_input_dir, recursive = T)
  if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = T)
  
  # filter and trim
  cli::cli_alert_info('Filtering and trimming')
  
  # Users will input the path up to the trimmomatic binary. I need to go up to one level wrt that, in order to find the adapters folder
  trimmomatic_maindir = dirname(trimmomatic_path)
  
  track_reads = data.frame(input = numeric(), trimmed = numeric(), merged = numeric())
  for (i in seq(1, length(fnFs))) {
    sample = sample.names[i]
    
    track_reads[map_file_sample[sample, 'sample'], 'input'] = length(Biostrings::readDNAStringSet(fnFs[i], format = 'fastq'))
    
    bashCallTrim = paste0("java -jar '", trimmomatic_path, "' PE '",
                          fnFs[i], 
                          "' '", fnRs[i], "' '", trimmed_dir, sample.names[i], "_R1_paired.fastq' '", 
                          trimmed_dir, sample.names[i], "_R1_unpaired.fastq' '", trimmed_dir, sample.names[i], "_R2_paired.fastq' '", trimmed_dir, sample.names[i], 
                          "_R2_unpaired.fastq' LEADING:10 TRAILING:10 MINLEN:20 SLIDINGWINDOW:5:10 ILLUMINACLIP:'", 
                          trimmomatic_maindir, "/adapters/TruSeq3-PE.fa':2:30:12")
    

    bashCallFlash = paste0("'", flash_path, "'", 
                           " --min-overlap 30 --max-overlap 250 --max-mismatch-density 0.05 --output-directory=", 
                           "'", flash_input_dir, sample,
                           "' '", trimmed_dir, sample.names[i],
                           "_R1_paired.fastq' '", trimmed_dir, sample.names[i], "_R2_paired.fastq'")

    
    system(bashCallTrim, wait = TRUE)
    track_reads[map_file_sample[sample, 'sample'], 'trimmed'] = length(Biostrings::readDNAStringSet(paste0(trimmed_dir, sample.names[i], "_R1_paired.fastq"), format = 'fastq'))
    system(bashCallFlash, wait = TRUE)
  }
  
  asvs_df = list()
  for (i in c(1:length(sample.names)) ){
    cat(sample.names[i], "\n")
    
    tmp = Biostrings::readDNAStringSet(paste0(flash_input_dir, sample.names[i], "/out.extendedFrags.fastq"), format = 'fastq')
    tmp_df = data.frame(tmp, row.names = NULL)
    colnames(tmp_df) = 'seq'
    
    track_reads[map_file_sample[sample.names[i], 'sample'], 'merged'] = nrow(tmp_df)
    
    
    tmp_df =tmp_df  %>% group_by(seq) %>% summarise(count = n())
    tmp_df$sample = map_file_sample[sample.names[i], 'sample']
    
    cat( nrow(tmp_df), " reads remain after merging and collapsing\n")
    #originalSequences[[sample.names[i]]] = nrow(tmp_df)
    
    asvs_df[[map_file_sample[sample.names[i], 'sample']]] = tmp_df
    
  }
  
    df = bind_rows(asvs_df) %>% # filter(count > 2) %>%
      tidyr::pivot_wider(names_from = sample, values_from = count) %>%
      mutate_if(is.numeric , tidyr::replace_na, replace = 0)
    
    
  
  df$seq_names = paste0("SEQ", formatC(c(1:(nrow(df))), 
                                       width = nchar(trunc(nrow(df))), 
                                       format = "d", flag = "0")) # -1 to start 00 with no changes sequence
  
  df = df %>% filter(!is.na(seq) & !stringr::str_detect(seq, 'N'))
  
  rownames(df) = df$seq
  
  seq_lengths = data.frame(seq = df$seq,
                           length = nchar(df$seq), stringsAsFactors = F)
  hist_seq_count = 
    ggplot2::ggplot(data=seq_lengths, aes(x = length)) + 
    geom_histogram(binwidth = 1, fill = '#B484A9') +
    scale_x_continuous(labels=scales::comma, 
                       breaks=c(1, seq(26, 520, 26)),
                       limits=c(0, 520)) +
    xlab("Sequence Length") +
    ylab("Sequence count") +
    theme(panel.border=element_blank(), axis.line = element_line()) +
    #lemon::coord_capped_cart(left="both", bottom="both") +
    geom_vline(xintercept=260, linetype="dotted", size=0.25, col="#84B48F") +
    barplot_nowaklab_theme() 
  
  ggsave(filename=file.path(figure_dir, 'asv_length_freq.pdf'), 
         plot=hist_seq_count, 
         #device=grDevices::cairo_pdf, 
         width=25, height=5, units = "cm") #17.5 for 4x
  
  # Save as Data csv
  utils::write.csv(track_reads, file.path(fastq_dir, "fastq_summary.csv"), quote = FALSE)
  
  return(list(seqtab=df, track = track_reads))
}


# This function checks whether there is correspondence between input files and organ/day code provided.
# In case users do not provide a map between samples and organ/day code, then there can be two options:
# On the one hand, if users provide the fastq files input directory, then this function assumes that fastq filenames
# are in the correct format and prints out the samples and corresponding organ/day codes identified.

# When users provide a map between samples and organ/day code, this function checks that the sample names
# provided in the list correspond to either the filenames of fastqs.
check_input = function(sample.names = NULL,
                       map_file_sample = NULL){
  if (! is.null(map_file_sample) & class(map_file_sample) == 'list'){
    # make sure mapping object is in the correct format
    map_file_sample = as.data.frame(t(as.data.frame(map_file_sample)))
    colnames(map_file_sample) = 'sample'
  } 
  if (!is.null(map_file_sample)) {
    # Alignment needs to be performed or it was already provided and a mapping between samples and organ/day is given
    # Thus, need to check whether the samples given match the ones found in the mapping
    if (!sum(rownames(map_file_sample) %in% sample.names) != length(sample.names)) {
      cli::format_error("File names from the data do not match the ones found in mapping\n")
      cat(" The following file names are found: ", sample.names, "\n")
      cat("While mapping information is given for the following files: ", rownames(map_file_sample), "\n")
      stop("Mapping error")
    } 
  } else {
    # No mapping information was provided -> need to create the mapping
    sample_code = gsub("^(?:[^_]+_){1}([^_]+).*", "\\1", sample.names)
    map_file_sample = data.frame(sample = sample_code, 
                                 row.names = sample.names, stringsAsFactors = FALSE)
  }
  
  
  # At this point in the code, map_file_sample exists either because it was provided
  # as a parameter, or because it was not provided and created from sample names.
  cli::cli_alert_info('The following files, mapped to the corresponding sample, were found: ')
  lapply(rownames(map_file_sample),  function(x) {cat(x, " : ", map_file_sample[x, 'sample'], "\n")})
  return(map_file_sample)
}



# Adjust Data Frame
adjust_seqtab = function(seqtab.nochim, map_file_sample) {
  # Transpose to Get Sequences that Are Now Rows with removed chimeras
  seqtab_df = data.frame(t(seqtab.nochim))
  colnames(seqtab_df) = map_file_sample[colnames(seqtab_df),'sample']
  # create a column with seq from renames
  seqtab_df$seq = rownames(seqtab_df)
  # create rownames,  sequence variant: SEQ###
  seqtab_df$seq_names = paste0("SEQ", formatC(c(1:(dim(seqtab_df)[1])), width = nchar(trunc(dim(seqtab_df)[1])), format = "d", flag = "0")) # -1 to start 00 with no changes sequence
  rownames(seqtab_df) = seqtab_df$seq
  

  return(tibble::tibble(seqtab_df))
}


