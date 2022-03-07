# This function actually performs all dada2 steps.
dada2_alignment = function(fnFs,
                           fnRs,
                           output_dir,
                           map_file_sample,
                           sample.names,
                           output_dir_files,
                           output_dir_dada2 = NULL, 
                           output_figures = TRUE,
                           multithread = TRUE,
                           dada2_pooled_analysis = FALSE,
                           dada2_chimeras_minFoldParentOverAbundance = 8,
                           verbose = TRUE,
                           dada2_errorRate_estimation = 'default',
                           ...) {
  dots = list(...)
  # Make directories for the filtered fastqs and (optionally) for figures
  if (is.null(output_dir_dada2)) {
    filt_input_dir = file.path(output_dir, "filtered_fastq")
  } else {
    filt_input_dir = output_dir_dada2
  }
  
  if (!dir.exists(filt_input_dir)) dir.create(filt_input_dir, recursive = T)
  
  if (output_figures) {
    figure_dir = file.path(output_dir, "dada2_figures")
    if (!dir.exists(figure_dir)) dir.create(figure_dir)
  } else {
    figure_dir = NULL
  }
  
  
  filtFs = file.path(filt_input_dir, paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs = file.path(filt_input_dir, paste0(sample.names, "_R_filt.fastq.gz"))
  
  if (output_figures) {
    cli::cli_alert_info('Creating quality profiles figures')
    # Inspect Read Quality Profiles
    # Visualize the quality profile of the forward reads:
    ggplot2::ggsave(filename = file.path(figure_dir, "quality_profile_forwardReads.pdf"), 
                    plot = dada2::plotQualityProfile(fnFs[1:length(fnFs)]), 
                    #device = grDevices::cairo_pdf, 
                    width = 10*length(fnFs), 
                    height = 10,#length(fnFs)+2, 
                    units = "cm")
    
    # Visualize the quality profile of the reverse reads:
    ggplot2::ggsave(filename = file.path(figure_dir, "quality_profile_reverseReads.pdf"), 
                    plot = dada2::plotQualityProfile(fnRs[1:length(fnFs)]), 
                    #device = grDevices::cairo_pdf, 
                    width = 10*length(fnFs), 
                    height = 10,#length(fnFs)+2, 
                    units = "cm")
  }
  
  # dada2 filter and trim
  cli::cli_alert_info('Filtering and trimming')
  out = do.call(dada2::filterAndTrim, 
                              c(list(fwd = fnFs, filt = filtFs, 
                                     rev=fnRs, filt.rev = filtRs, 
                                     verbose=verbose),
                                     #truncQ=0, minLen = 2),
                                get_args_from_dots(dots, dada2::filterAndTrim)))
  
  # dada2 Learn errors
  cli::cli_alert_info('Learning errors for forward reads')
  
  error_estimation = ifelse(dada2_errorRate_estimation == 'default', 'loessErrfun', dada2_errorRate_estimation)
  
  errF = do.call(dada2::learnErrors, c(list(fls = filtFs, multithread = multithread, errorEstimationFunction = get(error_estimation)), 
                                                     get_args_from_dots(dots, dada2::learnErrors)))
  cli::cli_alert_info('Learning errors for backward reads')
  errR = do.call(dada2::learnErrors, c(list(fls = filtRs, multithread = multithread, errorEstimationFunction = get(error_estimation)), 
                                                     get_args_from_dots(dots, dada2::learnErrors)))
  
  if (output_figures) {
    # Visualize the estimated error rates of the forward reads:
    ggsave(filename = file.path(figure_dir, "quality_errors_errF.pdf"), 
           plot = dada2::plotErrors(errF, nominalQ=TRUE), 
           #device = grDevices::cairo_pdf, 
           width = 15,#5+5*length(sample.names), 
           height = 15,#5+5*length(sample.names), 
           units = "cm")
    # Visualize the estimated error rates of the reverse reads:
    ggsave(filename = file.path(figure_dir, "quality_errors_errR.pdf"), 
           plot = dada2::plotErrors(errR, nominalQ=TRUE), 
           #device = grDevices::cairo_pdf, 
           width=15,#5+5*length(sample.names), 
           height=15,#5+5*length(sample.names), 
           units = "cm")
  }
  
  # dada2 De-replication of the filtered fastq files:
  cli::cli_alert_info('Dereplicating fastqs')
  derepFs = do.call(dada2::derepFastq, c(list(fls = filtFs, verbose=verbose), 
                                                       get_args_from_dots(dots, dada2::derepFastq)))
  derepRs = do.call(dada2::derepFastq, c(list(fls = filtRs, verbose=verbose), 
                                                       get_args_from_dots(dots, dada2::derepFastq)))
  
  # Name the derep-class objects by the sample names
  if (length(sample.names) > 1) {
    names(derepFs) = sample.names
    names(derepRs) = sample.names 
  }
  
  # dada2 Sample Inference
  cli::cli_alert_info('Inferring samples')
  dadaFs = do.call(dada2::dada, 
                                 c(list(derep = derepFs, err = errF, multithread = multithread, pool = dada2_pooled_analysis), 
                                   get_args_from_dots(dots, dada2::dada)))
  # (1) pool=FALSE; sample-by-sample inference; (2) pool=TRUE if conected samples i.e. mouse organs, days follow-ups 
  dadaRs = do.call(dada2::dada, 
                                 c(list(derep = derepRs, err = errR, multithread = multithread, pool = dada2_pooled_analysis), 
                                   get_args_from_dots(dots, dada2::dada)))
  
  # dada2 Merge Paired Reads
  mergers = do.call(dada2::mergePairs, 
                                  c(list(dadaF = dadaFs, derepF = derepFs, 
                                         dadaR = dadaRs, derepR = derepRs, 
                                         verbose=verbose),
                                    get_args_from_dots(dots, dada2::mergePairs)))
  

  # dada2 Construct sequence table for Forward and Revers Reads
  if (length(sample.names) == 1) {
    mergers = list(mergers)
    names(mergers) = sample.names
  }
  seqtab = do.call(dada2::makeSequenceTable, 
                                 c(list(samples = mergers), 
                                   get_args_from_dots(dots, dada2::makeSequenceTable)))
  
  if (output_figures) {
    seq_lengths = data.frame(seq = dada2::getSequences(seqtab),
                             length = nchar(dada2::getSequences(seqtab)))
    # freq_lengths = plyr::count(seq_lengths, 'length')
    # rownames(freq_lengths) = freq_lengths$length
    # freq_lengths = freq_lengths[,-1, drop=F]
    # pdf('prova.pdf', width = nrow(freq_lengths) / 2, height = 3)
    # pheatmap(t(freq_lengths), cluster_rows = F, cluster_cols = F,
    #          fontsize = 20,
    #          display_numbers = t(freq_lengths), 
    #          labels_row = '',
    #          main = 'Frequency of sequence length',
    #          angle_col = 45)
    # dev.off()
    
    hist_seq_count = 
      ggplot(data=seq_lengths, aes(x = length)) + 
      geom_histogram(binwidth = 1, fill = '#B484A9') +
      scale_x_continuous(labels=scales::comma, 
                         breaks=c(1, seq(52, 520, 26)),
                         limits=c(0, 520)) +
      xlab("Sequence Length") +
      ylab("Sequence count") +
      lemon::coord_capped_cart(left="both", bottom="both") +
      geom_vline(xintercept=260, linetype="dotted", size=0.25, col="#84B48F") +
      barplot_nowaklab_theme() 
    
    ggsave(filename=file.path(figure_dir, 'sequence_length.pdf'), 
           plot=hist_seq_count, 
           #device=grDevices::cairo_pdf, 
           width=25, height=5, units = "cm") #17.5 for 4x
  }
  
  # dada2 Remove Chimeric Sequences 
  cli::cli_alert_info('Removing bimeras')
  seqtab.nochim = do.call(dada2::removeBimeraDenovo, c(list(unqs = seqtab, 
                                                            minFoldParentOverAbundance = dada2_chimeras_minFoldParentOverAbundance, #changed wrt default 2
                                                            multithread = multithread, 
                                                            verbose=verbose),
                                                       get_args_from_dots(dots, dada2::removeBimeraDenovo)
  ))
  
  # seqtab vs 
  bimera_perc = sum(seqtab.nochim)/sum(seqtab)*100
  
  getN = function(x) sum(getUniques(x))
  if (class(dadaFs) != 'list') {
    # Only one sample 
    track = cbind(out, getN(dadaFs), getN(dadaRs), getN(mergers[[1]]), rowSums(seqtab.nochim))
  } else {
    track = cbind(out,
                  sapply(dadaFs, getN),
                  sapply(dadaRs, getN),
                  sapply(mergers, getN),
                  rowSums(seqtab.nochim))
  }
  
  colnames(track) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) = sample.names
  
  # change rows names: Names of interest between 1st "_" and 2nd "_" == {1} "\\1"
  rownames(track) = map_file_sample[rownames(track),,drop=F]$sample

  # Save as Data csv
  utils::write.csv(track, file.path(output_dir_files, "quality_track_reads.csv"))
  
  return(list(seqtab.nochim=seqtab.nochim, track = track, bimera_perc = bimera_perc, 
              nSequences_with_chimeras = dim(seqtab)[2]))
}

# This function checks whether there is correspondence between input files and organ/day code provided.
# In case users do not provide a map between samples and organ/day code, then there can be two options:
# On the one hand, if users provide the fastq files input directory, then this function assumes that fastq filenames
# are in the correct format and prints out the samples and corresponding organ/day codes identified.
# On the other hand, in case users provide the path to dada2 output, it is assumed that the file contains a dataframe
# where row names have the same format as the one expected for fastq files (i.e. SAMPLE_ORGAN|DAY_BARCODEVERSION_XXXXX).
# In both this two cases, this function will print the samples and organ/day codes idenfied.
# 
# When users provide a map between samples and organ/day code, this function checks that the sample names
# provided in the list correspond to either the rownames of dada2 output (when this is passed directly by the user) or
# to the filenames of fastqs.
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
adjust_seqtab = function(seqtab.nochim, map_file_sample, output_dir_files) {
  # Transpose to Get Sequences that Are Now Rows with removed chimeras
  seqtab_df = data.frame(t(seqtab.nochim))
  colnames(seqtab_df) = map_file_sample[colnames(seqtab_df),'sample']
  # create a column with seq from renames
  seqtab_df$seq = rownames(seqtab_df)
  # create rownames,  sequence variant: SEQ###
  seqtab_df$seq_names = paste0("SEQ", formatC(c(1:(dim(seqtab_df)[1])), width = nchar(trunc(dim(seqtab_df)[1])), format = "d", flag = "0")) # -1 to start 00 with no changes sequence
  rownames(seqtab_df) = seqtab_df$seq
  
  # save as data csv
  utils::write.csv(seqtab_df, file.path(output_dir_files, "dada2_asv_prefilter.csv"),
            row.names = FALSE)
  return(tibble::tibble(seqtab_df))
}
