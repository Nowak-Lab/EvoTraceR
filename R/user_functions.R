#' This function loads reticulate and needs to be run after the installation of the package
#' 
#' @title setup_reticulate
#' 
#' @export setup_reticulate
#' @import reticulate
setup_reticulate = function() {
  reticulate::py_install("git+https://github.com/YosefLab/Cassiopeia@master#egg=cassiopeia-lineage", pip = TRUE) 
  reticulate::py_install('umi_tools', pip = TRUE) 
}

#' This function initializes the EvoTraceR object, by computing the set of Amplicon Sequence Variants.
#' It calls \href{http://www.usadellab.org/cms/?page=trimmomatic}{trimmomatic} and \href{https://ccb.jhu.edu/software/FLASH/}{flash} to perform adapters trimming, discard low quality bases and merge
#' forward and reverse reads.
#' 
#' @title initialize_EvoTraceR
#' 
#' @examples 
#' \dontrun{
#' input_dir = system.file("extdata", "input", package = "EvoTraceR")
#' output_dir = system.file("extdata", "output", package = "EvoTraceR")
#' EvoTraceR_object = initialize_EvoTraceR(input_dir = input_dir, output_dir = output_dir, trimmomatic_path = '/path/to/trimmomatic.jar', flash_path = '/path/to/flash_bin_directory')
#' }
#' 
#' @param output_dir (Required). Path to the directory where all output files will be stored. 
#' This function will output  the file \code{quality_track_reads.csv}, containing a track of the number of sequences during the 
#' different steps.
#' 
#' @param input_dir Path to the directory containing \code{.fastq} files for forward and reverse reads.
#' This folder should contain the fastq files (2 for each sample) with the following name pattern:
#' FILEPREFIX_SAMPLE_BARCODEVERSION_R1.fastq FILEPREFIX_SAMPLE_BARCODEVERSION_R1.fastq. SAMPLE refers to either an organ (in case multiple organs were sequenced)
#' or timepoint (if longitudinal data are provided). Note that EvoTraceR does not support mixed sample types (i.e. samples must be either all from organs or all from timepoints).
#' @param trimmomatic_path (Required). Local path to the executable of \code{Trimmomatic}. 
#' For details about the download please see \href{http://www.usadellab.org/cms/?page=trimmomatic}{here}.
#' @param flash_path (Required). Local path to the executable of \code{Flash}.
#' For details about the download please see \href{https://ccb.jhu.edu/software/FLASH/}{here}.
#' @param map_file_sample (Optional). In case fastq files names are not in the format \cr FILEPREFIX_SAMPLE_FILESUFFIX_RX.fastq),
#' then users should provide a list that associates each filename (without the suffix _R1 and _R2) to the corresponding organ/day.
#' (e.g., if we have forward and reverse files named file1_R1.fastq and file1_R2.fastq that correspond to organ PRL 
#' (code for Prostate Left), than the parameter should be set as: \code{map_file_sample = c("file1" = "PRL")}).
#' @param sample_order (Optional). Vector containing the order in which the user wants samples to appear in all plots. If \code{NULL} the alphabetical order will be set.
#' 
#' @return An object of type EvoTraceR, which is a list that will contain the following fields: 
#' \itemize{
#' \item \code{fastq_directory}: directory where the input fastq files are located.
#' \item \code{output_directory}: directory where all the output files are being stored.
#' \item \code{map_file_sample}: dataframe has as many rows as the input datasets, and for each input stores the sample (e.g. organ or day for longitudinal data)
#' to which it is associated.
#' \item \code{asv_prefilter}: dataframe that stores all sequences detected after these preliminary steps. Note that
#' these sequences still need to be filtered (see also \code{\link{asv_analysis}}).
#' }. 
#' This function also saves the \code{.csv} file \code{quality_track_reads.csv}: track of the number of sequences during \code{trimming and merging} steps.
#' 
#' @export initialize_EvoTraceR
#' 
#' @rawNamespace import(ggplot2, except = c(element_render, CoordCartesian))
#' @rawNamespace import(dplyr, except = count)
#' @importFrom cli cli_alert_info
#' @importFrom utils write.csv read.csv unzip untar
#' @importFrom stringr str_remove_all str_replace
initialize_EvoTraceR = function(output_dir,
                                trimmomatic_path,
                                flash_path,
                                input_dir = NULL,
                                map_file_sample = NULL,
                                sample_order = 'alphabetical') {
  # Check that the user inserted the correct parameters
  if (is.null(input_dir) ) {
    stop('Please provide a directory containing fastqs')
  }
  if (is.null(output_dir)) {
    stop('Please provide a path where output files will be stored.')
  }
  if (is.null(trimmomatic_path)) {
    stop('Please provide a valid path to the jar for trimmomatic')
  } else if (!file.exists(file.path(trimmomatic_path))) {
    stop('Please provide a valid path to the jar for trimmomatic, file not found')
  }
  
  if (is.null(flash_path)) {
    stop('Please provide a valid path to flash')
  } else if (!file.exists(file.path(trimmomatic_path))) {
    stop('Please provide a valid path to flash, file not found')
  }
  
  EvoTraceR_object = list(fastq_directory = input_dir, output_directory = output_dir)
  
  class(EvoTraceR_object) = 'EvoTraceR'
  EvoTraceR_object$preprocessing = list()
  
  # Check if input files are compressed
  zip_files = list.files(input_dir, pattern = '.zip$', full.names = TRUE)
  tar_files = list.files(input_dir, pattern = '.tar$', full.names = TRUE)
  
  if (length(zip_files) > 0) {
    cli::cli_alert_info("Found zipped fastq files, extracting")
    for (zf in zip_files) {
      utils::unzip(zf, exdir = stringr::str_replace(input_dir, pattern = "/$", replacement=''))
    }
    cli::cli_alert_info("Done extracting")
  } else if (length(tar_files) > 0) {
    cli::cli_alert_info("Found tar fastq files, extracting")
    for (tf in tar_files) {
      utils::untar(tf, exdir = stringr::str_replace(input_dir, pattern = "/$", replacement=''))
    }
    cli::cli_alert_info("Done extracting")
  }
  fastqs = list.files(input_dir, pattern = '.fastq$')
  if (length(fastqs) == 0) {
    stop('No fastq files found in the input directory provided, stopping.')
  } else {
    cat("Found", length(fastqs), "fastq files\n")
  }
  fastqs = sort(fastqs) 
  fnFs = fastqs[grepl("_R1", fastqs)] 
  fnRs = fastqs[grepl("_R2", fastqs)] 
  # Get sample names, assuming files named as so: SAMPLENAME_XXX.fastq
  sample.names = stringr::str_remove_all(fnFs, "_R1.fastq")
  map_file_sample = check_input(sample.names = sample.names, 
                                map_file_sample = map_file_sample)
  # Specify the full path to the fnFs and fnRs
  fnFs = file.path(input_dir, fnFs)
  fnRs = file.path(input_dir, fnRs)
  
  align_output = alignment_pipeline(fnFs = fnFs,
                                    fnRs = fnRs,
                                    map_file_sample = map_file_sample,
                                    sample.names = sample.names,
                                    output_dir = output_dir,
                                    trimmomatic_path = trimmomatic_path,
                                    flash_path = flash_path
                                    )
  
  #seqtab.nochim = align_output$seqtab.nochim
  EvoTraceR_object$preprocessing$track = align_output$track
  
  
  EvoTraceR_object$map_file_sample = map_file_sample
  EvoTraceR_object$asv_prefilter = align_output$seqtab
  
  if (is.character(sample_order) & length(sample_order) == 1) {
    if (sample_order != 'alphabetical') {
      stop('The ordering provided is not valid. Either use the default alphabetical or give a list of sample names.')
    }
    sample_order = sort(unique(map_file_sample$sample))
    if ('PRL' %in% sample_order)
      sample_order = c('PRL', setdiff(sample_order, 'PRL'))
  } else {
    if (length(intersect(sample_order, EvoTraceR_object$map_file_sample$sample)) != length(sample_order)) {
      cli::cli_alert_danger('The samples found in variable sample_order do not match the samples contained in the dataset.')
      cat("The following samples need to be provided: ", unique(EvoTraceR_object$map_file_sample$sample))
      stop("Exiting")
    }
  }
  EvoTraceR_object$sample_order = sample_order
  return(EvoTraceR_object)
  
}

#' This function performs the analysis on ASV sequences identified by the previous steps and aligns them to the reference sequence.
#' First, it pools together those sequences characterized by a Hamming distance equals or lower than 2, summing their counts. To perform this step, EvoTraceR uses the
#' sequence clusering agorithm impemented in the python package \href{https://umi-tools.readthedocs.io/en/latest/QUICK_START.html}{UMI-tools}
#' Then it performs pairwise alignment using Needleman-Wunsch global alignment algorithm implemented in function \code{pairwiseAlignment}
#' in package \code{Biostrings}, aligning each sequence to the original barcode considered in the analysis.
#' (See the Biostrings documentation \href{https://www.rdocumentation.org/packages/Biostrings/versions/2.40.2/topics/pairwiseAlignment}{here} for more details).
#' After identifying all indels (insertions are identified by their start position and number of nucleoides inserted, while deletions are identified by their start and end position),
#' it removes the ones that are too small and don't span any cut site: to perform this filter, it exapnds the start and end position by the number of bases specified in the parameter \code{cleaning_window},
#' it counts the number of cut sites spanned after the expantion and it removes those that span 0 sites.
#' Next, it pools together those sequences that differ by one another only by substitutions, summing their counts, and then it removes
#' those sequences showing a frequency lower than the threshold specified in parameter \code{asv_count_cutoff}.
#' Next, it performs filtering of the sequences based on their flanking sequences and finally it computes the normalized counts of each ASV in each sample,
#' dividing each count by the total number of sequences for each sample and multiplying by 1e6, yelding Counts Per Million (CPM).
#'
#' Then, using the CPM it computes different statistics for the final ASVS, storing:
#'  the relative frequency of all ASVs in each sample. 
#'  the relative frequency of each ASV in the samples.
#'  the counts for each ASV normalized to the counts of the sample with maximum frequency.
#'  the frequency of the different ASVs in each sample.
#'  
#' @title ASV_analysis
#' 
#' @examples
#' \dontrun{
#' data(EvoTraceR_object)
#' EvoTraceR_object = asv_analysis(EvoTraceR_object = EvoTraceR_object)
#' }
#'
#' @param EvoTraceR_object (Required). Object of class EvoTraceR, result of the function \code{initialize_EvoTraceR}
#' @param ref_name String indicating the ID of the reference sequence used in the experiment. Default is 'BC10v0',
#' @param ref_seq String indicating the reference sequence used in the experimenti. Default is 'TCTACACGCGCGTTCAACCGAGGAAAACTACACACACGTTCAACCACGGTTTTTTACACACGCATTCAACCACGGACTGCTACACACGCACTCAACCGTGGATATTTACATACTCGTTCAACCGTGGATTGTTACACCCGCGTTCAACCAGGGTCAGATACACCCACGTTCAACCGTGGTACTATACTCGGGCATTCAACCGCGGCTTTCTGCACACGCCTACAACCGCGGAACTATACACGTGCATTCACCCGTGGATC',
#' @param ref_flank_left String indicating the first nucleotides of the reference sequence that never mutate over the
#' course of the experiment. Default is "^TCTAC",
#' @param ref_flank_right String indicating the first nucleotides of the reference sequence that never mutate over the
#' course of the experiment. Default is "CCCGTGGATC$",
#' @param flanking_filtering Which among the flaning regions to use to filter out contaminated sequences.
#' Must be one of c('left', 'right', 'both', 'either'), default is 'both'
#' @param ref_cut_sites Positions in the reference sequence of the cutting sites. Default is c(17, 42, 68, 94, 120, 146, 171, 198, 224, 251),
#' @param ref_border_sites c(26, 52, 78, 104, 130, 156, 182, 208, 234).
#' @param output_figures (Optional). Default TRUE: Boolean indicating whether a user whishes to store a figure indicating the number of ASV tracked during the different steps of the analysis.
#' @param asv_count_cutoff (Optional). Default to 2. Minimum number of counts for an ASV to be considered in the statistics.
#' @param pwa_gapOpening (Optional). Default is -25. Parameter \code{gapOpening} passed to \code{pairwiseAlignment} from \code{Biostrings} (See description).
#' @param pwa_gapExtension (Optional). Default is 0. Parameter \code{gapExtension} passed to \code{pairwiseAlignment} from \code{Biostrings} (See description). Default is 0.
#' @param pwa_mismatch (Optional). Default is -4. Parameter indicating the penalty for mismatch events during pairwise alignment with \code{Biostrings}.
#' @param pwa_match (Optional). Default is 15. Parameter indicating the score for matches during pairwise alignment with \code{Biostrings}. This parameter,
#' together with the previous one, are used to construct the substitution matrix used by the function \code{pairwiseAlignment}.
#' @param pwa_type (Optional). Parameter indicating the type of pairwise alignment. Must be one of One of "global", "local", "overlap", "global-local", and "local-global".
#' For more details see \href{https://www.rdocumentation.org/packages/Biostrings/versions/2.40.2/topics/pairwiseAlignment}{original documentation} 
#' @param cleaning_window (Optional). Default is c(3,3). Vector containing the number of nucleotides that we use for extending respectively the start and end position of each indel to determine the ones that don't span any cut sites and thus get removed. (See description for more information).
#' @param batch_size (Optional). Default is 100. Number of batches to split reads into for parallel execution of \code{pairwiseAlignment}. This parameter can be tuned together with the cores parameter to optimize the speed of alignment.
#' @param cores (Optional). Default is parallel::detectCores(). Number of cores to use for pairwise alignment.
#'
#' @return  The EvoTraceR object passed as a parameter with the following new fields:
#' \itemize{
#' \item \code{clean_asv_dataframe}: ASV sequences identified post-filtering (contamination removed,
#' sequences with a similarity higher than \code{pid_cutoff_nmbc} to the original barcode
#' aggregated to it and ASVs named in increasing order (ASV01, ASV02, etc.) according
#' to their total counts. 
#' \item \code{reference}: Info about the reference sequence used for the current analysis.
#'
#'  \item \code{statistics}: another list with the following sub-fields: 
#' \itemize{
#' 
#' \item \code{asv_df_percentages}: dataframe with six columns. \code{asv_names} is the name of the ASV.
#' \code{sample} is the sample identifier (e.g. ID of an organ or, in case of longitudinal data, of the timepoint);
#' \code{count}: total counts per million for a specific ASV in a specific sample;
#' \code{perc_in_sample}: CPM normalized to the total counts in the corresponding sample;
#' \code{perc_asv}: CPM normalized to the total counts for the corresponding ASV;
#' \code{perc_fold_to_max}: CPM normalized to the maximum counts observed for the corresponding ASV in a sample.
#' 
#' \item \code{asv_totalCounts}: for each ASV, total counts and number of samples in which it was detected.
#' \item \code{sample_totalcounts}: for each sample, total counts and number of distinct ASVs detected.
#' \item \code{asv_diversity_persample}: measures of clonal richness and measures of heterogeneity computed for each sample based on the ASVs detected.
#' \item \code{asv_persample_frequency}: counts for each ASV in each sample.
#' \item \code{asv_persample_detection}: binary matrix indicating whether a sequence has been detected in the corresponding sample.
#' \item \code{asv_toBarcode_similarity}: edit distance, percentage similarity and alignment score of each ASV compared to the original barcode.
#' \item \code{all_asv_statistics}: all the statistics computed on each ASV grouped together in the same tibble.
#' }
#' \item \code{alignment}, another list with the following fields:
#' \itemize{
#' \item \code{Binary_mutation_matrix}: binary matrix encoding the presence/absence of a mutation in an ASV.
#' \item \code{asv_barcode_alignment}: tibble where each line corresponds to a position in a ASV, and the columns encode the following information:
#' \itemize{
#' \item asv_names: name of the ASV
#' \item sample: sample identifier
#' \item position_bc260: position of the alteration in the original barcode. Note that insertions
#' are assigned to the position that coincides with their beginning.
#' \item alt: type of alteration. wt = Wild Type (i.e. non-mutated position). sub = substitution. del = deletion. ins = insertion.
#' \item{ref_asv}, \item{read_asv}: respectively, the reference nucleotide observed in the original barcode and the one observed on the sequence.
#' }
#' \item \code{mutations_coordinates}: dataframe containing a list of all mutations, with their start and end position.
#' \item \code{ASV_alterations_width}: dataframe containing the number of nucleotides affected by each type of mutation in each ASV.
#' }
#' }
#' This function saves the figure \code{asv_filtering_freq.pdf}, which is a barplot indicating the change in the number of sequences 
#' during all the processing steps.
#'   
#' @export asv_analysis
#'
#' @importFrom Biostrings pairwiseAlignment pid nedit score nucleotideSubstitutionMatrix DNAStringSet
#' @importFrom benthos total_abundance species_richness margalef rygg simpson hpie hill1 hill2 shannon
# @importFrom lemon coord_capped_cart facet_rep_grid
#' @importFrom scales comma
#' @importFrom utils write.csv
#' @importFrom stringr str_detect
#' @importFrom tidyr gather pivot_wider
#' @importFrom forcats fct_relevel
#' @importFrom data.table nafill
# @importFrom grDevices cairo_pdf
#' @rawNamespace import(dplyr, except = count)
#' @rawNamespace import(ggplot2, except = c(element_render, CoordCartesian))
#' @import tibble
#' @import foreach
#' @importFrom purrr map
#' @import parallel
#' @importFrom doParallel registerDoParallel 
asv_analysis = function(EvoTraceR_object,
                        #barcode = 'BC10v0',
                        ref_name = 'BC10v0',
                        ref_seq = 'TCTACACGCGCGTTCAACCGAGGAAAACTACACACACGTTCAACCACGGTTTTTTACACACGCATTCAACCACGGACTGCTACACACGCACTCAACCGTGGATATTTACATACTCGTTCAACCGTGGATTGTTACACCCGCGTTCAACCAGGGTCAGATACACCCACGTTCAACCGTGGTACTATACTCGGGCATTCAACCGCGGCTTTCTGCACACGCCTACAACCGCGGAACTATACACGTGCATTCACCCGTGGATC',
                        ref_flank_left = "^TCTAC",
                        ref_flank_right = "CCCGTGGATC$",
                        flanking_filtering = 'right',
                        ref_cut_sites = c(17, 43, 69, 95, 121, 147, 173, 199, 225, 251), # cut sites for Cas9 based on "ref_seq" 
                        ref_border_sites = c(1, 26, 52, 78, 104, 130, 156, 182, 208, 234),
                        output_figures = TRUE,
                        asv_count_cutoff = 2,
                        pwa_gapOpening = -25,
                        pwa_gapExtension = 0,
                        pwa_match = 15,
                        pwa_mismatch = -4,
                        pwa_type = 'global',
                        cleaning_window = c(3,3),
                        batch_size = 100,
			cores = parallel::detectCores()
) {
  
  barcodes_info = list(
    ref_name = ref_name,
    ref_seq = ref_seq,
    ref_flank_left = ref_flank_left, # 5x nts
    ref_flank_right = ref_flank_right,
    ref_cut_sites = ref_cut_sites,
    ref_border_sites = ref_border_sites,
    stringsAsFactors = FALSE)
  
  seqtab_df = EvoTraceR_object$asv_prefilter
  
  EvoTraceR_object$reference = barcodes_info
  
  # Store the original number of sequences
  orgseq <- nrow(seqtab_df)
  seqtab_df_original = seqtab_df
  
  # get organ list
  sample_columns = setdiff(colnames(seqtab_df), c("seq_names", "seq"))
  
  cli::cli_alert_info('Merging sequences based on Hamming distance')
  seqtab_df = merge_hamming(seqtab_df, sample_columns, cores)
  
  hamming_filter = nrow(seqtab_df)

  # for single organ, rowSums is not needed and will throw error
  if (length(sample_columns) > 1) {
    seqtab_df$totalCounts = rowSums(seqtab_df[,sample_columns])
  } else {
    seqtab_df$totalCounts = seqtab_df[,sample_columns]
  }
  #counts_filtering = nrow(seqtab_df)
  
  mx_crispr <- Biostrings::nucleotideSubstitutionMatrix(match = pwa_match, mismatch = pwa_mismatch, baseOnly = TRUE)
  ###### COLLAPSING
  # Collapse sequences that differ only by substitutions.
  collapse_result = asv_collapsing(seqtab_df, 
                                   barcodes_info$ref_seq,
                                   pwa_match,
                                   pwa_mismatch,
                                   pwa_gapOpening,
                                   pwa_gapExtension, 
                                   sample_columns,
                                   barcodes_info$ref_name,
                                   cut_sites = barcodes_info$ref_cut_sites,
                                   cleaning_window,
                                   batch_size,
				   cores)
  seqtab_df = collapse_result$seqtab_df
  tidy_alignment = collapse_result$tidy_alignment
  endseq_filter <- nrow(seqtab_df)
  
  cleaned_coordinate_matrix = collapse_result$mutations_coordinates
  binary_mutation_matrix = collapse_result$binary_matrix
  
  EvoTraceR_object$seqtab_dataframe_nonnorm = seqtab_df
  norm_seqtab_df = seqtab_df
  norm_seqtab_df[,sample_columns] = sapply(sweep(norm_seqtab_df[, sample_columns], 2, 
                                                EvoTraceR_object$preprocessing$track[sample_columns,'input'], '/') * 1e6, as.integer)  

  EvoTraceR_object$seqtab_dataframe_norm = norm_seqtab_df

  norm_seqtab_df = norm_seqtab_df %>% rowwise() %>% 
    mutate_if(is.numeric, function (x) if (x <= asv_count_cutoff) return(as.integer(0)) else (return(x))) %>%
    ungroup()
  norm_seqtab_df = norm_seqtab_df[rowSums(norm_seqtab_df[,sample_columns]) > 0,]
  norm_seqtab_df$totalCounts = rowSums(norm_seqtab_df[,sample_columns])
  counts_filtering = nrow(norm_seqtab_df)
  
  norm_seqtab_df = perform_flanking_filtering(barcodes_info = barcodes_info, norm_seqtab_df = norm_seqtab_df, flanking_filtering = flanking_filtering)
  tidy_alignment = tidy_alignment %>% filter(seq_names %in% norm_seqtab_df$seq_names)
  cleaned_coordinate_matrix = cleaned_coordinate_matrix %>% filter(seq_names %in% norm_seqtab_df$seq_names)
  binary_mutation_matrix = binary_mutation_matrix %>% filter(seq_names %in% norm_seqtab_df$seq_names) #[norm_seqtab_df$seq_names,]
  binary_mutation_matrix = binary_mutation_matrix %>% select_if(~ !is.numeric(.) || sum(.) != 0)
  flanking_filtering = nrow(norm_seqtab_df)
  
  # Replace the seq-name for those sequences that match exactly one of the original barcodes.
  # In case no original barcode is found then insert it with 0 counts
  if (sum(norm_seqtab_df$seq == barcodes_info$ref_seq) == 0){
    
    norm_seqtab_df = norm_seqtab_df %>% dplyr::add_row(seq_names = barcodes_info$ref_name, 
                                             seq = barcodes_info$ref_seq)
    norm_seqtab_df[norm_seqtab_df$seq_names == barcodes_info$ref_name, sample_columns] = 0 
    
  } else {
    barcode_seqname = as.character(norm_seqtab_df[norm_seqtab_df$seq == barcodes_info$ref_seq, "seq_names"])
    norm_seqtab_df[norm_seqtab_df$seq_names == barcode_seqname, "seq_names"] = paste0(barcodes_info$ref_name, '')#".NMBC")
    tidy_alignment[tidy_alignment$seq_names == barcode_seqname, 'seq_names'] = barcodes_info$ref_name
  }
  
  
  # Now sort ASV by total frequency and assign an ASV ID
  norm_seqtab_df_clean_asv <-
    tibble(norm_seqtab_df) %>% mutate(asv_total_freq = rowSums(across(where(is.numeric)))) %>% arrange(-asv_total_freq) %>%
    mutate(asv_names = seq_names)
  # Find Row for NMBC
  norm_seqtab_df_clean_nmbc <- norm_seqtab_df_clean_asv %>% filter(asv_names == barcodes_info$ref_name)  
  
  # skip NMBC
  norm_seqtab_df_clean_asv <-
    norm_seqtab_df_clean_asv %>% filter(seq_names != barcodes_info$ref_name)

  # create ASV count
  norm_seqtab_df_clean_asv$asv_names <- paste0("ASV", 
                                          formatC(c(1:nrow(norm_seqtab_df_clean_asv)), 
                                                  width = nchar(trunc(nrow(norm_seqtab_df_clean_asv))), 
                                                  format = "d", flag = "0")) # -1 to start 00 with no changes sequence
  
  # add back row with "norm_seqtab_df_perf_match"  
  norm_seqtab_df_clean_asv <-
    norm_seqtab_df_clean_asv %>%
    add_row(norm_seqtab_df_clean_nmbc) %>%
    arrange(-asv_total_freq) #%>%

  # Recompute tidy alignment matrix, mutations coordinates and binary mutation matrix.
  tidy_alignment_clean = tidy_alignment %>% ungroup() %>% left_join(norm_seqtab_df_clean_asv %>% select(seq_names, asv_names), by = 'seq_names') %>% 
    select(-c(seq_names))
  cleaned_coordinate_matrix = cleaned_coordinate_matrix  %>% left_join(norm_seqtab_df_clean_asv %>% select(seq_names, asv_names), by = 'seq_names') %>%
    select(-c(seq_names))
  binary_mutation_matrix = binary_mutation_matrix %>% left_join(norm_seqtab_df_clean_asv %>% select(seq_names, asv_names), by = 'seq_names') %>%
    tibble::column_to_rownames("asv_names") %>% select(-c(seq_names))
  
  bc_mut = as.list(rep(0, ncol(binary_mutation_matrix)))
  names(bc_mut) = colnames(binary_mutation_matrix)
  binary_mutation_matrix = dplyr::bind_rows(data.frame(bc_mut, row.names = barcodes_info$ref_name),
                                            binary_mutation_matrix)
  
  norm_seqtab_df_clean_asv = norm_seqtab_df_clean_asv %>% dplyr::select(-c("seq_names", "asv_total_freq")) %>% relocate(asv_names)
  # final number of ASVs for analysis
  clean_asv <- nrow(norm_seqtab_df_clean_asv)
  
  EvoTraceR_object$alignment$binary_mutation_matrix = binary_mutation_matrix
  
  cleaned_coordinate_matrix <- tibble::tibble(cleaned_coordinate_matrix) %>%
    dplyr::add_row(asv_names  = EvoTraceR_object$reference$ref_name, mutation_type = 'w', n_nucleotides = 0)
  
  EvoTraceR_object$alignment$asv_barcode_alignment = tidy_alignment_clean
  EvoTraceR_object$alignment$mutations_coordinates = cleaned_coordinate_matrix
  
  
  #norm_seqtab_df_clean_asv = seqtab_df_clean_asv
  #EvoTraceR_object$clean_asv_dataframe_nonnorm = seqtab_df_clean_asv
  
  #norm_seqtab_df_clean_asv[,sample_columns] = sapply(sweep(norm_seqtab_df_clean_asv[, sample_columns], 2, 
   #                                                        EvoTraceR_object$preprocessing$track[sample_columns,'input'], '/') * 1e6, as.integer)  
  
  
  #EvoTraceR_object$clean_asv_dataframe = norm_seqtab_df_clean_asv
  
  EvoTraceR_object$clean_asv_dataframe = norm_seqtab_df_clean_asv

  EvoTraceR_object = asv_statistics(EvoTraceR_object, 
                                    sample_columns, 
                                    asv_count_cutoff,
                                    nmbc = barcodes_info$ref_name)
  
  EvoTraceR_object$preprocessing$seq_filters = data.frame(name=as.factor(c("Starting ASVs", "Hamming Merging", "Substitutions Merging",  "Frequency Filter", "Flanking Seq. Filter", "Final ASVs")), 
                                                  num=c(orgseq, hamming_filter, endseq_filter, counts_filtering, flanking_filtering, clean_asv))
  
  seq_filtering_plot(EvoTraceR_object)

  return(EvoTraceR_object)
  
}

#' This function considers the pairwise alignment output performed with \code{asv_analysis} and computes the overall frequency in each 
#' sample of each mutation in every barcode position.
#' 
#' @title analyse_mutations
#' 
#' @examples 
#' \dontrun{
#' data(EvoTraceR_analyzed)
#' EvoTraceR_object = analyse_mutations(EvoTraceR_object)
#' }
#' 
#' @param EvoTraceR_object EvoTraceR object.
#' 
#' @return EvoTraceR object with the field \code{alignment}, updated with 
#' \code{mutations_df}: dataframe containing for each position in each ASV the corresponding mutation state and the frequency of that mutation in all sequences.
#'
#' @export analyse_mutations
#' @rawNamespace import(dplyr, except = count)
#' @rawNamespace import(ggplot2, except = c(element_render, CoordCartesian))
#' @import tibble
#' @rawNamespace import(data.table, except = c(last, first, between))
# @importFrom muscle muscle
#' @importFrom Biostrings DNAStringSet writeXStringSet DNAMultipleAlignment
# @importFrom ggmsa tidy_msa 
# @importFrom lemon coord_capped_cart facet_rep_grid
#' @importFrom utils write.csv
#' @importFrom forcats fct_relevel
# @importFrom grDevices cairo_pdf
#' @importFrom stringr str_replace str_detect
#' @importFrom methods as
#' @importFrom tidyr pivot_wider
#' @importFrom pheatmap pheatmap 
#' @importFrom plyr count
analyse_mutations = function(EvoTraceR_object) {#, cleaning_window = c(3,3)) {
  
  EvoTraceR_object = count_alterations(EvoTraceR_object)
  
  return(EvoTraceR_object)
  
}

#' This function reconstructs the phylogenetic tree using the \href{https://github.com/YosefLab/Cassiopeia}{Cassiopeia} suite. 
#' It uses the greedy algorithm, that iteratively splits the set of ASVs in clusters
#' based on the most common mutation at each iteration. Note that this requires installing cassiopeia in the anaconda/virtual enviroment used by reticulate.
#' 
#' @title infer_phylogeny
#' 
#' @examples 
#' \dontrun{
#' EvoTraceR_object = infer_phylogeny(EvoTraceR_object)
#' }
#' 
#' @param EvoTraceR_object (Required)
#' @param mutations_use (optional). Default = 'del_ins' A string indicating what type of mutations to use for phylogeny reconstruction.
#' Can be one of c('del', 'del_ins'). The first
#' corresponds to using only cleaned deletions. The second refers to the case where cleaned deletions
#' and insertions are considered.
#' 
#' @return EvoTraceR object with a new field called phylogeny, that stores:
#'  \itemize{
#'  \item \code{mutations_in_phylogeny}: string indicating which mutations were used for phylogeny reconstruction.
#'  \item \code{tree}: dataframe storing the tree returned b y Cassiopeia. This dataframe has a column "group" that indicates
#'  the clonal population to which the node belongs to (a clonal population is a cluster returned by Cassiopeia).
#'  \item \code{tree_phylo}: \code{phylo} object containing the tree returned by Cassiopeia.
#'  \item \code{tree_uncollapsed}: \code{phylo} object containing the tree returned by Cassiopeia prior to the collapsing of mutationless edges.
#'  }
#' 
#' @export infer_phylogeny
#' @import reticulate
#' @importFrom scales comma
#' @importFrom ggtree fortify
# @importFrom dynamicTreeCut cutreeDynamic
# @importFrom ape cophenetic.phylo bind.tree drop.tip
#' @rawNamespace import(dplyr, except = count)
#' @rawNamespace import(ggplot2, except = c(element_render, CoordCartesian))
# @import lemon
infer_phylogeny = function(EvoTraceR_object, mutations_use = 'del_ins') {
  
  if (! (mutations_use %in% c('del', 'del_ins'))) {
    cli::cli_alert_danger("Error, muations use must be one of 'del_ins', 'del'")
    stop('Exiting')
  }
  EvoTraceR_object$phylogeny$mutations_in_phylogeny = mutations_use
  output_dir = file.path(EvoTraceR_object$output_directory, paste0("phylogeny_analysis/phylogeny_", mutations_use))
  
  if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}
  
  asv_bin_var = EvoTraceR_object$alignment$binary_mutation_matrix
  barcode_var = asv_bin_var[EvoTraceR_object$reference$ref_name,]
  if (class(barcode_var) == "numeric") {
    first_name <- colnames(asv_bin_var)[1]
    barcode_var <- setNames(data.frame(barcode_var), first_name)
    rownames(barcode_var) <- EvoTraceR_object$reference$ref_name
  }
  if (mutations_use == 'del_ins') { 
    
    asv_bin_var = asv_bin_var %>%
      dplyr::select(starts_with('i_') | starts_with('d_')) %>% 
      filter(rowSums(dplyr::across(dplyr::everything())) > 0)
    
  } else {
    asv_bin_var = asv_bin_var %>% 
      dplyr::select(starts_with('del_')) %>%
      filter(rowSums(dplyr::across(dplyr::everything())) > 0)
  }
  
  asv_bin_var = dplyr::bind_rows(barcode_var, asv_bin_var)
  phyl_result = compute_tree_cassiopeia(asv_bin_var, EvoTraceR_object$reference$ref_name)
  ape::write.tree(phyl_result$tree_collapsed_phylo, #phyl_result$tree_uncollapsed, 
                  file = paste0(output_dir, "/tree_all_clones.newick"),
                  append = FALSE,
                  digits = 10, tree.names = FALSE)
  EvoTraceR_object$phylogeny$tree = phyl_result$tree_collapsed_df
  EvoTraceR_object$phylogeny$tree_uncollapsed = phyl_result$tree_uncollapsed
  EvoTraceR_object$phylogeny$tree_phylo = phyl_result$tree_collapsed_phylo
  write.csv(EvoTraceR_object$phylogeny$tree, file.path(output_dir, "tree_all_clones.csv"))
  
  return(EvoTraceR_object)
  
}

#' This function combines in one dataframe all the data that have been computed for the phylogenetic analysis.
#' 
#' @title create_df_summary
#' 
#' @examples
#' \dontrun{
#' EvoTraceR_object = create_df_summary(EvoTraceR_object)
#' }
#' 
#' @param EvoTraceR_object (Required).
#' 
#' @return EvoTraceR_object with a field named plot_summary containing dataframe df_to_plot_final.
#' This dataframe containins all the information computed within the package for each ASV in each sample.
#' 
#' @export create_df_summary
#' @import tibble
#' @import ggtree
#' @rawNamespace import(dplyr, except = count)
create_df_summary = function(EvoTraceR_object) {
  mut_in_phyl = EvoTraceR_object$phylogeny$mutations_in_phylogeny
  
  df_to_plot_perf_match = EvoTraceR_object$statistics$all_asv_statistics

  # merge with df_to_plot_perf_match
  sample_order = EvoTraceR_object$sample_order
  df_to_plot_perf_match$sample = factor(df_to_plot_perf_match$sample, levels = sample_order)
  
  df_to_plot_final <- tibble(merge(df_to_plot_perf_match, 
                                   EvoTraceR_object$alignment$ASV_alterations_width, 
                                   by.x ="asv_names", by.y="asv_names"))
  
  
  df_to_plot_final$asv_names <- as.factor(df_to_plot_final$asv_names)
  
  df_to_plot_final = tibble(merge(df_to_plot_final, 
                                  EvoTraceR_object$alignment$mutations_coordinates,
                                  by.x="asv_names", by.y="asv_names")) %>% dplyr::select(-c(spanned_cutSites))
  df_to_plot_final = tibble(dplyr::right_join(df_to_plot_final,
                                              as.data.frame(EvoTraceR_object$phylogeny$tree) %>% filter(isTip) %>% dplyr::select(label, group),
                                              by=c("asv_names" = "label")))
  output_dir = file.path(EvoTraceR_object$output_directory, paste0("phylogeny_analysis/phylogeny_", mut_in_phyl))
  
  if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}  
  
  write.csv(df_to_plot_final, file.path(output_dir, "asv_stat.csv"), quote = F, row.names = F)
  EvoTraceR_object$plot_summary$df_to_plot_final = df_to_plot_final

  return(EvoTraceR_object)
}
