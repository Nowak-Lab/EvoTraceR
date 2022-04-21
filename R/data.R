#' @name revo_initialized
#' @title example data obtained from running function \code{initialize_EvoTraceR}, where dada2 is run on the fastqs provided in input. 
#' @description example data obtained from running function \code{initialize_EvoTraceR}, where dada2 is run on the fastqs provided in input. 
#' @docType data
#' @usage data(revo_initialized)
#' @format An object of class EvoTraceR, which is a list with the following fields:
#' \describe{
#' \item{fastq_directory}{directory where the input fastq files are located.}
#' \item{output_directory}{directory where all the output files are being stored.}
#' \item{map_file_sample}{dataframe has as many rows as the input datasets, and for each input stores the sample (e.g. organ or day for longitudinal data)
#' to which it is associated.}
#' \item{dada2_asv_prefilter}{dataframe that stores all sequences detected by \code{dada2}. Note that
#' these sequences still need to be filtered.}
#' \item{dada2}{list which contains the percentage of chimeras found by \code{dada2} and a dataframe
#' that tracks the number of sequences during all \code{dada2} steps.}
#' }. 
#' An object of class EvoTraceR
NULL

#' @name revo_analyzed
#' @title example data obtained from running function \code{asv_analys}, where the original barcode is identified and all remaining ASVs are aligned to it. 
#' @description example data obtained from running function \code{asv_analys}, where the original barcode is identified and all remaining ASVs are aligned to it. 
#' @docType data
#' @usage data(revo_analyzed)
#' @format EvoTraceR object updated with the following fields: 
#' #' \describe{
#' \item{clean_asv_dataframe}{ASV sequences identified post-filtering (contamination removed,
#' sequences with a similarity higher than \code{pid_cutoff_nmbc} to the original barcode
#' aggregated to it and ASVs named in increasing order (ASV01, ASV02, etc.) according
#' to their total counts. }
#' \item{barcode}{Info about the barcode selected for the current analysis.}
#'
#'  \item{statistics}{another list with the following sub-fileds: 
#' \describe{
#' 
#' \item{asv_df_percentages}{dataframe with six columns. \code{asv_names} is the name of the ASV.
#' \code{day_organ} is the sample identifier (e.g. ID of an organ or, in case of longitudinal data, of the timepoint);
#' \code{count}: total counts for a specific ASV in a specific sample;
#' \code{perc_in_sample}: counts normalized to the total counts in the corresponding sample;
#' \code{perc_asv}: conuts normalized to the total counts for the corresponding ASV;
#' \code{perc_fold_to_max}: counts normalized to the maximum counts observed for the corresponding ASV in a sample.}
#' 
#' \item{asv_totalCounts}{for each ASV, total counts and number of samples in which it was detected.}
#' \item{sample_totalCounts}{for each sample, total counts and number of distinct ASVs detected.}
#' \item{asv_diversity_perSample}{measures of clonal richness and measures of heterogeneity computed for each sample based on the ASVs detected.}
#' \item{asv_persample_frequency}{counts for each ASV in each sample.}
#' \item{asv_persample_detection}{binary matrix indicating whether a sequence has been detected in the corresponding sample.}
#' \item{asv_toBarcode_similarity}{edit distance, percentage similarity and alignment score of each ASV compared to the original barcode.}
#' }}
#' }  
#' @return EvoTraceR object updated with statistics about the ASVs.
NULL

#' @name revo_msa
#' @title example data obtained from running function \code{compute_msa}, where alterations are identified in all the ASvs. 
#' @description example data obtained from running function \code{compute_msa}, where alterations are identified in all the ASvs.  
#' @docType data
#' @usage data(revo_msa)
#' @format EvoTraceR object with a new field named \code{alignment}, which is a list with the following fields:
#' \describe{
#' \item{msa_stringset}{output of MSA performed with MUSCLE}
#' \item{mutations_df}{tibble where each line corresponds to a position in a ASV, and the columns encode the name of the ASV, the sample
#' the position of the alteration in the original barcode, the type of alteration and the percentages of sequences that map to a sample and display the alteration.
#' }
#' } 
NULL


#' @name revo_phyl
#' @title example data obtained from running function \code{infer_phylogeny}, using only smoothed deletions for phylogeny reconstruction.
#' @description example data obtained from running function \code{infer_phylogeny}, using only smoothed deletions for phylogeny reconstruction.  
#' @docType data
#' @usage data(revo_msa)
#' @format EvoTraceR object with a new field named \code{phylogeny}, which is a list with the following fields:
#' \describe{
#' \item{mutations_in_phylogeny}{string indiating which mutations were used for phylogeny recontruction}
#' \item{tree}{Phylogenetic tree reconstructed by Rmix}
#' } 
NULL
