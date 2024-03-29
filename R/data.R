#' @name EvoTraceR_object
#' @title example data obtained from running function \code{asv_analys}, where the original barcode is identified and all ASVs are aligned to it. 
#' @description example data obtained from running function \code{asv_analys}, where the original barcode is identified and all ASVs are aligned to it. 
#' @docType data
#' @usage data(EvoTraceR_object)
#' @format EvoTraceR object updated with the following fields: 
#' \describe{
#' \item{\code{clean_asv_dataframe}:}{ASV sequences identified post-filtering (contamination removed,
#' sequences with a similarity higher than \code{pid_cutoff_nmbc} to the original barcode
#' aggregated to it and ASVs named in increasing order (ASV01, ASV02, etc.) according
#' to their total counts.}
#' \item{\code{reference}:}{Info about the reference sequence used for the current analysis.}
#'  \item{\code{statistics}:}{another list with the following sub-fields:}
#' \describe{
#' \item{\code{asv_df_percentages}:}{dataframe with six columns. \code{asv_names} is the name of the ASV.
#' \code{sample} is the sample identifier (e.g. ID of an organ or, in case of longitudinal data, of the timepoint);
#' \code{count}: total counts per million for a specific ASV in a specific sample;
#' \code{perc_in_sample}: CPM normalized to the total counts in the corresponding sample;
#' \code{perc_asv}: CPM normalized to the total counts for the corresponding ASV;
#' \code{perc_fold_to_max}: CPM normalized to the maximum counts observed for the corresponding ASV in a sample.}
#' \item{\code{asv_totalCounts}:}{for each ASV, total counts and number of samples in which it was detected.}
#' \item{\code{sample_totalcounts}:}{for each sample, total counts and number of distinct ASVs detected.}
#' \item{\code{asv_diversity_persample}:}{measures of clonal richness and measures of heterogeneity computed for each sample based on the ASVs detected.}
#' \item{\code{asv_persample_frequency}:}{counts for each ASV in each sample.}
#' \item{\code{asv_persample_detection}:}{binary matrix indicating whether a sequence has been detected in the corresponding sample.}
#' \item{\code{asv_toBarcode_similarity}:}{edit distance, percentage similarity and alignment score of each ASV compared to the original barcode.}
#' \item{\code{all_asv_statistics}:}{all the statistics computed on each ASV grouped together in the same tibble.}
#' }
#' \item{\code{alignment},}{another list with the following fields:}
#' \describe{
#' \item{\code{Binary_mutation_matrix}:}{binary matrix encoding the presence/absence of a mutation in an ASV.}
#' \item{\code{asv_barcode_alignment}:}{tibble where each line corresponds to a position in a ASV, and the columns encode the following information:}
#' \describe{
#' \item{asv_names:}{name of the ASV}
#' \item{sample:}{sample identifier}
#' \item{position_bc260:}{position of the alteration in the original barcode. Note that insertions
#' are assigned to the position that coincides with their beginning.}
#' \item{alt:}{type of alteration. wt = Wild Type (i.e. non-mutated position). sub = substitution. del = deletion. ins = insertion.}
#' \item{ref_asv and read_asv}{: respectively, the reference nucleotide observed in the original barcode and the one observed on the sequence.}
#' }
#' \item{\code{mutations_coordinates}:}{dataframe containing a list of all mutations, with their start and end position.}
#' \item{\code{ASV_alterations_width}:}{dataframe containing the number of nucleotides affected by each type of mutation in each ASV.}
#' }
#' }
NULL  

