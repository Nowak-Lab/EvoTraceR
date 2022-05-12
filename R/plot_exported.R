#' This function plots the phylogenetic tree and the frequency bubble plot for each cluster identified in the phylogenetic tree.
#' 
#' @title plot_cluster_summary
#' 
#' 
#' @examples
#' \dontrun{
#' data(revo_phyl)
#' output_dir = system.file("extdata", "output", package = "EvoTraceR")
#' revo_phyl$output_directory = output_dir
#' summary_plot = plot_cluster_summary(revo_phyl)
#' }
#' 
#' @param EvoTraceR_object (Required).
#' @param cluster_list (Optional). Subset of clusters for which one wants to create the visualization. 
#' 
#' @return NULL. This function stores in the phylogeny output directory a pdf file for each cluster.
plot_cluster_summary = function(EvoTraceR_object, cluster_list = NULL) {
  
  mut_in_phyl = EvoTraceR_object$phylogeny$mutations_in_phylogeny
  
  tree_mp_df = EvoTraceR_object$phylogeny$tree
  
  barcode_tip = tree_mp_df %>% filter(label == EvoTraceR_object$reference$ref_name)
  
  # Cassiopeia puts first the sequences that are not assigned to any cluster: it puts them at the bottom of the tree
  # So, if the barcode is not at the bottom I should put it there, swapping it with another sequence
  if (nrow(barcode_tip) > 0){
    first_tip = tree_mp_df %>% filter(y == 1)
    if (barcode_tip$y != 1) {
      current_barcode_y = barcode_tip$y
      tree_mp_df = tree_mp_df %>% 
        dplyr::mutate(label = ifelse(y == 1, barcode_tip$label, label)) %>%
        dplyr::mutate(label = ifelse(y == current_barcode_y, first_tip$label, label)) 
    }
    
  } 
  
  output_dir = file.path(EvoTraceR_object$output_directory, paste0("phylogeny_", mut_in_phyl))
  
  if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = T)}
  # 
  # write.csv(df_to_plot_final, file.path(output_dir, "df_to_plot_final.csv"), quote = F, row.names = F)
  # EvoTraceR_object$plot_summary$df_to_plot_final = df_to_plot_final
  
  sample_columns = sort(setdiff(colnames(EvoTraceR_object$clean_asv_dataframe), c("asv_names", "seq")))
  
  if (is.null(cluster_list)){
    clusters = levels(tree_mp_df$group)
  } else {
    clusters = cluster_list
    if (length(intersect(cluster_list, levels(tree_mp_df$group))) == 0) {
      stop("None of the clusters provided are found in the phylogenetic tree")
    }
  }
  cli::cli_progress_bar("Plotting clusters", total = length(clusters))
  
  for (c in setdiff(clusters, "0")) {
    
    cluster_tips = tree_mp_df %>% filter(group == c & isTip) %>% pull(label)
    
    toplot_tree_phylo = ape::drop.tip(EvoTraceR_object$phylogeny$tree_phylo,
                                      setdiff(tree_mp_df %>% filter(isTip) %>% pull(label), cluster_tips))
    
    toplot_tree_phylo = ape::drop.tip(return_list$tree_collapsed_phylo,
                                      setdiff(tree_mp_df %>% filter(isTip) %>% pull(label), cluster_tips))
    
    toplot_tree_df = ggtree::fortify(toplot_tree_phylo, ladderize = T, right=T)
    
    tip_order = toplot_tree_df %>% filter(!is.na(label)) %>% arrange(desc(y)) %>% pull(label)
    toplot_tree_phylo = ape::rotateConstr(toplot_tree_phylo, constraint = tip_order)
    
    toplot_tree = ggtree::fortify(toplot_tree_phylo)
    
    sequence_order = tip_order#toplot_tree %>% arrange(y) %>% filter(!is.na(label)) %>% pull(label)
    
    tree_cluster = plot_phylogenetic_tree(toplot_tree)
    tree_cluster = tree_cluster + scale_x_continuous(expand = expand_scale(0,0.6))
      #scale_x_continuous(breaks = seq(-5, length(cluster_tips) + 0.6, by=1))
    
    msa_cna_bc = plot_evotracer(EvoTraceR_object, what = 'msa', 
                                cleaned_deletions = FALSE, 
                                subset_asvs = sequence_order)
    msa_cna_bc = msa_cna_bc + 
      scale_y_discrete(expand = expand_scale(0,0.6))
      scale_y_discrete(breaks = seq(0, length(cluster_tips) - 1, by=1))
    bubble = plot_evotracer(EvoTraceR_object, what = 'frequency',  subset_asvs = sequence_order)
    
    msa.bubble = aplot::insert_right(msa_cna_bc, bubble, width = 0.2)
    
    tree.msa.bubble <- aplot::insert_left(msa.bubble, tree_cluster, width = 1)
    
    prova = gridExtra::grid.arrange(tree_cluster, msa_cna_bc, nrow=1)
    
    ggsave(filename=file.path(output_dir, paste0("summary_mutations_cluster", c, ".pdf")), 
           plot=tree.msa.bubble, width=50,
           height=dim(toplot_tree)[1]*0.6, units = "cm", limitsize = FALSE)
    #cli::cli_progress_update()
  }
  cli::cli_progress_done()
  
}

#' This function plots different aspects of the data.
#' 
#' @title plot_cluster_summary
#' 
#' 
#' @examples
#' \dontrun{
#' data(revo_phyl)
#' output_dir = system.file("extdata", "output", package = "EvoTraceR")
#' revo_phyl$output_directory = output_dir
#' summary_plot = plot_cluster_summary(revo_phyl)
#' }
#' 
#' @param EvoTraceR_object (Required).
#' @param what (Optional). Default to 'phylogeny'. Can be one of \code{c('phylogeny', 'frequency', 'pid', 'asv_length', 'alterations width', 'msa')}
#' 
#' @return NULL. This function stores in the phylogeny output directory a pdf file for each cluster.
#' 
#' @export seq_filtering_plot
#' 
seq_filtering_plot = function(EvoTraceR_object, figure_dir = file.path(EvoTraceR_object$output_directory, "asv_analysis_figures")) {
  track_data = EvoTraceR_object$dada2$seq_filters

  # # assemble data  with all number for each step of filtering
  # track_data <- data.frame(name=as.factor(c("Starting ASVs", "Frequency Filter", "Substitutions Filter", "Flanking Seq. Filter", "Final ASVs")), 
  #                           num=c(orgseq, counts_filtering, endseq_filter, flanking_filtering, clean_asv))
  # set order of columns DGN
  track_data <- dplyr::mutate(track_data, name = fct_relevel(name, c("Starting ASVs", "Frequency Filter", "Substitutions Filter", "Flanking Seq. Filter",  "Final ASVs")))
  track_data$num_names <- paste0(track_data$num, " x ASVs") # numbers of ASV included on the top of bar  
  
  # calculate percent change from previous filter
  track_data <-
    track_data %>%
    mutate(track_data, diff_perc = ceiling(x=(num/lag(num)-1)*100)) %>%
    mutate(diff_perc = paste0(diff_perc, " %")) %>%
    mutate(diff_perc = if_else(name == "Starting ASVs" | name == "Final ASVs", "", diff_perc)) %>%
    mutate(num_names_sf = if_else(name == "Starting ASVs" | name == "Final ASVs", num_names, "")) %>%
    mutate(num_names_ins = if_else(name == "Starting ASVs" | name == "Final ASVs", "", num_names))
  
  # start graph plotting 
  seqtab_df_clean_track <-
    ggplot(data=track_data) +
    geom_bar(aes(x=name, y=num, fill=name), position = "dodge", stat = "identity", width=0.8, size=0.2, show.legend = FALSE) +
    geom_text(aes(x=name, y=num, label=num_names_sf), vjust=-0.25, size=3) + # change order to have up whatever you choose, opposite to order
    geom_text(aes(x=name, y=num, label=num_names_ins), vjust=-1.75, size=3) + # change order to have up whatever you choose, opposite to order
    geom_text(aes(x=name, y=num, label=diff_perc), vjust=-0.25, size=3, col="blue") + # change order to have up whatever you choose, opposite to order
    scale_y_continuous(expand = c(0, 0), 
                        limits= c(0, plyr::round_any(max(track_data$num), 100, f = ceiling)+100/4), 
                        breaks = seq(0, (plyr::round_any(max(track_data$num), 100, f = ceiling)), 100)) +
    scale_fill_manual(values=c("#444c5c", "#aaaaaa", "#aaaaaa", "#aaaaaa", "#78a5a3")) +
    labs(x = "ASVs Filtering Steps", y = "Number of ASVs") + 
    lemon::coord_capped_cart(left="both") + # axis with lemon
    barplot_nowaklab_theme() + # add theme 
    theme(plot.margin = unit(c(0, 0, 0, 0), "mm"), # update theme specifically 
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # , hjust = 1, vjust = 1
          axis.line.x = element_blank(), # disable x axis lines
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank()) # disable x axis ticks lines
  # save pdf
  ggsave(filename=file.path(figure_dir, "track_asv_number.pdf"), plot=seqtab_df_clean_track, width=15, height=15, units = "cm")
  # save csv
  # write.csv(track_data, file.path(figure_dir, "/track_asv_number_data.csv"), row.names = FALSE, quote = FALSE)
}

#' This function plots different aspects of the data.
#' 
#' @title plot_cluster_summary
#' 
#' 
#' @examples
#' \dontrun{
#' data(revo_phyl)
#' output_dir = system.file("extdata", "output", package = "EvoTraceR")
#' revo_phyl$output_directory = output_dir
#' summary_plot = plot_cluster_summary(revo_phyl)
#' }
#' 
#' @param EvoTraceR_object (Required).
#' @param what (Optional). Default to 'phylogeny'. Can be one of \code{c('phylogeny', 'frequency', 'pid', 'asv_length', 'alterations width', 'msa')}
#' 
#' @return NULL. This function stores in the phylogeny output directory a pdf file for each cluster.
#' 
#' @export plot_evotracer
#' 
plot_evotracer = function(EvoTraceR_object, what = 'phylogeny', ...) {
  df = EvoTraceR_object$statistics$all_asv_statistics
  
  sample_order = EvoTraceR_object$sample_order
  df$sample = factor(df$sample, levels = sample_order)
  
  if (what == 'phylogeny') {
    tree_mp_df = EvoTraceR_object$phylogeny$tree
    
    barcode_tip = tree_mp_df %>% filter(label == EvoTraceR_object$reference$ref_name)
    
    # Cassiopeia puts first the sequences that are not assigned to any cluster: it puts them at the bottom of the tree
    # So, if the barcode is not at the bottom I should put it there, swapping it with another sequence
    if (nrow(barcode_tip) > 0){
      first_tip = tree_mp_df %>% filter(y == 1)
      if (barcode_tip$y != 1) {
        current_barcode_y = barcode_tip$y
        tree_mp_df = tree_mp_df %>% 
          dplyr::mutate(label = ifelse(y == 1, barcode_tip$label, label)) %>%
          dplyr::mutate(label = ifelse(y == current_barcode_y, first_tip$label, label)) 
      }
      
    } 
    return(plot_phylogenetic_tree(tree_mp_df, ...))
  } else if (what == 'frequency') {
    return(plot_percentage_asv_sample(df, ...))
  } else if (what == 'pid') {
    return(plot_similarity(df))
  } else if (what == 'asv length') {
    return(plot_asv_length(df))
  } else if (what == 'alteration width') {
    df = EvoTraceR_object$alignment$ASV_alterations_width
    df$sample = factor(df$sample, levels = sample_order)
    return(plot_mutations_width(df))
  } else if (what == 'msa') {
    return(plot_msa(EvoTraceR_object, ...))
  } else if (what == 'Edit frequency') {
    return(plot_alt_hist(EvoTraceR_object))
  } else {
    stop("Error, what must be one of c('phylogeny', 'frequency', 'pid', 'asv_length', 'alterations width', 'msa', 'Edit frequency')")
  }
  
}










