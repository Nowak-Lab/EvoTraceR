cut_phyl_dendogram = function(tree_mp, asv_bin_var) {
  # di default il barcode deve stare da solo o dalle sequenze non mutate
  set.seed(1)
  asv_sub = asv_bin_var[setdiff(rownames(asv_bin_var), c('BC10v0')),]
  dist_j = as.matrix(ade4::dist.binary(as.matrix(asv_bin_var), method=1, diag=F, upper=F))
  #hclust_avg <- hclust(dist_j, method = 'complete')
  dend = phylogram::as.dendrogram(tree_mp)
  set.seed(1)

  # dist jaccard
  
  clust_dist <- dynamicTreeCut::cutreeDynamic(as.hclust(dend), distM=dist_j, minClusterSize = 1)
  # clust_dist <- dynamicTreeCut::cutreeDynamic(hclust_avg, 
  #                                             distM=as.matrix(dist_j), 
  #                                             deepSplit = 3,
  #                                             minClusterSize = 2)
  
  df_clusters = tibble(asv_names = tree_mp$tip.label, cluster = clust_dist)
  
  return(df_clusters)
}

compute_phylogenetic_tree = function(asv_bin_var, phylip_package_path, barcode) {
  # Compute phylogeny
  ancestral <<- as.character(rep(0, dim(asv_bin_var)[2]))
  # asv_sub = asv_bin_var[setdiff(rownames(asv_bin_var), c('BC10v0')),]
  # Run phylip
  set.seed(1980)
  
  tree_mp_all_rmix <- Rphylip::Rmix(X = as.matrix(asv_bin_var),
                                    method = "Camin-Sokal", 
                                    ancestral = ancestral,
                                    path = phylip_package_path, 
                                    cleanup = T)
  
  #rm(ancestral)
  # Rename It

  
  ### Estimating the MRP (matrix representation parsimony) - SuperTree from a set of input trees (Baum 1992; Ragan 1992)) ------------------------------------------------------ 
  set.seed(1980)
  # if (class(tree_mp_all_rmix) != 'phylo') {
  #   tree_mp <- phytools::mrp.supertree(trees=tree_mp_all_rmix, 
  #                                      start="NJ", 
  #                                      method = "optim.parsimony", 
  #                                      root=TRUE)
  # } else {
  #   tree_mp = tree_mp_all_rmix
  # }

  # Option 2: Choose What Tree Number to Plot (Arbitrary Choice) -> need to argue why this one, the best score as 1?
  #tree_mp <- tree_mp_all[100] # plot single tree -> #1 chosen arbitrarily
  
  tree_mp = tree_mp_all_rmix[[1]]
  ### Fortify tree to data frame
  tree_mp_df <- ggtree::fortify(tree_mp)

  tree_mp_df$parent = tree_mp_df$parent + 1
  tree_mp_df$node = tree_mp_df$node + 1
  tree_mp_df$y = tree_mp_df$y + 1
  
  node1_data = tree_mp_df %>% filter(node == 2) 
  root_node = tree_mp_df %>% filter(x == 0) %>% pull(node)
  
  tree_mp_df = tree_mp_df %>% dplyr::add_row(parent = root_node, 
                                             node = 1, 
                                             label = barcode, 
                                             isTip = TRUE, 
                                             x = node1_data$x, 
                                             y = 1, 
                                             branch = node1_data$branch, 
                                             angle = node1_data$angle)
  

  return(list(tree_mp = tree_mp, tree_mp_df = tree_mp_df))
}

cluster_sequences_jaccard = function(asv_bin_var) {
  # The function from package ade4 dist.binary computes the distance matrix for binary data,
  # using a similarity index (parameter method)
  prova2 = as.matrix(ade4::dist.binary(as.matrix(asv_bin_var), method=1, diag=F, upper=F))
}



