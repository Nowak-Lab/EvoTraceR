cut_phyl_dendogram = function(tree_mp) {
  # di default il barcode deve stare da solo o dalle sequenze non mutate
  set.seed(1)
  dend = phylogram::as.dendrogram(tree_mp)
  set.seed(1)
  #dend = ape::as.hclust.phylo(tree_mp)
  dist = ape::cophenetic.phylo(tree_mp)
  clust_dist <- dynamicTreeCut::cutreeDynamic(as.hclust(dend), distM=dist, minClusterSize = 1)
  
  df_clusters = tibble(asv_names = tree_mp$tip.label, cluster = clust_dist)
  
  return(df_clusters)
}

compute_phylogenetic_tree = function(asv_bin_var, phylip_package_path, barcode) {
  # Compute phylogeny
  ancestral <<- as.character(rep(0, dim(asv_bin_var)[2]))
  
  # Run phylip
  set.seed(1980)
  
  tree_mp_all_rmix <- Rphylip::Rmix(X=as.matrix(asv_bin_var),
                                    method = "Camin-Sokal", 
                                    ancestral = ancestral,
                                    path = phylip_package_path, 
                                    cleanup = T)
  
  #rm(ancestral)
  # Rename It

  
  ### Estimating the MRP (matrix representation parsimony) - SuperTree from a set of input trees (Baum 1992; Ragan 1992)) ------------------------------------------------------ 
  set.seed(1980)
  if (class(tree_mp_all_rmix) != 'phylo') {
    tree_mp <- phytools::mrp.supertree(trees=tree_mp_all_rmix, 
                                       start="NJ", 
                                       method = "optim.parsimony", 
                                       root=FALSE)
  } else {
    tree_mp = tree_mp_all_rmix
  }

  # Option 2: Choose What Tree Number to Plot (Arbitrary Choice) -> need to argue why this one, the best score as 1?
  #tree_mp <- tree_mp_all[100] # plot single tree -> #1 chosen arbitrarily
  
  
  ### Fortify tree to data frame
  tree_mp_df <- ggtree::fortify(tree_mp)
  # Find node to root in the original sequence in "bc10_org"
  tree_root_mp <-
    tree_mp_df %>%
    filter(label == barcode) %>%
    pull(node)
  ### Re-root in the Original Sequence
  tree_mp <- TreeTools::RootTree(tree_mp, tree_root_mp)
  return(tree_mp)
}

# cluster_sequences_jaccard = function(asv_bin_var) {
#   # The function from package ade4 dist.binary computes the distance matrix for binary data,
#   # using a similarity index (parameter method)
#   prova2 = as.matrix(ade4::dist.binary(as.matrix(asv_bin_var), method=1, diag=F, upper=F))
# }



