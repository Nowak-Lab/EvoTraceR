

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
  if (class(tree_mp_all_rmix) != 'phylo') {
    tree_mp <- tree_mp_all_rmix[[1]]
  } else {
    tree_mp = tree_mp_all_rmix
  }

  # Option 2: Choose What Tree Number to Plot (Arbitrary Choice) -> need to argue why this one, the best score as 1?
  #tree_mp <- tree_mp_all[100] # plot single tree -> #1 chosen arbitrarily
  
  #tree_mp = tree_mp_all_rmix[[1]]
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


# compute_tree_cassiopeia = function(asv_bin_var, barcode) {
#   # Compute phylogeny
#   cas = reticulate::import('cassiopeia')
#   pd = reticulate::import('pandas')
#   np = reticulate::import('numpy')
#   
#   df = reticulate::py_to_r(asv_bin_var)
#   
#   cas_tree = cas$data$CassiopeiaTree(character_matrix=df)
#   
#   vanilla_greedy = cas$solver$VanillaGreedySolver()
#   
#   vanilla_greedy$solve(cas_tree, collapse_mutationless_edges=F)
#   
#   cas_tree$get_newick(record_branch_lengths = T)
#   
#   tree_uncollapsed = ape::read.tree(text=cas_tree$get_newick(record_branch_lengths = T))
#   
#   cas_tree$collapse_mutationless_edges(infer_ancestral_characters=T)
#   
#   tree_collapsed = ape::read.tree(text=cas_tree$get_newick(record_branch_lengths = T))
#   
#   ### Fortify tree to data frame
#   tree_mp_df <- ggtree::fortify(tree_collapsed)
#   
#   return(list(tree_uncollapsed = tree_uncollapsed, 
#               #tree_collapsed = tree_collapsed, 
#               tree_collapsed_df = tree_mp_df))
# }
# 
# 
# 
# 
# 
