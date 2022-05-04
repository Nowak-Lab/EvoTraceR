

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


compute_tree_cassiopeia = function(asv_bin_var, barcode) {
  # Compute phylogeny
  cas = reticulate::import('cassiopeia')
  pd = reticulate::import('pandas')
  
  np = reticulate::import('numpy')
  
  df = reticulate::r_to_py(asv_bin_var)
  
  cas_tree = cas$data$CassiopeiaTree(character_matrix=df)#, root_sample_name='BC10v0')
  
  vanilla_greedy = cas$solver$VanillaGreedySolver()
  
  vanilla_greedy$solve(cas_tree, collapse_mutationless_edges=F)
  
  return_list = list(tree_uncollapsed = ape::read.tree(text=cas_tree$get_newick(record_branch_lengths = T)))
  
  cas_tree$collapse_mutationless_edges(infer_ancestral_characters=T)
  
  tree = ape::read.tree(text=cas_tree$get_newick(record_branch_lengths = F))
  
  tree_df = ggtree::fortify(tree)
  # Find node to root in the original sequence in "bc10_org"
  tree_root_mp <-
    tree_df %>%
    filter(label == 'BC10v0') %>%
    pull(node)
  
  
  ### Re-root in the Original Sequence ------------------------------------------------------ 
  tree <- TreeTools::RootTree(tree, tree_root_mp)
  
  tree = ggtree::fortify(tree)
  # Take the parent node with the minimum label, which should be the the one to which all subtrees corresponding to clusters
  # are attached to.
  mock_root = min(tree$parent) + 1
  clusters_roots = tree %>% filter(parent == mock_root) %>% pull(node)
  tree_df = ggtree::groupClade(tree, .node=clusters_roots)
  #table(tree_df %>% filter(!is.na(label)) %>% pull(group))
  
  tree_df = tree_df %>% mutate(group = as.numeric(group))
  
  small_groups = tree_df %>% 
    filter(isTip) %>%
    plyr::count('group') %>% filter(freq == 1) %>% pull(group)
  
  # Riassegno lo 0 ai gruppi singoli
  tree_df$group[tree_df$group %in% small_groups] = 0
  
  # Re-compute te phylo object so that tips can be easily removed with ape::drop.tip
  
  tree = ape::read.tree(text=cas_tree$get_newick(record_branch_lengths = F))
  
  # Find node to root in the original sequence in "bc10_org"
  tree_root_mp <-
    ggtree::fortify(tree) %>%
    filter(label == 'BC10v0') %>%
    pull(node)
  
  
  ### Re-root in the Original Sequence ------------------------------------------------------ 
  tree <- TreeTools::RootTree(tree, tree_root_mp)
  
  tree_df$group = as.factor(tree_df$group)
  
  return_list$tree_collapsed_df = tree_df
  return_list$tree_collapsed_phylo = tree

  return(return_list)
}





