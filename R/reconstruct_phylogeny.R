
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
  
  #########################
  edge_dict = cas_tree[["_CassiopeiaTree__network"]][["_adj"]]
  
  for (e in names(edge_dict)) {
    if (length(edge_dict[[e]]) > 0) {
      subdict = edge_dict[[e]]
      nmut_parent = sum(cas_tree$get_character_states(e) == 1)
      for (child in names(subdict)) {
        nmut_child = sum(cas_tree$get_character_states(child) == 1)
        cas_tree$set_branch_length(parent = e, child = child, 
                                   length = (nmut_child - nmut_parent))
      }
    }
  }
  
  tree = ape::read.tree(text=cas_tree$get_newick(record_branch_lengths = T))
  
  tree_df = ggtree::fortify(tree, ladderize = T, right=T)
  
  tip_order = tree_df %>% filter(!is.na(label)) %>% arrange(desc(y)) %>% pull(label)
  tree = ape::rotateConstr(tree, constraint = tip_order)
  
  tree = ggtree::fortify(tree)
  # Take the parent node with the minimum label, which should be the the one to which all subtrees corresponding to clusters
  # are attached to, and then group the clades to assign the clonal populations label.
  mock_root = min(tree$parent) #+ 1
  clusters_roots = tree %>% filter(parent == mock_root) %>% pull(node)
  tree_df = ggtree::groupClade(tree, .node=clusters_roots)
  
  tree_df = tree_df %>% mutate(group = as.numeric(group) - 1)
  
  tree_df$group = paste0("CP", formatC(tree_df$group,
                                       width = nchar(trunc(max(tree_df$group))),
                                       format = "d", flag = "0"))
  
  # Re-compute te phylo object so that tips can be easily removed with ape::drop.tip
  
  tree = ape::read.tree(text=cas_tree$get_newick(record_branch_lengths = T))
  
  tree_df$group = as.factor(tree_df$group)
  
  return_list$tree_collapsed_df = tree_df
  return_list$tree_collapsed_phylo = tree
  
  return(return_list)
}





