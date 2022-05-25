plot_phylogenetic_tree_ggraph = function(toplot_tree) {# , sample_columns) {
  toplot_tree = toplot_tree %>% arrange(y)
  
  nodes <- data.frame(name=toplot_tree$node,
                      label=toplot_tree$label)
  relations <- data.frame(from=toplot_tree$parent,
                          to=toplot_tree$node,
                          length = toplot_tree$branch.length)
  
  g <- graph_from_data_frame(relations, directed=T, vertices=nodes)
  #prova_edge_list = toplot_tree %>% select(parent, node) %>% rename(from = parent, to = node)
  
  ggtree_mp = ggraph(g, layout = 'dendrogram', circular = FALSE) +#, length = length) +
    geom_node_text(aes(label=label), nudge_y = 0.5 ) +
    geom_edge_diagonal() +
    geom_node_point() +
    #scale_x_continuous(c(0, sum(!is.na(toplot_tree$label)))) +
    theme_void() +
    coord_flip() + 
    scale_y_reverse() +
    barplot_nowaklab_theme() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
          axis.text.y = element_blank(), # disable y axis text
          #axis.title.x = element_text(size=8, angle=0),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank()) +
    # scale_y_continuous(expand = c(0.05, 0.05),
    #                    limits=c(0, 1.15*max(toplot_tree$x)),
    #                    breaks=sort(c(0, 10, max(toplot_tree$x)))) +
    #lim_tree(1.1*max(toplot_tree$x)) +
    ylab("Phylogenetic Tree \n Cassiopeia Greedy") +
    xlab('') +
    theme(panel.border=element_blank(), axis.line = element_line()) #+
  #lemon::coord_capped_cart(bottom="both")
  
  return(ggtree_mp)
}



plot_phylogenetic_tree = function(tree_mp_df) {# , sample_columns) {
  
  # if (!('group' %in% colnames(tree_mp_df))) {
  #   tree_mp_df$group = 1
  #   tree_mp_df$group = as.factor(tree_mp_df$group)
  # }
  
  ggtree_mp <- 
    ggtree::ggtree(tree_mp_df) + #, layout="ellipse") + #%<+%
    # perc_max_tip_colors + # add data for labelling tips
    # geom_tippoint(aes(color = group), size=3) +
    #geom_text2(aes(subset=!isTip, label=node), hjust=-.3)  + 
    ggtree::geom_tiplab(#aes(fill=group, alpha = 0.5, color=NULL), 
                        geom = "text", 
                        align=TRUE, linesize=0.5, linetype="dotted", size=5) +
    #scale_colour_manual(values = sample_col[sample_columns], guide=guide_legend(keywidth=0.5, keyheight=0.5, order=4)) +
    scale_x_continuous(expand = c(0.05, 0.05), 
                       limits=c(0, 1.15*max(tree_mp_df$x)), 
                       breaks=sort(c(0, 10, max(tree_mp_df$x)))) +
    xlim_tree(1.1*max(tree_mp_df$x)) +
    xlab("Phylogenetic Tree \n Cassiopeia Greedy") +
    theme(panel.border=element_blank(), axis.line = element_line()) +
    lemon::coord_capped_cart(bottom="both") + # axis with lemon +
    scale_fill_manual(values=sample(rainbow(n = length(unique(tree_mp_df$group)))))
  
  set.seed(1)
  
  
  if (('group' %in% colnames(tree_mp_df))) {
    colors = sample(rainbow(n = length(unique(tree_mp_df$group))))
    names(colors) = unique(tree_mp_df$group)
    colors[['1']] = 'grey'
    for (c in unique(tree_mp_df$group)) {
      cluster_nodes = tree_mp_df %>% filter(group == c & isTip) %>% arrange(y) %>% pull(node)
      ggtree_mp = ggtree_mp +
        geom_strip(cluster_nodes[1], cluster_nodes[length(cluster_nodes)], barsize=8, color=colors[[c]], 
                   label= c, offset.text=0.25, offset = 1, angle = 90, align = T) 
      # geom_hilight(mapping=aes(subset = node %in% cluster_nodes, 
      #                               fill = group),
      #                   type = "gradient", gradient.direction = 'rt',
      #                   alpha = .8)
  }
  
  }
  
  
  ggtree_mp <-
    ggtree_mp +
    barplot_nowaklab_theme() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
          axis.text.y = element_blank(), # disable y axis text
          #axis.title.x = element_text(size=8, angle=0),
          axis.ticks.x = element_line(colour="black", size=0.5),
          axis.line.x = element_line(colour="black", size=0.5))
  return(ggtree_mp)
}

plot_msa = function(EvoTraceR_object, cleaned_deletions = FALSE, subset_asvs = NULL) {
  
  if (cleaned_deletions == 'del') {
    to_plot_df = EvoTraceR_object$alignment$asv_barcode_alignment %>% 
      mutate(alt = ifelse(alt == 'd', alt, 'w'))
  } else if (cleaned_deletions == 'del_ins') {
    to_plot_df = EvoTraceR_object$alignment$asv_barcode_alignment %>% 
      mutate(alt = ifelse(alt %in% c('d','i'), alt, 'w'))
  } else {
    to_plot_df = EvoTraceR_object$alignment$asv_barcode_alignment
  }
  
  if (! is.null(subset_asvs)) {
    to_plot_df = to_plot_df %>% filter(asv_names %in% subset_asvs)
    if (nrow(to_plot_df) == 0) {
      stop('None of the asv names provided in subset_asvs matches any of the ASVs in the EvoTraceR object.\nPlease select valid ASV names.')
    }
  } else {
    subset_asvs = unique(to_plot_df$asv_names)
  }
  # else if (smoothed_deletions == 'sub_smooth_del_ins') {
  #   to_plot_df = EvoTraceR_object$cleaned_deletions_insertions$asv_barcode_alignment
  # } else {
  #   to_plot_df = EvoTraceR_object$alignment$asv_barcode_alignment
  # }
  
  # if (smoothed_deletions != FALSE) {
  #   # expand deletions of 1nt in a window +-5nts (for visualization purpouses)
  #   del_to_expand = to_plot_df %>% dplyr::group_by(asv_names) %>% 
  #     mutate(next_alt = lead(alt), previous_alt = lag(alt)) 
  #   del_to_expand_left = del_to_expand %>%
  #     filter(alt == 'del' & alt != previous_alt) #& alt != next_alt)
  
  #   del_to_expand_right = del_to_expand %>%
  #     filter(alt == 'del' & alt != next_alt)
  
  #   asvs_with_del = unique(c(as.character(del_to_expand_left$asv_names), as.character(del_to_expand_right$asv_names)))
  
  #   for (asv in asvs_with_del) {
  #     to_plot_sub = to_plot_df %>% filter(asv_names == asv)
  #     to_plot_df = to_plot_df %>% filter(asv_names != asv)
  
  #     to_plot_sub$alt <- ifelse(sapply(to_plot_sub$position_bc260, function(p) 
  #       any((del_to_expand_left$position_bc260 - 5 <= p & p <= del_to_expand_left$position_bc260 & del_to_expand_left$asv_names == asv)|
  #           (del_to_expand_right$position_bc260 + 5 >= p & p >= del_to_expand_left$position_bc260 & del_to_expand_right$asv_names == asv))),
  #       "del",to_plot_sub$alt)
  
  #     to_plot_df = to_plot_df %>% bind_rows(to_plot_sub)
  #   }
  # }
  # Put insertions after wt for visualization
  to_plot_df = to_plot_df %>% dplyr::arrange(asv_names, position_bc260, desc(alt))
  
  
  # add full names for labeling of plots
  to_plot_df <-
    to_plot_df %>%
    mutate(alt_long_names = ifelse(alt == "i", "Insertion", 
                                   ifelse(alt == "d", "Deletion", 
                                          ifelse(alt == "s", "Substitutions", "No Edits"))))
  
  # Position of PAM in guides
  pam_pos <- EvoTraceR_object$reference$ref_cut_sites
  bc_len = nchar(EvoTraceR_object$reference$ref_seq)
  # Adjust for height of tiles -> "sub" smaller
  
  to_plot_df$tile_height <- ifelse(to_plot_df$alt == "s", 0.3, 0.75)
  # frames around msa apart ORG
  msa_frame <- data.frame(xmin= 1, xmax =bc_len, 
                          ymin = seq(from = 1.6, to = nlevels(as.factor(to_plot_df$asv_names)), by=1), 
                          ymax = seq(from = 2.4, to = nlevels(as.factor(to_plot_df$asv_names))+1, by=1))
  
  ### "msa Plot" - Barcode; Scale (1-260)  ------------------------------------------------------
  #to_plot_df$alt = factor(x = to_plot_df$alt, levels = c('sub', 'del', 'wt', 'ins'))
  to_plot_df$alt_long_names = factor(x = to_plot_df$alt_long_names)#, levels = c('Deletion', 'w', 'i'))
  to_plot_df$asv_names = factor(x = to_plot_df$asv_names, levels = subset_asvs)
  msa_cna_bc <- 
    ggplot(data=to_plot_df, aes(x=position_bc260, y=asv_names)) +
    geom_tile(aes(fill=alt_long_names, width=0.75, height=tile_height), colour = NA) +
    scale_fill_manual(values=c("Deletion"="#3366FF", "Insertion"="#FF0033", "No Edits"="#f2f2f2"), 
                      breaks=c("Insertion", "No Edits", "Deletion")) + 
    geom_vline(xintercept=pam_pos, linetype="solid", size=0.3, col="grey50") + # lines for guide targets
    geom_vline(xintercept=pam_pos, linetype="dashed", size=0.4, col="#ff8300") + # Cas9 Cleavage
    scale_x_continuous(labels=scales::comma, breaks=c(1, seq(26, 260, 26)), expand = c(0.014, 0.014)) +
    geom_rect(data=msa_frame, mapping=aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax"), colour = "grey50", fill=NA, inherit.aes = F, size=0.6) +
    geom_rect(xmin=1, xmax=260, ymin=0.6, ymax=1.4, colour = "#65A7F3", fill=NA, size=0.6) +
    xlab("Barcode Nucleotides \n (1-260)") +
    labs(fill = "Type of Editing") +
    lemon::coord_capped_cart(bottom="both") # axis with lemon
  #"sub"="#329932",
  # add theme
  msa_cna_bc <- 
    msa_cna_bc +
    barplot_nowaklab_theme() +
    scale_y_discrete(c(0,length(subset_asvs))) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
          axis.line.y = element_blank(), # disable y axis lines
          axis.ticks.y = element_blank(), # disable y axis ticks lines
          axis.title.y = element_blank(), # disable y axis lines
          axis.text.y = element_blank()) # disable y axis text
  return(msa_cna_bc)
}
# The input must be a tibble with the following columns: asv_names, width_total_del, width_total_ins, width_total_sub
# This tibble is found in EvoTraceR_object$alignment$asv_alterations_width
plot_mutations_width = function(df_to_plot_final) {
  bar_ins_del_sub_width <- 
    ggplot(data=dplyr::select(df_to_plot_final, asv_names, width_total_d, width_total_i) %>% #, width_total_s) %>% 
             unique() %>% 
             gather(key="alt", "width_total_d", "width_total_i",  value="width_sum") %>% #"width_total_s",
             filter(!str_detect(asv_names, "NMBC"))) +
    geom_bar(aes(x=width_sum, y=asv_names, fill=alt), position="stack", stat="identity", width=0.8, size=0.2) +
    scale_fill_manual(values=c("width_total_d" = "#3366FF", "width_total_i" = "#FF0033"), #"width_total_s" = "#329932" 
                      breaks=c("width_total_d", "width_total_i")) + #, "width_total_s"
    geom_vline(xintercept=0, linetype="solid", size=0.5, col="black") + # no indels
    scale_x_continuous(labels=scales::comma, expand = c(0.01, 0.01), breaks=c(0, 130, 260))+#, limits=c(0, 286)) +
    geom_vline(xintercept=130, linetype="dotted", size=0.5, col="gray50") +
    geom_vline(xintercept=260, linetype="dotted", size=0.5, col="gray50") +
    labs(x = "Cummulative \n Widths of Markings", fill = 'Width of Marking') +
    theme(panel.border=element_blank(), axis.line = element_line()) + #From: https://github.com/stefanedwards/lemon/issues/26
    lemon::coord_capped_cart(bottom="both") # axis with lemon
  
  # Add Theme  
  bar_ins_del_sub_width <- 
    bar_ins_del_sub_width + 
    barplot_nowaklab_theme() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
          axis.line.y = element_blank(), # disable y axis lines
          axis.ticks.y = element_blank(), # disable y axis ticks lines
          axis.title.y = element_blank(), # disable y axis lines
          axis.text.y = element_blank()) # disable y axis text
  
  return(bar_ins_del_sub_width)
  
}
# Input must be a tibble with the following columns: asv_names, seq_n (length of the ASV)
plot_asv_length = function(df_to_plot_final) {
  bar_seq_n <- 
    ggplot(data=dplyr::select(df_to_plot_final, asv_names, seq_n) %>% 
             unique() %>%
             mutate(seq_n_col=ifelse(seq_n == 260, "260", ifelse(seq_n < 260, "< 260", "> 260"))) %>%
             filter(str_detect(asv_names, "ASV"))) +
    geom_bar(aes(x=seq_n, y=asv_names, fill=seq_n_col), position = "dodge", stat = "identity", width=0.8, size=0.2) +
    geom_text(aes(x=1, y=asv_names, label = seq_n), col="white", size=4, hjust="left") + #, position=position_dodge(width=0.9)
    scale_fill_manual(values=c("260" = "#84B48F", "< 260" = "#377EB8", "> 260" = "#F39B7FFF"), breaks=c("260", "< 260", "> 260")) +
    geom_vline(xintercept=0, linetype="dotted", size=0.5, col="#377EB8") + # 0 
    geom_vline(xintercept=260, linetype="dotted", size=0.5, col="#84B48F") + # 260 bp = original length
    geom_vline(xintercept=520, linetype="dotted", size=0.5, col="#F39B7FFF") + # 0 
    scale_x_continuous(labels=scales::comma, expand = c(0.01, 0.01), limits=c(0, 572), breaks=c(0, 130, 260, 390, 520)) + # 572 = 1.1 * 520 -> nice separation between graphs
    labs(x = "ASV\nLength", fill = "ASV Length") +
    lemon::coord_capped_cart(bottom="both") # axis with lemon
  
  # Add Theme  
  bar_seq_n <- 
    bar_seq_n + 
    barplot_nowaklab_theme() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
          axis.line.y = element_blank(), # disable y axis lines
          axis.ticks.y = element_blank(), # disable y axis ticks lines
          axis.title.y = element_blank(), # disable y axis lines
          axis.text.y = element_blank()) # disable y axis text
  
  return(bar_seq_n)
}

# Plot with the percentage of sequence similarity between each ASV and the original barcode.
# Input must be a tibble that contains two columns: asv_names and pid.
plot_similarity = function(df_to_plot_final) {
  ### "Slider Plot" - Percent (%) Sequence Similarity ------------------------------------------------------ 
  # it shows similarity using Pairwise Alignment to original BC10.ORG and i.e. 100% means max similarity to the BC10.ORG
  bar_pid <- 
    ggplot(data=dplyr::select(df_to_plot_final, asv_names, pid) %>% unique() %>%
             filter(str_detect(asv_names, "ASV")) %>%
             dplyr::mutate_if(is.numeric, round, 0), 
           aes(x=pid, y=asv_names)) +
    geom_segment(aes(y=asv_names, yend=asv_names, x=0, xend=100), size=2, color="gray50") +
    geom_vline(xintercept=c(0, 25, 50, 75, 100), linetype="solid", size=0.5, col="white") + # 100 similar
    geom_point(size=9, fill="#ce5a57", col="white", shape=21, stroke=0.5) +
    geom_text(aes(label=pid), size=4, col="white") + 
    scale_x_continuous(labels=scales::comma, expand = c(0.01, 0.01), limits=c(0, 110), breaks=c(0, 25, 50, 75, 100)) +
    xlab("Sequence \n Identity (%)") +
    lemon::coord_capped_cart(bottom="both") # axis with lemon
  # Add Theme  
  bar_pid <- 
    bar_pid + 
    barplot_nowaklab_theme() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
          axis.line.y = element_blank(), # disable y axis lines
          axis.ticks.y = element_blank(), # disable y axis ticks lines
          axis.title.y = element_blank(), # disable y axis lines
          axis.text.y = element_blank()) # disable y axis text
  return(bar_pid)
}

# Dotplot where the size indicate the normalized asv count with respect to the total
# count in each sample, and the color indicates the counts of each ASV normalized with respect to
# to the counts for the same ASV in the sample with the highest abundance.
# input must be a tibble with this columns: asv_names, sample, perc_in_sample, perc_fold_to_max
plot_percentage_asv_sample = function(df_to_plot_final, subset_asvs = NULL) {
  
  
  
  if (! is.null(subset_asvs)) {
    df_to_plot_final = df_to_plot_final %>% filter(asv_names %in% subset_asvs)
    if (nrow(df_to_plot_final) == 0) {
      stop('None of the asv names provided in subset_asvs matches any of the ASVs in the EvoTraceR object.\nPlease select valid ASV names.')
    }
  } else {
    subset_asvs = unique(df_to_plot_final$asv_names)
  }
  
  df_to_plot_final$asv_names = factor(df_to_plot_final$asv_names, levels = subset_asvs)
  
  scale_bubble <- ceiling(max(df_to_plot_final %>% 
                                dplyr::pull(perc_in_sample), na.rm = T))
  
  scale_bubble_nonmbc <- ceiling(max(df_to_plot_final %>% 
                                       filter(str_detect(asv_names, "ASV")) %>% 
                                       dplyr::pull(perc_in_sample), na.rm = T))
  
  # df_to_plot_final = df_to_plot_final %>% filter(str_detect(asv_names, "ASV") | str_detect(asv_names, '.NMBC'))
  # df_to_plot_final$asv_names = stringr::str_replace_all(string = df_to_plot_final$asv_names, pattern = '.NMBC', '')
  # plot
  bubble <- 
    ggplot(data=df_to_plot_final, #%>% 
           #filter(str_detect(asv_names, "ASV")),
           aes(x=sample, y=asv_names, scale="globalminmax")) + # %>% darop_na() -> for BC10.ORG
    geom_point(aes(size=perc_in_sample, fill=perc_fold_to_max), shape=21, stroke=0.5, col="black") +
    #scale_fill_gradient2(low = "#417dd4", mid = "#f2f2f2", high = "#DC0000FF") + # for log2
    colorspace::scale_fill_continuous_sequential(palette = "Plasma", limits=c(0, 100), rev = F, na.value = "grey") +
    scale_size_area(#max_size = scale_bubble, 
      limits = c(0, scale_bubble), 
      trans = 'log1p',
      breaks = c(1, scale_bubble_nonmbc/2, scale_bubble_nonmbc, scale_bubble)) + # manual
    xlab("Analyzed \n Samples") +
    labs(x = "Analyzed \n Samples", size = "Frequency in\nSample", fill = "Frequency\nNormalized to Max") +
    lemon::coord_capped_cart(bottom="both") + # axis with lemon
    scale_y_discrete(c(0,length(subset_asvs))) +
    # Add theme
    barplot_nowaklab_theme() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
          legend.title.align = 0.5,
          axis.line.y = element_blank(), # disable y axis lines
          axis.ticks.y = element_blank(), # disable y axis ticks lines
          axis.title.y = element_blank(), # disable y axis lines
          axis.text.y = element_blank(), # disable y axis text
          #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.grid.major.y = element_line(colour="#b4b4b4FF", size=0.4, "dotted"), # y grid line 
          panel.grid.major.x = element_line(colour="#b4b4b4FF", size=0.4, "dotted")) # disable lines in grid on X-axis
  return(bubble)
  
}
