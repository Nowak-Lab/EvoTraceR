plot_phylogenetic_tree = function(tree_mp_df, sample_columns) {
  ggtree_mp <- 
    ggtree::ggtree(tree_mp_df) + #%<+%
    # perc_max_tip_colors + # add data for labelling tips
    #geom_tippoint(aes(color = cluster), size=3) +
    ggtree::geom_tiplab(#aes(fill=NA, alpha = 0.5), geom = "label", 
                        align=TRUE, linesize=0.5, linetype="dotted", size=6) +
    #scale_colour_manual(values = sample_col[sample_columns], guide=guide_legend(keywidth=0.5, keyheight=0.5, order=4)) +
    scale_x_continuous(expand = c(0.05, 0.05), limits=c(0, 1.15*max(tree_mp_df$x)), breaks=sort(c(0, 10, max(tree_mp_df$x)))) +
    xlim_tree(1.1*max(tree_mp_df$x)) +
    xlab("Phylogenetic Tree \n Maximum Parsimony Camin-Sokal") +
    theme(panel.border=element_blank(), axis.line = element_line()) +
    lemon::coord_capped_cart(bottom="both") # axis with lemon
  
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

plot_msa = function(REvoBC_object, smoothed_deletions = FALSE) {
  
  if (smoothed_deletions == 'smooth_del') {
    to_plot_df = REvoBC_object$smoothed_deletions_insertions$asv_barcode_alignment %>% 
      mutate(alt = ifelse(alt == 'del', alt, 'wt'))
  } else if (smoothed_deletions == 'smooth_del_ins'){
    to_plot_df = REvoBC_object$smoothed_deletions_insertions$asv_barcode_alignment %>% 
      mutate(alt = ifelse(alt %in% c('del','ins'), alt, 'wt'))
  } else if (smoothed_deletions == 'sub_smooth_del_ins') {
    to_plot_df = REvoBC_object$smoothed_deletions_insertions$asv_barcode_alignment
  } else {
    to_plot_df = REvoBC_object$alignment$asv_barcode_alignment
  }
  
  if (smoothed_deletions != FALSE) {
    # expand deletions of 1nt in a window +-5nts (for visualization purpouses)
    del_to_expand = to_plot_df %>% dplyr::group_by(asv_names) %>% 
      mutate(next_alt = lead(alt), previous_alt = lag(alt)) 
    del_to_expand_left = del_to_expand %>%
      filter(alt == 'del' & alt != previous_alt) #& alt != next_alt)
    
    del_to_expand_right = del_to_expand %>%
      filter(alt == 'del' & alt != next_alt)
    
    asvs_with_del = unique(c(as.character(del_to_expand_left$asv_names), as.character(del_to_expand_right$asv_names)))
    
    for (asv in asvs_with_del) {
      to_plot_sub = to_plot_df %>% filter(asv_names == asv)
      to_plot_df = to_plot_df %>% filter(asv_names != asv)
      
      to_plot_sub$alt <- ifelse(sapply(to_plot_sub$position_bc260, function(p) 
        any((del_to_expand_left$position_bc260 - 5 <= p & p <= del_to_expand_left$position_bc260 & del_to_expand_left$asv_names == asv)|
            (del_to_expand_right$position_bc260 + 5 >= p & p >= del_to_expand_left$position_bc260 & del_to_expand_right$asv_names == asv))),
        "del",to_plot_sub$alt)
      
      to_plot_df = to_plot_df %>% bind_rows(to_plot_sub)
    }
  }
  # Put insertions after wt for visualization
  to_plot_df = to_plot_df %>% dplyr::arrange(asv_names, position_bc260, desc(alt))
  
  # Position of PAM in guides
  pam_pos <- REvoBC_object$reference$ref_cut_sites
  bc_len = nchar(REvoBC_object$reference$ref_seq)
  # Adjust for height of tiles -> "sub" smaller
  
  to_plot_df$tile_height <- ifelse(to_plot_df$alt == "sub", 0.3, 0.75)
  # frames around msa apart ORG
  msa_frame <- data.frame(xmin= 1, xmax =bc_len, 
                          ymin = seq(from = 1.6, to = nlevels(as.factor(to_plot_df$asv_names)), by=1), 
                          ymax = seq(from = 2.4, to = nlevels(as.factor(to_plot_df$asv_names))+1, by=1))
  
  ### "msa Plot" - Barcode; Scale (1-260)  ------------------------------------------------------
  to_plot_df$alt = factor(x = to_plot_df$alt, levels = c('sub', 'del', 'wt', 'ins'))
  msa_cna_bc <- 
    ggplot(data=to_plot_df, aes(x=position_bc260, y=asv_names)) +
    geom_tile(aes(fill=alt, width=0.75, height=tile_height), colour = NA) +
    scale_fill_manual(values=c("del"="#3366FF", "sub"="#329932", "ins"="#FF0033", "wt"="#f2f2f2"), breaks=c("ins", "wt", "del", "sub")) +
    geom_vline(xintercept=REvoBC_object$reference$ref_cut_sites, linetype="solid", size=0.3, col="grey50") + # lines for guide targets
    geom_vline(xintercept=pam_pos, linetype="dashed", size=0.4, col="#ff8300") + # Cas9 Cleavage
    scale_x_continuous(labels=scales::comma, breaks=c(1, seq(26, 260, 26)), expand = c(0.014, 0.014)) +
    geom_rect(data=msa_frame, mapping=aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax"), colour = "grey50", fill=NA, inherit.aes = F, size=0.6) +
    geom_rect(xmin=1, xmax=260, ymin=0.6, ymax=1.4, colour = "#65A7F3", fill=NA, size=0.6) +
    xlab("Barcode Nucleotides \n (1-260)") #+
    #lemon::coord_capped_cart(bottom="both") # axis with lemon
  
  # add theme
  msa_cna_bc <- 
    msa_cna_bc +
    barplot_nowaklab_theme() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
          axis.line.y = element_blank(), # disable y axis lines
          axis.ticks.y = element_blank(), # disable y axis ticks lines
          axis.title.y = element_blank(), # disable y axis lines
          axis.text.y = element_blank()) # disable y axis text
  return(msa_cna_bc)
}
# The input must be a tibble with the following columns: asv_names, width_total_del, width_total_ins, width_total_sub
# This tibble is found in REvoBC_object$alignment$asv_alterations_width
plot_mutations_width = function(df_to_plot_final) {
  bar_ins_del_sub_width <- 
    ggplot(data=dplyr::select(df_to_plot_final, asv_names, width_total_del, width_total_ins, width_total_sub) %>% 
             unique() %>% 
             gather(key="alt", "width_total_del", "width_total_ins", "width_total_sub", value="width_sum") %>%
             filter(!str_detect(asv_names, "NMBC"))) +
    geom_bar(aes(x=width_sum, y=asv_names, fill=alt), position="stack", stat="identity", width=0.8, size=0.2) +
    scale_fill_manual(values=c("width_total_del" = "#3366FF", "width_total_ins" = "#FF0033", "width_total_sub" = "#329932"), 
                      breaks=c("width_total_del", "width_total_ins", "width_total_sub")) +
    geom_vline(xintercept=0, linetype="solid", size=0.5, col="black") + # no indels
    scale_x_continuous(labels=scales::comma, expand = c(0.01, 0.01), breaks=c(0, 130, 260))+#, limits=c(0, 286)) +
    geom_vline(xintercept=130, linetype="dotted", size=0.5, col="gray50") +
    geom_vline(xintercept=260, linetype="dotted", size=0.5, col="gray50") +
    xlab("Cummulative \n Widths of Alterations") +
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
             filter(!str_detect(asv_names, "NMBC"))) +
    geom_bar(aes(x=seq_n, y=asv_names, fill=seq_n_col), position = "dodge", stat = "identity", width=0.8, size=0.2) +
    geom_text(aes(x=1, y=asv_names, label = seq_n), col="white", size=4, hjust="left") + #, position=position_dodge(width=0.9)
    scale_fill_manual(values=c("260" = "#84B48F", "< 260" = "#377EB8", "> 260" = "#F39B7FFF"), breaks=c("260", "< 260", "> 260")) +
    geom_vline(xintercept=0, linetype="dotted", size=0.5, col="#377EB8") + # 0 
    geom_vline(xintercept=260, linetype="dotted", size=0.5, col="#84B48F") + # 260 bp = original length
    geom_vline(xintercept=520, linetype="dotted", size=0.5, col="#F39B7FFF") + # 0 
    scale_x_continuous(labels=scales::comma, expand = c(0.01, 0.01), limits=c(0, 572), breaks=c(0, 130, 260, 390, 520)) + # 572 = 1.1 * 520 -> nice separation between graphs
    xlab("Length of  \n ASV") +
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
             filter(!str_detect(asv_names, paste0(c("NMBC"), collapse = "|"))) %>%
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
# input must be a tibble with this columns: asv_names, perc_in_sample, perc_fold_to_max
plot_percentage_asv_sample = function(df_to_plot_final) {
  ### "Bubble Plot" - ASV % of colony  ------------------------------------------------------ 
  # bubble graph size coresponds to percentage of colony in i.e. Days or Organ
  scale_bubble <- ceiling(max(df_to_plot_final %>% 
                                dplyr::pull(perc_in_sample), na.rm = T))
  
  scale_bubble_nonmbc <- ceiling(max(df_to_plot_final %>% 
                        filter(str_detect(asv_names, "ASV")) %>% 
                        dplyr::pull(perc_in_sample), na.rm = T))
  
  df_to_plot_final = df_to_plot_final %>% filter(str_detect(asv_names, "ASV") | str_detect(asv_names, '.NMBC'))
  df_to_plot_final$asv_names = stringr::str_replace_all(string = df_to_plot_final$asv_names, pattern = '.NMBC', '')
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
    lemon::coord_capped_cart(bottom="both") + # axis with lemon
    # Add theme
    barplot_nowaklab_theme() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
          axis.line.y = element_blank(), # disable y axis lines
          axis.ticks.y = element_blank(), # disable y axis ticks lines
          axis.title.y = element_blank(), # disable y axis lines
          axis.text.y = element_blank(), # disable y axis text
          #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.grid.major.y = element_line(colour="#b4b4b4FF", size=0.4, "dotted"), # y grid line 
          panel.grid.major.x = element_line(colour="#b4b4b4FF", size=0.4, "dotted")) # disable lines in grid on X-axis
  return(bubble)
  
}
