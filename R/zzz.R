.onLoad <- function (libname, pkgname)
{
  utils::globalVariables("where")
  utils::globalVariables(c("ASV_Richness_day_organ", "ASV_Richness_sample", "H_ShannonsIndex", "alt",
                           "asv_names", "asv_total_freq", "cons_bin", "day_organ", "name", "perc_in_sample",
                           "position", "position_bc260", "read_asv", "ref_asv", "seq_names", "width",
                           "alt_bin", "condition", "seq_end", "seq_n", "seq_start", "sum_perc", "total_counts"))
}
