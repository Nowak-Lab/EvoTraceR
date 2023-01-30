data("EvoTraceR_object")

phylogeny_output_path = "../../output/phylogeny_analysis/phylogeny_del_ins/"
asv_stats_file = "asv_stat.csv"
tree_clones_csv = "tree_all_clones.csv"
tree_clones_newick = "tree_all_clones.newick"

test_that("EvoTraceR vignette produces correct output", {
    announce_snapshot_file(name = asv_stats_file)
    expect_snapshot_file(paste(phylogeny_output_path, asv_stats_file, sep=''), name=asv_stats_file)
    announce_snapshot_file(name = tree_clones_csv)
    expect_snapshot_file(paste(phylogeny_output_path, tree_clones_csv, sep=''), name=tree_clones_csv)
    announce_snapshot_file(name = tree_clones_newick)
    expect_snapshot_file(paste(phylogeny_output_path, tree_clones_newick, sep=''), name=tree_clones_newick)
    expect_snapshot(EvoTraceR_object$preprocessing)
    expect_snapshot(EvoTraceR_object$alignment)
})