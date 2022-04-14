context("REvoBC msa")

data("revo_msa")
test_that("REvoBC msa produces correct output", {
    expect_equal(names(revo_msa),c("fastq_directory","output_directory","dada2","map_file_sample","asv_prefilter","reference","clean_asv_dataframe","statistics","alignment", "cleaned_deletions_insertions"))
})
