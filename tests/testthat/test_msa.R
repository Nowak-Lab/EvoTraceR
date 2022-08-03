context("EvoTraceR msa")

data("revo_msa")
test_that("EvoTraceR msa produces correct output", {
    expect_equal(names(revo_msa),c("fastq_directory","output_directory","preprocessing","map_file_sample","asv_prefilter","reference","clean_asv_dataframe","statistics","alignment", "cleaned_deletions_insertions"))
})
