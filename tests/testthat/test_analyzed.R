context("EvoTraceR analyzed")

data("revo_analyzed")
test_that("EvoTraceR analyzed produces correct output", {
    expect_equal(names(revo_analyzed),c("fastq_directory","output_directory","dada2","map_file_sample","asv_prefilter","reference","clean_asv_dataframe","clean_asv_dataframe_countnorm","statistics","alignment"))
})
