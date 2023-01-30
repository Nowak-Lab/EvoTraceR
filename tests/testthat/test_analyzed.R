data("EvoTraceR_object")
test_that("EvoTraceR analyzed produces correct output", {
    expect_equal(names(EvoTraceR_object),c("fastq_directory","output_directory","preprocessing","map_file_sample","asv_prefilter","sample_order", "reference","alignment","clean_asv_dataframe","clean_asv_dataframe_nonnorm","statistics"))
})
