context("EvoTraceR initialized")

data("revo_initialized")
test_that("EvoTraceR initialized produces correct output", {
    expect_equal(names(revo_initialized),c("fastq_directory","output_directory","dada2","map_file_sample","asv_prefilter"))
})
