context("REvoBC initialized")

data("revo_initialized")
test_that("REvoBC initialized produces correct output", {
    expect_equal(names(revo_initialized),c("fastq_directory","output_directory","dada2","map_file_sample","dada2_asv_prefilter"))
})
