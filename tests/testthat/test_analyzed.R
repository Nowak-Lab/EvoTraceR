context("REvoBC analyzed")

data("revo_analyzed")
test_that("REvoBC analyzed produces correct output", {
    expect_equal(names(revo_analyzed),c("fastq_directory","output_directory","dada2","map_file_sample","dada2_asv_prefilter","reference","clean_asv_dataframe","statistics","alignment"))
})
