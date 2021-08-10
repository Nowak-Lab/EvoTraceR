context("REvoBC msa")

data("revo_msa")
test_that("REvoBC msa produces correct output", {
    expect_equal(names(revo_msa),c("fastq_directory","output_directory","dada2","map_file_sample","dada2_asv_prefilter","barcode","clean_asv_dataframe","statistics","alignment"))
})
