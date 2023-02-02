data("EvoTraceR_object")

test_that("EvoTraceR vignette produces expected starting ASVs", {
    expect_equal(EvoTraceR_object$preprocessing$seq_filters$num[1], 12241)
})


test_that("EvoTraceR vignette produces expected final ASVs", {
    expect_true(EvoTraceR_object$preprocessing$seq_filters$num[6] >= 707 & EvoTraceR_object$preprocessing$seq_filters$num[6] <= 709)
})