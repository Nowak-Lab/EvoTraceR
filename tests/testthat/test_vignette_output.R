data("EvoTraceR_object")

test_that("EvoTraceR vignette produces expected starting ASVs", {
    expect_equal(EvoTraceR_object$preprocessing$seq_filters$num[1], 9434)
})


test_that("EvoTraceR vignette produces expected final ASVs", {
    expect_true(EvoTraceR_object$preprocessing$seq_filters$num[6] >= 19 & EvoTraceR_object$preprocessing$seq_filters$num[6] <= 21)
})
