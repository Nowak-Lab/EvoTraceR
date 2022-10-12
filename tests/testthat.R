Sys.setenv("R_TESTS" = "")

library("testthat")
library("EvoTraceR")

test_check("EvoTraceR")
