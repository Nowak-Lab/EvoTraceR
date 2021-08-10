Sys.setenv("R_TESTS" = "")

library("testthat")
library("REvoBC")

test_check("REvoBC")
