test_that("error for same sign", {
  x0 <- 1
  x1 <- 2
  f <- function(x) x
  expect_error(bisection_method(f, x0, x1))
})
