test_that("error for same sign", {
  x0 <- 1
  x1 <- 2
  f <- function(x) x
  expect_error(bisection_method(f, x0, x1))
})

test_that("error for wrong endpoints", {
  x0 <- 1
  x1 <- -1
  f <- function(x) x
  expect_error(bisection_method(f, x0, x1))
})

test_that("sample rootfinding", {
  expect_equal(bisection_method(sin, -pi, pi), 0)
})
