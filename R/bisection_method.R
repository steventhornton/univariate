#' Bisection method
#'
#' The bisection method is an iterative root-finding method with linear
#' convergence that does not use derivatives.
#'
#' @details The bisection method finds the root of a univariate function \eqn{f}
#' given two initial guesses \eqn{x_0} and \eqn{x_1} where \eqn{x_0 < x_1} and
#' \eqn{sign(f(x_0)) = -sign(f(x_1))} by sequentially bisecting the interval
#' (\eqn{x_2 = (x_0 + x_1) / 2}).
#'
#' Convergence is guaranteed as long as f is continuous in the interval
#' \eqn{[x_0, x_1]} and \eqn{f(x_0)} and \eqn{f(x_1)} have opposite signs
#'
#'
#' The algorithm terminates when:
#'   \itemize{
#'     \item the algorithm exceeds 1000 iterations,
#'     \item the value of \eqn{f} is non-finite for an iterate (\eqn{x_n}),
#'     \item the algorithm converges and \eqn{|f(x_{n}) - f(x_{n + 1})| < tol}.
#'   }
#'
#' @param f Univariate function to find root of
#'
#' @param x0 A point
#'
#' @param x1 A point larger than \code{x0}
#'
#' @param tol Tolerance for convergence.
#'
#' @return A root of \code{f}. If the algorithm does not converge, \code{NA}
#'   is returned.
#'
#' @examples
#' bisection_method(cos, 0, pi)
#' bisection_method(function(x) x ^ 3 - x - 2, 1, 2)
#'
#' @export
bisection_method <- function(f, x0, x1, tol = 1e-8) {

  if (x0 >= x1) {
    stop("x0 must be less than x1")
  }

  if (sign(f(x0)) == sign(f(x1))) {
    stop("Function must be of opposite sign at initial points")
  }

  f0 <- f(x0)
  f1 <- f(x1)

  for (i in 1:1000) {

    xc <- (x1 + x0) / 2

    if (abs(f1 - f0) < tol) {
      return(xc)
    }

    fc <- f(xc)

    if (sign(fc) == sign(f0)) {
      x0 <- xc
      f0 <- fc
    } else {
      x1 <- xc
      f1 <- fc
    }

  }

  if (abs(f1 - f0) > tol) {
    return(NA_real_)
  }

  xc

}
