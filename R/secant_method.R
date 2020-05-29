#' Secant method
#'
#' The secant method is an iterative root-finding method with order of
#' convergence phi (1.618...) that does not use derivatives.
#'
#' @details The secant method finds the root of a univariate function \eqn{f}
#' given two initial guesses \eqn{x_0} and \eqn{x_1} by the
#' iteration:
#' \deqn{x_{n + 1} = x_n - f(x_n) \frac{x_n - x_{n - 1}}{f(x_n) - f(x_{n - 1})}}.
#'
#' The algorithm terminates when:
#'   \itemize{
#'     \item the algorithm exceeds 1000 iterations,
#'     \item the value of \eqn{f} is non-finite for an iterate (\eqn{x_n}),
#'     \item the iterate (\eqn{x_n}) becomes non-finite, or
#'     \item the algorithm converges and \eqn{|f(x_{n}) - f(x_{n + 1})| < tol}.
#'   }
#'
#' @param f Univariate function to find root of
#'
#' @param x0 A point close to the root of \code{f}
#'
#' @param x1 A point close to the root of \code{f} not equal to \code{x0}
#'
#' @param tol Tolerance for convergence.
#'
#' @return A root of \code{f}. If the algorithm does not converge, \code{NA}
#'   is returned.
#'
#' @examples
#' secant_method(cos, 0, 1)
#' secant_method(function(x) x ^ 3 - x - 2, 1, 2)
#'
#' @export
secant_method <- function(f, x0, x1, tol = 1e-8) {

  f0 <- f(x0)
  f1 <- f(x1)

  for (i in 1:1000) {

    if (abs(f1 - f0) < tol) {
      return(x1)
    }

    x2 <- x1 - f1 * (x1 - x0) / (f1 - f0)

    x0 <- x1
    x1 <- x2

    f0 <- f1
    f1 <- f(x1)

  }

  if (abs(f1 - f0) > tol) {
    return(NA_real_)
  }

  x1

}
