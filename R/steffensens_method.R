#' Steffensen's method
#'
#' Steffensen's method is an iterative root-finding method with quadratic
#' convergence that does not use derivatives.
#'
#' @details Steffensen's method finds the root of a univariate function \eqn{f}
#' given an initial guess \eqn{x_0} by the iteration:
#' \deqn{x_{n + 1} = x_n - \frac{f(x_n)}{g(x_n)}}
#' where
#' \deqn{g(x_n) = \frac{f(x_n + f(x_n))}{f(x_n)} - 1}.
#'
#' The algorithm terminates when:
#'   \itemize{
#'     \item the algorithm exceeds 1000 iterations,
#'     \item the value of \code{f} is non-finite for an iterate (\eqn{x_n}),
#'     \item the iterate (\eqn{x_n}) becomes non-finite, or
#'     \item the algorithm converges and \eqn{|f(x_{n}) - f(x_{n + 1})| < tol}.
#'   }
#'
#' @param f Univariate function to find root of
#'
#' @param x0 A point close to the root of \code{f}
#'
#' @param tol Tolerance for convergence.
#'
#' @return A root of \code{f} near \code{x0}. If the algorithm does not
#'   converge, \code{NA} is returned.
#'
#' @export
steffensens_method <- function(f, x0, tol = 1e-8) {

  xnp1 <- x0
  fxnp1 <- f(x0)

  for (i in 1:1000) {
    xn <- xnp1
    fxn <- fxnp1
    gxn <- f(xn + fxn) / fxn - 1
    xnp1 <- xn - fxn / gxn
    fxnp1 <- f(xnp1)

    if (abs(fxnp1 - fxn) < tol) {
      return(xnp1)
    }
  }

  if (abs(fxnp1 - fxn) > tol) {
    return(NA_real_)
  }

}
