#' Halley's method
#'
#' Halley's method is an iterative root-finding method with cubic
#' convergence that requires the first and second derivative.
#'
#' @details Halley's method finds the root of a univariate function \eqn{f}
#' with first derivative \eqn{f'} and second derivative \eqn{f''} given an
#' initial guess \eqn{x_0} by the iteration:
#' \deqn{x_{n + 1} = x_n - \frac{2f(x_n)f'(x_n)}{2(f'(x_n))^2 - f(x_n)f''(x_n)}}.
#'
#' The algorithm terminates when:
#'   \itemize{
#'     \item the algorithm exceeds 1000 iterations,
#'     \item the value of \eqn{f}, \eqn{f'}, or \eqn{f''} is non-finite for an
#'     iterate (\eqn{x_n}),
#'     \item the iterate (\eqn{x_n}) becomes non-finite, or
#'     \item the algorithm converges and \eqn{|f(x_{n}) - f(x_{n + 1})| < tol}.
#'   }
#'
#' @param f Univariate function to find root of
#'
#' @param fp First derivative of \code{f}
#'
#' @param fpp Second derivative of \code{f}
#'
#' @param x0 A point close to the root of \code{f}
#'
#' @param tol Tolerance for convergence.
#'
#' @return A root of \code{f} near \code{x0}. If the algorithm does not
#'   converge, \code{NA} is returned.
#'
#' @export
halleys_method <- function(f, fp, fpp, x0, tol = 1e-8) {

  xnp1 <- x0

  fxn <- f(x0)
  fpxn <- fp(x0)
  fppxn <- fpp(x0)

  for (i in 1:1000) {

    xn <- xnp1

    xnp1 <- xn - (2 * fxn * fpxn) / (2 * (fpxn) ^ 2 - fxn * fppxn)

    fxnp1 <- fxn
    fxn <- f(xnp1)
    fpxn <- fp(xnp1)
    fppxn <- fpp(xnp1)

    if (abs(fxnp1 - fxn) < tol) {
      return(xnp1)
    }

  }

  if (abs(fxnp1 - fxn) > tol) {
    return(NA_real_)
  }

  return(xnp1)

}
