#' Newtons's method
#'
#' Newtons's method is an iterative root-finding method with quadratic
#' convergence that requires the first derivative.
#'
#' @details Newtons's method finds the root of a univariate function \eqn{f}
#' with first derivative \eqn{f'} given an initial guess \eqn{x_0} by the
#' iteration:
#' \deqn{x_{n + 1} = x_n - \frac{f(x_n)}{f'(x_n)}}.
#'
#' The algorithm terminates when:
#'   \itemize{
#'     \item the algorithm exceeds 1000 iterations,
#'     \item the value of \eqn{f} or \eqn{f'} is non-finite for an iterate
#'     (\eqn{x_n}),
#'     \item the iterate (\eqn{x_n}) becomes non-finite, or
#'     \item the algorithm converges and \eqn{|f(x_{n}) - f(x_{n + 1})| < tol}.
#'   }
#'
#' @param f Univariate function to find root of
#'
#' @param fp First derivative of \code{f}
#'
#' @param x0 A point close to the root of \code{f}
#'
#' @param tol Tolerance for convergence.
#'
#' @return A root of \code{f} near \code{x0}. If the algorithm does not
#'   converge, \code{NA} is returned.
#'
#' @examples
#' newtons_method(cos,
#'                function(x) -sin(x),
#'                0.5)
#' newtons_method(function(x) x ^ 3 - x - 2,
#'                function(x) 3 * x ^ 2 - 1,
#'                1)
#'
#' @export
newtons_method <- function(f, fp, x0, tol = 1e-8) {

  xnp1 <- x0
  fxnp1 <- f(x0)
  fpxnp1 <- fp(x0)

  for (i in 1:1000) {
    xn <- xnp1

    fxn <- fxnp1
    fpxn <- fpxnp1

    xnp1 <- xn - fxn / fpxn

    fxnp1 <- f(xnp1)
    fpxnp1 <- fp(xnp1)

    if (abs(fxnp1 - fxn) < tol) {
      return(xnp1)
    }
  }

  if (abs(fxnp1 - fxn) > tol) {
    return(NA_real_)
  }

}
