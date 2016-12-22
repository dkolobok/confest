#' fits a model to data
#'
#' a wrapper over optim function
#' @param par initial vector of parameters
#' @param left left margin
#' @param right left margin
#' @param cm an object of class custommodel
#' @param method method for optim
#' @param control control list
#' @keywords likelihood
#' @export
#' @examples
#' m2lL(par, cm, ode_solver = cm_ode_c, ...)
#'
#'
#

fit <- function(par = NULL, left = NULL, right = NULL, cm,
                method = 'L-BFGS-B',
                control = list(maxit = 1e4, factr = 10, trace = 0)) {
  if (is.null(par)) par <- cm$fit_start
  if (is.null(left)) left <- cm$fit_start * 0.5
  if (is.null(right)) right <- cm$fit_start * 2
  fn1 <- function(th) m2lL(th, cm = cm)

  optim(par = par, fn = fn1,
         #                gr = function(x) grad(fn1, x, method = 'simple', method.args = list(eps = 1e-8)),
         method = method,
         lower = left, upper = right,
         control = control
  )
}
