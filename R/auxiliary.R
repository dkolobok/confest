#character or factor to numeric preseving the names
ch2n <- function(x) {
  y <- as.numeric(as.character(x))
  names(y) <- names(x)
  y
}

# auxiliary function: optimization over all parameters but chosen

#par = parameter vector (length = d)
#fun = function to be optimized
#arg = index of profile argument (optimization is conducted over all parameters but this)
optim.fix <- function (par, fun, arg, cm, lower = NULL, upper = NULL,
                       control = NULL, ...) {
  fun.fix <- function(par.fix) {
    all_par <- par
    all_par[names(par.fix)] <- par.fix
    r <- fun(all_par)
    if (!is.finite(r)) return(1e6)
    r
  }
  if (is.null(lower)) lower <- cm$parbox[-arg, 1]
  if (is.null(upper)) upper <- cm$parbox[-arg, 2]
  if (is.null(names(lower))) names(lower) <- names(par[-arg])
  if (is.null(names(upper))) names(upper) <- names(par[-arg])
  if (is.null(control)) control <- list(maxit = 1e4, factr = 10, trace = 0,
                                        parscale = cm$scalevec[-arg])
  res <- optim(par = par[-arg],
               gr = function(x) grad(fun.fix, x, method = 'simple', method.args = list(eps = 1e-8)),
               fn = fun.fix,
               method = 'L-BFGS-B', lower = lower, upper = upper,
               control = control)
  res
}
