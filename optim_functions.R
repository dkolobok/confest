profile2D <- function(v) #v = new vector, par.ind = their indices, theta.opt = MLE
{
  #browser()
  theta.help <- theta.opt
  theta.help[par.ind[1]] <- v[1] #theta.help = start fitting value
  theta.help[par.ind[2]] <- v[2] #theta.help = start fitting value
  #optimization over all parameters but par.ind-th
  return(optim.fix(par = theta.help, fun = RSS1, arg = par.ind)$value)
}

project2D <- function(v) #v = new vector, par.ind = their indices, theta.opt = MLE
{
  theta.help <- theta.opt
  theta.help[par.ind[1]] <- v[1] #theta.help = start fitting value
  theta.help[par.ind[2]] <- v[2] #theta.help = start fitting value
  #optimization over all parameters but par.ind-th
  return(RSS1(theta.help))
}

#optimization function
RSS <- function (k, dots, weights=1,fun=f)
  sum(weights*(fun(dots$x,k)-dots$y)^2)

#auxiliary
RSS1<- function (k) RSS(k, dots = dots, fun=f)

# auxiliary function: optimization over all parameters but chosen

#par = parameter vector (length = d)
#fun = function to be optimized
#arg = index of profile argument (optimization is conducted over all parameters but this)
optim.fix <- function (par, fun, arg, cm, lower = NULL, upper = NULL,
                       control = NULL, ...) {
  fun.fix <- function(par.fix) {
    all_par <- par
    all_par[names(par.fix)] <- par.fix
    return(fun(all_par))
  }
  if (is.null(lower)) lower <- par[-arg] * 0.2
  if (is.null(upper)) upper <- par[-arg] * 5
  if (is.null(control)) control <- list(maxit = 1e4, factr = 1e12, trace = 0,
       parscale = cm$scalevec[-arg])
  res <- optim(par = par[-arg],
               gr = function(x) grad(fun.fix, x, method = 'simple', method.args = list(eps = 1e-8)),
               fn = fun.fix,
               method = 'L-BFGS-B', lower = lower, upper = upper,
               control = control)
  res
}



