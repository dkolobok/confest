###########################
### PL confidence bands ###
###########################
CB <- function(xseq, fname, cm) {
  lapply(xseq, CB_point, fname, cm)
}

CB_point <- function(x, fname, cm) {
  if (fname %in% names(cm$implfunlist))
    q1 <- implicit_fun(x = x, par = cm$fit$par, fname = fname, cm = cm) else
      q1 <- cm$explfunlist[[fname]](x, cm$fit$par)

  seq <- seq(q1 * 0.95, q1 * 1.05, length.out = 9)
  crit <- c(cm$fit$value + qchisq(0.99, df = 1))

  fn1 <- function(th) m2lL(th, cm)
  midind <- ceiling(length(seq) / 2)
  plvalues <- rep(NA, length(seq))
  heqvec <- rep(NA, length(seq))
  p <- length(cm$fit$par)
  plvalues[midind] <- fn1(cm$fit$par)
  heqvec[midind] <- 0
  res <- list(pars = cm$fit$par)
  resmid <- res #save the middle
  #res <- resmid

  # right half
  ## start from previouly found values
  #  browser()
  for (i in (midind + 1):length(seq)) {
    #    browser()
    heq <- function(theta) implicit_fun(x, theta, fname, cm = cm) - seq[i]
    system.time(res <- solnp(pars = res$pars,
                             fun = fn1,
                             eqfun = heq,
                             eqB = 0,
                             LB = res$pars - abs(cm$fit$par) * 0.2,
                             UB = res$pars + abs(cm$fit$par) * 0.2
    ))
    heqvec[i] <- heq(res$pars)
    plvalues[i] <- fn1(res$pars)
    if (class(plvalues[i]) != 'numeric') browser()
    #    temp <<- i
    if (plvalues[i] > crit) break
  }

  # left half
  res <- resmid
  ## start from previouly found values
  for (i in 1:(midind - 1)) {
    #    browser()
    heq <- function(theta) implicit_fun(x, theta, fname, cm = cm) - seq[midind - i]
    system.time(res <- solnp(pars = res$pars,
                             fun = fn1,
                             eqfun = heq,
                             eqB = 0,
                             LB = res$pars - abs(cm$fit$par) * 0.2,
                             UB = res$pars + abs(cm$fit$par) * 0.2
    ))
    heqvec[midind - i] <- heq(res$pars)
    plvalues[midind - i] <- fn1(res$pars)
    if (class(plvalues[midind - i]) != 'numeric') browser()
    if (plvalues[midind - i] > crit) break
  }

  return(list(plvalues = plvalues, heqvalues = heqvec, xseq = seq))
}
