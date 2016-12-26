###########################
### PL confidence bands ###
###########################
CB <- function(xseq, fname, cm, qrelseq = c(0.5, 0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.2, 1.5)) {
  lapply(xseq, CB_point, fname, cm, qrelseq)
}

CB_point <- function(x, fname, cm, qrelseq = c(0.5, 0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.2, 1.5)) {
  if (fname %in% names(cm$implfunlist)) {
    q1 <- implicit_fun(x = x, par = cm$fit$par, fname = fname, cm = cm)
    fn <- function(theta) implicit_fun(x, theta, fname, cm = cm)
  }
  else {
    q1 <- cm$explfunlist[[fname]](x, cm$fit$par)
    fn <- function(theta) cm$explfunlist[[fname]](x, theta)
  }


  seq <- q1 * qrelseq
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
    heq <- function(theta) fn(theta) - seq[i]
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
    heq <- function(theta) fn(theta) - seq[midind - i]
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


band_point <- function(prof_res, conf) {
  prof <- approxfun(x = prof_res$xseq[!is.na(prof_res$plvalues)], y = prof_res$plvalues[!is.na(prof_res$plvalues)])
  crit <- min(prof_res$plvalues, na.rm = TRUE) + qchisq(conf, df = 1)

  l <- try(uniroot(function(x) prof(x) - crit, lower = min(prof_res$xseq[!is.na(prof_res$plvalues)]),
                   upper = prof_res$xseq[which(prof_res$plvalues == min(prof_res$plvalues, na.rm = TRUE))], tol = 1e-64)$root)
  r <- try(uniroot(function(x) prof(x) - crit, upper = max(prof_res$xseq[!is.na(prof_res$plvalues)]),
                   lower = prof_res$xseq[which(prof_res$plvalues == min(prof_res$plvalues, na.rm = TRUE))], tol = 1e-64)$root)
  if (class(l) != 'numeric') l <- NA
  if (class(r) != 'numeric') r <- NA
  c(l = l, r = r)
}

band_constructor <- function(band_profile_list, conf) {
  df1 <- sapply(band_profile_list, band_point, conf) #for each x
  data.frame(x = as.numeric(colnames(df1)), l = df1[1,], r = df1[2,])
}
