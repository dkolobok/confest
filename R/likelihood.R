#' -2lnL
#'
#' gets parameter and state values, return -2lnL
#' @param par a vector of parameters
#' @param cm an object of class custommodel
#' @param ode_solver ode solver function (defaults to cm_ode_c)
#' @keywords likelihood
#' @export
#' @examples
#' m2lL(par, cm, ode_solver = cm_ode_c, ...)
#'
#'
#
m2lL <- function(par, cm, ode_solver = cm_ode_c, ...) {
  # ode experimental loop (process each experiment =
  # points with the same experimental conditions and parameter values)
#  browser()
  if (class(par) == 'matrix') par <- par[, 1]
  if (is.null(names(par))) names(par) <- names(cm$fit_start)
    l <- list()
#  ep <- list()
#  enoise <- list()
  eps <- list() #residuals

  #substitute parameters
  tab <- cm$data_ode
  tab <- as.data.frame(apply(tab, c(1, 2), function(x) ifelse(x %in% names(par), return(as.character(par[names(par) == x])), return(x))))
  sep_data <- tab %>% separate_data

  for (i in 1:length(sep_data)) {
    # conditions for the ith ode and explicit list element and
    # for ith sep_data list element are exactly the same
    assigned <- assign_all(sep_data[[i]], cm)
#    ode <- cm_ode(sep_data[[i]], state = assigned$state, parms = assigned$parms, time = assigned$time)
    xall <- c(assigned$state, assigned$parms)
    names(xall) <- c(paste('init', names(assigned$state), sep = '_'), names(assigned$parms))
    ode <- ode_solver(sep_data[[i]], state = assigned$state, parms = xall[cm$ode_parnames], time = assigned$time, cm = cm)[, 1:(length(cm$ode_names) + 1)]

    state <- ch2n(assigned$state)
    parms <- ch2n(assigned$parms)
    shared_names <- intersect(names(parms), names(cm$scalevec))
    parms[shared_names] <- parms[shared_names] * cm$scalevec[shared_names]
    snames <- names(cm$scalevec)[grep('init_', names(cm$scalevec))]
    state[gsub('init_', '', snames)] <- state[gsub('init_', '', snames)] * cm$scalevec[snames]

    attr(ode, 'state') <- state
    attr(ode, 'parms') <- parms
    explicit_ode <- ode_explicit(ode)

    l[[i]] <- rep(NA, nrow(sep_data[[i]]))
#    ep[[i]] <- rep(NA, nrow(sep_data[[i]]))
#    enoise[[i]] <- rep(NA, nrow(sep_data[[i]]))
    eps[[i]] <- rep(NA, nrow(sep_data[[i]]))

    for (j in 1:nrow(sep_data[[i]])) {
#      browser()
      af <- approxfun(x = explicit_ode$time, y = explicit_ode[, as.character(sep_data[[i]][j, 'yname'])])
      eps[[i]][j] <- af(as.numeric(as.character(sep_data[[i]][j, 'time']))) - as.numeric(as.character(sep_data[[i]][j, 'yvalue']))
#      browser()
      yn <- gsub('_cpm', '', sep_data[[i]][j,]$yname)
      yn <- paste('sd_', yn, sep = '')
      noise <- as.numeric(as.character(sep_data[[i]][j,]$noise)) * cm$scalevec[yn]
      l[[i]][j] <- eval_likelihood(eps[[i]][j], noise)
#            ep[[i]][j] <- eval_pi(eps[[i]][j], as.numeric(as.character(sep_data[[i]][j,]$noise)))
#            enoise[[i]][j] <- eval_noise(eps[[i]][j], as.numeric(as.character(sep_data[[i]][j,]$noise)))
    }
  }
  res <- sum(sapply(l, sum))
#  res_pi <- sum(sapply(ep, sum))
#  res_noise <- sum(sapply(enoise, sum))

#  browser()
  l_extra <- list()
  eps_extra <- list()

  for (i in 1:length(cm$data_expl_list)) {
    l_extra[[i]] <- c(NA)
    eps_extra[[i]] <- c(NA)
    range <- range(cm$data_expl_list[[i]][, 1])
    xseq <- seq(range[1], range[2], length.out = 100)
    explicit_extra <- data.frame(xseq, yseq = cm$explfunlist[[i]](xseq, par))
    af <- approxfun(x = explicit_extra$xseq, y = explicit_extra$yseq)
    for (j in 1:nrow(cm$data_expl_list[[i]])) {
      eps_extra[[i]][j] <- af(cm$data_expl_list[[i]][j, 1]) - cm$data_expl_list[[i]][j, 2]
      #          browser()
      if (is.na(suppressWarnings(sd <- as.numeric(as.character(cm$data_expl_list[[i]][j, 'noise']))) * cm$scalevec['sd_Epo_bound'])) sd <- as.numeric(par[as.character(cm$data_expl_list[[i]][j, 'noise'])]) * cm$scalevec['sd_Epo_bound']
      l_extra[[i]][j] <- eval_likelihood(eps_extra[[i]][j], sd)
      #    ep_extra[j] <- eval_pi(eps_extra[j], sd)
      #    enoise_extra[j] <- eval_noise(eps_extra[j], sd)
    }

  }

  res_extra <- sum(sapply(l_extra, sum))
#  res_pi_extra <- sum(ep_extra)
#  res_enoise_extra <- sum(enoise_extra)


#browser()
  2 * (res + res_extra)
#  2 * (res_pi + res_pi_extra)
#  2 * (res_noise + res_enoise_extra)
}

#gets data table, extracts character inputs (to be fitted)
get_char <- function(dat) {
  dat <- dat %>% select(-time, -yvalue, -yname)
  dat <- as.character(unlist(dat))
  res <- c()
  for (i in 1:length(dat)) if (is.na(try(as.numeric(dat[i])))) res <- c(res, dat[i])
  unique(res)
}


#sigma1 = 1
#eval(parse(text = l[[1]][1]))

##################
#takes fitted parameter name and sequence
#returns profile (as a function)

profilefun <- function(parname, seq) {
  par <- cm$fit$par
  parind <- which(names(par) == parname)
  profileseq <- sapply(seq, function(x) {
    par[parind] <- x
    try(optim.fix(par, m2lL, parind, control = list(abstol = 1e-16))$value)
  })
  approxfun(x = seq, y = profileseq)
}

PL <- function(x, parname) {
  par <- cm$fit$par
  parind <- which(names(par) == parname)
  par[parind] <- x
  try(optim.fix(par, m2lL, parind, control = list(abstol = 1e-16))$value)
}

PLexpl <- function(x, q, fname) {
#  feq <- function(x) fexpl(x) - q
  p <- length(cm$fit$par)
  heq <- function(theta) explicit_fun(x, theta, fname) - q
  lL <- function (k) m2lL(k[1:p]) + k[p + 1] * heq(k[1:p])
  lL1 <- function (k) explicit_fun(x, k[1:p], fname) - q
  grad_exp <- function(fn, x, ...) grad(fn, c(exp(x[1:p]), x[p + 1]), ...)
  browser()
  system.time(res <- sane(par = c(log(cm$fit$par), 0),
                            fn = function (k) grad_exp(lL, k, method = 'simple',
                                                                  method.args = list(eps = 1e-8)),
                            method=2, control=list(tol = 0.002)))

  system.time(res <- nleqslv(x = c(log(cm$fit$par), 0),
                             fn = function (k) grad_exp(lL, k, method = 'simple',
                                                    method.args = list(eps = 1e-8)),
                             control = list(ftol = 5e-3)))
  res$x[1:p] <- exp(res$x)[1:p]
  #  system.time(res <- cobyla(x = cm$fit$par, fn = m2lL,
#                            #upper = cm$fit$par*1.1, lower = cm$fit$par*0.9,
#                hin = function(theta) -abs(heq(theta)),
#                control = list(maxeval = 1e5, ftol_abs = 1e-4)))

#  system.time(res <- alabama::auglag(par = cm$fit$par, fn = m2lL, heq = heq))

#  res <- nloptr::auglag(x0 = cm$fit$par, fn = m2lL, heq = heq)
  #solnp(pars = cm$fit$par, fun = m2lL, eqfun = feq)
  return(m2lL(res$x[1:p]))
}

error_list <- vector('list')
explicit_fun_debug <- function(x, k, fname) {
  tryCatch(explicit_fun(x, k, fname), error = function(e) error_list <<- list(error_list, k))
}

PLexpl2 <- function(x, q, fname) {
  feq <- function(k) explicit_fun_debug(x, k, fname) - q
  auglag(par = cm$fit$par, fn = m2lL_debug, heq = feq)
}

#gets residual and likelihood model
#returns likelihood value in a point
eval_likelihood <- function(eps, noise)
  -log(exp(-(eps) ^ 2 / 2 / (noise ^ 2)) / noise / sqrt(2 * pi))
#-log(exp(-(eps) ^ 2 / 2 / (noise ^ 2)) / 2 / noise / sqrt(2 * pi))

