#takes a data table with the same experiment and solves an ode for it
#parameters & states are considered equal for all data points

#' ODE solver
#'
#' Takes a data table with the same experiment and solves an ode for it. Parameters & states are considered equal for all data points
#' @param dat a data.frame with ode data
#' @param state a vector of initial states
#' @param parms a vector of parameters
#' @param parms a vector of times
#' @param cm an object of class custommodel
#' @param maxsteps maxsteps parameter, defaults to 1e6
#' @param hmax hmax parameter, defaults to 0
#' @keywords ode
#' @export
#' @examples
#' cm_ode_c(dat, state, parms, time, cm)
#'
#'

cm_ode_c <- function(dat, state, parms, time, cm, maxsteps = 1e6, hmax = 0, ...) {
  #  browser()
  #  sol <- ode(y = ch2n(state), time = seq(0, max(as.numeric(as.character(unlist(dat$time)))) * 1.1, length.out = 100), func = becker,
  #             parms = ch2n(parms))
  #rescaling
  oldstate <- state
  oldparms <- parms
  state <- ch2n(state)
  parms <- ch2n(parms)
  shared_names <- intersect(names(parms), names(cm$scalevec))
  parms[shared_names] <- parms[shared_names] * cm$scalevec[shared_names]
  snames <- names(cm$scalevec)[grep('init_', names(cm$scalevec))]
  state[gsub('init_', '', snames)] <- state[gsub('init_', '', snames)] * cm$scalevec[snames]

  #  sol <- ode(y = ch2n(state), time = seq(0, max(as.numeric(as.character(unlist(dat$time)))) * 1.1, length.out = 100), func = becker,
  #             parms = ch2n(parms))
  sol <- ode(y = state[cm$ode_names], times = seq(0, max(as.numeric(as.character(unlist(dat$time)))) * 1.1, length.out = 100),
             func = "derivs",
             parms = parms,
             jacfunc = "jac",
             jactype = 'fullusr',
             dllname = cm$cname,
             initfunc = "initmod", nout = length(cm$ode_names),
             maxsteps = maxsteps,
             hmax = hmax, ...#,
             #atol = 1e-10,
             #rtol = 0
             #atol = 1e-128
             #      , method = 'lsode'
             #      ,hmax = 0.02
  )
  #  browser()
  #  attr(sol, 'state') <- state
  #  attr(sol, 'parms') <- parms
  sol
}

cm_ode_c_nj <- function(dat, state, parms, time) {
  #  browser()
  #  sol <- ode(y = ch2n(state), time = seq(0, max(as.numeric(as.character(unlist(dat$time)))) * 1.1, length.out = 100), func = becker,
  #             parms = ch2n(parms))
  #rescaling
  oldstate <- state
  oldparms <- parms
  state <- ch2n(state)
  parms <- ch2n(parms)
  shared_names <- intersect(names(parms), names(cm$scalevec))
  parms[shared_names] <- parms[shared_names] * cm$scalevec[shared_names]
  snames <- names(cm$scalevec)[grep('init_', names(cm$scalevec))]
  state[gsub('init_', '', snames)] <- state[gsub('init_', '', snames)] * cm$scalevec[snames]

  #  sol <- ode(y = ch2n(state), time = seq(0, max(as.numeric(as.character(unlist(dat$time)))) * 1.1, length.out = 100), func = becker,
  #             parms = ch2n(parms))
  sol <- ode(y = state[cm$ode_state_names], times = seq(0, max(as.numeric(as.character(unlist(dat$time)))) * 1.1, length.out = 100),
             func = "derivs",
             parms = parms,
#             jacfunc = "jac",
             dllname = "becker2010",
             initfunc = "initmod", nout = 6,
             maxsteps = 1e6,
             hmax = 0,
             atol = 1e-128
             #      , method = 'lsode'
             #      ,hmax = 0.02
  )
  #  browser()
  #  attr(sol, 'state') <- state
  #  attr(sol, 'parms') <- parms
  sol
}

cm_ode_r <- function(dat, state, parms, time) {
  #  browser()
  #  sol <- ode(y = ch2n(state), time = seq(0, max(as.numeric(as.character(unlist(dat$time)))) * 1.1, length.out = 100), func = becker,
  #             parms = ch2n(parms))
  #rescaling
  oldstate <- state
  oldparms <- parms
  state <- ch2n(state)
  parms <- ch2n(parms)
  shared_names <- intersect(names(parms), names(cm$scalevec))
  parms[shared_names] <- parms[shared_names] * cm$scalevec[shared_names]
  snames <- names(cm$scalevec)[grep('init_', names(cm$scalevec))]
  state[gsub('init_', '', snames)] <- state[gsub('init_', '', snames)] * cm$scalevec[snames]

  #  sol <- ode(y = ch2n(state), time = seq(0, max(as.numeric(as.character(unlist(dat$time)))) * 1.1, length.out = 100), func = becker,
  #             parms = ch2n(parms))
  sol <- ode(y = state[cm$ode_state_names], times = seq(0, max(as.numeric(as.character(unlist(dat$time)))) * 1.1, length.out = 100),
             #func = "derivs",
             func = becker,
             parms = parms,
                         jacfunc = jac,
             jactype = 'fullusr',
             #             dllname = "becker2010",
             #      initfunc = "initmod", nout = 6,
             maxsteps = 1e6,
             hmax = 0,
             atol = 1e-128
             #      , method = 'lsode'
             #      ,hmax = 0.02
  )
  #  browser()
  #  attr(sol, 'state') <- state
  #  attr(sol, 'parms') <- parms
  sol
}

cm_ode_r_nj <- function(dat, state, parms, time) {
  #  browser()
  #  sol <- ode(y = ch2n(state), time = seq(0, max(as.numeric(as.character(unlist(dat$time)))) * 1.1, length.out = 100), func = becker,
  #             parms = ch2n(parms))
  #rescaling
  oldstate <- state
  oldparms <- parms
  state <- ch2n(state)
  parms <- ch2n(parms)
  shared_names <- intersect(names(parms), names(cm$scalevec))
  parms[shared_names] <- parms[shared_names] * cm$scalevec[shared_names]
  snames <- names(cm$scalevec)[grep('init_', names(cm$scalevec))]
  state[gsub('init_', '', snames)] <- state[gsub('init_', '', snames)] * cm$scalevec[snames]

  #  sol <- ode(y = ch2n(state), time = seq(0, max(as.numeric(as.character(unlist(dat$time)))) * 1.1, length.out = 100), func = becker,
  #             parms = ch2n(parms))
  sol <- ode(y = state[cm$ode_state_names], times = seq(0, max(as.numeric(as.character(unlist(dat$time)))) * 1.1, length.out = 100),
             #func = "derivs",
             func = becker,
             parms = parms,
             #            jacfunc = "jac",
             #             dllname = "becker2010",
             #      initfunc = "initmod", nout = 6,
             maxsteps = 1e6,
             hmax = 0,
             atol = 1e-128
             #      , method = 'lsode'
             #      ,hmax = 0.02
  )
  #  browser()
  #  attr(sol, 'state') <- state
  #  attr(sol, 'parms') <- parms
  sol
}
