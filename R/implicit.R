#' fits a model to data
#'
#' a wrapper over optim function
#' @param x an independent variable (one point or a sequence)
#' @param par parameter vector
#' @param fname implicit function name (should be in cm$implfunlist)
#' @param cm a custommodel object
#' @keywords implicit
#' @export
#' @examples
#' implicit_fun(x, par, fname, cm)
#'
#'
#

implicit_fun <- function(x, par, fname, cm) {
  # ode experimental loop (process each experiment =
  # points with the same experimental conditions and parameter values)
  #  browser()
  if (is.null(names(par))) names(par) <- names(cm$fit$par)
  l <- list()
  eps <- list() #residuals

  #substitute parameters
  tab <- cm$data_ode
  tab <- as.data.frame(apply(tab, c(1, 2), function(x) ifelse(x %in% names(par), return(as.character(par[names(par) == x])), return(x))))
  sep_data <- tab %>% separate_data

  # conditions for the ith ode and explicit list element and
  # for ith sep_data list element are exactly the same
  assigned <- assign_all(sep_data[[1]], cm = cm)
  xall <- c(assigned$state, assigned$parms)
  names(xall) <- c(paste('init', names(assigned$state), sep = '_'), names(assigned$parms))
  xall[intersect(names(xall), names(par))] <- par[intersect(names(xall), names(par))]
  edstates <- assigned$state
  edstates[intersect(names(edstates), gsub('init_', '', names(par)))] <- par[paste('init_', intersect(names(edstates), gsub('init_', '', names(par))), sep = '')]
  ode <- cm_ode_c(sep_data[[1]], state = assigned$state, parms = xall[cm$ode_parnames], time = assigned$time, cm = cm)[, 1:(length(cm$ode_names) + 1)]
  attr(ode, 'state') <- assigned$state
  attr(ode, 'parms') <- xall
#  explicit_ode <- ode_explicit(ode)
  #browser()
#  if (fname %in% names(explicit_ode)) fn <- approxfun(explicit_ode$time, explicit_ode[, fname])
#  else {
#    explicit_extra <- data.frame(epo_seq, epo_bound = epo_binding(epo_seq, par))
#    fn <- approxfun(x = explicit_extra[, 1], y = explicit_extra[, fname])
#  }
  fn <- cm$implfunlist[[fname]](ode, cm)
  sapply(x, fn)
}

#q <- implicit_fun(x = 150, par = cm$fit$par, fname = 'Epo_ext_cpm')
