source('R/packages.R')
source('R/create_cm.R')
source('R/likelihood.R')
source('R/ode_solver.R')
source('R/ode_explicit.R')
source('R/separate_data.R')
source('R/auxiliary.R')
source('R/fit.R')
source('R/implicit.R')
source('R/confidence_band.R')

fit_start <- 10^(c(3.12952554975111, 2.695617, 2.15189998283210, -1.92124706391398, -2.89528621922769, -1.25601966555971, -3.23806881827911, -1.09320136599830, -0.820036349822232, -1.79554039332794, -4.99999999999996, -0.00884943952210257, -1.40070145472325, -2.07455830963986, -1.25688223082166, -1.32032555968201))
names(fit_start) <- c('init_Epo', 'init_EpoR', 'kD', 'kde', 'kdi', 'ke', 'kex', 'koff', 'kon', 'kt', 'offset', 'scale', 'sd_Epo_bound', 'sd_Epo_ext', 'sd_Epo_int', 'sd_Epo_mem')
#fit_start['kD'] <- fit_start['koff'] / fit_start['kon']
fit_start['kon'] <- fit_start['kon'] / fit_start['init_Epo']
fit_start['scale'] <- fit_start['scale'] / fit_start['init_Epo']
fit_start['sd_Epo_bound'] <- fit_start['sd_Epo_bound'] * 1500

#initialize explicit functions

epo_binding <- function(x, par) {
  arg <- par * scalevec
  for (i in 1:length(arg)) assign(names(arg)[i], as.numeric(arg[i]))
  init_EpoR / 4 * (10 ^ x) / (10 ^ kD + (10 ^ x))
}

explfunlist <- list(epo_binding = epo_binding)

#initialize implicit functions

Epo_ext_cpm <- function(sol, cm) {
  res <- ch2n(attr(sol, 'parms')['offset']) * cm$scalevec['offset'] +
    ch2n(attr(sol, 'parms')['scale']) * cm$scalevec['scale'] *
    (sol[, 'Epo'] + sol[, 'dEpo_e'])
  approxfun(x = sol[, 'time'], y = res)
}

Epo_mem_cpm <- function(sol, cm) {
  res <- ch2n(attr(sol, 'parms')['offset']) * cm$scalevec['offset'] +
    ch2n(attr(sol, 'parms')['scale']) * cm$scalevec['scale'] *
    sol[, 'Epo_EpoR']
  approxfun(x = sol[, 'time'], y = res)
}

Epo_int_cpm <- function(sol, cm) {
  res <- ch2n(attr(sol, 'parms')['offset']) * cm$scalevec['offset'] +
    ch2n(attr(sol, 'parms')['scale']) * cm$scalevec['scale'] *
    (sol[, 'Epo_EpoR_i'] + sol[, 'dEpo_i'])
  approxfun(x = sol[, 'time'], y = res)
}

implfunlist <- list(Epo_ext_cpm = Epo_ext_cpm, Epo_mem_cpm = Epo_mem_cpm, Epo_int_cpm = Epo_int_cpm)

#create model

cm1 <- create_cm(data_ode = 'data/data_ode_fit.xlsx',
                 data_expl_list = 'data/data_expl.xls',
                 cname = 'becker2010',
                 ode_names = c('EpoR', 'Epo', 'Epo_EpoR', 'Epo_EpoR_i', 'dEpo_i', 'dEpo_e'),
                 ode_parnames = c('kon', 'koff', 'kt', 'init_EpoR', 'kex', 'ke', 'kdi', 'kde'),
                 explfunlist = explfunlist,
                 implfunlist = implfunlist,
                 fit_start = fit_start)
load('scalevec.RData')
cm1$scalevec <- scalevec

load('cm_true.RData')
#cm_true$data_expl_list = 'data/data_expl.xls'
#cm_true$cname <- 'becker2010'
#cm_true$scalevec <- scalevec
#cm_true$ode_names = c('EpoR', 'Epo', 'Epo_EpoR', 'Epo_EpoR_i', 'dEpo_i', 'dEpo_e')
#cm_true$ode_parnames = c('kon', 'koff', 'kt', 'init_EpoR', 'kex', 'ke', 'kdi', 'kde')
#cm_true$explfun <- list(epo_binding)


m2lL(cm_true$fit$par, cm = cm1)

system.time(res1 <- fit(cm = cm1, par = cm_true$fit$par,
                        left = cm_true$fit$par * 0.2,
                        right = cm_true$fit$par * 5))
cm1$fit <- res1

implicit_fun(c(50, 100, 150), cm_true$fit$par, 'Epo_ext_cpm', cm1)

system.time(temp1 <- CB(xseq = 50, fname = 'Epo_ext_cpm', cm = cm1))
