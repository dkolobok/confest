#gets an ode solution matrix and returns an explicit functions df
ode_explicit <- function(sol) {
  #  browser()
  explicit <- data.frame(time = sol[, 'time'],
                         Epo_ext = sol[, 'Epo'] + sol[, 'dEpo_e'],
                         Epo_int = sol[, 'Epo_EpoR_i'] + sol[, 'dEpo_i']
  )
  explicit$Epo_ext_cpm <- ch2n(attr(sol, 'parms')['offset']) + ch2n(attr(sol, 'parms')['scale']) * explicit$Epo_ext
  explicit$Epo_mem_cpm <- ch2n(attr(sol, 'parms')['offset']) + ch2n(attr(sol, 'parms')['scale']) * sol[, 'Epo_EpoR']
  explicit$Epo_int_cpm <- ch2n(attr(sol, 'parms')['offset']) + ch2n(attr(sol, 'parms')['scale']) * explicit$Epo_int
  attr(explicit, 'parms') <- attr(sol, 'parms')
  attr(explicit, 'state') <- attr(sol, 'state')
  explicit
}


#######################
Epo_ext_cpm <- function(sol) {
  res <- ch2n(attr(sol, 'parms')['offset']) + ch2n(attr(sol, 'parms')['scale']) *
    (sol[, 'Epo'] + sol[, 'dEpo_e'])
  approxfun(x = sol[, 'time'], y = res)
}

Epo_mem_cpm <- function(sol) {
  res <- ch2n(attr(sol, 'parms')['offset']) + ch2n(attr(sol, 'parms')['scale']) * sol[, 'Epo_EpoR']
  approxfun(x = sol[, 'time'], y = res)
}

Epo_int_cpm <- function(sol) {
  res <- ch2n(attr(sol, 'parms')['offset']) + ch2n(attr(sol, 'parms')['scale']) *
    (sol[, 'Epo_EpoR_i'] + sol[, 'dEpo_i'])
  approxfun(x = sol[, 'time'], y = res)
}


