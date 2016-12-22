source('R/packages.R')
source('R/create_cm.R')
source('R/parbox.R')
source('R/likelihood.R')
source('R/ode_solver.R')
source('R/ode_explicit.R')
source('R/separate_data.R')
source('R/auxiliary.R')
source('R/fit.R')
source('R/implicit.R')
source('R/profile.R')
source('R/confidence_band.R')
source('R/design.R')
dyn.load(paste('data/becker2010', .Platform$dynlib.ext, sep = ""))

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

cm_becker <- create_cm(data_ode = 'data/data_ode_fit.xlsx',
                 data_expl_list = 'data/data_expl.xls',
                 cname = 'becker2010',
                 ode_names = c('EpoR', 'Epo', 'Epo_EpoR', 'Epo_EpoR_i', 'dEpo_i', 'dEpo_e'),
                 ode_parnames = c('kon', 'koff', 'kt', 'init_EpoR', 'kex', 'ke', 'kdi', 'kde'),
                 explfunlist = explfunlist,
                 implfunlist = implfunlist,
                 fit_start = fit_start)
load('scalevec.RData')
cm_becker$scalevec <- scalevec
cm_becker$parbox <- cm_ped_list[[1]]$parbox

load('cm_true.RData')
#cm_true$data_expl_list = 'data/data_expl.xls'
#cm_true$cname <- 'becker2010'
#cm_true$scalevec <- scalevec
#cm_true$ode_names = c('EpoR', 'Epo', 'Epo_EpoR', 'Epo_EpoR_i', 'dEpo_i', 'dEpo_e')
#cm_true$ode_parnames = c('kon', 'koff', 'kt', 'init_EpoR', 'kex', 'ke', 'kdi', 'kde')
#cm_true$explfun <- list(epo_binding)


m2lL(cm_true$fit$par, cm = cm_becker)

system.time(res1 <- fit(cm = cm_becker, par = cm_true$fit$par,
                        left = cm_true$fit$par * 0.2,
                        right = cm_true$fit$par * 5))
cm_becker$fit <- res1

implicit_fun(c(50, 100, 150), cm_true$fit$par, 'Epo_ext_cpm', cm_becker)

system.time(temp1 <- CB(xseq = 50, fname = 'Epo_ext_cpm', cm = cm_becker))
band_point(temp1[[1]], 0.95)



cm2 <- des1(x = 3, cm = cm_becker)


#function to optimize
fn1 <- function(x) crit1(des1(x, cm_becker))
system.time(res1 <- optim(par = mean(range(cm_becker$data_expl_list[[1]][, 1])), fn1,
      lower = min(cm_becker$data_expl_list[[1]][, 1]),
      upper = max(cm_becker$data_expl_list[[1]][, 1]),
      method = 'Brent'))
save(res1, file = 'design_res.RData')

#visualize
pdf('design.pdf')
ggplot(cm_becker$data_expl_list[[1]], aes(x = Log_Epo_free, y = Epo_bound)) +
  geom_point() +
  stat_function(fun = function(x) cm_becker$explfunlist[[1]](x, cm_becker$fit$par),
                geom = "line") +
  geom_point(data = data.frame(Log_Epo_free = res1$par,
                        Epo_bound = cm_becker$explfunlist[[1]](res1$par,
                                                               cm_becker$fit$par)),
             mapping = aes(x = Log_Epo_free, y = Epo_bound,
                 color = 'red')) +
  guides(color = FALSE)
dev.off()

temp1 <- fn1(mean(range(cm_becker$data_expl_list[[1]][, 1])))
temp2 <- fn1(min(range(cm_becker$data_expl_list[[1]][, 1])))
temp3 <- fn1(max(range(cm_becker$data_expl_list[[1]][, 1])))


############
# profiles #
############

cm_becker$parbox <- parbox(cm_becker)
cm_becker$parbox[1, 2] <-

system.time(cm_becker$profile <- mclapply(names(cm_becker$fit$par), function(x) {
  ind <- which(!is.na(temp1[[x]]$profileseq)) #take not na indeces
  if (min(ind) > 1) ind <- c(min(ind) - 1, ind) #add extra element to the left
  if (max(ind) < length(temp1[[x]]$profileseq)) ind <- c(ind, max(ind) + 1) #to the right
  seqran <- range(temp1[[x]]$parseq[ind])
  res <- try(profile(parname = x, cm = cm_becker,
                     seq = seq(seqran[1], seqran[2],
                     l = 21)))
  save(res, file = paste('RData/', x, '_profile.RData', sep = ''))
}))

cm_becker$profile <- list()
for (i in names(cm_becker$fit$par)) {
  load(paste('RData/', i, '_profile.RData', sep = ''))
  cm_becker$profile[[i]] <- res
}

temp1 <- profile_vis(names(cm_becker$fit$par), cm_becker, plotfile = 'profiles.pdf')




profseqlist <- lapply(cm_becker$profile, function(x) x$parseq)
profseqlist1 <- list()

profseqlist1$init_Epo <- seq(min(profseqlist$init_Epo), 1.3, l = 41)
profseqlist1$init_EpoR <- seq(min(profseqlist$init_Epo), max(profseqlist$init_Epo), l = 41)
profseqlist1$kdi <- seq(0.3, 1.7, l = 41)
profseqlist1$kex <- seq(-5, 7, l = 21)
profseqlist1$koff <- seq(0.01, 2.5, l = 20)
profseqlist1$offset <- seq(0.01, 2.5, l = 20)
profseqlist1$sd_Epo_bound <- seq(min(profseqlist$sd_Epo_bound), 2.5, l = 20)
profseqlist1$sd_Epo_ext <- seq(min(profseqlist$sd_Epo_ext), 2.5, l = 20)
profseqlist1$sd_Epo_mem <- seq(min(profseqlist$sd_Epo_mem), 2.5, l = 20)
profseqlist1$sd_Epo_int <- seq(min(profseqlist$sd_Epo_int), 2.5, l = 20)


system.time(temp3 <- mclapply(names(profseqlist1), function(x) {
  res <- try(profile(parname = x, cm = cm_becker,
                     seq = profseqlist1[[x]]))
  save(res, file = paste('RData2/', x, '_profile.RData', sep = ''))
}))
system.time(temp3 <- lapply(names(profseqlist1)[4], function(x) {
  res <- try(profile(parname = x, cm = cm_becker,
                     seq = profseqlist1[[x]]))
  save(res, file = paste('RData2/', x, '_profile.RData', sep = ''))
}))

for (i in names(profseqlist1)[4]) {
  load(paste('RData2/', i, '_profile.RData', sep = ''))
  cm_becker$profile[[i]] <- res
}

temp2 <- profile_vis(names(cm_becker$fit$par), cm_becker, plotfile = 'profiles_new.pdf')

########################
### Confidence bands ###
########################
t1 <- system.time(band_Epo_ext <- mclapply(seq(0, 300, by = 10),
                                           function(x)
                                             try(CB_point(x = x, fname = 'Epo_ext_cpm', cm = cm_becker))))
names(band_Epo_ext) <- seq(0, 300, by = 10)
b1 <- band_constructor(band_Epo_ext, 0.95)
