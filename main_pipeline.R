# 1. fitting 

fit <- function(cm) {
#  cm1 <- ped_creator(cm$fit$par)
  fn1 <- function(th) m2lL(th, cm = cm)
  
#  kl <- c(rep(0.8, 12), rep(0.5, 4))
#  names(kl) <- names(cm$fit$par)
#  ku <- c(rep(1.2, 12), rep(5, 4))
#  names(ku) <- names(cm$fit$par)
#  old <- cm$fit
#  succu <- FALSE
#  succl <- FALSE
  
#  while (!(succu & succl)) {
#    succl <- is.numeric(try(fn1(cm$fit$par * kl)))
#    succu <- is.numeric(try(fn1(cm$fit$par * ku)))
#    if (!succu) ku[1:12] <- ku[1:12] - 0.1
#    if (!succl) kl[1:12] <- kl[1:12] + 0.1
#  }
  
  optimx(par = cm$fit$par, fn = fn1,
         #                gr = function(x) grad(fn1, x, method = 'simple', method.args = list(eps = 1e-8)),
         method = 'L-BFGS-B',
         lower = cm_true$fit$par * 0.2, upper = cm_true$fit$par * 5,
         control = list(maxit = 1e4, factr = 10, trace = 0,
                        parscale = sc1)
  )
}

fit2 <- function(cm) {
  fn1 <- function(th) m2lL(th, cm = cm)
  
  exc_l <- cm$fit$par - cm_true$fit$par * 0.2 < 1e-2
  exc_u <- -cm$fit$par + cm_true$fit$par * 5 < 1e-2
  if (all(exc_l == FALSE) & all(exc_u == FALSE)) return(cm$fit)
  left_m <- cm_true$fit$par * 0.2
  right_m <- cm_true$fit$par * 5
  left_m[exc_l] <- cm$parbox[exc_l, 'l']
  right_m[exc_u] <- cm$parbox[exc_u, 'r']
  optimx(par = cm$fit$par, fn = fn1,
         #                gr = function(x) grad(fn1, x, method = 'simple', method.args = list(eps = 1e-8)),
         method = 'L-BFGS-B',
         lower = left_m, upper = right_m,
         control = list(maxit = 1e4, factr = 10, trace = 0,
                        parscale = sc1)
  )
}

fit_bs <- function(cm_bs, start, cm) {
  cm1 <- cm
  cm1$data_ode <- cm_bs$data_ode
  cm1$data_expl <- cm_bs$data_expl
  
  fn1 <- function(th) m2lL(th, cm = cm1)
  
  res <- optimx(par = start, fn = fn1,
         #                gr = function(x) grad(fn1, x, method = 'simple', method.args = list(eps = 1e-8)),
         method = 'L-BFGS-B',
         lower = start * 0.1, upper = start * 10,
         control = list(maxit = 1e4, factr = 10, trace = 0,
                        parscale = sc1)
  )
  gr <- grad(fn1, res[1:16], method = 'simple', method.args = list(eps = 1e-8))
  list(res, gr)
}


# 2. visualization

#visualize an ode solution - returns a ggplot object
ode_gg <- function(sol) {
  l <- list()
  for (i in 1:(ncol(sol) - 1)) 
    l[[i]] <- ggplot(data = as.data.frame(sol), aes_string(x = 'time', y = colnames(sol)[i + 1])) + 
      geom_line()
  l
}

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

# 3. profiles
# надо считать профили с адаптивным шагом, не считать выше 99% инт.

profile <- function(parname, cm, seq_length = 20) {
  fn1 <- function(th) m2lL(th, cm)
  parfit <- as.numeric(cm$fit$par[parname])
  crit <- c(cm$fit$value + qchisq(0.99, df = 1))
  
  # MLE is a center of the profiling interval
  #  ifelse(parfit < (cm$parbox[parname, 1] + cm$parbox[parname, 2]) / 2,
  #         seq <- seq(cm$parbox[parname, 1], 2 * parfit - cm$parbox[parname, 1], length.out = seq_length),
  #         seq <- seq(2 * parfit - cm$parbox[parname, 2], cm$parbox[parname, 2], length.out = seq_length)) 
  
  #determine if log-scale should be used
  cent <- (cm$parbox[parname, 1] + cm$parbox[parname, 2]) / 2
  logcent <- 10 ^ ((log(cm$parbox[parname, 1], base = 10) + log(cm$parbox[parname, 2], base = 10)) / 2)
#  logflag <- ifelse(abs(logcent - as.numeric(cm$fit$par[parname])) < abs(cent - as.numeric(cm$fit$par[parname])), 1, 0)
#  if (is.na(logflag)) logflag <- 0
  
  ifelse(logflag[parname] == 1, 
         seq <- 10 ^ seq(log(cm$parbox[parname, 1], base = 10), log(cm$parbox[parname, 2], base = 10), length.out = seq_length),
         seq <- seq(cm$parbox[parname, 1], cm$parbox[parname, 2], length.out = seq_length)
  )
  #  browser()
  seq <- sort(c(seq, as.numeric(cm$fit$par[parname])))
  par <- cm$fit$par
  midind <- which(seq %in% par[parname])
  parind <- which(names(par) == parname)
  profileseq <- rep(NA, length(seq))
  profileseq[midind] <- cm$fit$value
  pl_fit_list <- as.list(rep(NA, length(seq)))
  
  #right
  for (i in (midind + 1):length(seq)) {
    par[parind] <- seq[i]
    res <- try(optim.fix(par, fn1, parind))
    pl_fit_list[[i]] <- res
    if (class(res) == 'list') {
      par[-parind] <- res$par
      profileseq[i] <- res$value
      if (res$value > crit) break
    }
  }
  
  #left
  par <- cm$fit$par
  for (i in 1:(midind - 1)) {
    par[parind] <- seq[midind - i]
    res <- try(optim.fix(par, fn1, parind))
    pl_fit_list[[midind - i]] <- res
    if (class(res) == 'list') {
      par[-parind] <- res$par
      profileseq[midind - i] <- res$value
      if (res$value > crit) break
    }
  }
  
  lseq <- as.numeric(unlist(sapply(seq, function(x) {
    par1 <- cm$fit$par
    par1[parname] <- x
    try(fn1(par1))})))
  
  list(parseq = seq, profileseq = profileseq, projectseq = lseq, pl_fit_list = pl_fit_list)
}

# profile visualization
profile_vis <- function(parname, cm) {
  fn1 <- function(th) m2lL(th, cm)
  parfit <- as.numeric(cm$fit$par[parname])
  
  df <- data.frame(x = cm$profiles[[parname]]$parseq, y = cm$profiles[[parname]]$profileseq)
  df$y <- try(as.numeric(as.character(df$y)))
  keepind <- which(!is.na(df$y))
  df <- df[keepind,]
  profile <- approxfun(x = df$x, y = df$y)
  df$l <- cm$profiles[[parname]]$projectseq[keepind]
  
  conf <- c(0.5, 0.8, 0.9, 0.95, 0.99)
  crit <- cm$fit$value + qchisq(conf, df = 1)
  
  #  browser()
  plot <- ggplot(data = df) +
    geom_line(aes(x = x, y = y)) +
    #    geom_line(aes(x = x, y = l), color = 'blue', linetype = "dotted") +
    geom_hline(data = data.frame(y = crit), aes(yintercept = y), color = 'red') +
    geom_vline(data = data.frame(x = cm_true$fit$par[parname]), aes(xintercept = x), color = 'green') +
    geom_point(data = data.frame(x = parfit, y = cm$fit$value), aes(x = x, y = y), col = 'red') +
    xlab(parname) +
    ylab('-2lnL') +
    coord_cartesian(xlim = c(min(df$x), max(df$x)))
  
  plot
}

# 4. confidence intervals
CI_cm <- function(parname, cm, conf = c(0.5, 0.8, 0.9, 0.95, 0.99), seq_length = 100) {
  range <- range(cm$profiles[[parname]]$parseq[!is.na(cm$profiles[[parname]]$profileseq)])
  prof <- approxfun(x = cm$profiles[[parname]]$parseq, y = cm$profiles[[parname]]$profileseq)
  #  crit <- cm$fit$value + qchisq(conf, df = 1)
  #  t2<-u.crit(theta_MLE[j],t1,crit)
  #  l.LRT[j,]<-t2[,1]
  #  r.LRT[j,]<-t2[,2]
  conf <- c(conf)
  crit <- c(cm$fit$value + qchisq(conf, df = 1))
  
  ci_df <- data.frame(emp = c(l = 1, r = 1))[, -1]
  for (i in 1:length(conf)) {
    l <- try(uniroot(function(x) prof(x) - crit[i], lower = range[1], upper = as.numeric(cm$fit$par[parname]), tol = 1e-64)$root)
    r <- try(uniroot(function(x) prof(x) - crit[i], lower = as.numeric(cm$fit$par[parname]), upper = range[2], tol = 1e-64)$root)
    ci_df <- data.frame(ci_df, v = c(as.numeric(l), as.numeric(r)))
    names(ci_df)[ncol(ci_df)] <- as.character(conf[i])
  }
  
  ci_df
}

CI_cm_mod <- function(parname, cm, conf = c(0.5, 0.8, 0.9, 0.95, 0.99), seq_length = 100) {
  range <- range(cm$profiles[[parname]]$parseq[!is.na(cm$profiles[[parname]]$profileseq)])
  prof <- approxfun(x = cm$profiles[[parname]]$parseq, y = cm$profiles[[parname]]$profileseq)
  #  crit <- cm$fit$value + qchisq(conf, df = 1)
  #  t2<-u.crit(theta_MLE[j],t1,crit)
  #  l.LRT[j,]<-t2[,1]
  #  r.LRT[j,]<-t2[,2]
  conf <- c(conf)
  n <- nrow(cm$data_expl) + nrow(cm$data_ode) 
  crit <- cm$fit$value + n / (n - length(cm$fit$par)) * qchisq(conf, df = 1)
  
  ci_df <- data.frame(emp = c(l = 1, r = 1))[, -1]
  for (i in 1:length(conf)) {
    l <- try(uniroot(function(x) prof(x) - crit[i], lower = range[1], upper = as.numeric(cm$fit$par[parname]), tol = 1e-64)$root)
    r <- try(uniroot(function(x) prof(x) - crit[i], lower = as.numeric(cm$fit$par[parname]), upper = range[2], tol = 1e-64)$root)
    ci_df <- data.frame(ci_df, v = c(as.numeric(l), as.numeric(r)))
    names(ci_df)[ncol(ci_df)] <- as.character(conf[i])
  }
  
  ci_df
}

#determines parbox
parbox <- function(cm) {
  fn1 <- function(th) m2lL(th, cm)
  df <- data.frame(l = cm$fit$par, r = cm$fit$par)
  for (i in 1:length(cm$fit$par)) {
    paru <- cm$fit$par
    parl <- cm$fit$par
    ifelse(logflag[i], {
      paru[i] <- cm$fit$par[i] * 10
      parl[i] <- cm$fit$par[i] * 0.1
    }, {
      paru[i] <- cm$fit$par[i] * 5
      parl[i] <- -cm$fit$par[i] * 4
    })
    
    succu <- FALSE
    succl <- FALSE
    #browser()
    while (!(succu & succl)) {
      succl <- is.finite(try(fn1(parl)))
      succu <- is.finite(try(fn1(paru)))
      if (!succu) ifelse(logflag[i], 
                         paru[i] <- paru[i] * 0.9,
                         paru[i] <- paru[i] - cm$fit$par[i] * 0.1)
      if (!succl) ifelse(logflag[i], 
                         parl[i] <- parl[i] * 1.1,
                         parl[i] <- parl[i] + cm$fit$par[i] * 0.1)
    }
    df[i,] <- c(parl[i], paru[i])
  }
  df
}

#############################
### 5. Confidence regions ###
#############################

CR_isin <- function(cm, conf = c(0.5, 0.8, 0.9, 0.95, 0.99)) {
  m2lL(cm_true$fit$par, cm) - cm$fit$value < qchisq(conf, df = length(cm$fit$par))
}

CR_mod_isin <- function(cm, conf = c(0.5, 0.8, 0.9, 0.95, 0.99)) {
  n <- nrow(cm$data_ode) + nrow(cm$data_expl)
  m2lL(cm_true$fit$par, cm) - cm$fit$value < n / (n - length(cm$fit$par)) * qchisq(conf, df = length(cm$fit$par))
}

##########################
### Visualize accuracy ###
##########################

# takes a data.frame returns a ggplot
ci_acc_vis <- function(df) {
  df1 <- data.frame(p = as.numeric(as.character(colnames(df))), pi = df['T',] / (df['T',] + df['F',]))
  ggplot(data = df1, mapping = aes(x = p, y = pi, group = 1)) +
    geom_point() +
    geom_segment(aes(x = 0.5, y = 0.5, xend = 0.99, yend = 0.99), col = 'green')
}


###########################
### PL confidence bands ###
###########################
CB <- function(xseq, fname, cm) {
  lapply(xseq, CB_point, fname, cm)
}

CB_point <- function(x, fname, cm) {
  q1 <- explicit_fun(x = x, par = cm$fit$par, fname = fname, cm = cm)
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
    heq <- function(theta) explicit_fun(x, theta, fname, cm = cm) - seq[i]
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
    heq <- function(theta) explicit_fun(x, theta, fname, cm = cm) - seq[midind - i]
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

CB_point_naive <- function(x, fname, cm, conf = 0.95) {
  fn1 <- function(th) m2lL(th, cm)

  hobj <- function(theta) explicit_fun(x, theta, fname, cm = cm)
  mhobj <- function(theta) (-1) * explicit_fun(x, theta, fname, cm = cm)
  
  resu <- list(pars = cm$fit$par)
  system.time(resu <- solnp(pars = cm$fit$par, 
                            fun = hobj,
                            ineqfun = fn1,
                            ineqUB = cm$fit$value + qchisq(conf, df = length(cm$fit$par)),
                            ineqLB = cm$fit$value,
                            LB = cm$fit$par - abs(cm$fit$par) * 0.2,
                            UB = cm$fit$par + abs(cm$fit$par) * 0.2
  ))
  resl <- list(pars = cm$fit$par)
  system.time(resl <- solnp(pars = cm$fit$par, 
                            fun = mhobj,
                            ineqfun = fn1,
                            ineqUB = cm$fit$value + qchisq(conf, df = length(cm$fit$par)),
                            ineqLB = cm$fit$value,
                            LB = cm$fit$par - abs(cm$fit$par) * 0.2,
                            UB = cm$fit$par + abs(cm$fit$par) * 0.2
  ))
  list(l = resl, u = resu)
}

band_point <- function(prof_res, conf) {
  prof <- approxfun(x = prof_res$xseq[!is.na(prof_res$plvalues)], y = prof_res$plvalues[!is.na(prof_res$plvalues)])
  crit <- min(prof_res$plvalues, na.rm = TRUE) + qchisq(conf, df = 1)
  
  l <- try(uniroot(function(x) prof(x) - crit, lower = min(prof_res$xseq[!is.na(prof_res$plvalues)]), 
                   upper = prof_res$xseq[which(prof_res$plvalues == min(prof_res$plvalues, na.rm = TRUE))], tol = 1e-64)$root)
  r <- try(uniroot(function(x) prof(x) - crit, upper = max(prof_res$xseq[!is.na(prof_res$plvalues)]), 
                   lower = prof_res$xseq[which(prof_res$plvalues == min(prof_res$plvalues, na.rm = TRUE))], tol = 1e-64)$root)
  c(l = l, r = r)
}

#check band accuracy
band_acc <- function(num = 50, xseq = c(50, 150, 250), expl_name) {
  ifelse(expl_name %in% names(cm_true$expl_fun$ode),
         df1 <- cm_true$expl_fun$ode[, c('time', expl_name)],
         df1 <- cm_true$expl_fun$extra)
  truefun <- approxfun(x = df1[, 1], y = df1[, 2])
  checklist <- list()
  for (k in c('0.5', '0.8', '0.9', '0.95', '0.99')) {
    checklist[[k]] <- matrix(NA, nrow = num, ncol = length(xseq))
    for (i in 1:num) {
      for (j in 1:length(xseq)) {
        tr <- explicit_fun(x = xseq[j], par = cm_true$fit$par, fname = expl_name, cm = cm_true)
        b <- band_constructor(cm_ped_list[[i]]$band[[expl_name]], as.numeric(k))
        checklist[[k]][i, j] <- tr > b['l', as.character(xseq[j])] &
          tr < b['r', as.character(xseq[j])]
      }
    }
    
  }
  checklist
}

#####################
### bootstrapping ###
#####################

#take a model and create a single non-parametric bootstrap models from it
bootstrap_model <- function(cm) {
  bs <- list()
  y1_ind <- which(cm$data_ode$yname == 'Epo_ext_cpm')
  y2_ind <- which(cm$data_ode$yname == 'Epo_mem_cpm')
  y3_ind <- which(cm$data_ode$yname == 'Epo_int_cpm')
  new_i <- c(sample(y1_ind, length(y1_ind), replace = TRUE), 
             sample(y2_ind, length(y2_ind), replace = TRUE), 
             sample(y3_ind, length(y3_ind), replace = TRUE))
  bs$data_ode <- cm$data_ode[new_i,]
  bs$data_expl <- cm$data_expl[sample(1:nrow(cm$data_expl), nrow(cm$data_expl), replace = TRUE),]
  bs
}

