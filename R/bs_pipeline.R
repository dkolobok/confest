#####################
### bootstrapping ###
#####################

#take a model and create a single non-parametric bootstrap models from it
bootstrap_model <- function(cm) {
  bs <- list()
  bs$fit_start <- cm$fit$par
  bs$ode_parnames <- cm$ode_parnames
  bs$ode_names <- cm$ode_names
  bs$scalevec <- cm$scalevec
  bs$explfunlist <- cm$explfunlist
  bs$cname <- cm$cname
  y1_ind <- which(cm$data_ode$yname == 'Epo_ext_cpm')
  y2_ind <- which(cm$data_ode$yname == 'Epo_mem_cpm')
  y3_ind <- which(cm$data_ode$yname == 'Epo_int_cpm')
  new_i <- c(sample(y1_ind, length(y1_ind), replace = TRUE),
             sample(y2_ind, length(y2_ind), replace = TRUE),
             sample(y3_ind, length(y3_ind), replace = TRUE))
  bs$data_ode <- cm$data_ode[new_i,]
  bs$data_expl_list <- list()
  for (i in 1:length(cm$data_expl_list))
    bs$data_expl_list[[i]] <- cm$data_expl_list[[i]][sample(1:nrow(cm$data_expl_list[[i]]), nrow(cm$data_expl_list[[i]]), replace = TRUE),]
  bs
}

fit_bs <- function(cm_bs, start, cm) {
#  cm1 <- cm
#  cm1$data_ode <- cm_bs$data_ode
#  cm1$data_expl <- cm_bs$data_expl

  fn1 <- function(th) m2lL(th, cm = cm_bs)

  res <- optim(par = start, fn = fn1,
                #                gr = function(x) grad(fn1, x, method = 'simple', method.args = list(eps = 1e-8)),
                method = 'L-BFGS-B',
                lower = start * 0.9, upper = start * 1.1,
                control = list(maxit = 1e4, factr = 10, trace = 0,
                               parscale = cm$scalevec)
  )
  gr <- grad(fn1, res$par, method = 'simple', method.args = list(eps = 1e-8))
  list(fit = res, gr = gr)
}

ci_bs <- function(parname, cm, p = c(0.99, 0.95, 0.9, 0.8, 0.5)) {
  v <- sort(sapply(cm$bs_models, function(bsm) bsm$fit$par[parname]))
  df <- data.frame(l = quantile(v, probs = (1 - p) / 2) * cm_becker$scalevec[parname],
                   r = quantile(v, probs = (1 + p) / 2) * cm_becker$scalevec[parname])
  rownames(df) <- p
  df
}

cb_bs <- function(fname, xseq, cm, p = c(0.99, 0.95, 0.9, 0.8, 0.5)) {
  lapply(xseq, function(x) {
    v <- sort(sapply(cm$bs_models, function(bsm) {
      if (fname %in% names(cm$implfunlist)) implicit_fun(x = x, par = bsm$fit$par, fname = fname, cm = cm)
      else cm$explfunlist[[fname]](x, bsm$fit$par)
    }))
    df <- data.frame(l = quantile(v, probs = (1 - p) / 2), r = quantile(v, probs = (1 + p) / 2))
    rownames(df) <- p
    df
  })
}

