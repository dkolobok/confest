profile <- function(parname, cm, seq,
                    logflag = 0) {
  fn1 <- function(th) m2lL(th, cm)
  parfit <- as.numeric(cm$fit$par[parname])
  crit <- c(cm$fit$value + qchisq(0.995, df = 1))
  #  browser()
  seq <- unique(sort(c(seq, as.numeric(cm$fit$par[parname]))))
  par <- cm$fit$par
  midind <- which(seq %in% par[parname])
  parind <- which(names(par) == parname)
  profileseq <- rep(NA, length(seq))
  profileseq[midind] <- cm$fit$value
  pl_fit_list <- as.list(rep(NA, length(seq)))

  #right
  for (i in (midind + 1):length(seq)) {
    par[parind] <- seq[i]
    res <- try(optim.fix(par, fn1, parind, cm = cm))
    pl_fit_list[[i]] <- res
    if (class(res) == 'list') {
      par[-parind] <- res$par
      profileseq[i] <- res$value
      if (res$value > crit) break
    }
  }

  #left
#  browser()
  par <- cm$fit$par
  for (i in 1:(midind - 1)) {
    par[parind] <- seq[midind - i]
    res <- try(optim.fix(par, fn1, parind, cm = cm))
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
profile_vis <- function(parname, cm, plotfile = NULL) {
  pl <- lapply(parname, profile_vis_ind, cm)
  if (!is.null(plotfile)) {
    number_ticks <- function(n) {function(limits) pretty(limits, n)}
    pl <- lapply(pl, function(x) x + scale_x_continuous(breaks = number_ticks(3)))
    pdf(plotfile)
      grid.arrange(grobs = pl, cols = floor(sqrt(length(pl))))
    dev.off()
  }
  pl
}

profile_vis_ind <- function(parname, cm) {
  fn1 <- function(th) m2lL(th, cm)
  parfit <- as.numeric(cm$fit$par[parname])

  df <- data.frame(x = cm$profile[[parname]]$parseq, y = cm$profile[[parname]]$profileseq)
  df$y <- try(as.numeric(as.character(df$y)))
  keepind <- which(!is.na(df$y))
  df <- df[keepind,]
  profile <- approxfun(x = df$x, y = df$y)
  df$l <- cm$profile[[parname]]$projectseq[keepind]

  conf <- c(0.5, 0.8, 0.9, 0.95, 0.99)
  crit <- cm$fit$value + qchisq(conf, df = 1)

  #  browser()
  plot <- ggplot(data = df) +
    geom_line(aes(x = x, y = y)) +
    #    geom_line(aes(x = x, y = l), color = 'blue', linetype = "dotted") +
    geom_hline(data = data.frame(y = crit), aes(yintercept = y), color = 'red') +
#    geom_vline(data = data.frame(x = cm_true$fit$par[parname]), aes(xintercept = x), color = 'green') +
    geom_point(data = data.frame(x = parfit, y = cm$fit$value), aes(x = x, y = y), col = 'red') +
    xlab(parname) +
    ylab('-2lnL') +
    coord_cartesian(xlim = c(min(df$x), max(df$x)))

  plot
}
