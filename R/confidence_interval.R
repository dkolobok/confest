CI <- function(parname, cm, conf = c(0.5, 0.8, 0.9, 0.95, 0.99), seq_length = 100) {
  range <- range(cm$profile[[parname]]$parseq[!is.na(cm$profile[[parname]]$profileseq)])
  prof <- approxfun(x = cm$profile[[parname]]$parseq, y = cm$profile[[parname]]$profileseq)
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
  #  names(ci_df)[ncol(ci_df)] <- as.character(conf[i])
  }
  names(ci_df) <- conf
  t(ci_df)
}
