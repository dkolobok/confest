#determines parbox
parbox <- function(cm) {
  fn1 <- function(th) m2lL(th, cm)
  df <- data.frame(l = cm$fit$par, r = cm$fit$par)
  for (i in 1:length(cm$fit$par)) {
    paru <- cm$fit$par
    parl <- cm$fit$par
    ifelse(cm$logflag[i], {
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
      if (!succu) ifelse(cm$logflag[i],
                         paru[i] <- paru[i] * 0.9,
                         paru[i] <- paru[i] - cm$fit$par[i] * 0.1)
      if (!succl) ifelse(cm$logflag[i],
                         parl[i] <- parl[i] * 1.1,
                         parl[i] <- parl[i] + cm$fit$par[i] * 0.1)
    }
    df[i,] <- c(parl[i], paru[i])
  }
  df
}
