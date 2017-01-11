#CIdf <- as.data.frame(Reduce(cbind, CI_list))
#names(CIdf) <- paste(rep(names(cm_becker$fit$par), each = 2),
#                     c('l', 'r'), sep = '_')
cm_becker$ci <- CI_list
CIdf <- t(sapply(CI_list, function(x) {
  #  browser()
  z <- c(rev(x[, 1]), x[, 2])
  names(z) <- c(paste(rev(rownames(x)), 'l', sep = '_'),
                paste(rownames(x), 'r', sep = '_'))
  z
})
)
BSCIdf <- t(sapply(cm_becker$ci_bs, function(x) {
  #  browser()
  z <- c(x[, 1], rev(x[, 2]))
  names(z) <- c(paste(rownames(x), 'l', sep = '_'),
                paste(rev(rownames(x)), 'r', sep = '_'))
  z
})
)

require(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
xtable(CIdf, display = rep('e', ncol(CIdf)))
xtable(BSCIdf, display = 'e')

#other output
xtable(cm_becker$data_expl_list[[1]])
xtable(cm_becker$data_ode, floating = TRUE, floating.environment = "sidewaystable")
