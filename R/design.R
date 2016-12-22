crit1 <- function(cm) {
  #refit
  system.time({
  cm$fit_start <- cm$fit$par
  cm$fit <- fit(cm = cm)

  #evaluate criterion
  temp1 <- CB(xseq = 50, fname = 'Epo_ext_cpm', cm = cm)
  temp2 <- band_point(temp1[[1]], 0.95)})
  as.numeric(temp2['r'] - temp2['l'])
}


# takes new experimental points, returns a new custommodel object
# f(x, theta_true) + eps is presumed equal to f(x, theta_mle)?
des1 <- function(x, cm) {
  newcm <- cm
  #ensure that the first two columns remain numeric
  newcm$data_expl_list[[1]][nrow(newcm$data_expl_list[[1]]) + 1,] <-
    c(x, newcm$explfunlist[[1]](x, par = newcm$fit$par), newcm$data_expl_list[[1]]$noise[1])
  #ensure that the first two columns remain numeric
  newcm$data_expl_list[[1]][, 1] <- as.numeric(newcm$data_expl_list[[1]][, 1])
  newcm$data_expl_list[[1]][, 2] <- as.numeric(newcm$data_expl_list[[1]][, 2])

  newcm
}
