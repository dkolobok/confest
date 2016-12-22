source('packages.R')
dyn.load(paste("becker2010", .Platform$dynlib.ext, sep = ""))

# create a list of ped models
cm_ped_list <- list()
for (i in 1:1000) {
  cm_ped_list[[i]] <- ped_creator(cm$fit$par)
  cm_ped_list[[i]]$m2lL_true <- m2lL(cm$fit$par, cm = cm_ped_list[[i]])
}

#library(foreach)
#library(doMC)
registerDoMC(4)
#t1 <- system.time(fit_list <- foreach (i = 1:1000) %dopar% try(fit(cm_ped_list[[i]])))
t1 <- system.time(fit_list <- mclapply(cm_ped_list, function(x) try(fit(x))))
#check fit distribution
fitvalues <- sapply(fit_list, '[[', 'value')
pdf('fitval_dist.pdf')
plot(density(fitvalues))#, xlim = range(fitvalues))
abline(v = cm$fit$value)
dev.off()

#############
### refit ###
#############
registerDoMC(4)
#t1 <- system.time(fit_list <- foreach (i = 1:1000) %dopar% try(fit(cm_ped_list[[i]])))
t2 <- system.time(fit_list2 <- mclapply(cm_ped_list[1:100], function(x) try(fit2(x))))


registerDoMC(4)
t00 <- system.time(fit_list00 <- foreach (i = 1:100) %dopar% try(fit2(cm_ped_list[[i]])))
#t1 <- system.time(fit_list <- foreach (i = 1:1000) %dopar% try(fit(cm_ped_list[[i]])))
t00 <- system.time(fit_list00 <- mclapply(cm_ped_list[1:4], function(x) try(fit2(x))))
registerDoMC(4)
#t1 <- system.time(fit_list <- foreach (i = 1:1000) %dopar% try(fit(cm_ped_list[[i]])))
t3 <- system.time(fit_list3 <- mclapply(cm_ped_list[101:1000], function(x) try(fit2(x))))
fit_list <- c(fit_list2, fit_list3)
for (i in 1:length(cm_ped_list)) cm_ped_list[[i]]$fit <- fit_list[[i]]

# add to cm_ped_list
for (i in 1:length(fit_list2)) {
  ifelse(length(fit_list2[[i]]) == 2, cm_ped_list[[i]]$fit <- fit_list2[[i]],{
    cm_ped_list[[i]]$fit <- list()
    cm_ped_list[[i]]$fit$par <- as.numeric(fit_list2[[i]][1:16])
    names(cm_ped_list[[i]]$fit$par) <- names(cm$fit$par)
    cm_ped_list[[i]]$fit$value <- as.numeric(fit_list2[[i]]['value'])
  })
}



t3 <- system.time(
  res2 <- optimx(par = cm$fit$par, fn = fn1,
                 #                gr = function(x) grad(fn1, x, method = 'simple', method.args = list(eps = 1e-8)),
                 method = 'L-BFGS-B',
                 lower = cm$fit$par * 0.2, upper = cm$fit$par * 5,
                 control = list(maxit = 1e4, factr = 10, trace = 0,
                                parscale = sc1)
  )
)
respar <- as.numeric(res2[1, 1:16])
names(respar) <- names(cm$fit$par)
grad(fn1, respar, method = 'simple', method.args = list(eps = 1e-8))

t4 <- system.time(
  res3 <- optimx(par = respar, fn = fn1,
                 #                gr = function(x) grad(fn1, x, method = 'simple', method.args = list(eps = 1e-8)),
                 method = 'L-BFGS-B',
                 lower = respar * 0.8, upper = respar * 1.2,
                 control = list(maxit = 1e4, factr = 10, trace = 0,
                                parscale = sc1)
  )
)
respar1 <- as.numeric(res3[1, 1:16])
names(respar1) <- names(cm$fit$par)
grad(fn1, respar1, method = 'simple', method.args = list(eps = 1e-8))


# check grad
temp <- sapply(cm_ped_list[1:100], function(x) max(abs(grad(m2lL, x$fit$par, method = 'simple', method.args = list(eps = 1e-8), cm = x))))

####################################
# determine parboxes for profiling #
####################################

logflag <- c()
for (i in 1:length(cm$fit$par)) {
  par1 <- cm$fit$par
  par1[i] <- -par1[i] * 0.1
  logflag[i] <- !is.finite(try(m2lL(par1, cm)))
}
logflag[2] <- TRUE
names(logflag) <- names(cm$fit$par)

#cm$parbox <- parbox(cm)
#for (i in 1:length(cm_ped_list)) system.time(cm_ped_list[[i]]$parbox <- parbox(cm_ped_list[[i]]))
system.time(parbox_list <- foreach (i = 1:length(cm_ped_list[1:100])) %dopar% parbox(cm_ped_list[[i]]))
for (i in 1:length(cm_ped_list[1:100])) cm_ped_list[[i]]$parbox <- parbox_list[[i]]
################
### profiles ###
################
#system.time(temp1 <- try(lapply(names(cm$fit$par), profile, cm = cm_ped_list[[1]])))
#cm_ped_list[[1]]$profiles <- list()
#for (j in 1:length(cm$fit$par)) cm_ped_list[[1]]$profiles[[j]] <- temp1[[j]]
#names(cm_ped_list[[1]]$profiles) <- names(cm$fit$par)
prof_list <- list()
registerDoMC(4)
t1 <- system.time(prof_list <- foreach (i = 1:100) %dopar% 
                    try(lapply(names(cm$fit$par), profile, cm = cm_ped_list[[i]])))
t2 <- system.time(prof_list2 <- foreach (i = 101:300) %dopar% 
                    try(lapply(names(cm$fit$par), profile, cm = cm_ped_list[[i]])))
debseq <- which(sapply(prof_list, length) == 1)
prof_list1 <- foreach (i = 1:length(debseq)) %dopar% 
  try(lapply(names(cm$fit$par), profile, cm = cm_ped_list[[debseq[i]]]))

prof_list[debseq] <- prof_list1

# add information to model list
for (i in 1:length(prof_list)) {
  cm_ped_list[[i]]$profiles <- list()
  for (j in 1:length(cm$fit$par)) cm_ped_list[[i]]$profiles[[j]] <- prof_list[[i]][[j]]
  names(cm_ped_list[[i]]$profiles) <- names(cm$fit$par)
}

prof_list2 <- foreach (i = 100:length(debseq)) %dopar% 
  try(lapply(names(cm$fit$par), profile, cm = cm_ped_list[[debseq[i]]]))

t2 <- system.time(prof_list1 <- foreach (i = 101:160) %dopar% 
                    try(lapply(names(cm$fit$par), profile, cm = cm_ped_list[[i]])))

for (i in 1:length(prof_list)) {
  cm_ped_list[[i]]$profiles <- prof_list[[i]]
  names(cm_ped_list[[i]]$profiles) <- names(cm$fit$par)
}

############################
### Confidence intervals ###
############################

system.time(CI_list <- lapply(cm_ped_list[1:100], 
                              function(cm) lapply(names(cm$fit$par), CI_cm, cm)))

for (i in 1:100) {
  cm_ped_list[[i]]$ci <- CI_list[[i]]
  names(cm_ped_list[[i]]$ci) <- names(cm$fit$par)
}

for (i in 1:100) {
  cm_ped_list[[i]]$ci_isin <- list()
  for (j in 1:length(cm$fit$par)) cm_ped_list[[i]]$ci_isin[[j]] <- sapply(cm_ped_list[[i]]$ci[[j]], function(x) cm_true$fit$par[j] > x[1] & cm_true$fit$par[j] < x[2])
  names(cm_ped_list[[i]]$ci_isin) <- names(cm$fit$par)
}

#calculate summary
rearr <- list()
for (i in names(cm$fit$par)) {
  rearr[[i]] <- list()
  for (j in names(cm_ped_list[[1]]$ci_isin[[i]])) {
    rearr[[i]][[j]] <- sapply(cm_ped_list[1:100], function(x) as.logical(x$ci_isin[[i]][j]))
  }
  names(rearr[[i]]) <- c('0.5', '0.8', '0.9', '0.95', '0.99')
}

summ <- list()
for (i in names(cm$fit$par)) {
  summ[[i]] <- list()
  for (j in names(rearr[[i]])) {
    summ[[i]][[j]] <- table(rearr[[i]][[j]], useNA = 'always')
  }
  #names(rearr[[i]]) <- c('0.5', '0.8', '0.9', '0.95', '0.99')
}

summ <- lapply(summ, function(x) sapply(x, function(y) {
  v <- vector(mode = "numeric", length = 3) 
  names(v) <- c('F', 'T', 'NA')
  ifelse('FALSE' %in% names(y), v['F'] <- y['FALSE'], v['F'] <- 0)
  ifelse('TRUE' %in% names(y), v['T'] <- y['TRUE'], v['T'] <- 0)
  ifelse('NA' %in% names(y), v['NA'] <- y['NA'], v['NA'] <- 0)
  v
}))


#visualize
ci_acc_vis_list <- lapply(summ, ci_acc_vis)
# add titles
for (i in 1:length(cm$fit$par)) 
  ci_acc_vis_list[[i]] <- ci_acc_vis_list[[i]] + ggtitle(names(cm$fit$par)[i])
pdf('ci_acc.pdf')
multiplot(ci_acc_vis_list, cols = 4)
dev.off()

temp <- list()
for (j in 1:length(cm$fit$par)) temp[[j]] <- CI_cm(names(cm$fit$par)[j], cm_ped_list[[1]])





# debug
deb <- list()
for (j in 1:length(names(cm$fit$par)))
  deb[[j]] <-  profile(names(cm$fit$par)[j], cm = cm_ped_list[[2]])
cm_ped_list[[2]]$profiles <- deb
names(cm_ped_list[[2]]$profiles) <- names(cm_true$fit$par)
for (j in 1:length(names(cm$fit$par)))
  deb[[j]] <-  profile_vis(names(cm$fit$par)[j], cm = cm_ped_list[[2]])
pdf('deb_prof.pdf')
deb
dev.off()


prof_res <- list()
system.time(for (i in 1:length(names(cm$fit$par))) {
  system.time(prof_res[[i]] <- profile(names(cm$fit$par)[i], cm))
})

###################################################
### determine approximate 0.99 CI ratio margins ###
###################################################

CIratio <- list()
for (i in 1:100) {
  CIratio[[i]] <- list()
  for (j in 1:length(cm_true$fit$par)) {
    CIratio[[i]][[j]] <- cm_ped_list[[i]]$ci[[j]][, '0.99'] / cm_ped_list[[i]]$fit$par[j]
  }
  names(CIratio[[i]]) <- names(cm_true$fit$par)
}

lmarg <- sapply(names(cm_true$fit$par), 
                function(x) min(sapply(CIratio, function(y) y[[x]][1]), na.rm = TRUE))
rmarg <- sapply(names(cm_true$fit$par), 
                function(x) max(sapply(CIratio, function(y) y[[x]][2]), na.rm = TRUE))
marg_df <- data.frame(l = lmarg, r = rmarg)
rownames(marg_df) <- names(cm_true$fit$par)

##########################
### Confidence regions ###
##########################

res_CR <- t(as.data.frame(lapply(cm_ped_list[1:100], CR_isin)))
rownames(res_CR) <- c()
df_CR_acc <- data.frame(p = c(0.5, 0.8, 0.9, 0.95, 0.99), pi = sapply(as.data.frame(res_CR), function(x) length(x[x == TRUE]) / length(x)))

res_CR_mod <- t(as.data.frame(lapply(cm_ped_list[1:100], CR_mod_isin)))
rownames(res_CR_mod) <- c()
df_CR_mod_acc <- data.frame(p = c(0.5, 0.8, 0.9, 0.95, 0.99), pi = sapply(as.data.frame(res_CR_mod), function(x) length(x[x == TRUE]) / length(x)))

########################
### Confidence bands ###
########################

## "naive" bands ##
t1 <- system.time(temp <- CB_point_naive(x = 50, fname = 'Epo_ext', cm = cm_ped_list[[1]], conf = 0.99))

## "correct" bands ##

registerDoMC(4)
t1 <- system.time(band_Epo_ext <- foreach (i = 1:50) %dopar% 
                    CB(xseq = c(50, 150, 250), fname = 'Epo_ext_cpm', cm = cm_ped_list[[i]]))

#51-100
registerDoMC(4)
t2 <- system.time(band_Epo_ext2 <- foreach (i = 1:50) %dopar% 
                    CB(xseq = c(50, 150, 250), fname = 'Epo_ext_cpm', cm = cm_ped_list[[i + 50]]))

registerDoMC(4)
t1 <- system.time(band_Epo_mem <- foreach (i = 1:100) %dopar% 
                    CB(xseq = c(50, 150, 250), fname = 'Epo_mem_cpm', cm = cm_ped_list[[i]]))
save.image()
registerDoMC(4)
t2 <- system.time(band_Epo_int <- foreach (i = 1:100) %dopar% 
                    CB(xseq = c(50, 150, 250), fname = 'Epo_int_cpm', cm = cm_ped_list[[i]]))

registerDoMC(4)
t1 <- system.time(band_Epo_bound <- foreach(i = 1:100) %dopar%
                    CB(xseq = 1:3, fname = 'epo_bound', cm = cm_ped_list[[i]]))

#write to cm_ped_list
for (i in 1:length(band_Epo_ext)) {
  cm_ped_list[[i]]$band <- list(NA, NA, NA, NA)
  names(cm_ped_list[[i]]$band) <- c('Epo_ext_cpm', 'Epo_mem_cpm', 'Epo_int_cpm', 'Epo_bound')
  cm_ped_list[[i]]$band$Epo_ext_cpm <- band_Epo_ext[i][[1]]
  names(cm_ped_list[[i]]$band$Epo_ext_cpm) <- c('50', '150', '250')  
}

for (i in 1:length(band_Epo_ext2)) {
  cm_ped_list[[i + 50]]$band <- list(NA, NA, NA, NA)
  names(cm_ped_list[[i + 50]]$band) <- c('Epo_ext_cpm', 'Epo_mem_cpm', 'Epo_int_cpm', 'Epo_bound')
  cm_ped_list[[i + 50]]$band$Epo_ext_cpm <- band_Epo_ext[i][[1]]
  names(cm_ped_list[[i + 50]]$band$Epo_ext_cpm) <- c('50', '150', '250')  
}

for (i in 1:length(band_Epo_mem)) {
  cm_ped_list[[i]]$band$Epo_mem_cpm <- band_Epo_mem[i][[1]]
  names(cm_ped_list[[i]]$band$Epo_mem_cpm) <- c('50', '150', '250')  
}

for (i in 1:length(band_Epo_int)) {
  cm_ped_list[[i]]$band$Epo_int_cpm <- band_Epo_int[i][[1]]
  names(cm_ped_list[[i]]$band$Epo_int_cpm) <- c('50', '150', '250')  
}

for (i in 1:length(band_Epo_bound)) {
  cm_ped_list[[i]]$band$Epo_bound <- band_Epo_bound[i][[1]]
  names(cm_ped_list[[i]]$band$Epo_bound) <- c('1', '2', '3')  
}



# now calculate confidence bands
band_constructor(band_Epo_ext[[1]], 0.95)
#get cm_true expl. funs
cm_true$expl_fun <- list()

explicit <- lapply(ode, ode_explicit)

cm_true$expl_fun$ode <- explicit[[1]]
cm_true$expl_fun$extra <- explicit[[2]]



#check accuracy
isin <- band_acc(num = 100, xseq = c(50, 150, 250), expl_name = 'Epo_ext_cpm')
  acc <- list()
  acc$band <- list()
  isin <- band_acc(num = 100, xseq = c(50, 150, 250), expl_name = 'Epo_ext_cpm')
  acc$band$Epo_ext_cpm <- sapply(isin, function(x) apply(x, 2, function(y) {print(y); sum(y) / length(y)}))
  isin <- band_acc(num = 100, xseq = c(50, 150, 250), expl_name = 'Epo_mem_cpm')
  acc$band$Epo_mem_cpm <- sapply(isin, function(x) apply(x, 2, function(y) {print(y); sum(y) / length(y)}))
  isin <- band_acc(num = 100, xseq = c(50, 150, 250), expl_name = 'Epo_int_cpm')
  acc$band$Epo_int_cpm <- sapply(isin, function(x) apply(x, 2, function(y) {print(y); sum(y) / length(y)}))
  isin <- band_acc(num = 100, xseq = c(1, 2, 3), expl_name = 'Epo_bound')
  acc$band$Epo_bound <- sapply(isin, function(x) apply(x, 2, function(y) {print(y); sum(y) / length(y)}))

#####################
### bootstrapping ###
#####################

#create bs models

  system.time(cm_ped_list <- lapply(cm_ped_list, function(x) {
  x$bs_models <- replicate(100, bootstrap_model(x), simplify = FALSE)
  x
}))

#fit bs models
system.time(temp <- fit_bs(cm_ped_list[[1]]$bs_models[[1]], cm_ped_list[[1]]$fit$par, cm_ped_list[[1]]))

system.time(temp2 <- fit_bs(cm_ped_list[[1]]$bs_models[[2]], cm_ped_list[[1]]$fit$par, cm_ped_list[[1]]))

system.time(temp <- lapply(cm_ped_list[[1]]$bs_models[1:5], fit_bs, cm_ped_list[[1]]$fit$par, cm_ped_list[[1]]))

