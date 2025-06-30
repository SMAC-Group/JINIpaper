cc <- c(2.8, 2.6, 2.2)
settings <- 1:3
eps <- c(0,0,1)
fn <- c(0,2,2)
theta0 <- c(.3 ,-2, -4, 0)
methods <- c("JINI", "consistent", "robBY")

results <- matrix(nrow = 3 * length(settings) * length(eps), ncol = 4)
colnames(results) <- c("beta0", "beta1", "beta2", "beta00")
rownames(results) <- paste0(rep(methods, times = length(settings) * length(eps)), 
                            "_eps_", rep(eps, times = length(methods) * length(settings)),
                            "_fn_", rep(fn, times = length(methods) * length(settings)),
                            "_setting_", rep(settings, each = length(methods) * length(eps)))
k <- 1L
for(i in seq_along(settings)){
  for(j in 1:length(eps)){
    load(paste0(".simu/data/rm_logistic_setting_",settings[i],"_cc_",cc[i],"_eps_",eps[j],"_fn_",fn[j],"_new.rds"))
    theta <- rep(0, ncol(res$jimi))
    for(t in 1:3) theta[t] <- theta0[t]
    if(i==1){bias <- apply(res$jimi,2,median,na.rm = TRUE) - theta}else{bias <- colMeans(res$jimi,na.rm = TRUE) - theta}
    results[k,] <- c(bias[1:3], mean(bias[4:ncol(res$jimi)], na.rm = TRUE))
    # bias <- colMeans(res$consistent,na.rm = TRUE) - theta
    if(i==1){bias <- apply(res$consistent,2,median,na.rm = TRUE) - theta}else{bias <- colMeans(res$consistent,na.rm = TRUE) - theta}
    results[k+1,] <- c(bias[1:3], mean(bias[4:ncol(res$jimi)], na.rm = TRUE))
    # bias <- colMeans(res$robBY,na.rm = TRUE) - theta
    if(i==1){bias <- apply(res$robBY,2,median,na.rm = TRUE) - theta}else{bias <- colMeans(res$robBY,na.rm = TRUE) - theta}
    results[k+2,] <- c(bias[1:3], mean(bias[4:ncol(res$jimi)], na.rm = TRUE))
    k <- k + 3L
  }
}

load(".simu/data/rm_logistic_setting_3_cc_2.2_eps_0_fn_0_new.rds")

# simu 1: 
# small p large n JINI robust, JINI mle, mle, BY
