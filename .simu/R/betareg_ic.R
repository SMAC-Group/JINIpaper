# for verifying a confidence interval
check_ci <- function(ci,theta){
  if(any(is.na(ci))) return(NA_real_)
  if(!is.matrix(ci)) ci <- as.matrix(ci)
  if(ncol(ci)!=2) ci <- t(ci)
  if(nrow(ci)!=length(theta)) stop("verify ci and theta dimensions match")
  apply(cbind(ci,theta), 1, function(x){
    if(x[3] <= min(x[1:2])){"left"}
    else if(x[3] >= max(x[1:2])){"right"}else{"center"}
  })
}

# -----------------
# Inconsistent case
# -----------------
load(".simu/data/betareg_1_n_1000_p_70.rds")
load(".simu/data/betareg_1_n_200_p_50.rds")
load(".simu/data/betareg_1_n_500_p_60.rds")
n_methods <- length(res$setting$method) - 1 # minus initial
p <- as.integer(res$setting$p)
MC <- as.integer(res$setting$MC)

eval(str2expression(res$setting$beta))
theta <- c(beta, as.numeric(res$setting$gamma))

alpha <- .05
CI_center <- matrix(nr=2,nc=length(theta))
CI_left <- matrix(nr=2,nc=length(theta))
CI_right <- matrix(nr=2,nc=length(theta))
CI_length <- matrix(nr=2,nc=length(theta))
for(i in seq_along(theta)){
  which_beta <- i
  CI <- list(
    jimi = res$jimi[,which_beta] + matrix(qnorm(c(alpha/2,1-alpha/2)),nr=MC,nc=2,byr=T) * cbind(res$jimi_sd[,which_beta],res$jimi_sd[,which_beta]),#matrix(apply(res$jimi,2,sd,na.rm=T)[which_beta],nr=MC,nc=2),#cbind(res$jimi_sd[,which_beta],res$jimi_sd[,which_beta]),
    consistent = res$consistent[,which_beta] + matrix(qnorm(c(alpha/2,1-alpha/2)),nr=MC,nc=2,byr=T) * cbind(res$consistent_sd[,which_beta],res$consistent_sd[,which_beta])
  )
  CI_center[,i] <- sapply(CI,FUN=function(x) mean(apply(x,1,function(y){check_ci(ci=y,theta=theta[which_beta])})=="center",na.rm=T))
  CI_left[,i] <- sapply(CI,FUN=function(x) mean(apply(x,1,function(y){check_ci(ci=y,theta=theta[which_beta])})=="left",na.rm=T))
  CI_right[,i] <- sapply(CI,FUN=function(x) mean(apply(x,1,function(y){check_ci(ci=y,theta=theta[which_beta])})=="right",na.rm=T))
  CI_length[,i] <- sapply(CI,FUN=function(x) mean(apply(x,1,diff)))
}

plot(x=seq_along(theta),y=CI_center[1,],type="b",pch=19,ylim=c(.9,.97))
lines(x=seq_along(theta), y=CI_center[2,],pch=19,type="b", col="green")
ci <- binom.test(950,1000,p=.95)$conf.int
abline(h=ci)

plot(x=seq_along(theta),y=CI_right[1,],type="b",pch=19,ylim=c(0,.05))
lines(x=seq_along(theta), y=CI_right[2,],pch=19,type="b", col="green")
ci <- binom.test(25,1000,p=.95)$conf.int
abline(h=ci)

plot(x=seq_along(theta),y=CI_left[1,],type="b",pch=19,ylim=c(0,.05))
lines(x=seq_along(theta), y=CI_left[2,],pch=19,type="b", col="green")
ci <- binom.test(25,1000,p=.95)$conf.int
abline(h=ci)


# lapply(CI,FUN=function(x) mean(apply(x,1,function(y){check_ci(ci=y,theta=theta[which_beta])})=="left",na.rm=T))
# lapply(CI,FUN=function(x) mean(apply(x,1,function(y){check_ci(ci=y,theta=theta[which_beta])})=="right",na.rm=T))
# lapply(CI,FUN=function(x) median(apply(x,1,diff),na.rm=T))
