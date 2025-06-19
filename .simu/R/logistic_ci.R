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
# Consistent case
# -----------------
library(JINIpaper)
load(".simu/data/logistic_2_n_200_p_25.rds")
n_methods <- 4
p <- as.integer(res$setting$p)
n <- as.integer(res$setting$n)
H <- as.integer(res$setting$H)
MC <- as.integer(res$setting$MC)
lab_beta <- c("$\\beta_0$","$\\beta_{1}$","$\\beta_{2}$",
              "$\\beta_{4:p}$")
lab_methods <- c("JIMI","BBC","MLE")
lab_main <- paste0("$n=",res$setting$n,"$")

eval(str2expression(res$setting$beta))
alpha <- .05

# Sigma <- str2expression(res$setting$Sigma)
# eval(Sigma)
# design <- str2expression(res$setting$Design)
# 
# set.seed(as.integer(res$setting$seed))
# seed <- vector("list",3)
# seed$process <- sample.int(1e7,MC)
# seed$design <- sample.int(1e7,MC)
# seed$sc <- sample.int(1e7,MC)
# 
# res$mle_sd <- matrix(ncol=p+1,nrow=MC)
# res$jimi_sd <- matrix(ncol=p+1,nrow=MC)
# res$bbc_sd <- matrix(ncol=p+1,nrow=MC)
# res$br_sd <- matrix(ncol=p+1,nrow=MC)
# 
# for(i in seq_len(MC)) {
#   # set the seed
#   set.seed(seed$process[i])
#   # simulate the design
#   eval(design)
#   
#   # compute covariance matrices
#   # MLE
#   if(!any(is.na(res$jimi[i,]))){
#     logistic_object <- make_logistic(x, res$mle[i,], robust=FALSE)
#     probs <- predict(logistic_object, type = "response")
#     W <- diag(probs * (1.0 - probs))
#     cov_hat <- NULL
#     try(cov_hat <- solve(t(cbind(1,x)) %*% W %*% cbind(1,x)),silent=T)
#     if(!is.null(cov_hat))  res$mle_sd[i,] <- sqrt(diag(cov_hat))
#     
#     # BR
#     logistic_object <- make_logistic(x, res$br[i,], robust=FALSE)
#     probs <- predict(logistic_object, type = "response")
#     W <- diag(probs * (1.0 - probs))
#     cov_hat <- NULL
#     try(cov_hat <- solve(t(cbind(1,x)) %*% W %*% cbind(1,x)),silent=T)
#     if(!is.null(cov_hat))  res$br_sd[i,] <- sqrt(diag(cov_hat))
#     
#     # JIMI
#     logistic_object <- make_logistic(x, res$jimi[i,], robust=FALSE)
#     probs <- predict(logistic_object, type = "response")
#     W <- diag(probs * (1.0 - probs))
#     cov_hat <- NULL
#     try(cov_hat <- solve(t(cbind(1,x)) %*% W %*% cbind(1,x)),silent=T)
#     if(!is.null(cov_hat))  res$jimi_sd[i,] <- sqrt(diag((1+1/H) * cov_hat))
#     
#     # BBC
#     logistic_object <- make_logistic(x, res$bbc[i,], robust=FALSE)
#     probs <- predict(logistic_object, type = "response")
#     W <- diag(probs * (1.0 - probs))
#     cov_hat <- NULL
#     try(cov_hat <- solve(t(cbind(1,x)) %*% W %*% cbind(1,x)),silent=T)
#     if(!is.null(cov_hat))  
#       res$bbc_sd[i,] <- sqrt(diag((1+1/H) * cov_hat))}
#   
#   cat(i,"\n")
# }
# 
# save(res,file=".simu/data/logistic_2_n_200_p_25.rds")

CI_center <- matrix(nr=n_methods,nc=length(beta))
CI_left <- matrix(nr=n_methods,nc=length(beta))
CI_right <- matrix(nr=n_methods,nc=length(beta))
CI_length <- matrix(nr=n_methods,nc=length(beta))
for(i in seq_along(beta)){
  which_beta <- i
  CI <- list(
    jimi = res$jimi[,which_beta] + matrix(qnorm(c(alpha/2,1-alpha/2)),nr=MC,nc=2,byr=T) * cbind(res$jimi_sd[,which_beta],res$jimi_sd[,which_beta]),
    mle = res$mle[,which_beta] + matrix(qnorm(c(alpha/2,1-alpha/2)),nr=MC,nc=2,byr=T) * cbind(res$mle_sd[,which_beta],res$mle_sd[,which_beta]),
    bbc = res$bbc[,which_beta] + matrix(qnorm(c(alpha/2,1-alpha/2)),nr=MC,nc=2,byr=T) * cbind(res$bbc_sd[,which_beta],res$bbc_sd[,which_beta]),
    br = res$br[,which_beta] + matrix(qnorm(c(alpha/2,1-alpha/2)),nr=MC,nc=2,byr=T) * cbind(res$br_sd[,which_beta],res$br_sd[,which_beta])
  )
  CI_center[,i] <- sapply(CI,FUN=function(x) mean(apply(x,1,function(y){check_ci(ci=y,theta=beta[which_beta])})=="center",na.rm=T))
  CI_left[,i] <- sapply(CI,FUN=function(x) mean(apply(x,1,function(y){check_ci(ci=y,theta=beta[which_beta])})=="left",na.rm=T))
  CI_right[,i] <- sapply(CI,FUN=function(x) mean(apply(x,1,function(y){check_ci(ci=y,theta=beta[which_beta])})=="right",na.rm=T))
  CI_length[,i] <- sapply(CI,FUN=function(x) median(apply(x,1,diff),na.rm=T))
}

plot(x=seq_along(beta),y=CI_center[1,],type="b",pch=19,ylim=c(.85,1), ylab="Coverage", xlab="Beta")
for(i in 2:n_methods) lines(x=seq_along(beta), y=CI_center[i,],pch=19,type="b", col=i)
ci <- binom.test(9500,1e4,p=.95)$conf.int
abline(h=ci)

id <- order(CI$mle[,1])
plot(CI$mle[id,1],type="l",ylim=c(-10,10))
lines(x=1:1e4,y=CI$mle[id,2])
abline(h=0)
lines(x=1:1e4,y=CI$jimi[id,1],col=3)
lines(x=1:1e4,y=CI$jimi[id,2],col=3)
plot(spline(1:1e4,y=CI$jimi[id,2]))

load(".simu/data/logistic_2_n_200_p_25.rds")
mle1 <- res$mle
mle1_sd <- res$mle_sd
load(".simu/data/roblogistic_2_n_200_p_20.rds")
mle2 <- res$mle
mle2_sd <- res$mle_sd
boxplot(mle1_sd[,5]-sd(mle1[,5],na.rm=T),mle2_sd[,5]-sd(mle1[,5],na.rm=T),outline=F,names=c("Asy. variance","Bootstrap"))

load(".simu/data/logistic_2_n_200_p_25.rds")
boxplot(res$mle_sd[,5],res$bbc_sd[,5],res$jimi_sd[,5],outline=F,names=c("MLE","BBC","JIMI"))
