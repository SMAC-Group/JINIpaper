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
load(".simu/data/roblogistic_2_n_200_p_20.rds")
n_methods <- length(res$setting$method) - 1 # minus initial
p <- as.integer(res$setting$p)
MC <- as.integer(res$setting$MC)
lab_beta <- c("$\\beta_0$","$\\beta_{1}$","$\\beta_{2}$",
              "$\\beta_{4:p}$")
lab_methods <- c("JIMI","Consistent","CR","CR2","BY","BY2","MLE","BR")
lab_main <- paste0("$n=",res$setting$n,"$")

n <- as.integer(res$setting$n)
beta <- eval(str2expression(res$setting$beta))
Sigma <- str2expression(res$setting$Sigma)
eval(Sigma)
design <- str2expression(res$setting$Design)
cc <- as.numeric(Sys.getenv("C"))

set.seed(as.integer(res$setting$seed))
seed <- vector("list",3)
seed$process <- sample.int(1e7,MC)
seed$design <- sample.int(1e7,MC)
seed$sc <- sample.int(1e7,MC)

res$robCR2_sd <- matrix(ncol=p+1,nrow=MC)
res$robBY2_sd <- matrix(ncol=p+1,nrow=MC)

library(JINIpaper)
library(stocapp)
library(robustbase)

for(i in seq_len(MC)) {
  # set the seed
  set.seed(seed$process[i])
  # simulate the design
  eval(design)
  logistic_object <- make_logistic(x, beta, robust=TRUE)
  # simulate logistic
  y <- simulation(logistic_object,
                  control = list(seed=seed$process[i]))
  
  ##------ Robust estimator (Cantoni-Ronchetti) ----------------
  fit_robCR <- NULL
  try(fit_robCR <- glmrob(y ~ x, family=binomial(link="logit"), control = glmrobMqle.control(maxit=200), method = "Mqle"), silent=TRUE)
  if(!is.null(fit_robCR)) res$robCR2_sd[i,] <- sqrt(diag(fit_robCR$cov))
  
  ##------ Robust estimator (Branco-Yohai) ----------------
  fit_robBY <- NULL
  try(fit_robBY <- glmrob(y ~ x, family=binomial(link="logit"), control = glmrobBY.control(maxit=200), method = "BY"), silent=T)
  if(!is.null(fit_robBY) && !is.null(fit_robBY$cov)) res$robBY2_sd[i,] <- sqrt(diag(fit_robBY$cov))
  
  cat(i,"\n")
}

save(res,file=".simu/data/roblogistic_2_n_200_p_20.rds")

alpha <- .05

CI_center <- matrix(nr=8,nc=length(beta))
CI_left <- matrix(nr=8,nc=length(beta))
CI_right <- matrix(nr=8,nc=length(beta))
for(i in seq_along(beta)){
  which_beta <- i
  CI <- list(
    jimi = res$jimi[,which_beta] + matrix(qnorm(c(alpha/2,1-alpha/2)),nr=MC,nc=2,byr=T) * matrix(apply(res$jimi,2,mad,na.rm=T,constant=1.75)[which_beta],nr=MC,nc=2),#cbind(res$jimi_sd[,which_beta],res$jimi_sd[,which_beta]),
    consistent = res$consistent[,which_beta] + matrix(qnorm(c(alpha/2,1-alpha/2)),nr=MC,nc=2,byr=T) * cbind(res$consistent_sd[,which_beta],res$consistent_sd[,which_beta]),
    robCR = res$robCR[,which_beta] + matrix(qnorm(c(alpha/2,1-alpha/2)),nr=MC,nc=2,byr=T) * cbind(res$robCR_sd[,which_beta],res$robCR_sd[,which_beta]),
    robCR2 = res$robCR[,which_beta] + matrix(qnorm(c(alpha/2,1-alpha/2)),nr=MC,nc=2,byr=T) * cbind(res$robCR2_sd[,which_beta],res$robCR2_sd[,which_beta]),
    robBY = res$robBY[,which_beta] + matrix(qnorm(c(alpha/2,1-alpha/2)),nr=MC,nc=2,byr=T) * cbind(res$robBY_sd[,which_beta],res$robBY_sd[,which_beta]),
    robBY2 = res$robBY[,which_beta] + matrix(qnorm(c(alpha/2,1-alpha/2)),nr=MC,nc=2,byr=T) * cbind(res$robBY2_sd[,which_beta],res$robBY2_sd[,which_beta]),
    mle = res$mle[,which_beta] + matrix(qnorm(c(alpha/2,1-alpha/2)),nr=MC,nc=2,byr=T) * cbind(res$mle_sd[,which_beta],res$mle_sd[,which_beta]),
    br = res$br[,which_beta] + matrix(qnorm(c(alpha/2,1-alpha/2)),nr=MC,nc=2,byr=T) * cbind(res$br_sd[,which_beta],res$br_sd[,which_beta])
  )
  CI_center[,i] <- sapply(CI,FUN=function(x) mean(apply(x,1,function(y){check_ci(ci=y,theta=beta[which_beta])})=="center",na.rm=T))
  CI_left[,i] <- sapply(CI,FUN=function(x) mean(apply(x,1,function(y){check_ci(ci=y,theta=beta[which_beta])})=="left",na.rm=T))
  CI_right[,i] <- sapply(CI,FUN=function(x) mean(apply(x,1,function(y){check_ci(ci=y,theta=beta[which_beta])})=="right",na.rm=T))
}

plot(NA,xlim=c(1,30),ylim=c(.88,1),xlab="beta",ylab="coverage")
for(i in 1:8) lines(x=seq_along(beta), y=CI_center[i,],pch=19,type="b", col=i)
ci <- binom.test(950,1000,p=.95)$conf.int
# abline(h=ci)
lines(x=seq_along(beta),y=rep(ci[1],length(beta)),lwd=2)
lines(x=seq_along(beta),y=rep(ci[2],length(beta)),lwd=2)
legend("topright",legend=lab_methods,col = 1:8,horiz = F,pch=19,bty = "n",cex = .8,xjust = 0.5)

par(mar = c(7.5,4,4,2)+.1)
plot(NA, xlim = c(0.5, length(MC) + 0.5), ylim = c(-1, 1), 
     main = paste0("Robust logistic, beta ", which_beta-1), xaxt="n", xlab = "Trial",
     ylab = "Confidence interval")
grid()
abline(h = 0)
beta_non_zero <- NULL
bnz <- NULL
confidence_intervals <- matrix(nrow = length(MC), ncol = 4)
for(i in seq_along(MC)){
  CI <- res$consistent[which_beta,i] + qnorm(c(alpha/2,1-alpha/2)) * res$consistent_sd[which_beta,i]
  confidence_intervals[i,1:2] <- CI
  if (CI[1]*CI[2] > 0){
    lines(c(i,i) - 0.1, CI, col = 1)
    points(i - 0.1, res$consistent[which_beta,i], pch = 16, cex = 0.5, col = 1)
    bnz <- c(bnz, i)
  }else{
    lines(c(i,i) - 0.1, CI, col = "grey80")
  }
  
  CI <- res$jimi[which_beta,i] + qnorm(c(alpha/2,1-alpha/2)) * res$jimi_sd[which_beta,i]
  confidence_intervals[i,3:4] <- CI
  if (CI[1]*CI[2] > 0){
    lines(c(i,i) + 0.1, CI, col = 1)
    points(i + 0.1, res$jimi[which_beta,i], pch = 16, cex = 0.5, col = 1)
    beta_non_zero <- c(beta_non_zero, i)
  }else{
    lines(c(i,i) + 0.1, CI, col = "grey80")
  }
}
axis(1, at = seq_along(FN), labels = paste0(FN,"%"), las = 1, col.axis = "black")


# -----------------
# Consistent case
# -----------------
library(JINIpaper)
load(".simu/data/logistic_2_n_800_p_45.rds")
n_methods <- 4#length(res$setting$method)
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

Sigma <- str2expression(res$setting$Sigma)
eval(Sigma)
design <- str2expression(res$setting$Design)

set.seed(as.integer(res$setting$seed))
seed <- vector("list",3)
seed$process <- sample.int(1e7,MC)
seed$design <- sample.int(1e7,MC)
seed$sc <- sample.int(1e7,MC)

res$mle_sd <- matrix(ncol=p+1,nrow=MC)
res$jimi_sd <- matrix(ncol=p+1,nrow=MC)
res$bbc_sd <- matrix(ncol=p+1,nrow=MC)
res$br_sd <- matrix(ncol=p+1,nrow=MC)

for(i in seq_len(MC)) {
  # set the seed
  set.seed(seed$process[i])
  # simulate the design
  eval(design)
  
  # compute covariance matrices
  # MLE
  if(!any(is.na(res$jimi[i,]))){
  logistic_object <- make_logistic(x, res$mle[i,], robust=FALSE)
  probs <- predict(logistic_object, type = "response")
  W <- diag(probs * (1.0 - probs))
  cov_hat <- NULL
  try(cov_hat <- solve(t(cbind(1,x)) %*% W %*% cbind(1,x)),silent=T)
  if(!is.null(cov_hat))  res$mle_sd[i,] <- sqrt(diag(cov_hat))
  
  # BR
  logistic_object <- make_logistic(x, res$br[i,], robust=FALSE)
  probs <- predict(logistic_object, type = "response")
  W <- diag(probs * (1.0 - probs))
  cov_hat <- NULL
  try(cov_hat <- solve(t(cbind(1,x)) %*% W %*% cbind(1,x)),silent=T)
  if(!is.null(cov_hat))  res$br_sd[i,] <- sqrt(diag(cov_hat))
  
  # JIMI
  logistic_object <- make_logistic(x, res$jimi[i,], robust=FALSE)
  probs <- predict(logistic_object, type = "response")
  W <- diag(probs * (1.0 - probs))
  cov_hat <- NULL
  try(cov_hat <- solve(t(cbind(1,x)) %*% W %*% cbind(1,x)),silent=T)
  if(!is.null(cov_hat))  res$jimi_sd[i,] <- sqrt(diag((1+1/H) * cov_hat))
  
  # BBC
  logistic_object <- make_logistic(x, res$bbc[i,], robust=FALSE)
  probs <- predict(logistic_object, type = "response")
  W <- diag(probs * (1.0 - probs))
  cov_hat <- NULL
  try(cov_hat <- solve(t(cbind(1,x)) %*% W %*% cbind(1,x)),silent=T)
  if(!is.null(cov_hat))  
  res$bbc_sd[i,] <- sqrt(diag((1+1/H) * cov_hat))}
  
  cat(i,"\n")
}

CI_center <- matrix(nr=n_methods,nc=length(beta))
CI_left <- matrix(nr=n_methods,nc=length(beta))
CI_right <- matrix(nr=n_methods,nc=length(beta))
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
}

plot(x=seq_along(beta),y=CI_center[1,],type="b",pch=19,ylim=c(.85,1), ylab="Coverage", xlab="Beta")
for(i in 2:n_methods) lines(x=seq_along(beta), y=CI_center[i,],pch=19,type="b", col=i)
ci <- binom.test(9500,1e4,p=.95)$conf.int
abline(h=ci)

