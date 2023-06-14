# Plot point estimates with 1-alpha CI
source(".simu/R/figures.R")
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
n_setting <- 1
i=3
load(paste0(".simu/data/roblogistic_1_sd_setting_",i,".rds"))
n_methods <- length(res$setting$method) - 1 # minus initial
p <- as.integer(res$setting$p) +1
MC <- as.integer(res$setting$MC)
lab_beta <- c("$\\beta_0$","$\\beta_{1}$","$\\beta_{2}$",
              "$\\beta_{4:p}$")
lab_methods <- c("JIMI","Consistent","CR","BY","MLE","BR")
lab_main <- paste0("$n=",res$setting$n,"$")

eval(str2expression(res$setting$beta))
alpha <- .05

which_beta <- 3
CI <- list(
  jimi = res$jimi[,which_beta] + matrix(qnorm(c(alpha/2,1-alpha/2)),nr=MC,nc=2,byr=T) * cbind(res$jimi_sd[,which_beta],res$jimi_sd[,which_beta]),
  consistent = res$consistent[,which_beta] + matrix(qnorm(c(alpha/2,1-alpha/2)),nr=MC,nc=2,byr=T) * cbind(res$consistent_sd[,which_beta],res$consistent_sd[,which_beta]),
  robCR = res$robCR[,which_beta] + matrix(qnorm(c(alpha/2,1-alpha/2)),nr=MC,nc=2,byr=T) * cbind(res$robCR_sd[,which_beta],res$robCR_sd[,which_beta]),
  robBY = res$robBY[,which_beta] + matrix(qnorm(c(alpha/2,1-alpha/2)),nr=MC,nc=2,byr=T) * cbind(res$robBY_sd[,which_beta],res$robBY_sd[,which_beta]),
  mle = res$mle[,which_beta] + matrix(qnorm(c(alpha/2,1-alpha/2)),nr=MC,nc=2,byr=T) * cbind(res$mle_sd[,which_beta],res$mle_sd[,which_beta]),
  br = res$br[,which_beta] + matrix(qnorm(c(alpha/2,1-alpha/2)),nr=MC,nc=2,byr=T) * cbind(res$br_sd[,which_beta],res$br_sd[,which_beta])
)
lapply(CI,FUN=function(x) mean(apply(x,1,function(y){check_ci(ci=y,theta=beta[which_beta])})=="center",na.rm=T))
lapply(CI,FUN=function(x) mean(apply(x,1,function(y){check_ci(ci=y,theta=beta[which_beta])})=="left",na.rm=T))
lapply(CI,FUN=function(x) mean(apply(x,1,function(y){check_ci(ci=y,theta=beta[which_beta])})=="right",na.rm=T))
lapply(CI,FUN=function(x) median(apply(x,1,diff),na.rm=T))

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
