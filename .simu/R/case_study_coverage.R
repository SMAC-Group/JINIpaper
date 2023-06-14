#-------------
# CASE STUDY
#-------------
# packages
library(ib)
library(JINIpaper)

# for verifying a confidence interval
check_ci <- function(ci,theta){
  if(!is.matrix(ci)) ci <- as.matrix(ci)
  if(ncol(ci)!=2) ci <- t(ci)
  if(nrow(ci)!=length(theta)) stop("verify ci and theta dimensions match")
  apply(cbind(ci,theta), 1, function(x){
    if(x[3] <= min(x[1:2])){"left"}
    else if(x[3] >= max(x[1:2])){"right"}else{"center"}
  })
}

# for simulating logistic with misclassification
misclassification <- function(object, control, extra){
  ftd <- fitted(object)
  n <- length(ftd)
  N <- n * control$H
  u1 <- 0.0
  u2 <- FN
  mu_star <- u1 * (1 - rep(ftd,control$H)) + (1 - u2) * rep(ftd,control$H)
  matrix(ifelse(runif(N) > mu_star,0,1), ncol=control$H)
}

# setting
H <- 200
B <- 100
ib_control <- ib_control2 <- ibControl(H = H, sim = misclassification, maxit = 100L)
set.seed(432)
seeds <- sample.int(1e7,B)

#-------------
# Estimation for different FN rate
#-------------
fit_initial <- glm(alc ~ ., family = binomial(link = "logit"), data = student)
x <- as.matrix(student[,-1])
y <- student[,1]

p <- length(coef(fit_initial))
tmp <- matrix(nrow = p, ncol = 8)
colnames(tmp) <- paste0(3:10, "%")
consistent <- jimi <- consistent_sd <- jimi_sd <- tmp
rm(tmp)

i <- 1L
for(FN in 3:10 / 100){
  # Consistent MLE
  fit_mle <- logistic_misclassification_mle(x, y, fp = 0, fn = FN)
  consistent[,i] <- fit_mle
  consistent_sd[,i] <- sqrt(diag(inverse_FIM(x, fit_mle, 0.0, FN)))

  # Iterative bootstrap
  fit_ib <- ib(fit_initial, control = ib_control)
  beta_tilde <- coef(fit_ib)
  jimi[,i] <- beta_tilde

  # bootstrap estimator of the covariance matrix
  logistic_object <- make_logistic(x, beta_tilde)
  boot <- matrix(nrow = B, ncol = p)
  for(b in 1:B){
    ##------- simulate the process ---------
    # simulate logistic
    y <- simulation(logistic_object, control=list(seed=b, sim=misclassification))

    ##------ MLE estimation ----------------
    fit_tmp <- glm(y~x, family=binomial(link="logit"))

    ##------ IB bias correction ------------
    ib_control2$seed <- seeds[b]
    t1 <- Sys.time()
    fit_ib <- ib(fit_tmp, control = ib_control2)
    t2 <- Sys.time()
    boot[b,] <- getEst(fit_ib)

    cat(b,"Time:",difftime(t2,t1,units = "secs"),"\n")
  }

  jimi_sd[,i] <- sqrt(diag(cov(boot)))
  
  i <- i + 1L
  cat(FN,"\n")
}

res <- list()
res$consistent <- consistent
res$consistent_sd <- consistent_sd
res$jimi <- jimi
res$jimi_sd <- jimi_sd
save(res, file=".simu/data/cs_coverage.rds")

# Plot point estimates with 1-alpha CI
load(".simu/data/cs_coverage.rds")
FN <- 3:10
which_beta <- 42
alpha <- .05
par(mar = c(7.5,4,4,2)+.1)
plot(NA, xlim = c(0.5, length(FN) + 0.5), ylim = c(-1, 1), 
     main = paste0("Case Study, beta ", which_beta), xaxt="n", xlab = "False negative rate",
     ylab = "Confidence interval")
grid()
abline(h = 0)
beta_non_zero <- NULL
bnz <- NULL
confidence_intervals <- matrix(nrow = length(FN), ncol = 4)
for(i in seq_along(FN)){
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

# axis(1, at = beta_non_zero, labels = var_names[beta_non_zero], las = 2)
# axis(1, at = setdiff(bnz,beta_non_zero), labels = var_names[setdiff(bnz,beta_non_zero)], las = 2, col.axis = "red")

