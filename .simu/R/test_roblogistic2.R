# -----------
# Simulations
# -----------
# general setting

# packages
library(JINIpaper)
library(stocapp)
# library(glmrob2)
# library(brglm)
library(robustbase)

## simulation specifics
MC <- 100 # number of simulations

source(".simu/compiled/logistic_approxNR.cpp.R")
NR_rob <- function(y, x, c, start, max_iter=200, stopping=1e-6){
  x1 <- cbind(1,x)
  betas <- start
  p <- length(betas)
  
  for (i in 2:max_iter){
    betas_old <- betas
    betas <- make_rob_step_cpp(y=y, X=x1, betas=betas_old, crob=c)
    
    if(norm(betas_old-betas, "2")/p < stopping){
      return(list(betas=betas, iter=i, convergence=0))
    }
  }
  
  list(betas=betas, iter=i, convergence=1)
}

logistic_wmleNR_ib <- function(x, thetastart, c = 4.685061, H = 200, maxit=200, tol=1e-7, verbose=FALSE, seed=321){
  p <- length(thetastart)
  pi0 <- t0 <- thetastart
  
  # test diff between thetas
  test_theta <- tol + 1
  
  # iterator
  k <- 0L
  diff <- rep(NA_real_, maxit)
  
  # Iterative bootstrap algorithm:
  while(test_theta > tol && k < maxit){
    # update object for simulation
    
    # approximate
    tmp_pi <- matrix(NA_real_, nrow=p, ncol=H)
    for(h in seq_len(H)){
      seed1 <- seed + h
      sim <- r_logistic(t0, x, seed1)
      fit_tmp <- logistic_wmle(sim, x, c)
      iter <- 1L
      while(fit_tmp$convergence != 0 && iter < 10L){
        seed1 <- seed + H * h + iter
        sim <- r_logistic(t0, x, seed1)
        fit_tmp <- NR_rob(sim, x, c, start = t0)
        iter <- iter + 1L
      }
      if(fit_tmp$convergence != 0) next
      tmp_pi[,h] <- fit_tmp$betas
    }
    pi_star <- rowMeans(tmp_pi, na.rm=TRUE)
    
    # update value
    delta <- pi0 - pi_star
    t1 <- t0 + delta
    
    # test diff between thetas
    test_theta <- sum(delta^2)
    if(k>0) diff[k] <- test_theta
    
    # initialize test
    if(!k) tt_old <- test_theta+1
    
    # Alternative stopping criteria, early stop :
    # if(control$early_stop){
    #   if(tt_old <= test_theta){
    #     warning("Algorithm stopped because the objective function does not reduce")
    #     break
    #   }
    # }
    
    # Alternative stopping criteria, "statistically flat progress curve" :
    if(k > 10L){
      try1 <- diff[k:(k-10)]
      try2 <- k:(k-10)
      if(var(try1)<=1e-3) break
      mod <- lm(try1 ~ try2)
      if(summary(mod)$coefficients[2,4] > 0.2) break
    }
    
    # update increment
    k <- k + 1L
    
    # Print info
    if(verbose){
      cat("Iteration:",k,"Norm between theta_k and theta_(k-1):",test_theta,"\n")
    }
    
    # update theta
    t0 <- t1
    
    # update test
    tt_old <- test_theta
  }
  # warning for reaching max number of iterations
  if(k>=maxit) warning("maximum number of iteration reached")
  
  list(iteration = k,
       of = sqrt(drop(crossprod(delta))),
       estimate = t0,
       test_theta = test_theta,
       boot = tmp_pi)
}


# ----------------
# Case 1: logistic
# ----------------
model <- "roblogistic_2"

# TODO Setting 1
# TODO Setting 2
# TODO Setting 3

# +---------+------+-----+-------+
# | Setting |   n  |  p  | p/n   |
# +---------+------+-----+-------+
# |       1 | 2000 | 200 | 0.100 |
# |       2 | 4000 | 300 | 0.075 |
# |       3 | 8000 | 400 | 0.050 |
# +---------+------+-----+-------+

# Non-zero coefficients:
# intercept is 1
# beta = c(rep(3,5),rep(-5,5),rep(7,5))
# beta_nz <- c(-1,rep(2,5),rep(-3,5),rep(4,5))
# New:
# beta_nz <- c(1,-3)

# Design:
# rnorm(0,1/p)

##------------------ Specific functions -------
## Function to simulate responses with outliers
## see ?ib::simulation and ?ib::control for more details 
outlying_mechanism <- function(object, control, extra){
  ftd <- fitted(object)
  n <- length(ftd)
  N <- n * control$H
  n_out <- round(n * control$eps)
  y <- matrix(ifelse(runif(N) > rep(ftd,control$H), FALSE, TRUE), ncol=control$H)
  # adverserial approach:
  # select_outliers <- order(abs(fitted(object)-0.5), decreasing = T)[1:n_out]
  # adverserial approach 2:
  # select_outliers <- order(fitted(object), decreasing = T)[1:n_out]
  # totally at random approach:
  select_outliers <- sample.int(n, n_out)
  y[select_outliers] <- !y[select_outliers]
  as.integer(y)
}

# outlying_mechanism <- function(object, x, control){
#   ftd <- fitted(object)
#   n <- length(ftd)
#   N <- n
#   n_out <- round(n * control$eps)
#   y <- matrix(ifelse(runif(N) > ftd, FALSE, TRUE), ncol=1)
#   # # adverserial approach:
#   # select_outliers <- order(abs(ftd-0.5), decreasing = T)[1:n_out]
#   # 
#   # # random split: half outliers to y, half to x
#   # select_outliers <- select_outliers[sample.int(n_out,n_out)]
#   # n12 <- round(n_out / 2)
#   # y[select_outliers[1:n12]] <- !y[select_outliers[1:n12]]
#   # n_out_x <- n_out - n12
#   # n12 <- n12 + 1
#   # select_cols <- sample.int(ncol(x),n_out_x,replace=TRUE)
#   # x[select_outliers[n12:n_out], select_cols] <- x[select_outliers[n12:n_out], select_cols] * -10
#   # return(list(y = as.integer(y), x = x))
#   select_outliers <- sample.int(n, n_out)
#   y[select_outliers] <- !y[select_outliers]
#   as.integer(y)
# }

##------------------ Setting ----------------
# setting <- as.integer(Sys.getenv("SETTING"))
# n <- as.integer(Sys.getenv("N"))
n <- 200
# p <- as.integer(Sys.getenv("P"))
# p <- 50
# beta <- c(beta_nz,rep(0,p-15))
beta <- c(0.3,-2,-4,rep(0,17))
p <- length(beta)-1
# mat_design <- Sys.getenv("mat")
mat_design <- "toeplitz"
if(mat_design=="toeplitz") Sigma <- t(chol(toeplitz(0.8^(0:(p-1))))) else Sigma <- diag(p)
k <- 4

set.seed(895724)
seed <- vector("list",3)
seed$process <- sample.int(1e7,MC)
seed$design <- sample.int(1e7,MC)
seed$sc <- sample.int(1e7,MC)

##------------------ Regular logistic --------
jini_control <- stocappControl(maxit = 1e4, tol=1e-7) # control for stochastic approximation

res <- list(mle = matrix(ncol=p+1,nrow=MC),
            br = matrix(ncol=p+1,nrow=MC),
            initial = matrix(ncol=p+1,nrow=MC),
            initial1 = matrix(ncol=p+1,nrow=MC),
            jini = matrix(ncol=p+1,nrow=MC),
            consistent = matrix(ncol=p+1,nrow=MC),
            robBY = matrix(ncol=p+1,nrow=MC))




for(m in 1:MC){
  ##------- simulate the process ---------
  # set the seed
  set.seed(seed$process[m])
  # simulate the design
  # x <- t(Sigma %*% matrix(rnorm(n*p,sd= k/sqrt(n)), nr=p))
  x <- matrix(rnorm(n*p,sd= k/sqrt(n)), nr=n) %*% Sigma
  # x <- matrix(rnorm(n * p, sd = 4 / sqrt(n)), nr = n)
  logistic_object <- make_logistic(x, beta, robust=TRUE)
  # simulate logistic
  y <- simulation(logistic_object,
                  control = list(seed=seed$process[m]))
  
  
  ##------ MLE estimation ----------------
  fit_mle <- glm(y ~ x, family=binomial(link="logit"), control = glm.control(maxit=200))
  res$mle[m,] <- coef(fit_mle)
  
  ##------ Robust estimator (Cantoni-Ronchetti) ----------------
  # fit_robHub <- NULL
  # try(fit_robHub <- glmrob(y ~ x, family=binomial(link="logit"), control = glmrobMqle.control(maxit=200), method = "Mqle"), silent=T)
  # if(!is.null(fit_robHub))  res$consistent[m,] <- coef(fit_robHub)
  
  ##------ Robust estimator (Branco-Yohai) ----------------
  fit_robBY <- glmrob(y ~ x, family=binomial(link="logit"), control = glmrobBY.control(maxit=200), method = "BY")
  res$robBY[m,] <- coef(fit_robBY)
  # 
  # ##------ Initial estimator (inconsistent, Tukey bisquare) ----------------
  # cc <- find_tuning_constantWmle(x, beta, .3, 20)
  # cc <- 1.85
  cc <- 3
  fit_initial1 <- roblogisticWmle1(y, x, start = coef(fit_mle), c = cc)
  # res$initial1[m,] <- fit_initial1$coefficients
  # fit_initial2 <- logistic_wmle(y, x, c = cc)
  res$initial1[m,] <- fit_initial1$coefficients
  # fit_initial <- NR_rob(y, x, start = coef(fit_mle), c = cc)
  # res$initial[m,] <- fit_initial$betas
  
  # fit_initial2 <- NR_rob(y, x, cc, coef(fit_mle))
  
  # fit_jini <- logistic_wmle_ib(x, thetastart=fit_initial1$coefficients, c=cc, seed=seed$sc[m], verbose = T, maxit = 5, H = 20)
  # res$jini[m,] <- fit_jini$estimate
  # 
  # fit_jini <- roblogisticWmle1_ib(x, thetastart=fit_initial2$coefficients, c=cc, seed=seed$sc[m], verbose= T, maxit = 5, H=20)
  # res$jini[m,] <- fit_jini$estimate
  # 
  cat(m,"\n")
}
res_no_outliers <- res
res <- list(mle = matrix(ncol=p+1,nrow=MC),
            br = matrix(ncol=p+1,nrow=MC),
            initial = matrix(ncol=p+1,nrow=MC),
            initial1 = matrix(ncol=p+1,nrow=MC),
            jini = matrix(ncol=p+1,nrow=MC),
            consistent = matrix(ncol=p+1,nrow=MC),
            robBY = matrix(ncol=p+1,nrow=MC))

for(m in 1:MC){
  ##------- simulate the process ---------
  # set the seed
  set.seed(seed$process[m])
  # simulate the design
  # x <- t(Sigma %*% matrix(rnorm(n*p,sd= k / sqrt(sqrt(n))), nr=p))
  # x <- matrix(rnorm(n * p, sd = 4 / sqrt(n)), nr = n)
  x <- matrix(rnorm(n*p,sd= k/sqrt(n)), nr=n) %*% Sigma
  logistic_object <- make_logistic(x, beta, robust=TRUE)
  # simulate logistic with contaminated y
  y <- simulation(logistic_object,
                  control = list(seed=seed$process[m],
                                 sim=outlying_mechanism, eps=0.02))
  # jini_control$seed <- seed$process[m]
  # jini_control$eps <- 0.04
  # sim <- outlying_mechanism(logistic_object, x, control = jini_control)
  # y <- sim$y
  # x <- sim$x
  
  
  ##------ MLE estimation ----------------
  fit_mle <- glm(y ~ x, family=binomial(link="logit"), control = glm.control(maxit=200))
  res$mle[m,] <- coef(fit_mle)
  
  # ##------ Robust estimator (Cantoni-Ronchetti) ----------------
  # fit_robHub <- glmrob(y ~ x, family=binomial(link="logit"), control = glmrobMqle.control(maxit=200), method = "Mqle")
  # res$consistent[m,] <- coef(fit_robHub)
  # 
  # ##------ Robust estimator (Branco-Yohai) ----------------
  fit_robBY <- glmrob(y ~ x, family=binomial(link="logit"), control = glmrobBY.control(maxit=200), method = "BY")
  res$robBY[m,] <- coef(fit_robBY)
  # 
  # ##------ Initial estimator (inconsistent, Tukey bisquare) ----------------
  # cc <- find_tuning_constantWmle(x, beta, .3, 20)
  # cc <- 2.17
  cc <- 3
  fit_initial1 <- roblogisticWmle1(y, x, start = coef(fit_mle), c = cc)
  res$initial1[m,] <- fit_initial1$coefficients
  # fit_initial <- NR_rob(y, x, start = coef(fit_mle), c = cc)
  # res$initial[m,] <- fit_initial$betas
  # fit_initial1 <- logistic_wmle(y, x, c = cc)
  # res$initial1[m,] <- fit_initial1
  # 
  # fit_jini <- roblogisticWmle1_ib(x, thetastart=fit_initial1$coefficients, c=cc, seed=seed$sc[m])
  # res$jini[m,] <- fit_jini$estimate
  
  cat(m,"\n")
}
res_outliers <- res

ind <- 1:15
y_lim <- c(-5,5)
par(mfrow=c(3,2))
boxplot(res_no_outliers$mle-matrix(beta,nr=MC,nc=length(beta),byr=T),ylim=y_lim);abline(h=0)
boxplot(res_outliers$mle-matrix(beta,nr=MC,nc=length(beta),byr=T),ylim=y_lim);abline(h=0)

# par(mfrow=c(2,1))
# boxplot(res_no_outliers$initial-matrix(beta,nr=MC,nc=length(beta),byr=T),ylim=y_lim);abline(h=0)
# boxplot(res_outliers$initial-matrix(beta,nr=MC,nc=length(beta),byr=T),ylim=y_lim);abline(h=0)
boxplot(res_no_outliers$initial1-matrix(beta,nr=MC,nc=length(beta),byr=T),ylim=y_lim);abline(h=0)
boxplot(res_outliers$initial1-matrix(beta,nr=MC,nc=length(beta),byr=T),ylim=y_lim);abline(h=0)

# 
# par(mfrow=c(2,1))
# boxplot(res_no_outliers$jini-matrix(beta,nr=MC,nc=length(beta),byr=T),ylim=y_lim);abline(h=0)
# boxplot(res_outliers$jini-matrix(beta,nr=MC,nc=length(beta),byr=T),ylim=y_lim);abline(h=0)
# 
# par(mfrow=c(2,1))
# boxplot(res_no_outliers$consistent-matrix(beta,nr=MC,nc=length(beta),byr=T),ylim=y_lim);abline(h=0)
# boxplot(res_outliers$consistent-matrix(beta,nr=MC,nc=length(beta),byr=T),ylim=y_lim);abline(h=0)
# 
# par(mfrow=c(2,1))
boxplot(res_no_outliers$robBY - matrix(beta,nr=MC,nc=length(beta),byr=T),ylim=y_lim);abline(h=0)
boxplot(res_outliers$robBY - matrix(beta,nr=MC,nc=length(beta),byr=T),ylim=y_lim);abline(h=0)
