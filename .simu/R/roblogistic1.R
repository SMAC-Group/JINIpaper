# -----------
# Simulations
# -----------
# general setting

# packages
library(JINIpaper)
library(stocapp)
library(glmrob2)
library(brglm)
library(MASS)  # for mvrnorm()
library(Rcpp)
library(RcppArmadillo)
library(brglm)  # for brglm()
library(robustbase)  # for glmrob()
sourceCpp(".simu/src/logistic_approxNR.cpp")

simu_data = function(X, betas){
  ### simulate response y
  probs = expit_vec_cpp(X%*%betas)
  sapply(probs, FUN=rbinom, n=1, size=1)
}
NR_rob = function(y, X, crob, gamma, start, max_iter, stopping){
  p = ncol(X)
  max_iter = p*max_iter
  y = gamma*(1-y) + (1-gamma)*y
  
  if(sum(start=="glm") != 0){
    betas = glm(y ~ X, family=binomial(link="logit"))$coefficients
  }else{
    betas = start
  }
  
  for (i in 2:max_iter){
    betas_old = betas
    betas = make_rob_step_cpp(y=y, X=X, betas=betas_old, crob=crob)
    
    if(norm(betas_old-betas, "2")/p < stopping){
      return(list(betas=betas, iter=i, convergence=0))
    }
  }
  
  list(betas=betas, iter=i, convergence=1)
}
find_pi <- function(y, X, crob, gamma, start){
  NR_rob(y=y, X=X, crob=crob, gamma=gamma, start=start, 
         max_iter=200, stopping=10^(-6))
}
compute_bias = function(betas, X, crob, gamma, start, H){
  p = ncol(X)
  inter = rep(0,p)
  for (h in 1:H) {
    set.seed(h)
    y_sim = simu_data(X=X, betas=betas)
    estim = find_pi(y=y_sim, X=X, crob=crob, gamma=gamma, start=start)$betas
    inter = inter + estim/H
  }
  betas-inter
}
ib <- function(y, X, crob, gamma, start, H=100, iter_max=20, epsilon=10^(-6)){
  p = ncol(X)
  theta = matrix(NA, nrow=iter_max+1, ncol=p)
  theta[1,] = find_pi(y=y, X=X, crob=crob, gamma=gamma, start=start)$betas
  diff = norm(theta[1,], type="2")
  
  for (i in 1:iter_max) {
    theta[(i+1),] = theta[1,] + compute_bias(betas=theta[i,], X=X, crob=crob, gamma=gamma, start=start, H=H)
    diff_old = diff
    diff = norm(theta[(i+1),]-theta[i,], type="2")
    if (diff < epsilon || diff_old < diff){return(theta[(i+1),])}
  }
  
  return(theta[(i+1),])  
}

## simulation specifics
MC <- 100 # number of simulations

# ----------------
# Case 1: logistic
# ----------------
model <- "roblogistic_1"

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
# Old:
# intercept is 1
# beta = c(rep(3,5),rep(-5,5),rep(7,5))
# beta_nz <- c(1,rep(3,5),rep(-5,5),rep(7,5))
# New:
beta_nz <- c(1, 5)

# Design:
# Old:
# rnorm(0,1/p)
# New:

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
  select_outliers <- order(abs(fitted(object)-0.5), decreasing = T)[1:n_out]
  # adverserial approach 2:
  # select_outliers <- order(fitted(object), decreasing = F)[1:n_out]
  # totally at random approach:
  # select_outliers <- sample.int(n, n_out)
  y[select_outliers] <- !y[select_outliers]
  as.integer(y)
}

##------------------ Setting ----------------
# setting <- as.integer(Sys.getenv("SETTING"))
# n <- as.integer(Sys.getenv("N"))
n <- 100
# p <- as.integer(Sys.getenv("P"))
p <- 2
beta <- c(beta_nz,rep(0,p-1))
# mat_design <- Sys.getenv("mat")
mat_design <- "diag"
if(mat_design=="toeplitz") Sigma <- t(chol(toeplitz(0.8^(0:(p-1))))) else Sigma <- diag(p)

set.seed(895724)
seed <- vector("list",3)
seed$process <- sample.int(1e7,MC)
seed$design <- sample.int(1e7,MC)
seed$sc <- sample.int(1e7,MC)

##------------------ Regular logistic --------
jini_control_out <- stocappControl(maxit = 1e4, sim = outlying_mechanism, eps = 0.05) # control for stochastic approximation
jini_control <- stocappControl(maxit = 1e4) # control for stochastic approximation

res <- list(mle = matrix(ncol=p+1,nrow=MC),
            br = matrix(ncol=p+1,nrow=MC),
            initial = matrix(ncol=p+1,nrow=MC),
            jini = matrix(ncol=p+1,nrow=MC),
            consistent = matrix(ncol=p+1,nrow=MC),
            robBY = matrix(ncol=p+1,nrow=MC))

##------------------ Slurm specs --------------
n_array <- 1000
ind <- matrix(seq_len(MC),nr=n_array,byr=T)
id_slurm <- Sys.getenv("SLURM_ARRAY_TASK_ID")

# for(m in na.omit(ind[as.numeric(id_slurm),])){
for(m in 1:MC){
 ##------- simulate the process ---------
  # set the seed
  set.seed(seed$process[m])
  # simulate the design
  # x <- matrix(rnorm(n*p,sd=p^(-.5)), nr=n)
  x <- t(Sigma %*% matrix(rnorm(n*p,sd= 2 / sqrt(sqrt(n))), nr=p))
  # x_out <- x
  # x_out[sample.int(nrow(x),10,replace=T), sample.int(ncol(x),10 ,replace = T)] <- 30
  logistic_object <- make_logistic(x, beta, robust=TRUE)
  # simulate logistic
  y <- simulation(logistic_object,
                  control = list(seed=seed$process[m], sim=outlying_mechanism, eps=0.01))
  # y <- simulation(logistic_object, 
  #                 control = list(seed=seed$process[m]))
  
  
  ##------ MLE estimation ----------------
  fit_mle <- glm(y ~ x, family=binomial(link="logit"), control = glm.control(maxit=200))
  res$mle[m,] <- coef(fit_mle)

  # ##------ "brglm" bias reduction --------
  # fit_br <- brglm(y~x, family=binomial(link="logit"), control.glm = glm.control1(maxit=200))
  # res$br[m,] <- coef(fit_br)
  # 
  # ##------ Robust estimator (Cantoni-Ronchetti) ----------------
  fit_robHub <- glmrob(y ~ x, family=binomial(link="logit"), control = glmrobMqle.control(maxit=200), method = "Mqle")
  res$consistent[m,] <- coef(fit_robHub)
  # 
  # ##------ Robust estimator (Branco-Yohai) ----------------
  fit_robBY <- glmrob(y ~ x, family=binomial(link="logit"), control = glmrobBY.control(maxit=200), method = "BY")
  res$robBY[m,] <- coef(fit_robBY)
  # 
  # ##------ Initial estimator (inconsistent, Tukey bisquare) ----------------
  fit_initial <- glmrob2(y ~ x, family=binomial(link="logit"), control = glmrob2.control(maxit=200))
  res$initial[m,] <- coef(fit_initial)
  # 
  # ##------ Stochastic approximation bias correction ------------
  jini_control$seed <- seed$sc[m] # update the seed
  # t1 <- Sys.time()
  fit_jini <- stocapp(fit_initial, control = jini_control) # compute JINI using stochastic approximation, see ? stocapp::stocapp for more details
  # t2 <- Sys.time()
  res$jini[m,] <- getEst(fit_jini) # retrieve estimator
  # res$iteration[m] <- getIteration(fit_jini)
  # res$time_ib[m] <- difftime(t2,t1,units = "secs")
  # 
  # ##------ Stochastic approximation bias correction (oracle) ------------
  # # jini_control_out$seed <- seed$sc[m] # update the seed
  # # fit_jini2 <- stocapp(fit_initial, control = jini_control_out) # compute JINI using stochastic approximation, see ? stocapp::stocapp for more details
  # # res$jini_oracle[m,] <- getEst(fit_jini2) # retrieve estimator
  # 
  # # save results
  # save(res, file=paste0("tmp/",model,"_setting_",setting,"_id_",id_slurm,".rds"))
  cat(m,"\n")
}
res_outliers <- res
boxplot(res_outliers$mle-matrix(beta,nr=MC,nc=length(beta),byr=T));abline(h=0)
boxplot(res_outliers$robBY-matrix(beta,nr=MC,nc=length(beta),byr=T));abline(h=0)
boxplot(res_outliers$jini-matrix(beta,nr=MC,nc=length(beta),byr=T));abline(h=0)
boxplot(res_outliers$mle[,3]-beta[3],res_outliers$consistent[,3]-beta[3],
        res_outliers$robBY[,3]-beta[3],res_outliers$jini[,3]-beta[3],
        res_outliers$initial[,3]-beta[3]);abline(h=0)
boxplot(res_outliers$mle[,2]-beta[2],res_outliers$consistent[,2]-beta[2],
        res_outliers$robBY[,2]-beta[2],res_outliers$jini[,2]-beta[2],
        res_outliers$initial[,2]-beta[2]);abline(h=0)
