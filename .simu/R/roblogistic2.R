# -----------
# Simulations
# -----------
# general setting

# packages
library(JINIpaper)
library(stocapp)
library(glmrob2)
library(brglm)

## simulation specifics
MC <- 100 # number of simulations

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
# beta_nz <- c(1,rep(3,5),rep(-5,5),rep(7,5))
# New:
beta_nz <- c(1,5)

# Design:
# rnorm(0,1/p)

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
  logistic_object <- make_logistic(x, beta, robust=TRUE)
  # simulate logistic
  y <- simulation(logistic_object, 
                  control = list(seed=seed$process[m]))

  ##------ MLE estimation ----------------
  fit_mle <- glm(y ~ x, family=binomial(link="logit"), control = glm.control(maxit=200))
  res$mle[m,] <- coef(fit_mle)

  ##------ "brglm" bias reduction --------
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
  # # save results
  # save(res, file=paste0("tmp/",model,"_setting_",setting,"_id_",id_slurm,".rds"))
  cat(m,"\n")
}
res_no_outliers <- res
boxplot(res_no_outliers$mle-matrix(beta,nr=MC,nc=length(beta),byr=T));abline(h=0)

boxplot(res_no_outliers$consistent-matrix(beta,nr=MC,nc=length(beta),byr=T));abline(h=0)

par(mfrow=c(2,1))
boxplot(res_no_outliers$mle-matrix(beta,nr=MC,nc=length(beta),byr=T));abline(h=0)
boxplot(res_outliers$mle-matrix(beta,nr=MC,nc=length(beta),byr=T));abline(h=0)

par(mfrow=c(2,1))
boxplot(res_no_outliers$jini-matrix(beta,nr=MC,nc=length(beta),byr=T));abline(h=0)
boxplot(res_outliers$jini-matrix(beta,nr=MC,nc=length(beta),byr=T));abline(h=0)

par(mfrow=c(2,1))
boxplot(res_no_outliers$initial-matrix(beta,nr=MC,nc=length(beta),byr=T));abline(h=0)
boxplot(res_outliers$initial-matrix(beta,nr=MC,nc=length(beta),byr=T));abline(h=0)

par(mfrow=c(2,1))
boxplot(res_no_outliers$consistent-matrix(beta,nr=MC,nc=length(beta),byr=T),outline=F);abline(h=0)
boxplot(res_outliers$consistent-matrix(beta,nr=MC,nc=length(beta),byr=T),outline=F);abline(h=0)

par(mfrow=c(2,1))
boxplot(res_no_outliers$robBY-matrix(beta,nr=MC,nc=length(beta),byr=T));abline(h=0)
boxplot(res_outliers$robBY-matrix(beta,nr=MC,nc=length(beta),byr=T));abline(h=0)
