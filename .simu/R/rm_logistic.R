# -----------
# Simulations
# -----------
# general setting

# packages
library(JINIpaper)
library(stocapp)
library(robustbase)
library(brglm)

## simulation specifics
MC <- as.integer(Sys.getenv("MC")) # number of simulations
model0 <- as.character(Sys.getenv("MODEL"))
model <- paste0(model0,"_1")

##------------------ Specific functions -------
## Function to simulate responses with outliers
## see ?ib::simulation and ?ib::control for more details 
outlying_mechanism <- function(object, control, extra){
  ftd <- fitted(object)
  n <- length(ftd)
  N <- n * control$H
  n_out <- round(n * control$eps)

  # misclassification mechanism
  u1 <- 0.0
  u2 <- FN
  mu_star <- u1 * (1 - rep(ftd,control$H)) + (1 - u2) * rep(ftd,control$H)
  y <- matrix(ifelse(runif(N) > rep(ftd,control$H), FALSE, TRUE), ncol=control$H)

  # adverserial approach:
  select_outliers <- apply(matrix(mu_star, ncol = control$H), 2, function(x) order(abs(x - 0.5), decreasing = TRUE)[1:n_out])
  # select_outliers <- order(abs(fitted(object)-0.5), decreasing = T)[1:n_out]
  # adverserial approach 2:
  # select_outliers <- order(fitted(object), decreasing = T)[1:n_out]
  # totally at random approach:
  # select_outliers <- sample.int(n, n_out)
  y[select_outliers] <- !y[select_outliers]
  as.integer(y)
}

r_logistic_rm <- function(thetas, x, seed, FN, FP, eps){
  set.seed(seed)
  n <- nrow(x)
  p <- ncol(x)
  eta <- cbind(1,x) %*% thetas
  prob <- exp(eta) / (1 + exp(eta))
  mu_star <- FP * (1 - prob) + (1 - FN) * prob
  # y <- rbinom(n, size=1, prob=mu_star)
  y <- ifelse(runif(n) > mu_star, FALSE, TRUE)
  
  n_out <- round(n * eps)
  select_outliers <- order(abs(mu_star - 0.5), decreasing = TRUE)[1:n_out]
  y[select_outliers] <- !y[select_outliers]
  as.integer(y)
}
##------------------ Setting ----------------
setting <- as.integer(Sys.getenv("SETTING"))
n <- as.integer(Sys.getenv("N"))
p <- as.integer(Sys.getenv("P"))
beta <- str2expression(Sys.getenv("BETA"))
eval(beta)
Sigma <- str2expression(Sys.getenv("SIGMA"))
eval(Sigma)
design <- str2expression(Sys.getenv("DESIGN"))
cc <- as.numeric(Sys.getenv("C"))
FN <- as.numeric(Sys.getenv("FN"))

set.seed(as.integer(Sys.getenv("SEED")))
seed <- list()
seed$process <- sample.int(1e7,MC)
seed$design <- sample.int(1e7,MC)
seed$sc <- sample.int(1e7,MC)

##------------------ Regular logistic --------
res <- list(mle = matrix(ncol=p+1,nrow=MC),
            br = matrix(ncol=p+1,nrow=MC),
            initial = matrix(ncol=p+1,nrow=MC),
            jimi = matrix(ncol=p+1,nrow=MC),
	    time = matrix(ncol=7, nrow=MC), 
            consistent = matrix(ncol=p+1,nrow=MC),
            robCR = matrix(ncol=p+1,nrow=MC),
            robBY = matrix(ncol=p+1,nrow=MC))
colnames(res$time) <- c("mle", "br", "initial", "jimi", "consistent", "robCR", "robBY") 

##------------------ Slurm specs --------------
n_array <- as.integer(Sys.getenv("N_ARRAY"))
ind <- matrix(seq_len(MC),nr=n_array,byr=T)
id_slurm <- Sys.getenv("SLURM_ARRAY_TASK_ID")

for(m in na.omit(ind[as.numeric(id_slurm),])){
 ##------- simulate the process ---------
  # set the seed
  set.seed(seed$process[m])
  # simulate the design
  eval(design)
  # logistic_object <- make_logistic(x, beta, robust=TRUE)
  # simulate logistic
  # y <- simulation(logistic_object, 
  #                 control = list(seed=seed$process[m], sim=outlying_mechanism, eps=0.00)) 
  y <- r_logistic_rm(thetas=beta, x=x, seed=seed$process[m], FN=FN, FP=0.0, eps=0.00)
  
  ##------ MLE estimation ----------------
  t1 <- Sys.time()
  fit_mle <- glm(y ~ x, family=binomial(link="logit"), control = glm.control(maxit=200))
  t2 <- Sys.time()
  res$mle[m,] <- coef(fit_mle)
  res$time[,"mle"][m] <- difftime(t2,t1,units="secs")

  ##------ Consistent MLE  ----------------
  t1 <- Sys.time()
  fit_consistent <- logistic_misclassification_mle(x, y, fp = 0, fn = FN)
  t2 <- Sys.time()
  res$consistent[m,] <- fit_consistent
  res$time[,"consistent"][m] <- difftime(t2,t1,units="secs")

#  ##------ "brglm" bias reduction --------
#  t1 <- Sys.time()
#  fit_br <- brglm(y~x, family=binomial(link="logit"), control.glm = glm.control1(maxit=200))
#  t2 <- Sys.time()
#  res$br[m,] <- coef(fit_br)
#  res$time[,"br"][m] <- difftime(t2,t1,units="secs")
#  
#  ##------ Robust estimator (Cantoni-Ronchetti) ----------------
#  t1 <- Sys.time()
#  fit_robCR <- NULL
#  try(fit_robCR <- glmrob(y ~ x, family=binomial(link="logit"), control = glmrobMqle.control(maxit=200), method = "Mqle"), silent=TRUE)
#  t2 <- Sys.time()
#  if(!is.null(fit_robCR)) res$robCR[m,] <- coef(fit_robCR)
#  res$time[,"robCR"][m] <- difftime(t2,t1,units="secs")
#  
  ##------ Robust estimator (Branco-Yohai) ----------------
  t1 <- Sys.time()
  fit_robBY <- glmrob(y ~ x, family=binomial(link="logit"), control = glmrobBY.control(maxit=200), method = "BY")
  t2 <- Sys.time()
  res$robBY[m,] <- coef(fit_robBY)
  res$time[,"robBY"][m] <- difftime(t2,t1,units="secs")
  
  ##------ Initial estimator (inconsistent, Tukey's weights) ----------------
  t1 <- Sys.time()
  fit_initial <- roblogisticWmle1(y, x, start = coef(fit_mle), c = cc) 
  t2 <- Sys.time()
  res$initial[m,] <- fit_initial$coefficients
  res$time[,"initial"][m] <- difftime(t2,t1,units="secs")
#
#  ##------ Consistent WMLE (Tukey's weights) ----------------
#  t1 <- Sys.time()
#  fit_consistent <- roblogisticWmle(y, x, start = coef(fit_mle), c = cc) 
#  t2 <- Sys.time()
#  res$consistent[m,] <- fit_consistent$coefficients
#  res$time[,"consistent"][m] <- difftime(t2,t1,units="secs")
#  
  ##------ Iterative bootstrap bias correction ------------
  t1 <- Sys.time()
  #fit_jimi <- roblogisticWmle1_ib(x, thetastart=fit_initial$coefficients, c=cc, seed=seed$sc[m])
  fit_jimi <- robmisclogisticWmle1_ib(x, thetastart=fit_initial$coefficients, c=cc, seed=seed$sc[m], FN=FN)
  t2 <- Sys.time()
  if(!is.finite(fit_jimi$test_theta)) next
  res$jimi[m,] <- fit_jimi$estimate
  res$time[,"jimi"][m] <- difftime(t2,t1,units="secs")
#  
#  # save results
#  save(res, file=paste0("tmp/",model,"_setting_",setting,"_id_",id_slurm,".rds"))
  save(res, file=paste0("tmp/",model,"_setting_",setting,"_id_",id_slurm,".rds"))
  cat(m,"\n")
}
