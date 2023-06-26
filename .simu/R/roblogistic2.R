# -----------
# Simulations
# -----------
# general setting

# packages
library(JINIpaper)
library(stocapp)
library(robustbase)
library(brglm)

readRenviron(".simu/R/setting_roblogistic.sh")

## simulation specifics
MC <- as.integer(Sys.getenv("MC")) # number of simulations
model0 <- as.character(Sys.getenv("MODEL"))
model <- paste0(model0,"_2")

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

set.seed(as.integer(Sys.getenv("SEED")))
seed <- vector("list",3)
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
  logistic_object <- make_logistic(x, beta, robust=TRUE)
  # simulate logistic
  y <- simulation(logistic_object, 
                  control = list(seed=seed$process[m])) 

  ##------ MLE estimation ----------------
  t1 <- Sys.time()
  fit_mle <- glm(y ~ x, family=binomial(link="logit"), control = glm.control(maxit=200))
  t2 <- Sys.time()
  res$mle[m,] <- coef(fit_mle)
  res$time[,"mle"][m] <- difftime(t2,t1,units="secs")

  ##------ "brglm" bias reduction --------
  t1 <- Sys.time()
  fit_br <- brglm(y~x, family=binomial(link="logit"), control.glm = glm.control1(maxit=200))
  t2 <- Sys.time()
  res$br[m,] <- coef(fit_br)
  res$time[,"br"][m] <- difftime(t2,t1,units="secs")
  
  ##------ Robust estimator (Cantoni-Ronchetti) ----------------
  t1 <- Sys.time()
  fit_robCR <- NULL
  try(fit_robCR <- glmrob(y ~ x, family=binomial(link="logit"), control = glmrobMqle.control(maxit=200), method = "Mqle"), silent=TRUE)
  t2 <- Sys.time()
  if(!is.null(fit_robCR)) res$robCR[m,] <- coef(fit_robCR)
  res$time[,"robCR"][m] <- difftime(t2,t1,units="secs")
  
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

  ##------ Consistent WMLE (Tukey's weights) ----------------
  t1 <- Sys.time()
  fit_consistent <- roblogisticWmle(y, x, start = coef(fit_mle), c = cc) 
  t2 <- Sys.time()
  res$consistent[m,] <- fit_consistent$coefficients
  res$time[,"consistent"][m] <- difftime(t2,t1,units="secs")
  
  ##------ Iterative bootstrap bias correction ------------
  t1 <- Sys.time()
  fit_jimi <- roblogisticWmle1_ib(x, thetastart=fit_initial$coefficients, c=cc, seed=seed$sc[m])
  t2 <- Sys.time()
  res$jimi[m,] <- fit_jimi$estimate
  res$time[,"jimi"][m] <- difftime(t2,t1,units="secs")
  
  # save results
  save(res, file=paste0("tmp/",model,"_setting_",setting,"_id_",id_slurm,".rds"))
  cat(m,"\n")
}
