# -----------
# Simulations
# -----------
# general setting
readRenviron(".simu/R/setting.sh")

# packages
library(JINIpaper)
library(stocapp)
library(robustbase)
library(brglm)

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
B <- as.integer(Sys.getenv("B"))

set.seed(as.integer(Sys.getenv("SEED")))
seed <- vector("list",3)
seed$process <- sample.int(1e7,MC)
seed$design <- sample.int(1e7,MC)
seed$sc <- sample.int(1e7,MC)
seed$boot <- sample.int(1e7,MC)

##------------------ Regular logistic --------
load(paste0(".simu/data/",model,"_setting_",setting,".rds"))
mle_sd <- matrix(ncol=p+1,nrow=MC)
br_sd <- matrix(ncol=p+1,nrow=MC)
jimi_sd <- matrix(ncol=p+1,nrow=MC)
robBY_sd <- matrix(ncol=p+1,nrow=MC)
robCR_sd <- matrix(ncol=p+1,nrow=MC)
consistent_sd <- matrix(ncol=p+1,nrow=MC)
time_sd <- matrix(ncol=7, nrow=MC)
colnames(time_sd) <- c("mle", "br", "initial", "jimi", "consistent", "robCR", "robBY") 

##------------------ Slurm specs --------------
n_array <- as.integer(Sys.getenv("N_ARRAY"))
ind <- matrix(seq_len(MC),nr=n_array,byr=T)
id_slurm <- Sys.getenv("SLURM_ARRAY_TASK_ID")

load(".simu/data/roblogistic_2_sd_setting_1.rds")
which(is.na(res$mle_sd[,1]))
which(is.na(res$mle[,1]))
m <- 47

for(m in na.omit(ind[as.numeric(id_slurm),])){
 ##------- simulate the process ---------
  # set the seed
  set.seed(seed$process[m])
  # simulate the design
  eval(design)

  ##------ MLE estimation ----------------
  if(!any(is.na(res$mle[m,]))){
	boot <- matrix(ncol=p+1, nrow=B)
	t1 <- Sys.time()
	for(b in 1:B) {
		logistic_object <- make_logistic(x, res$mle[m,], robust=TRUE)
		y <- simulation(logistic_object, control = list(seed=seed$boot[m]+b))
		fit_mle <- glm(y ~ x, family=binomial(link="logit"), control = glm.control(maxit=200))
  		boot[b,] <- coef(fit_mle)
	}
  	t2 <- Sys.time()
 	mle_sd[m,] <- apply(boot,2,sd,na.rm=TRUE)	
  	time_sd[,"mle"][m] <- difftime(t2,t1,units="secs")
  }
  res$mle_sd <- mle_sd

  ##------ "brglm" bias reduction --------
  if(!any(is.na(res$br[m,]))){
	boot <- matrix(ncol=p+1, nrow=B)
	t1 <- Sys.time()
	for(b in 1:B) {
		logistic_object <- make_logistic(x, res$br[m,], robust=TRUE)
		y <- simulation(logistic_object, control = list(seed=seed$boot[m]+b))
  		fit_br <- brglm(y~x, family=binomial(link="logit"), control.glm = glm.control1(maxit=200))
  		boot[b,] <- coef(fit_br)
	}
  	t2 <- Sys.time()
 	br_sd[m,] <- apply(boot,2,sd,na.rm=TRUE)	
  	time_sd[,"br"][m] <- difftime(t2,t1,units="secs")
  }
  res$br_sd <- br_sd

  ##------ Robust estimator (Cantoni-Ronchetti) ----------------
  if(!any(is.na(res$robCR[m,]))){
	boot <- matrix(ncol=p+1, nrow=B)
	t1 <- Sys.time()
	for(b in 1:B) {
		logistic_object <- make_logistic(x, res$robCR[m,], robust=TRUE)
		y <- simulation(logistic_object, control = list(seed=seed$boot[m]+b))
  		fit_robCR <- NULL
		try(fit_robCR <- glmrob(y ~ x, family=binomial(link="logit"), control = glmrobMqle.control(maxit=200), method = "Mqle"), silent=TRUE)
		if(!is.null(fit_robCR)) boot[b,] <- coef(fit_robCR)
	}
  	t2 <- Sys.time()
 	robCR_sd[m,] <- apply(boot,2,sd,na.rm=TRUE)	
  	time_sd[,"robCR"][m] <- difftime(t2,t1,units="secs")
  }
  res$robCR_sd <- robCR_sd
  
  ##------ Robust estimator (Branco-Yohai) ----------------
   if(!any(is.na(res$robBY[m,]))){
	boot <- matrix(ncol=p+1, nrow=B)
	t1 <- Sys.time()
	for(b in 1:B) {
		logistic_object <- make_logistic(x, res$robBY[m,], robust=TRUE)
		y <- simulation(logistic_object, control = list(seed=seed$boot[m]+b))
  		try(fit_robBY <- glmrob(y ~ x, family=binomial(link="logit"), control = glmrobBY.control(maxit=200), method = "BY"), silent=TRUE)
  		if(!is.null(fit_robBY)) boot[b,] <- coef(fit_robBY)
	}
  	t2 <- Sys.time()
 	robBY_sd[m,] <- apply(boot,2,sd,na.rm=TRUE)	
  	time_sd[,"robBY"][m] <- difftime(t2,t1,units="secs")
  }
  res$robBY_sd <- robBY_sd

  ##------ Consistent WMLE (Tukey's weights) ----------------
  if(!any(is.na(res$consistent[m,]))){
	boot <- matrix(ncol=p+1, nrow=B)
	t1 <- Sys.time()
	for(b in 1:B) {
		logistic_object <- make_logistic(x, res$consistent[m,], robust=TRUE)
		y <- simulation(logistic_object, control = list(seed=seed$boot[m]+b))
  		fit_consistent <- roblogisticWmle(y, x, start = res$consistent[m,], c = cc) 
  		boot[b,] <- fit_consistent$coefficients
	}
  	t2 <- Sys.time()
 	consistent_sd[m,] <- apply(boot,2,sd,na.rm=TRUE)	
  	time_sd[,"consistent"][m] <- difftime(t2,t1,units="secs")
  }
  res$consistent_sd <- consistent_sd

  ##------ Iterative bootstrap bias correction ------------
  if(!any(is.na(res$jimi[m,]))){
	boot <- matrix(ncol=p+1, nrow=B)
	t1 <- Sys.time()
	for(b in 1:B) {
		logistic_object <- make_logistic(x, res$jimi[m,], robust=TRUE)
		y <- simulation(logistic_object, control = list(seed=seed$boot[m]+b))
  		fit_initial <- roblogisticWmle1(y, x, start = coef(fit_mle), c = cc) 
  		fit_jimi <- roblogisticWmle1_ib(x, thetastart=fit_initial$coefficients, c=cc, seed=seed$sc[m])
  		boot[b,] <- fit_jimi$estimate
	}
  	t2 <- Sys.time()
 	jimi_sd[m,] <- apply(boot,2,sd,na.rm=TRUE)	
  	time_sd[,"jimi"][m] <- difftime(t2,t1,units="secs")
  }
  res$jimi_sd <- jimi_sd
  
  # save results
  res$time_sd <-  time_sd
  save(res, file=paste0("tmp/",model,"_sd_setting_",setting,"_id_",id_slurm,".rds"))
  cat(m,"\n")
}
