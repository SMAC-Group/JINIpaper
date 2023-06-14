# -----------
# Simulations
# -----------
# general setting

# packages
library(JINIpaper)
library(ib)
library(betareg)

readRenviron(".simu/R/setting_betareg.sh")

# simulation specifics
MC <- as.integer(Sys.getenv("MC")) # number of simulations
model0 <- as.character(Sys.getenv("MODEL"))
model <- paste0(model0,"_1")

##------------------ Specific functions -------
# Function to simulate responses with rounding
rounding <- function(object, control, extra=NULL){
  simulate_betareg <- getFromNamespace("simulate_betareg", ns = "ib")
  y <- simulate_betareg(object, control$H)
  n <- length(y) / control$H
  y <- round(y,1)
  y <- (y*(n-1) + 0.5)/n
  matrix(y, ncol=control$H)
}

##------------------ Setting ----------------
setting <- as.integer(Sys.getenv("SETTING"))
n <- as.integer(Sys.getenv("N"))
p <- as.integer(Sys.getenv("P"))
beta <- str2expression(Sys.getenv("BETA"))
eval(beta)
gamma <- as.numeric(Sys.getenv("GAMMA"))
theta <- c(beta, log(gamma))
Sigma <- str2expression(Sys.getenv("SIGMA"))
eval(Sigma)
design <- str2expression(Sys.getenv("DESIGN"))
B <- as.integer(Sys.getenv("B"))

set.seed(as.integer(Sys.getenv("SEED")))
seed <- vector("list",3)
seed$process <- sample.int(1e7,MC)
seed$design <- sample.int(1e7,MC)
seed$sc <- sample.int(1e7,MC)

##------------------ Rounding --------
# control for IB
H <- as.integer(Sys.getenv("H"))
jimi_control <- ibControl(sim=rounding, H=H, constraint=FALSE)
load(paste0(".simu/data/",model,"_setting_",setting,".rds"))
jimi_sd <- matrix(ncol=p+2,nrow=MC)
initial_sd <- matrix(ncol=p+2,nrow=MC)
consistent_sd <- matrix(ncol=p+2,nrow=MC)
time_sd <- matrix(ncol=3,nrow=MC)
colnames(time_sd) <- c("initial", "jimi", "consistent") 

##------------------ Slurm specs --------------
n_array <- as.integer(Sys.getenv("N_ARRAY"))
ind <- matrix(seq_len(MC),nr=n_array,byr=T)
id_slurm <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

for(m in na.omit(ind[id_slurm,])){
  ##------- simulate rounding process ---------
  # set the seed
  set.seed(seed$design[m])
  # simulate the design
  eval(design)
  
  ##------ Initial (MLE) ----------------
# if(!any(is.na(res$initial[m,]))){
#	boot <- matrix(ncol=p+2, nrow=B)
#	beta <- res$initial[m,]
#	beta[p+2] <- log(beta[p+2])
#  	betareg_object <- make_betareg(x, beta)
#	t1 <- Sys.time()
#	for(b in 1:B) {
#  		y <- simulation(betareg_object, control = list(sim=rounding, seed=seed$process[m]+b))
#  		fit_mle <- NULL
#		try(fit_mle <- betareg(y ~ x), silent=T)
#  		if(is.null(fit_mle)) next
#  		boot[b,] <- coef(fit_mle)
#	}
#  	t2 <- Sys.time()
# 	initial_sd[m,] <- apply(boot,2,sd,na.rm=TRUE)	
#  	time_sd[,"initial"][m] <- difftime(t2,t1,units="secs")
#  }
#  res$initial_sd <- initial_sd

  ##------ EM estimation (consistent) ----------------
  if(!any(is.na(res$consistent[m,]))){
	boot <- matrix(ncol=p+2, nrow=B)
	beta <- res$consistent[m,]
	beta[p+2] <- log(beta[p+2])
  	betareg_object <- make_betareg(x, beta)
	t1 <- Sys.time()
	for(b in 1:B) {
  		y <- simulation(betareg_object, control = list(sim=rounding, seed=seed$process[m]+b))
  		fit_mle <- NULL
		try(fit_mle <- betareg(y ~ x), silent=T)
  		if(is.null(fit_mle)) next
  		fit_consistent <- NULL
  		sv <- coef(fit_mle)
		sv[p+2] <- log(sv[p+2])
		try(fit_consistent <- em(y, cbind(1,x), sv), silent=T)
  		if(is.null(fit_consistent)) next
  		boot[b,] <- fit_consistent$par
		boot[b,p+2] <- exp(fit_consistent$par[p+2])
	}
  	t2 <- Sys.time()
 	consistent_sd[m,] <- apply(boot,2,sd,na.rm=TRUE)	
  	time_sd[,"consistent"][m] <- difftime(t2,t1,units="secs")
  }
  res$consistent_sd <- consistent_sd
  
  ##------ IB bias correction ------------
  if(!any(is.na(res$jimi[m,]))){
	boot <- matrix(ncol=p+2, nrow=B)
	beta <- res$jimi[m,]
	beta[p+2] <- log(beta[p+2])
  	betareg_object <- make_betareg(x, beta)
	t1 <- Sys.time()
	for(b in 1:B) {
  		y <- simulation(betareg_object, control = list(sim=rounding, seed=seed$process[m]+b))
  		fit_mle <- NULL
		try(fit_mle <- betareg(y ~ x), silent=T)
  		if(is.null(fit_mle)) next
  		jimi_control$seed <- seed$sc[m]+b
  		fit_jimi <- NULL
		try(fit_jimi <- ib(fit_mle, control = jimi_control),silent=T)
  		if(is.null(fit_jimi)) next
  		boot[b,] <- getEst(fit_jimi)
	}
  	t2 <- Sys.time()
 	jimi_sd[m,] <- apply(boot,2,sd,na.rm=TRUE)	
  	time_sd[,"jimi"][m] <- difftime(t2,t1,units="secs")
  }
  res$jimi_sd <- jimi_sd
  
  # save results
  res$time_sd <- time_sd
  save(res, file=paste0("tmp/",model,"_sd_setting_",setting,"_id_",id_slurm,".rds"))
  cat(m,"\n")
}

