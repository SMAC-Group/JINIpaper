# -----------
# Simulations
# -----------
# general setting

# packages
library(JINIpaper)
library(stocapp)
library(ib)
library(betareg)

readRenviron(".simu/R/setting.sh")

# simulation specifics
MC <- as.integer(Sys.getenv("MC")) # number of simulations
model0 <- as.character(Sys.getenv("MODEL"))
model <- paste0(model0,"_2")

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
seed$ib <- sample.int(1e7,MC)

##------------------ Rounding --------
# control for IB
jimi_control <- stocappControl(constraint=FALSE, maxit=1e4)
H <- as.integer(Sys.getenv("H"))
bbc_control <- ibControl(H=H, maxit=1L, constraint=FALSE)
load(paste0(".simu/data/",model,"_setting_",setting,".rds"))
jimi_sd <- matrix(ncol=p+2,nrow=MC)
mle_sd <- matrix(ncol=p+2,nrow=MC)
bbc_sd <- matrix(ncol=p+2,nrow=MC)
time_sd <- matrix(ncol=3,nrow=MC)
colnames(time_sd) <- c("mle", "jimi", "bbc") 

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
 if(!any(is.na(res$mle[m,]))){
	boot <- matrix(ncol=p+2, nrow=B)
	t1 <- Sys.time()
	for(b in 1:B) {
	    beta <- res$mle[m,]
	    beta[p+2] <- log(beta[p+2])
  		betareg_object <- make_betareg(x, beta)
  		y <- simulation(betareg_object, control = list(seed=seed$process[m]+b))
  		fit_mle <- NULL
		try(fit_mle <- betareg(y ~ x), silent=T)
  		if(is.null(fit_mle)) next
  		boot[b,] <- coef(fit_mle)
	}
  	t2 <- Sys.time()
 	mle_sd[m,] <- apply(boot,2,sd,na.rm=TRUE)	
  	time_sd[,"mle"][m] <- difftime(t2,t1,units="secs")
  }
  res$mle_sd <- mle_sd

  ##------ Bootstrap bias correction -----
  if(!any(is.na(res$bbc[m,]))){
	boot <- matrix(ncol=p+2, nrow=B)
	t1 <- Sys.time()
	for(b in 1:B) {
  		betareg_object <- make_betareg(x, res$bbc[m,])
  		y <- simulation(betareg_object, control = list(seed=seed$process[m]+b))
  		fit_mle <- NULL
		try(fit_mle <- betareg(y ~ x), silent=T)
  		if(is.null(fit_mle)) next
  		bbc_control$seed <- seed$ib[m]+b
  		fit_bbc <- NULL
		try(fit_bbc <- ib(fit_mle, control=bbc_control), silent=T)
  		if(is.null(fit_bbc)) next
  		boot[b,] <- ib::getEst(fit_bbc)
	}
  	t2 <- Sys.time()
 	bbc_sd[m,] <- apply(boot,2,sd,na.rm=TRUE)	
  	time_sd[,"bbc"][m] <- difftime(t2,t1,units="secs")
  }
  res$bbc_sd <- bbc_sd
  
  ##------ IB bias correction ------------
  if(!any(is.na(res$jimi[m,]))){
	boot <- matrix(ncol=p+2, nrow=B)
	t1 <- Sys.time()
	for(b in 1:B) {
  		betareg_object <- make_betareg(x, res$jimi[m,])
  		y <- simulation(betareg_object, control = list(seed=seed$process[m]+b))
  		fit_mle <- NULL
		try(fit_mle <- betareg(y ~ x), silent=T)
  		if(is.null(fit_mle)) next
  		jimi_control$seed <- seed$sc[m]+b
  		fit_jimi <- NULL
		try(fit_jimi <- stocapp(fit_mle, control = jimi_control),silent=T)
  		if(is.null(fit_jimi)) next
  		boot[b,] <- stocapp::getEst(fit_jimi)
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

