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

set.seed(as.integer(Sys.getenv("SEED")))
seed <- vector("list",3)
seed$process <- sample.int(1e7,MC)
seed$design <- sample.int(1e7,MC)
seed$sc <- sample.int(1e7,MC)

##------------------ Rounding --------
# control for IB
H <- as.integer(Sys.getenv("H"))
jimi_control <- ibControl(sim=rounding, constraint=FALSE, maxit=1e2)

res <- list(jimi = matrix(ncol=p+2,nrow=MC),
            initial = matrix(ncol=p+2,nrow=MC),
            consistent = matrix(ncol=p+2,nrow=MC),
            time = matrix(ncol=3, nrow=MC)) 
colnames(res$time) <- c("initial", "jimi", "consistent") 

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
  betareg_object <- make_betareg(x, theta)
  y <- simulation(betareg_object, control = list(sim=rounding, seed=seed$process[m]))
  
  ##------ Initial (MLE) ----------------
  t1 <- Sys.time()
  fit_mle <- NULL
  try(fit_mle <- betareg(y ~ x), silent=T)
  t2 <- Sys.time()
  if(is.null(fit_mle)) next
  res$initial[m,] <- coef(fit_mle)
  res$time[,"initial"][m] <- difftime(t2,t1,units="secs")

  ##------ EM estimation (consistent) ----------------
  fit_em <- NULL
  sv <- res$initial[m,]
  sv[p+2] <- log(sv[p+2])
  t1 <- Sys.time()
  try(fit_em <- em(y,cbind(1,x),sv), silent=T)
  t2 <- Sys.time()
  if(!is.null(fit_em)){
    res$consistent[m,] <- fit_em$par
    res$consistent[m,p+2] <- exp(fit_em$par[p+2])
    }
  res$time[,"consistent"][m] <- difftime(t2,t1,units="secs")
  
  ##------ IB bias correction ------------
  jimi_control$seed <- seed$sc[m]
  t1 <- Sys.time()
  fit_jimi <- ib(fit_mle, control = jimi_control)
  t2 <- Sys.time()
  res$jimi[m,] <- getEst(fit_jimi)
  res$time[,"jimi"][m] <- difftime(t2,t1,units="secs")
  
  # save results
  save(res, file=paste0("tmp/",model,"_setting_",setting,"_id_",id_slurm,".rds"))
  cat(m,"\n")
}

