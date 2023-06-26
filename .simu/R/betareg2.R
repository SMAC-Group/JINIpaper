# -----------
# Simulations
# -----------
# general setting

# packages
library(JINIpaper)
library(stocapp)
library(ib)
library(betareg)

readRenviron(".simu/R/setting_betareg.sh")

# simulation specifics
MC <- as.integer(Sys.getenv("MC")) # number of simulations
model0 <- as.character(Sys.getenv("MODEL"))
model <- paste0(model0,"_2")

##------------------ Setting ----------------
setting <- as.integer(Sys.getenv("SETTING"))
n <- as.integer(Sys.getenv("N"))
p <- as.integer(Sys.getenv("P"))
# beta <- str2expression(Sys.getenv("BETA"))
beta <- str2expression("beta <- c(-.5, rep(1,5), rep(-1.5,5), rep(2,5), rep(c(-.1,.1),p-15)[1:(p-15)])")
eval(beta)
# gamma <- as.numeric(Sys.getenv("GAMMA"))
gamma = 5
theta <- c(beta, log(gamma))
Sigma <- str2expression(Sys.getenv("SIGMA"))
eval(Sigma)
# design <- str2expression(Sys.getenv("DESIGN"))
design <- str2expression("x <- matrix(rnorm(n*p,sd= 2/sqrt(n)), nr=n)")

set.seed(as.integer(Sys.getenv("SEED")))
seed <- vector("list",3)
seed$process <- sample.int(1e7,MC)
seed$design <- sample.int(1e7,MC)
seed$sc <- sample.int(1e7,MC)

##------------------ Regular beta regression --------
# control for IB
jimi_control <- stocappControl(constraint=FALSE, maxit=1e4)
H <- as.integer(Sys.getenv("H"))
bbc_control <- ibControl(H=H, maxit=1L, constraint=FALSE)

res <- list(jimi = matrix(ncol=p+2,nrow=MC),
            mle = matrix(ncol=p+2,nrow=MC),
            bbc = matrix(ncol=p+2,nrow=MC),
            time = matrix(ncol=3, nrow=MC)) 
colnames(res$time) <- c("mle", "jimi", "bbc") 

##------------------ Slurm specs --------------
n_array <- as.integer(Sys.getenv("N_ARRAY"))
ind <- matrix(seq_len(MC),nr=n_array,byr=T)
id_slurm <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

for(m in na.omit(ind[id_slurm,])){
  ##------- simulate the process ---------
  # set the seed
  set.seed(seed$design[m])
  # simulate the design
  eval(design)
  betareg_object <- make_betareg(x, theta)
  y <- simulation(betareg_object, control = list(seed=seed$process[m]))
  
  ##------ Initial (MLE) ----------------
  t1 <- Sys.time()
  fit_mle <- NULL
  try(fit_mle <- betareg(y ~ x), silent=T)
  t2 <- Sys.time()
  if(is.null(fit_mle)) next
  res$mle[m,] <- coef(fit_mle)
  res$time[,"mle"][m] <- difftime(t2,t1,units="secs")

  ##------ Bootstrap bias correction -----
  t1 <- Sys.time()
  bbc_control$seed <- seed$ib[m]
  t2 <- Sys.time()
  fit_bbc <- ib(fit_mle, control = bbc_control)
  res$bbc[m,] <- ib::getEst(fit_bbc)
  res$time[,"bbc"][m] <- difftime(t2,t1,units="secs")
 
  ##------ IB bias correction ------------
  jimi_control$seed <- seed$sc[m]
  t1 <- Sys.time()
  fit_jimi <- stocapp(fit_mle, control = jimi_control)
  t2 <- Sys.time()
  res$jimi[m,] <- stocapp::getEst(fit_jimi)
  res$time[,"jimi"][m] <- difftime(t2,t1,units="secs")
  
  # save results
  save(res, file=paste0("tmp/",model,"_setting_",setting,"_id_",id_slurm,".rds"))
  cat(m,"\n")
}

