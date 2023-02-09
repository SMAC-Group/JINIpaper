# -----------
# Simulations
# -----------
# packages
library(MASS)
library(JINIpaper)
library(stocapp)

# simulation specifics
MC <- 1000 # number of simulations

# IB specifics
H <- 200 # number of simulated estimators

# Function to simulate responses with rounding
random_censoring <- function(object, control, extra=NULL){
  simulate_negbin <- getFromNamespace("simulate_negbin", ns = "stocapp")
  y <- simulate_negbin(object, control$H)
  u <- rpois(length(y), lambda = 3)
  y[y<=u] <- u[y<=u]
  matrix(y, ncol=control$H)
}

# -------------------------
# Case 3: beta regression 
# -------------------------
model <- "negbin_1"

##------------------ Setting  ----------------
setting <- as.integer(Sys.getenv("SETTING"))
n <- as.integer(Sys.getenv("N"))
p <- as.integer(Sys.getenv("P"))
beta_nz <- c(2.0,1.0,-1.0)
beta <- c(beta_nz,rep(0,p-2))
alpha <- 0.7
set.seed(895724)
seed <- vector("list",3)
seed$process <- sample.int(1e7,MC)
seed$design <- sample.int(1e7,MC)
seed$sc <- sample.int(1e7,MC)

##------------------ Rounding --------
# control for IB
jini_control <- stocappControl(sim=random_censoring, maxit = 1e4)

res <- list(jini = matrix(ncol=p+2,nrow=MC),
            initial = matrix(ncol=p+2,nrow=MC),
            consistentn = matrix(ncol=p+2,nrow=MC),
            time = rep(NA_real_,MC),
            iteration = rep(NA_integer_,MC))

##------------------ Slurm specs --------------
n_array <- 1000
ind <- matrix(seq_len(MC),nr=n_array,byr=T)
id_slurm <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

for(m in na.omit(ind[id_slurm,])){
  ##------- simulate rounding process ---------
  # set the seed
  set.seed(seed$design[m])
  # simulate the design
  x <- matrix(rnorm(n*p, sd = 2 * sqrt(2) / sqrt(p)), nr=n)
  negbin_object <- make_negbin(x, beta, 1 / alpha)
  # simulate negative binomial responses with random censoring
  y <- simulation(negbin_object, control = list(sim=random_censoring, seed=seed$process[m]))
  
  ##------ MLE estimation ----------------
  fit_mle <- NULL
  try(fit_mle <- glm.nb(y~x), silent=T)
  if(is.null(fit_mle)) next
  res$initial[m,1:(p+1)] <- coef(fit_mle)
  res$initial[m,p+2] <- 1 / fit_mle$theta

  ##------ MLE estimation (consistent) ----------------
  fit_em <- NULL
  try(fit_em <- glm_negbin(y, x, lambda = 3, maxit = 1000L, trace = T), silent = TRUE)
  #sv <- res$mle[m,]
  #sv[p+2] <- log(sv[p+2])
  #try(fit_em <- optim(par = sv, fn = nll_max, method = "BFGS", y = y, x = cbind(1,x),
	#	      lambda = 3, cens = 3), silent=T)
  if(!is.null(fit_em)){
    res$consistent[m,] <- fit_em$par
  #  res$em[m,p+2] <- exp(fit_em$par[p+2])
  }
  

  ##------ IB bias correction ------------
  jini_control$seed <- seed$sc[m]
  t1 <- Sys.time()
  fit_jini <- stocapp(fit_mle, control = jini_control, extra_param = TRUE)
  t2 <- Sys.time()
  res$jini[m,] <- getEst(fit_jini)
  res$iteration[m] <- getIteration(fit_jini)
  res$time[m] <- difftime(t2,t1,units = "secs")
  
  # save results
  save(res, file=paste0("tmp/",model,"_setting_",setting,"_id_",id_slurm,".rds"))
  cat(m,"\n")
}

