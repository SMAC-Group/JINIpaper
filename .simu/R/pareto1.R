# -----------
# Simulations
# -----------
# general setting

# packages
library(JINIpaper)
library(VGAM)

## simulation specifics
MC <- as.integer(Sys.getenv("MC")) # number of simulations
model0 <- as.character(Sys.getenv("MODEL"))
model <- paste0(model0,"_1")

##------------------ Specific function ----------------
starting_values <- function(y, p){
  x <- log(y[-which.min(y)])
  stat <- boxplot.stats(x)
  c(exp(1/mean(x[!x%in%stat$out])), rep(0,p))
}

##------------------ Setting ----------------
setting <- as.integer(Sys.getenv("SETTING"))
n <- as.integer(Sys.getenv("N"))
p <- as.integer(Sys.getenv("P"))
beta <- str2expression(Sys.getenv("BETA"))
eval(beta)
scale <- as.numeric(Sys.getenv("SCALE"))
Sigma <- str2expression(Sys.getenv("SIGMA"))
eval(Sigma)
design <- str2expression(Sys.getenv("DESIGN"))
cc <- as.numeric(Sys.getenv("C"))
H <- as.integer(Sys.getenv("H"))

set.seed(as.integer(Sys.getenv("SEED")))
seed <- vector("list",3)
seed$process <- sample.int(1e7,MC)
seed$design <- sample.int(1e7,MC)
seed$sc <- sample.int(1e7,MC)

##------------------ Regular logistic --------
res <- list(mle = matrix(ncol=p+2,nrow=MC),
            initial = matrix(ncol=p+2,nrow=MC),
            jimi_rob = matrix(ncol=p+2,nrow=MC),
            jimi_mle = matrix(ncol=p+2,nrow=MC),
	    time = matrix(ncol=4, nrow=MC)) 
colnames(res$time) <- c("mle", "initial", "jimi_rob", "jimi_mle")

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
  # simulate logistic
  shapes <- exp(cbind(1,x) %*% beta)
  y <- rpareto(n, scale = scale, shape = shapes)
  n_out <- round(.04 * n)
  y[sample.int(n, n_out)] <- rpareto(n_out, scale = 8 * scale, shape = 3)

  ##------ Starting value ----------------
  sv <- starting_values(y, p)
  
  ##------ Initial estimator (inconsistent, Tukey's weights) ----------------
  t1 <- Sys.time()
  fit_initial <- paretoWmle1(y, x, start = sv, c = cc, maxit = 1e3) 
  t2 <- Sys.time()
  if(fit_initial$conv != 0) next
  if(any(abs(fit_initial$coefficients)>1e2)) next
  res$initial[m,] <- c(fit_initial$coefficients, fit_initial$scale)
  res$time[,"initial"][m] <- difftime(t2,t1,units="secs")

  ##------ MLE estimation ----------------
  t1 <- Sys.time()
  fit_mle <- paretoMle(y, x, start = sv, maxit = 1e3)
  t2 <- Sys.time()
  if(fit_mle$conv == 0) {
	  res$mle[m,] <- c(fit_mle$coefficients, fit_mle$scale)
	  res$time[,"mle"][m] <- difftime(t2,t1,units="secs")
  }

  ##------ Iterative bootstrap bias correction ------------
  if(fit_mle$conv == 0) {
  t1 <- Sys.time()
  fit_jimi_mle <- paretoWmle1_ib(x, thetastart=res$mle[m,], c=Inf, seed=seed$sc[m], H=H)
  t2 <- Sys.time()
  if(!is.finite(fit_jimi_mle$test_theta)) next
  res$jimi_mle[m,] <- fit_jimi_mle$estimate
  res$time[,"jimi_mle"][m] <- difftime(t2,t1,units="secs")
  }

  ##------ Iterative bootstrap bias correction ------------
  t1 <- Sys.time()
  fit_jimi_rob <- paretoWmle1_ib(x, thetastart=res$initial[m,], c=cc, seed=seed$sc[m], H=H)
  t2 <- Sys.time()
  if(!is.finite(fit_jimi_rob$test_theta)) next
  res$jimi_rob[m,] <- fit_jimi_rob$estimate
  res$time[,"jimi_rob"][m] <- difftime(t2,t1,units="secs")
  
  # save results
  save(res, file=paste0("tmp/",model,"_setting_",setting,"_id_",id_slurm,".rds"))
  cat(m,"\n")
}
