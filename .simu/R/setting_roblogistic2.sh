export SETTING=3		 # setting number (internal meaning)
export N=800			 # sample size
export P=40			 # number of regression coefficients 
export MODEL=roblogistic	 # model studied
export H=100			 # number of Monte Carlo replicates for IB
export B=100			 # number of bootstrap replicates
export MC=1000			 # number of simulated samples
export SEED=895724		 # seed for reproducibility
export SIGMA="Sigma <- t(chol(toeplitz(0.8^(0:(p-1)))))" 
				 # The way variance of design is computed
export DESIGN="x <- matrix(rnorm(n*p,sd= 4/sqrt(n)), nr=n) %*% Sigma"
				 # The way the design is generated
export BETA="beta <- c(c(0.3, -2, -4),rep(0,p-2))" 
export C=Inf		 # tuning parameter for robust estimator
export N_ARRAY=1000		 # number of array for HPC  
