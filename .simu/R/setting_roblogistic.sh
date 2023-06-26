SETTING=1		 
N=200			 
P=20			 
MODEL=roblogistic	 
H=100			 
B=100			 
MC=1000			 
SEED=895724		 
SIGMA="Sigma <- t(chol(toeplitz(0.8^(0:(p-1)))))" 
DESIGN="x <- matrix(rnorm(n*p,sd= 4/sqrt(n)), nr=n) %*% Sigma"
BETA="beta <- c(c(0.3, -2, -4),rep(0,p-2))" 
C=2.8			 
N_ARRAY=1000		 
