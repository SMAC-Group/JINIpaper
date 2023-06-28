SETTING=1		 
N=150			 
P=40			 
MODEL=pareto		 
B=100			 
H=100			 
MC=1000			 
SEED=895724		 
SIGMA="Sigma <- diag(p)" 
DESIGN="x <- matrix(rnorm(n*p,sd= 1/sqrt(n)), nr=n) %*% Sigma"
BETA="beta <- c(1.5, rep(-2,2), rep(2,2), rep(c(-.1,.1),p-4)[1:(p-4)])" 
SCALE=5
N_ARRAY=1000		 
C=13			 