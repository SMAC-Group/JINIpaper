SETTING=1
N=200	
P=50	
MODEL=betareg	
B=100	
H=100
MC=1000	
SEED=895724
SIGMA="Sigma <- diag(p)" 
DESIGN="x <- matrix(rnorm(n*p,sd= 2/sqrt(n)), nr=n) %*% Sigma"
BETA="beta <- c(-.5, rep(1,5), rep(-2,5), rep(3,5), rep(c(-.1,.1),p-15)[1:(p-15)])" 
GAMMA=10
N_ARRAY=1000	
