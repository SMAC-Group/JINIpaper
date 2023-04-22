###################################################################
### Simulation study on logistic regression (with approximated NR)
###################################################################

library(MASS)  # for mvrnorm()
library(Rcpp)
library(RcppArmadillo)
library(brglm)  # for brglm()
library(robustbase)  # for glmrob()

sourceCpp(".simu/src/logistic_approxNR.cpp")

#######################
### Simulate data
#######################

simu_X = function(n, p, cont_X){
  ### simulate X
  X = matrix(NA, nrow=n, ncol=p)
  
  # Bernoulli ones
  X[,1] = rbinom(n=n, size=1, p=0.5)
  
  # Normal ones
  Sigma = matrix(NA, nrow=(p-1), ncol=(p-1))
  for (i in 1:(p-1)) {
    for (j in 1:(p-1)) {
      Sigma[i,j] = 4/sqrt(n)*0.8^(abs(i-j))
    }
  }
  X[,2:p] = mvrnorm(n = n, mu = rep(0,(p-1)), Sigma = Sigma)
  X_nocont = X
  
  ### simulate contaminated X
  cont_individual = sample(1:n, ceiling(cont_X*n))
  cont_variable = sample(2:p, length(cont_individual), replace = TRUE)
  for (i in 1:length(cont_individual)) {
    myrow = cont_individual[i]
    mycol = cont_variable[i]
    X[myrow, mycol] = X[myrow, mycol] * 10
  }
  X_cont = X
  
  list(uncont = X_nocont, cont = X_cont)  
}



simu_data = function(X, betas, cont_y){
  ### simulate response y
  probs = expit_vec_cpp(X%*%betas)
  y = sapply(probs, FUN=rbinom, n=1, size=1)
  y_uncont = y
  
  ### contaminate observed y
  # compute probs of all observed data
  one_minus_probs = 1 - probs
  
  # create r = min(prob, 1-prob)
  n = length(y)
  r = rep(NA, n)
  for (i in 1:n) {
    r[i] = min(probs[i], one_minus_probs[i])
  }
  
  # pick con_y % of the smallest r to contaminate
  n_cont = ceiling(n*cont_y)
  cont_index = which(r <= sort(r, decreasing=FALSE)[n_cont])
  for (i in 1:n_cont) {
    index = cont_index[i]
    if(y[index]==1){ y[index]=0 }else{ y[index]=1 }
  }
  y_cont = y
  
  list(uncont = y_uncont, cont = y_cont)  
}


ggplot_like_colors = function(n, alpha = 1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}

#####################################################################
### Newton-Raphson to find pi (WMLE + IB, with weights on X all 1)
#####################################################################

NR_rob = function(y, X, crob, gamma, start, max_iter, stopping){
  p = ncol(X)
  max_iter = p*max_iter
  y = gamma*(1-y) + (1-gamma)*y
  
  if(sum(start=="glm") != 0){
    betas = glm(y ~ X, family=binomial(link="logit"))$coefficients
  }else{
    betas = start
  }
  
  for (i in 2:max_iter){
    betas_old = betas
    betas = make_rob_step_cpp(y=y, X=X, betas=betas_old, crob=crob)
    
    if(norm(betas_old-betas, "2")/p < stopping){
      return(list(betas=betas, iter=i, convergence=0))
    }
  }
  
  list(betas=betas, iter=i, convergence=1)
}

#########################################
### WMLE + IB (with weights on X all 1)
#########################################

find_pi = function(y, X, crob, gamma, start){
  NR_rob(y=y, X=X, crob=crob, gamma=gamma, start=start, 
         max_iter=200, stopping=10^(-6))
}

compute_bias = function(betas, X, crob, gamma, start, H){
  p = ncol(X)
  inter = rep(0,p)
  for (h in 1:H) {
    set.seed(h)
    y_sim = simu_data(X=X, betas=betas, cont_y=0.02)$uncont
    estim = find_pi(y=y_sim, X=X, crob=crob, gamma=gamma, start=start)$betas
    inter = inter + estim/H
  }
  betas-inter
}

ib = function(y, X, crob, gamma, start, H=100, iter_max=20, epsilon=10^(-6)){
  p = ncol(X)
  theta = matrix(NA, nrow=iter_max+1, ncol=p)
  theta[1,] = find_pi(y=y, X=X, crob=crob, gamma=gamma, start=start)$betas
  diff = norm(theta[1,], type="2")
  
  for (i in 1:iter_max) {
    theta[(i+1),] = theta[1,] + compute_bias(betas=theta[i,], X=X, crob=crob, gamma=gamma, start=start, H=H)
    diff_old = diff
    diff = norm(theta[(i+1),]-theta[i,], type="2")
    if (diff < epsilon || diff_old < diff){return(theta[(i+1),])}
  }
  
  return(theta[(i+1),])  
}

###############################################################
### Newton-Raphson to find pi (WMLE + IB, with weights on X)
###############################################################

NR_rob2 = function(y, X, crob, weightX, gamma, start, max_iter, stopping){
  p = ncol(X)
  max_iter = p*max_iter
  y = gamma*(1-y) + (1-gamma)*y
  
  if(sum(start=="glm") != 0){
    betas = glm(y ~ X-1, family=binomial(link="logit"))$coefficients
  }else{
    betas = start
  }
  
  for (i in 2:max_iter){
    betas_old = betas
    betas = make_rob_step2_cpp(y=y, X=X, betas=betas_old, crob=crob, weightX=weightX)
    
    if(norm(betas_old-betas, "2")/p < stopping){
      return(list(betas=betas, iter=i, convergence=0))
    }
  }
  
  list(betas=betas, iter=i, convergence=1)
}


#########################################
### WMLE + IB (with weights on X)
#########################################

find_pi2 = function(y, X, crob, weightX, gamma, start){
  NR_rob2(y=y, X=X, crob=crob, weightX=weightX, gamma=gamma, start=start, 
          max_iter=200, stopping=10^(-6))
}


compute_bias2 = function(betas, X, crob, weightX, gamma, start, H){
  p = ncol(X)
  inter = rep(0,p)
  for (h in 1:H) {
    set.seed(h)
    y_sim = simu_data(X=X, betas=betas, cont_y=0.02)$uncont
    estim = find_pi2(y=y_sim, X=X, crob=crob, weightX=weightX, gamma=gamma, start=start)$betas
    inter = inter + estim/H
  }
  betas-inter
}

ib2 = function(y, X, crob, weightX, gamma, start, H=100, iter_max=20, epsilon=10^(-6)){
  p = ncol(X)
  theta = matrix(NA, nrow=iter_max+1, ncol=p)
  theta[1,] = find_pi2(y=y, X=X, crob=crob, weightX=weightX, gamma=gamma, start=start)$betas
  diff = norm(theta[1,], type="2")
  
  for (i in 1:iter_max) {
    theta[(i+1),] = theta[1,] + compute_bias2(betas=theta[i,], X=X, crob=crob, weightX=weightX, gamma=gamma, start=start, H=H)
    diff_old = diff
    diff = norm(theta[(i+1),]-theta[i,], type="2")
    if (diff < epsilon || diff_old < diff){return(theta[(i+1),])}
  }
  
  return(theta[(i+1),])  
}

######################################
### Simulations
######################################

n = 200
betas_true = c(0,0.3,-2,-4,0,0)

c_glmrob = 2.12
c_ib = 2.21

crob = 3.88  # default value of c for IB

B = 200
gamma = 0.005

mle_nocont = mle_contX = mle_contY = mle_contXY = matrix(NA, B, length(betas_true))
brglm_nocont = brglm_contX = brglm_contY = brglm_contXY = matrix(NA, B, length(betas_true))
glmrob_nocont = glmrob_contX = glmrob_contY = glmrob_contXY = matrix(NA, B, length(betas_true))
glmrob2_nocont = glmrob2_contX = glmrob2_contY = glmrob2_contXY = matrix(NA, B, length(betas_true))
glmrob3_nocont = glmrob3_contX = glmrob3_contY = glmrob3_contXY = matrix(NA, B, length(betas_true))
glmrob4_nocont = glmrob4_contX = glmrob4_contY = glmrob4_contXY = matrix(NA, B, length(betas_true))
ib_nocont = ib_contX = ib_contY = ib_contXY = matrix(NA, B, length(betas_true))
ib2_nocont = ib2_contX = ib2_contY = ib2_contXY = matrix(NA, B, length(betas_true))

for (i in 1:B) {
  set.seed(525+i)
  
  X_list = simu_X(n = n, p = length(betas_true)-1, cont_X = 0.02)
  y_list = simu_data(X = cbind(1,X_list$uncont), betas = betas_true, cont_y = 0.02) # add intercept
  
  # ----- uncont
  X = X_list$uncont
  weightX = get_weightX_cpp(X)
  y = y_list$uncont
  
  # MLE
  mle_nocont[i,] = glm(y ~ X, family=binomial(link="logit"))$coefficients
  
  # brglm
  brglm_nocont[i,] = brglm(y ~ X, family=binomial)$coefficients
  
  # glmrob (with weight on X, default c)
  glmrob_try = try(glmrob(y ~ X, family=binomial(link="logit"), weights.on.x = "hat")$coefficients)
  if(!("try-error" %in% class(glmrob_try))){
    glmrob_nocont[i,] = glmrob_try
  }
  remove(glmrob_try)
  
  # glmrob (with weight on X, chosen c)
  glmrob_try = try(glmrob(y ~ X, family=binomial(link="logit"), weights.on.x = "hat",
                          control = glmrobMqle.control(tcc = c_glmrob))$coefficients)
  if(!("try-error" %in% class(glmrob_try))){
    glmrob2_nocont[i,] = glmrob_try
  }
  remove(glmrob_try)
  
  # glmrob (Bianco-Yohai method)
  glmrob_try = try(glmrob(y ~ X, method="BY", family=binomial, weights.on.x = "none")$coefficients)
  if(!("try-error" %in% class(glmrob_try))){
    glmrob3_nocont[i,] = glmrob_try
  }
  remove(glmrob_try)
  
  # glmrob (Bianco-Yohai method with weights)
  glmrob_try = try(glmrob(y ~ X, method="WBY", family=binomial, weights.on.x = "none")$coefficients)
  if(!("try-error" %in% class(glmrob_try))){
    glmrob4_nocont[i,] = glmrob_try
  }
  remove(glmrob_try)
  
  # WMLE + IB (with weight on X, default c)
  # ib_try = try(ib(y = y, X = X, crob = crob, gamma = gamma, start = betas_true))
  # if(!("try-error" %in% class(ib_try))){
  #   ib_nocont[i,] = ib_try
  # }
  # remove(ib_try)
  
  # WMLE + IB (with weight on X, chosen c)
  # ib_try = try(ib(y = y, X = X, crob = c_ib, gamma = gamma, start = betas_true))
  # if(!("try-error" %in% class(ib_try))){
  #   ib2_nocont[i,] = ib_try
  # }
  # remove(ib_try)
  
  # ----- cont X
  X = X_list$cont
  weightX = get_weightX_cpp(X)
  y = y_list$uncont
  
  # MLE
  mle_contX[i,] = glm(y ~ X, family=binomial(link="logit"))$coefficients
  
  # brglm
  brglm_contX[i,] = brglm(y ~ X, family=binomial)$coefficients
  
  # glmrob (with weight on X, default c)
  glmrob_try = try(glmrob(y ~ X, family=binomial(link="logit"), weights.on.x = "hat")$coefficients)
  if(!("try-error" %in% class(glmrob_try))){
    glmrob_contX[i,] = glmrob_try
  }
  remove(glmrob_try)
  
  # glmrob (with weight on X, chosen c)
  glmrob_try = try(glmrob(y ~ X, family=binomial(link="logit"), weights.on.x = "hat",
                          control = glmrobMqle.control(tcc = c_glmrob))$coefficients)
  if(!("try-error" %in% class(glmrob_try))){
    glmrob2_contX[i,] = glmrob_try
  }
  remove(glmrob_try)
  
  # glmrob (Bianco-Yohai method)
  glmrob_try = try(glmrob(y ~ X, method="BY", family=binomial, weights.on.x = "none")$coefficients)
  if(!("try-error" %in% class(glmrob_try))){
    glmrob3_contX[i,] = glmrob_try
  }
  remove(glmrob_try)
  
  # glmrob (Bianco-Yohai method with weights)
  glmrob_try = try(glmrob(y ~ X, method="WBY", family=binomial, weights.on.x = "none")$coefficients)
  if(!("try-error" %in% class(glmrob_try))){
    glmrob4_contX[i,] = glmrob_try
  }
  remove(glmrob_try)
  
  # WMLE + IB (with weights on X, default c)
  # ib_try = try(ib(y = y, X = X, crob = crob, gamma = gamma, start = betas_true))
  # if(!("try-error" %in% class(ib_try))){
  #   ib_contX[i,] = ib_try
  # }
  # remove(ib_try)
  
  # WMLE + IB (with weights on X, chosen c)
  # ib_try = try(ib(y = y, X = X, crob = c_ib, gamma = gamma, start = betas_true))
  # if(!("try-error" %in% class(ib_try))){
  #   ib2_contX[i,] = ib_try
  # }
  # remove(ib_try)
  
  
  # ----- cont Y
  X = X_list$uncont
  weightX = get_weightX_cpp(X)
  y = y_list$cont
  
  # MLE
  mle_contY[i,] = glm(y ~ X, family=binomial(link="logit"))$coefficients
  
  # brglm
  brglm_contY[i,] = brglm(y ~ X, family=binomial)$coefficients
  
  # glmrob (with weight on X, default c)
  glmrob_try = try(glmrob(y ~ X, family=binomial(link="logit"), weights.on.x = "hat")$coefficients)
  if(!("try-error" %in% class(glmrob_try))){
    glmrob_contY[i,] = glmrob_try
  }
  remove(glmrob_try)
  
  # glmrob (with weight on X, chosen c)
  glmrob_try = try(glmrob(y ~ X, family=binomial(link="logit"), weights.on.x = "hat",
                          control = glmrobMqle.control(tcc = c_glmrob))$coefficients)
  if(!("try-error" %in% class(glmrob_try))){
    glmrob2_contY[i,] = glmrob_try
  }
  remove(glmrob_try)
  
  # glmrob (Bianco-Yohai method)
  glmrob_try = try(glmrob(y ~ X, method="BY", family=binomial, weights.on.x = "none")$coefficients)
  if(!("try-error" %in% class(glmrob_try))){
    glmrob3_contY[i,] = glmrob_try
  }
  remove(glmrob_try)
  
  # glmrob (Bianco-Yohai method with weights)
  glmrob_try = try(glmrob(y ~ X, method="WBY", family=binomial, weights.on.x = "none")$coefficients)
  if(!("try-error" %in% class(glmrob_try))){
    glmrob4_contY[i,] = glmrob_try
  }
  remove(glmrob_try)
  
  # WMLE + IB (with weights on X, default c)
  # ib_try = try(ib(y = y, X = X, crob = crob, gamma = gamma, start = betas_true))
  # if(!("try-error" %in% class(ib_try))){
  #   ib_contY[i,] = ib_try
  # }
  # remove(ib_try)
  
  # WMLE + IB (with weights on X, chosen c)
  # ib_try = try(ib(y = y, X = X, crob = c_ib, gamma = gamma, start = betas_true))
  # if(!("try-error" %in% class(ib_try))){
  #   ib2_contY[i,] = ib_try
  # }
  # remove(ib_try)
  
  
  # ----- cont XY
  X = X_list$cont
  weightX = get_weightX_cpp(X)
  y = y_list$cont
  
  # MLE
  mle_contXY[i,] = glm(y ~ X, family=binomial(link="logit"))$coefficients
  
  # brglm
  brglm_contXY[i,] = brglm(y ~ X, family=binomial)$coefficients
  
  # glmrob (with weight on X, default c)
  glmrob_try = try(glmrob(y ~ X, family=binomial(link="logit"), weights.on.x = "hat")$coefficients)
  if(!("try-error" %in% class(glmrob_try))){
    glmrob_contXY[i,] = glmrob_try
  }
  remove(glmrob_try)
  
  # glmrob (with weight on X, chosen c)
  glmrob_try = try(glmrob(y ~ X, family=binomial(link="logit"), weights.on.x = "hat",
                          control = glmrobMqle.control(tcc = c_glmrob))$coefficients)
  if(!("try-error" %in% class(glmrob_try))){
    glmrob2_contXY[i,] = glmrob_try
  }
  remove(glmrob_try)
  
  # glmrob (Bianco-Yohai method)
  glmrob_try = try(glmrob(y ~ X, method="BY", family=binomial, weights.on.x = "none")$coefficients)
  if(!("try-error" %in% class(glmrob_try))){
    glmrob3_contXY[i,] = glmrob_try
  }
  remove(glmrob_try)
  
  # glmrob (Bianco-Yohai method with weights)
  glmrob_try = try(glmrob(y ~ X, method="WBY", family=binomial, weights.on.x = "none")$coefficients)
  if(!("try-error" %in% class(glmrob_try))){
    glmrob4_contXY[i,] = glmrob_try
  }
  remove(glmrob_try)
  
  # WMLE + IB (with weights on X, default c)
  # ib_try = try(ib(y = y, X = X, crob = crob, gamma = gamma, start = betas_true))
  # if(!("try-error" %in% class(ib_try))){
  #   ib_contXY[i,] = ib_try
  # }
  # remove(ib_try)
  
  # WMLE + IB (with weights on X, chosen c)
  # ib_try = try(ib(y = y, X = X, crob = c_ib, gamma = gamma, start = betas_true))
  # if(!("try-error" %in% class(ib_try))){
  #   ib2_contXY[i,] = ib_try
  # }
  # remove(ib_try)
  
  print(i)
  
  
  # ----- boxplot
  if (i %% 10 == 0){
    mycol = ggplot_like_colors(7)
    par(mfrow = c(1, 6), oma = c(5,4,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)
    for (i in 1:6){
      boxplot(mle_nocont[,i],
              mle_contX[,i],
              mle_contY[,i],
              mle_contXY[,i],
              brglm_nocont[,i],
              brglm_contX[,i],
              brglm_contY[,i],
              brglm_contXY[,i],
              glmrob_nocont[,i],
              glmrob_contX[,i],
              glmrob_contY[,i],
              glmrob_contXY[,i],
              glmrob2_nocont[,i],
              glmrob2_contX[,i],
              glmrob2_contY[,i],
              glmrob2_contXY[,i],
              glmrob3_nocont[,i],
              glmrob3_contX[,i],
              glmrob3_contY[,i],
              glmrob3_contXY[,i],
              glmrob4_nocont[,i],
              glmrob4_contX[,i],
              glmrob4_contY[,i],
              glmrob4_contXY[,i],
              # ib_nocont[,i],
              # ib_contX[,i],
              # ib_contY[,i],
              # ib_contXY[,i],
              # ib2_nocont[,i],
              # ib2_contX[,i],
              # ib2_contY[,i],
              # ib2_contXY[,i],
              outline = F, col = rep(mycol, each = 4))
      abline(h = betas_true[i], col = 2, lwd = 2)
    }
  }
}









