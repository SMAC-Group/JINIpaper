`.sourceCpp_1_DLLInfo` <- dyn.load('/home/orsos/Github/SMAC-Group/JINIpaper/.simu/compiled/sourceCpp_2.so')

expit_ele_cpp <- Rcpp:::sourceCppFunction(function(x) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_expit_ele_cpp')
expit_vec_cpp <- Rcpp:::sourceCppFunction(function(x) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_expit_vec_cpp')
get_W_cpp <- Rcpp:::sourceCppFunction(function(y, X, betas, crob) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_get_W_cpp')
get_gradient_cpp <- Rcpp:::sourceCppFunction(function(y, X, betas, crob) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_get_gradient_cpp')
get_D_cpp <- Rcpp:::sourceCppFunction(function(X, betas) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_get_D_cpp')
get_hessian_cpp <- Rcpp:::sourceCppFunction(function(y, X, betas, crob) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_get_hessian_cpp')
make_rob_step_cpp <- Rcpp:::sourceCppFunction(function(y, X, betas, crob) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_make_rob_step_cpp')
get_weightX_cpp <- Rcpp:::sourceCppFunction(function(X) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_get_weightX_cpp')
get_gradient2_cpp <- Rcpp:::sourceCppFunction(function(y, X, betas, crob, weightX) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_get_gradient2_cpp')
get_hessian2_cpp <- Rcpp:::sourceCppFunction(function(y, X, betas, crob, weightX) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_get_hessian2_cpp')
make_rob_step2_cpp <- Rcpp:::sourceCppFunction(function(y, X, betas, crob, weightX) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_make_rob_step2_cpp')

rm(`.sourceCpp_1_DLLInfo`)
