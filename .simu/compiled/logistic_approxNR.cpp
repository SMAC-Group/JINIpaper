#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// ------------------
// --- Helpers
// ------------------


// [[Rcpp::export]]
double expit_ele_cpp(double x) {
  double result = exp(x) / (exp(x) + 1);
  return(result);
}


// [[Rcpp::export]]
vec expit_vec_cpp(vec x) {
  int n=x.n_elem;
  vec result(n);
  for ( int i = 0; i < n; i++ ){
    result[i] = expit_ele_cpp(x[i]);
  }
  return(result);
}


// ------------------
// --- NR algo
// ------------------

// [[Rcpp::export]]
mat get_W_cpp(vec y, mat X, vec betas, double crob){
  int n = y.n_elem;
  
  mat W = zeros(n,n) ;
  vec mu = expit_vec_cpp(X * betas);
  for ( int i = 0; i < n; i++ ){
    vec s = trans(X.row(i)) * (y[i]-mu[i]);
    double s_norm = norm(s,2);
    if(s_norm < crob){
      double a = pow(s_norm/crob,2);
      W(i,i) = pow(1-a, 2);
    }else{
      W(i,i) = 0;
    }
  }
  
  return(W);
}


// [[Rcpp::export]]
vec get_gradient_cpp(vec y, mat X, vec betas, double crob){
  mat W = get_W_cpp(y, X, betas, crob);
  vec mu = expit_vec_cpp(X * betas);
  vec grad = trans(X) * W * (y-mu);
  return(grad);
}


// [[Rcpp::export]]
mat get_D_cpp(mat X, vec betas){
  int n = X.n_rows;
  mat D = zeros(n,n);
  for ( int i = 0; i < n; i++ ){
    vec tmp = X.row(i) * betas;
    double z = tmp[0];
    D(i,i) = exp(z) / pow(1+exp(z), 2);
  }
  return(D);
}


// [[Rcpp::export]]
mat get_hessian_cpp(vec y, mat X, vec betas, double crob){
  mat W = get_W_cpp(y, X, betas, crob);
  mat D = get_D_cpp(X, betas);
  mat result = - trans(X) * W * D * X;
  return(result);
}
  

// [[Rcpp::export]]  
vec make_rob_step_cpp(vec y, mat X, vec betas, double crob){
  vec G = get_gradient_cpp(y, X, betas, crob);
  mat H = get_hessian_cpp(y, X, betas, crob);
  vec betas_new = betas + inv(-H) * G;
  return(betas_new);
}


// ------------------
// --- NR algo 2
// ------------------

// [[Rcpp::export]]
mat get_weightX_cpp(mat X){
  int n = X.n_rows;
  mat H = X * inv(trans(X)*X) * trans(X);
  mat weightX = zeros(n,n);
  for ( int i = 0; i < n; i++ ){
    weightX(i,i) = sqrt(1-H(i,i));
  }
  return(weightX);
}


// [[Rcpp::export]]
vec get_gradient2_cpp(vec y, mat X, vec betas, double crob, mat weightX){
  mat W = get_W_cpp(y, X, betas, crob);
  vec mu = expit_vec_cpp(X * betas);
  vec grad = trans(X) * weightX * W * (y-mu);
  return(grad);
}


// [[Rcpp::export]]
mat get_hessian2_cpp(vec y, mat X, vec betas, double crob, mat weightX){
  mat W = get_W_cpp(y, X, betas, crob);
  mat D = get_D_cpp(X, betas);
  mat result = - trans(X) * weightX * W * D * X;
  return(result);
}


// [[Rcpp::export]]  
vec make_rob_step2_cpp(vec y, mat X, vec betas, double crob, mat weightX){
  vec G = get_gradient2_cpp(y, X, betas, crob, weightX);
  mat H = get_hessian2_cpp(y, X, betas, crob, weightX);
  vec betas_new = betas + inv(-H) * G;
  return(betas_new);
}














#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// expit_ele_cpp
double expit_ele_cpp(double x);
RcppExport SEXP sourceCpp_1_expit_ele_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(expit_ele_cpp(x));
    return rcpp_result_gen;
END_RCPP
}
// expit_vec_cpp
vec expit_vec_cpp(vec x);
RcppExport SEXP sourceCpp_1_expit_vec_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(expit_vec_cpp(x));
    return rcpp_result_gen;
END_RCPP
}
// get_W_cpp
mat get_W_cpp(vec y, mat X, vec betas, double crob);
RcppExport SEXP sourceCpp_1_get_W_cpp(SEXP ySEXP, SEXP XSEXP, SEXP betasSEXP, SEXP crobSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< vec >::type betas(betasSEXP);
    Rcpp::traits::input_parameter< double >::type crob(crobSEXP);
    rcpp_result_gen = Rcpp::wrap(get_W_cpp(y, X, betas, crob));
    return rcpp_result_gen;
END_RCPP
}
// get_gradient_cpp
vec get_gradient_cpp(vec y, mat X, vec betas, double crob);
RcppExport SEXP sourceCpp_1_get_gradient_cpp(SEXP ySEXP, SEXP XSEXP, SEXP betasSEXP, SEXP crobSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< vec >::type betas(betasSEXP);
    Rcpp::traits::input_parameter< double >::type crob(crobSEXP);
    rcpp_result_gen = Rcpp::wrap(get_gradient_cpp(y, X, betas, crob));
    return rcpp_result_gen;
END_RCPP
}
// get_D_cpp
mat get_D_cpp(mat X, vec betas);
RcppExport SEXP sourceCpp_1_get_D_cpp(SEXP XSEXP, SEXP betasSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< vec >::type betas(betasSEXP);
    rcpp_result_gen = Rcpp::wrap(get_D_cpp(X, betas));
    return rcpp_result_gen;
END_RCPP
}
// get_hessian_cpp
mat get_hessian_cpp(vec y, mat X, vec betas, double crob);
RcppExport SEXP sourceCpp_1_get_hessian_cpp(SEXP ySEXP, SEXP XSEXP, SEXP betasSEXP, SEXP crobSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< vec >::type betas(betasSEXP);
    Rcpp::traits::input_parameter< double >::type crob(crobSEXP);
    rcpp_result_gen = Rcpp::wrap(get_hessian_cpp(y, X, betas, crob));
    return rcpp_result_gen;
END_RCPP
}
// make_rob_step_cpp
vec make_rob_step_cpp(vec y, mat X, vec betas, double crob);
RcppExport SEXP sourceCpp_1_make_rob_step_cpp(SEXP ySEXP, SEXP XSEXP, SEXP betasSEXP, SEXP crobSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< vec >::type betas(betasSEXP);
    Rcpp::traits::input_parameter< double >::type crob(crobSEXP);
    rcpp_result_gen = Rcpp::wrap(make_rob_step_cpp(y, X, betas, crob));
    return rcpp_result_gen;
END_RCPP
}
// get_weightX_cpp
mat get_weightX_cpp(mat X);
RcppExport SEXP sourceCpp_1_get_weightX_cpp(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(get_weightX_cpp(X));
    return rcpp_result_gen;
END_RCPP
}
// get_gradient2_cpp
vec get_gradient2_cpp(vec y, mat X, vec betas, double crob, mat weightX);
RcppExport SEXP sourceCpp_1_get_gradient2_cpp(SEXP ySEXP, SEXP XSEXP, SEXP betasSEXP, SEXP crobSEXP, SEXP weightXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< vec >::type betas(betasSEXP);
    Rcpp::traits::input_parameter< double >::type crob(crobSEXP);
    Rcpp::traits::input_parameter< mat >::type weightX(weightXSEXP);
    rcpp_result_gen = Rcpp::wrap(get_gradient2_cpp(y, X, betas, crob, weightX));
    return rcpp_result_gen;
END_RCPP
}
// get_hessian2_cpp
mat get_hessian2_cpp(vec y, mat X, vec betas, double crob, mat weightX);
RcppExport SEXP sourceCpp_1_get_hessian2_cpp(SEXP ySEXP, SEXP XSEXP, SEXP betasSEXP, SEXP crobSEXP, SEXP weightXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< vec >::type betas(betasSEXP);
    Rcpp::traits::input_parameter< double >::type crob(crobSEXP);
    Rcpp::traits::input_parameter< mat >::type weightX(weightXSEXP);
    rcpp_result_gen = Rcpp::wrap(get_hessian2_cpp(y, X, betas, crob, weightX));
    return rcpp_result_gen;
END_RCPP
}
// make_rob_step2_cpp
vec make_rob_step2_cpp(vec y, mat X, vec betas, double crob, mat weightX);
RcppExport SEXP sourceCpp_1_make_rob_step2_cpp(SEXP ySEXP, SEXP XSEXP, SEXP betasSEXP, SEXP crobSEXP, SEXP weightXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< vec >::type betas(betasSEXP);
    Rcpp::traits::input_parameter< double >::type crob(crobSEXP);
    Rcpp::traits::input_parameter< mat >::type weightX(weightXSEXP);
    rcpp_result_gen = Rcpp::wrap(make_rob_step2_cpp(y, X, betas, crob, weightX));
    return rcpp_result_gen;
END_RCPP
}
