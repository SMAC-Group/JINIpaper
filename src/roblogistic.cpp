// Robust logistic regression with Tukey's weight (inconsistent)
#include "common.h"

// ------------------
// Implementation
// ------------------
// --------------
// Functions for logistic
// --------------
double g(double x){return std::log(x / (0.1e1 - x));}
// double g1(double x){return 0.1e1 / (x - x * x);}
// double g2(double x){return (0.2e1 * x - 0.1e1) / (x*x - x*x*x*x);}
double V(double x){return x - x * x;}
double V1(double x){return 0.1e1 - 0.2e1 * x;}
// Tukey's weight:
double psi(double x, double c){return (std::abs(x) <= c) ? x * (x*x/c/c-0.1e1) * (x*x/c/c-0.1e1) : 0.0;}
double wc(double x, double c){return (std::abs(x) <= c) ? 0.1e1 - 0.2e1 * x*x/c/c + x*x*x*x/c/c/c/c : 0.0;}
double wc1(double x, double c){return (std::abs(x) <= c) ? - 0.4e1 * x/c/c + 0.4e1 * x*x*x/c/c/c/c : 0.0;}
double psi1(double x, double c){return (std::abs(x) <= c) ? 0.1e1 - 6*x*x/c/c + 5*x*x*x*x/c/c/c/c : 0.0;}
// Hubers' weight:
// double psi(double x, double c){return (std::abs(x) <= c) ? x : c;}
// double psi1(double x, double c){return (std::abs(x) <= c) ? 0.1e1 : 0.0;}

// Expectation of psi(r)
double Epsi(double mu, double c){
  // data storage
  double varmu, mu2, mu4, mu_2, mu_4, c2, c4, v12, v32, t1(0.0), t2(0.0);
  bool is_below, is_above, is_zero;
  
  // computation
  varmu = V(mu);
  v12 = std::sqrt(varmu);
  is_below = mu <= c * v12;
  is_above = mu > 0.1e1 - c * v12;
  is_zero = !is_below && !is_above;
  if(is_zero) {
    return 0.0;
  }
  c2 = c * c;
  c4 = c2 * c2;
  v32 = v12 * varmu; 
  if(is_below) {
    mu2 = mu * mu;
    mu4 = mu2 * mu2;
    t1 = -v12 + 0.2e1 * mu2 / v12 / c2 - mu4 / c4 / v32;
  }
  if(is_above) {
    mu_2 = 0.1e1 - 0.2e1 * mu + mu * mu;
    mu_4 = mu_2 * mu_2;
    t2 = v12 - 0.2e1 * mu_2 / v12 / c2 + mu_4 / c4 / v32;
  }
  return t1 + t2;
}

// First derivative of expectation of psi(r)
double dEpsi(double mu, double c){
  // data storage
  double varmu, varmup, mu2, mu3, mu4, mu_, mu_2, mu_3, mu_4, c2, c4, v12, v2, t1(0.0), t2(0.0);
  bool is_below, is_above, is_zero;
  
  // computation
  varmu = V(mu);
  v12 = std::sqrt(varmu);
    is_below = mu <= c * v12;
  is_above = mu > 0.1e1 - c * v12;
  is_zero = !is_below && !is_above;
  if(is_zero) {
    return 0.0;
  }
  c2 = c * c;
  c4 = c2 * c2;
  v2 = varmu * varmu; 
  varmup = V1(mu);
  if(is_below) {
    mu2 = mu * mu;
    mu3 = mu2 * mu;
    mu4 = mu2 * mu2;
    t1 = -varmup / 0.2e1 + 0.4e1 * mu / c2 - mu2 * varmup / varmu / c2 - 0.4e1 * mu3 / c4 / varmu + 0.3e1 * mu4 * varmup / 0.2e1 / c4 / v2;
  }
  if(is_above) {
    mu_ = 0.1e1 - mu;
    mu_2 = mu_ * mu_;
    mu_3 = mu_2 * mu_;
    mu_4 = mu_2 * mu_2;
    t2 = varmup / 0.2e1 + 0.4e1 * mu_ / c2 + mu_2 * varmup / varmu / c2 - 0.4e1 * mu_3 / c4 / varmu - 0.3e1 * mu_4 * varmup / 0.2e1 / c4 / v2;
  }
  return (t1 + t2) / v12;
}

// Expectation of first derivative of psi(r)
double Edpsi(double mu, double c){
  // data storage
  double varmu, mu3, mu_, mu_3, c2, c4, v12, t1(0.0), t2(0.0);
  bool is_below, is_above, is_zero;
  
  // computation
  varmu = V(mu);
  v12 = std::sqrt(varmu);
  is_below = mu <= c * v12;
  is_above = mu > 0.1e1 - c * v12;
  is_zero = !is_below && !is_above;
  if(is_zero) {
    return 0.0;
  }
  c2 = c * c;
  c4 = c2 * c2;
  mu_ = 0.1e1 - mu;
  if(is_below) {
    mu3 = mu * mu * mu;
    t1 = mu_ - 0.6e1 * mu / c2 + 0.5e1 * mu3 / c4 / varmu;
  }
  if(is_above) {
    mu_3 = mu_ * mu_ * mu_;
    t2 = mu - 0.6e1 * mu_ / c2 + 0.5e1 * mu_3 / c4 / varmu;
  }
  return t1 + t2;
}

// Expectation of first derivative of psi(r) times r
double Edpsir(double mu, double c){
  // data storage
  double varmu, mu2, mu4, mu_, mu_2, mu_4, c2, c4, v12, v32, t1(0.0), t2(0.0);
  bool is_below, is_above, is_zero;
  
  // computation
  varmu = V(mu);
  v12 = std::sqrt(varmu);
  is_below = mu <= c * v12;
  is_above = mu > 0.1e1 - c * v12;
  is_zero = !is_below && !is_above;
  if(is_zero) {
    return 0.0;
  }
  c2 = c * c;
  c4 = c2 * c2;
  v32 = varmu * v12; 
  if(is_below) {
    mu2 = mu * mu;
    mu4 = mu2 * mu2;
    t1 = -v12 + 0.6e1 * mu2 / c2 / v12 - 0.5e1 * mu4 / c4 / v32;
  }
  if(is_above) {
    mu_ = 0.1e1 - mu;
    mu_2 = mu_ * mu_;
    mu_4 = mu_2 * mu_2;
    t2 = v12 - 0.6e1 * mu_2 / c2 / v12 + 0.5e1 * mu_4 / c4 / v32;
  }
  return t1 + t2;
}

// Expectation of psi(r) squared
double Epsi2(double mu, double c){
  // data storage
  double varmu, mu2, mu3, mu5, mu7, mu9, mu_, mu_2, mu_3, mu_5, mu_7, mu_9, c2, c4, c6, c8, v12, v2, v3, v4, t1(0.0), t2(0.0);
  bool is_below, is_above, is_zero;
  
  // computation
  varmu = V(mu);
  v12 = std::sqrt(varmu);
  is_below = mu <= c * v12;
  is_above = mu > 0.1e1 - c * v12;
  is_zero = !is_below && !is_above;
  if(is_zero) {
    return 0.0;
  }
  c2 = c * c;
  c4 = c2 * c2;
  c6 = c4 * c2;
  c8 = c4 * c4;
  v2 = varmu * varmu; 
  v3 = v2 * varmu;
  v4 = v2 * v2;
  if(is_below) {
    mu2 = mu * mu;
    mu3 = mu2 * mu;
    mu5 = mu3 * mu2;
    mu7 = mu5 * mu2;
    mu9 = mu7 * mu2;
    t1 = mu - 0.4e1 * mu3 / varmu / c2 + 0.6e1 * mu5 / v2 / c4 - 0.4e1 * mu7 / v3 / c6 + mu9 / c8 / v4;
  }
  if(is_above) {
    mu_ = 0.1e1 - mu;
    mu_2 = mu_ * mu_;
    mu_3 = mu_2 * mu_;
    mu_5 = mu_3 * mu_2;
    mu_7 = mu_5 * mu_2;
    mu_9 = mu_7 * mu_2;
    t2 = mu_ - 0.4e1 * mu_3 / varmu / c2 + 0.6e1 * mu_5 / v2 / c4 - 0.4e1 * mu_7 / v3 / c6 + mu_9 / c8 / v4;
  }
  return t1 + t2;
}

// Expectation of Tukey's weight
// Note : xp stands for p-norm of vector x
double Ew(double mu, double c, double xp){
  // data storage
  double varmu, mu3, mu_, mu_3, xpc, xpc2, xpc4, t1(0.0), t2(0.0);
  bool is_below, is_above, is_zero;
  
  // computation
  xpc = xp / c; 
  is_below = mu <= 0.1e1 / xpc;
  is_above = mu > 0.1e1 - 0.1e1 / xpc;
  is_zero = !is_below && !is_above;
  if(is_zero) {
    return 0.0;
  }
  varmu = V(mu);
  xpc2 = xpc * xpc;
  xpc4 = xpc2 * xpc2;
  mu_ = 0.1e1 - mu;
  if(is_below) {
    mu3 = mu * mu * mu;
    t1 = mu_ - 0.2e1 * varmu * xpc2 * mu + varmu * mu3 * xpc4;
  }
  if(is_above) {
    mu_3 = mu_ * mu_ * mu_;
    t2 = mu - 0.2e1 * varmu * mu_ * xpc2 + varmu * mu_3 * xpc4;
  }
  return t1 + t2;
}

// Expectation of first derivative of Tukey's weight times -1^y(y-mu)
// Note : xp stands for p-norm of vector x
double Edwymu(double mu, double c, double xp){
  // data storage
  double varmu, mu3, mu_, mu_3, xpc, xpc3, k1, t1(0.0), t2(0.0);
  bool is_below, is_above, is_zero;
  
  // computation
  xpc = xp / c; 
  is_below = mu <= 0.1e1 / xpc;
  is_above = mu > 0.1e1 - 0.1e1 / xpc;
  is_zero = !is_below && !is_above;
  if(is_zero) {
    return 0.0;
  }
  varmu = V(mu);
  xpc3 = xpc * xpc * xpc;
  mu_ = 0.1e1 - mu;
  k1 = 0.4e1 * varmu / c;
  if(is_below) {
    mu3 = mu * mu * mu;
    t1 = mu * xpc - mu3 * xpc3;
  }
  if(is_above) {
    mu_3 = mu_ * mu_ * mu_;
    t2 = mu_ * xpc - mu_3 * xpc3;
  }
  return k1 * t1 + k1 * t2;
}

// Expectation of Tukey's weight times (y-mu)
// Note : xp stands for p-norm of vector x
double Ewymu(double mu, double c, double xp){
  // data storage
  double varmu, mu2, mu4, mu_, mu_2, mu_4, xpc, xpc2, xpc4, t1(0.0), t2(0.0);
  bool is_below, is_above, is_zero;
  
  // computation
  xpc = xp / c; 
  is_below = mu <= 0.1e1 / xpc;
  is_above = mu > 0.1e1 - 0.1e1 / xpc;
  is_zero = !is_below && !is_above;
  if(is_zero) {
    return 0.0;
  }
  varmu = V(mu);
  xpc2 = xpc * xpc;
  xpc4 = xpc2 * xpc2;
  mu_ = 0.1e1 - mu;
  if(is_below) {
    mu2 = mu * mu;
    mu4 = mu2 * mu2;
    t1 = -0.1e1 + 0.2e1 * mu2 * xpc2 - mu4 * xpc4;
  }
  if(is_above) {
    mu_2 = mu_ * mu_;
    mu_4 = mu_2 * mu_2;
    t2 = 0.1e1 - 0.2e1 * mu_2 * xpc2 + mu_4 * xpc4;
  }
  return varmu * t1 + varmu * t2;
}

// First derivative of expectation of Tukey's weight times (y-mu)
// Note : xp stands for p-norm of vector x
double dEwymu(double mu, double c, double xp){
  // data storage
  double varmu, varmup, mu2, mu3, mu4, mu_, mu_2, mu_3, mu_4, xpc, xpc2, xpc4, t1(0.0), t2(0.0);
  bool is_below, is_above, is_zero;
  
  // computation
  xpc = xp / c; 
  is_below = mu <= 0.1e1 / xpc;
  is_above = mu > 0.1e1 - 0.1e1 / xpc;
  is_zero = !is_below && !is_above;
  if(is_zero) {
    return 0.0;
  }
  varmu = V(mu);
  varmup = V1(mu);
  xpc2 = xpc * xpc;
  xpc4 = xpc2 * xpc2;
  mu_ = 0.1e1 - mu;
  if(is_below) {
    mu2 = mu * mu;
    mu3 = mu2 * mu;
    mu4 = mu2 * mu2;
    t1 = -varmup + 0.2e1 * xpc2 * (0.2e1 * mu * varmu + mu2 * varmup) - xpc4 * (0.4e1 * mu3 * varmu + mu4 * varmup);
  }
  if(is_above) {
    mu_2 = mu_ * mu_;
    mu_3 = mu_2 * mu_;
    mu_4 = mu_2 * mu_2;
    t2 = varmup - 0.2e1 * xpc2 * (mu_2 * varmup - 0.2e1 * mu_ * varmu) + xpc4 * (mu_4 * varmup - 0.2e1 * mu_* varmu);
  }
  return t1 + t2;
}

// Expectation of Tukey's weight times (y-mu) all squared
// Note : xp stands for p-norm of vector x
double Ewymu2(double mu, double c, double xp){
  // data storage
  double mu2, mu4, mu6, mu8, mu_, mu_2, mu_4, mu_6, mu_8, xpc, xpc2, xpc4, xpc6, xpc8, t1(0.0), t2(0.0);
  bool is_below, is_above, is_zero;
  
  // computation
  xpc = xp / c; 
  is_below = mu <= 0.1e1 / xpc;
  is_above = mu > 0.1e1 - 0.1e1 / xpc;
  is_zero = !is_below && !is_above;
  if(is_zero) {
    return 0.0;
  }
  xpc2 = xpc * xpc;
  xpc4 = xpc2 * xpc2;
  xpc6 = xpc4 * xpc2;
  xpc8 = xpc6 * xpc2;
  mu_ = 0.1e1 - mu;
  if(is_below) {
    mu2 = mu * mu;
    mu4 = mu2 * mu2;
    mu6 = mu4 * mu2;
    mu8 = mu6 * mu2;
    t1 = 0.1e1 - 0.4e1 * mu2 * xpc2 + 0.6e1 * mu4 * xpc4 - 0.4e1 * mu6 * xpc6 + mu8 * xpc8;
  }
  if(is_above) {
    mu_2 = mu_ * mu_;
    mu_4 = mu_2 * mu_2;
    mu_6 = mu_4 * mu_2;
    mu_8 = mu_6 * mu_2;
    t2 = 0.1e1 - 0.4e1 * mu_2 * xpc2 + 0.6e1 * mu_4 * xpc4 - 0.4e1 * mu_6 * xpc6 + mu_8 * xpc8;
  }
  return mu2 * mu_ * t1 + mu_2 * mu * t2;
}



// //' Robust logistic regression initial (inconsistent) estimator with Tukey's weights
// //'
// //' @param y a vector of responses
// //' @param x a n x p matrix of design
// //' @param beta a p-vector of parameter (starting values)
// //' @param c tuning parameter for Tukey's weights (default value is 4.685061)
// //' @param maxit max number of iteration for IRWLS
// //' @param tol tolerance for stopping criterion
// //' @param verbose print info
// //' @export
// // [[Rcpp::export]]
// Rcpp::List roblogisticMqle1(
//      Eigen::ArrayXd& y,
//      Eigen::MatrixXd& x,
//      Eigen::VectorXd& start,
//      double c = 4.685061,
//      unsigned int maxit=200,
//      double tol=1e-7,
//      bool verbose=false
//  ){
//    unsigned int n = y.size();
//    unsigned int p = x.cols();
//    unsigned int p1 = p+1;
// 
//    // data storage
//    Eigen::VectorXd mu(n),eta(n),z(n);
//    Eigen::MatrixXd x1(n,p+1);
//    Eigen::MatrixXd w = Eigen::MatrixXd::Zero(n,n);
//    double t1, t2, alpha, varmu, varmup, r, c2, c4;
//    // double psi_r, gg, ;
//    Eigen::MatrixXd A(p1,p1);
//    Eigen::VectorXd b(p1),beta(p1);
// 
//    // pre-computation
//    x1.rightCols(p) = x;
//    x1.col(0) = Eigen::VectorXd::Constant(n,1.0);
// 
//    // Iterative Re-Weighted Least Squares
//    unsigned int it = 0;
//    while(it < maxit){
//      // Compute w, z
//      eta = x1 * start;
//      mu = (eta).unaryExpr(Sigmoid());
//      for(unsigned int j=0;j<n;++j){
//        // Computation of alpha (Newton-Raphson IRLS)
//        // varmu = V(mu(j));
//        // t1 = std::sqrt(varmu);
//        // r = y(j) / t1 - mu(j) / t1; // Pearson's residual
//        // t2 = V1(mu(j)) / t1 / 0.2e1;
//        // psi_r = psi(r,c);
//        // gg = g1(mu(j));
//        // alpha = psi_r * t2 + psi_r * t1 * g2(mu(j)) / gg + psi1(r,c) * (0.1e1 + r * t2);
//        // z(j) = eta(j) + psi_r * gg * t1 / alpha;
//        // w(j,j) = alpha / varmu / gg / gg;
// 
//        // Computation of E[alpha] (expectation) for Fisher scoring IRLS
//        varmu = V(mu(j));
//        t1 = std::sqrt(varmu);
//        r = y(j) / t1 - mu(j) / t1; // Pearson's residual
//        if(std::abs(r) > c) {
//          z(j) = eta(j);
//          w(j,j) = 0;
//        } else {
//          c2 = c * c;
//          c4 = c2 * c2;
//          varmup = V1(mu(j));
//          t2 = varmup * varmup;
//          alpha = 0.1e1 - 1.5e1 / c4 + (0.5e1 / c4 + t2 / c2 + t2 / c4) / varmu - 0.6e1 * varmu / c2  - t2 / varmu / varmu / 0.2e1 / c4;
//          z(j) = eta(j) + psi(r,c) / alpha / t1;
//          w(j,j) = alpha * varmu;
//        }
//      }
// 
//      // IRWLS
//      A = x1.transpose() * w * x1;
//      b = x1.transpose() * w * z;
//      // if(method == "s") {
//      beta = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
//      // }
//      // if(method == "q") {
//      //   beta = A.fullPivHouseholderQr().solve(b);
//      // }
// 
//      // Check convergence
//      Eigen::VectorXd Dbeta = start - beta;
//      double relE = std::sqrt(Dbeta.squaredNorm() / std::max(1e-20, start.squaredNorm()));
//      ++it;
// 
//      if(verbose) {
//        Rcpp::Rcout << "Relative error " << relE << " at iteration " << it << std::endl;
//      }
// 
//      if(relE <= tol){
//        start = beta;
//        break;
//      }
//      else {
//        start = beta;
//      }
//    }
// 
//    return Rcpp::List::create(
//      Rcpp::Named("weights") = w.diagonal(),
//      Rcpp::Named("z") = z,
//      Rcpp::Named("coefficients") = beta,
//      Rcpp::Named("maxit") = it
//    );
//  }

// //' Robust logistic regression initial (inconsistent) estimator with Tukey's weights
// //'
// //' @param y a vector of responses
// //' @param x a n x p matrix of design
// //' @param beta a p-vector of parameter (starting values)
// //' @param c tuning parameter for Tukey's weight (default value is 4.685061)
// //' @param maxit max number of iteration for IRWLS
// //' @param tol tolerance for stopping criterion
// //' @param verbose print info
// //' @export
// // [[Rcpp::export]]
// Rcpp::List roblogisticMqle1(
//      Eigen::ArrayXd& y,
//      Eigen::MatrixXd& x,
//      Eigen::VectorXd& start,
//      double c = 4.685061,
//      unsigned int maxit=200,
//      double tol=1e-7,
//      bool verbose=false
//  ){
//    unsigned int n = y.size();
//    unsigned int p = x.cols();
//    unsigned int p1 = p+1;
// 
//    // data storage
//    Eigen::VectorXd mu(n),eta(n),z(n);
//    Eigen::MatrixXd x1(n,p+1);
//    Eigen::MatrixXd w = Eigen::MatrixXd::Zero(n,n);
//    double varmu, varmup, r, v12;
//    Eigen::MatrixXd A(p1,p1);
//    Eigen::VectorXd b(p1),beta(p1);
// 
//    // pre-computation
//    x1.rightCols(p) = x;
//    x1.col(0) = Eigen::VectorXd::Constant(n,1.0);
// 
//    // Iterative Re-Weighted Least Squares
//    unsigned int it = 0;
//    while(it < maxit){
//      // Compute w, z
//      eta = x1 * start;
//      mu = (eta).unaryExpr(Sigmoid());
//      for(unsigned int j=0;j<n;++j){
//        varmu = V(mu(j));
//        v12 = std::sqrt(varmu);
//        varmup = V1(mu(j));
//        r = y(j) / v12 - mu(j) / v12; // Pearson's residual
//        w(j,j) = varmu * Edpsi(mu(j),c) + v12 * varmup * Edpsir(mu(j),c) / 0.2e1 + v12 * varmu* dEpsi(mu(j),c);
//        z(j) = eta(j) + v12 * psi(r,c) / w(j,j); 
//      }
// 
//      // IRWLS
//      A = x1.transpose() * w * x1;
//      b = x1.transpose() * w * z;
//      beta = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
// 
//      // Check convergence
//      Eigen::VectorXd Dbeta = start - beta;
//      double relE = std::sqrt(Dbeta.squaredNorm() / std::max(1e-20, start.squaredNorm()));
//      ++it;
// 
//      if(verbose) {
//        Rcpp::Rcout << "Relative error " << relE << " at iteration " << it << std::endl;
//      }
// 
//      if(relE <= tol){
//        start = beta;
//        break;
//      }
//      else {
//        start = beta;
//      }
//    }
// 
//    return Rcpp::List::create(
//      Rcpp::Named("weights") = w.diagonal(),
//      Rcpp::Named("z") = z,
//      Rcpp::Named("coefficients") = beta,
//      Rcpp::Named("maxit") = it
//    );
//  }

//' Robust logistic regression initial (inconsistent) estimator with Tukey's weights
//'
//' @param y a vector of responses
//' @param x a n x p matrix of design
//' @param beta a p-vector of parameter (starting values)
//' @param c tuning parameter for Tukey's weight (default value is 4.685061)
//' @param maxit max number of iteration for IRWLS
//' @param tol tolerance for stopping criterion
//' @param verbose print info
//' @export
// [[Rcpp::export]]
Rcpp::List roblogisticMqle1(
     Eigen::ArrayXd& y,
     Eigen::MatrixXd& x,
     Eigen::VectorXd& start,
     double c = 4.685061,
     unsigned int maxit=200,
     double tol=1e-7,
     bool verbose=false
 ){
   unsigned int n = y.size();
   unsigned int p = x.cols();
   unsigned int p1 = p+1;
   
   // data storage
   Eigen::VectorXd mu(n),eta(n),z(n);
   Eigen::MatrixXd x1(n,p+1);
   Eigen::MatrixXd w = Eigen::MatrixXd::Zero(n,n);
   double varmu, varmup, r, v12;
   Eigen::MatrixXd A(p1,p1);
   Eigen::VectorXd b(p1),Dbeta(p1),beta(p1);
   
   // pre-computation
   x1.rightCols(p) = x;
   x1.col(0) = Eigen::VectorXd::Constant(n,1.0);
   
   // Fisher scoring
   unsigned int it = 0;
   while(it < maxit){
     // Compute w, z
     eta = x1 * start;
     mu = (eta).unaryExpr(Sigmoid());
     for(unsigned int j=0;j<n;++j){
       varmu = V(mu(j));
       v12 = std::sqrt(varmu);
       varmup = V1(mu(j));
       r = y(j) / v12 - mu(j) / v12; // Pearson's residual
       w(j,j) = varmu * Edpsi(mu(j),c) + v12 * varmup * Edpsir(mu(j),c) / 0.2e1 + v12 * varmu* dEpsi(mu(j),c);
       z(j) = v12 * psi(r,c); 
     }
     
     // Fisher scoring
     A = x1.transpose() * w * x1;
     b = x1.transpose() * z;
     Dbeta = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
     // Dbeta = A.colPivHouseholderQr().solve(b);
     beta = start + Dbeta;
     
     // Check convergence
     double relE = std::sqrt(Dbeta.squaredNorm() / std::max(1e-20, start.squaredNorm()));
     ++it;
     
     if(verbose) {
       Rcpp::Rcout << "Relative error " << relE << " at iteration " << it << std::endl;
     }
     
     if(relE <= tol){
       start = beta;
       break;
     }
     else {
       start = beta;
     }
   }
   
   return Rcpp::List::create(
     Rcpp::Named("weights") = w.diagonal(),
     Rcpp::Named("z") = z,
     Rcpp::Named("coefficients") = beta,
     Rcpp::Named("maxit") = it
   );
 }


//' Robust logistic regression estimator with Tukey's weight
//'
//' @param y a vector of responses
//' @param x a n x p matrix of design
//' @param beta a p-vector of parameter (starting values)
//' @param c tuning parameter for Tukey's weight (default value is 4.685061)
//' @param maxit max number of iteration for IRWLS
//' @param tol tolerance for stopping criterion
//' @param verbose print info
//' @export
// [[Rcpp::export]]
Rcpp::List roblogisticMqle(
     Eigen::ArrayXd& y,
     Eigen::MatrixXd& x,
     Eigen::VectorXd& start,
     double c = 4.685061,
     unsigned int maxit=200,
     double tol=1e-7,
     bool verbose=false
 ){
   unsigned int n = y.size();
   unsigned int p = x.cols();
   unsigned int p1 = p+1;

   // data storage
   Eigen::VectorXd mu(n),eta(n),z(n);
   Eigen::MatrixXd x1(n,p+1);
   Eigen::MatrixXd w = Eigen::MatrixXd::Zero(n,n);
   double varmu, varmup, r, v12;
   Eigen::MatrixXd A(p1,p1);
   Eigen::VectorXd b(p1),beta(p1);

   // pre-computation
   x1.rightCols(p) = x;
   x1.col(0) = Eigen::VectorXd::Constant(n,1.0);

   // Iterative Re-Weighted Least Squares
   unsigned int it = 0;
   while(it < maxit){
     // Compute w, z
     eta = x1 * start;
     mu = (eta).unaryExpr(Sigmoid());
     for(unsigned int j=0;j<n;++j){
       varmu = V(mu(j));
       v12 = std::sqrt(varmu);
       varmup = V1(mu(j));
       r = y(j) / v12 - mu(j) / v12; // Pearson's residual
       w(j,j) = varmu * Edpsi(mu(j),c) + v12 * varmup * Edpsir(mu(j),c) / 0.2e1 + v12 * varmu* dEpsi(mu(j),c);
       z(j) = eta(j) + v12 * psi(r,c) / w(j,j) - v12 * Epsi(mu(j),c) / w(j,j);
     }

     // IRWLS
     A = x1.transpose() * w * x1;
     b = x1.transpose() * w * z;
     beta = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

     // Check convergence
     Eigen::VectorXd Dbeta = start - beta;
     double relE = std::sqrt(Dbeta.squaredNorm() / std::max(1e-20, start.squaredNorm()));
     ++it;

     if(verbose) {
       Rcpp::Rcout << "Relative error " << relE << " at iteration " << it << std::endl;
     }

     if(relE <= tol){
       start = beta;
       break;
     }
     else {
       start = beta;
     }
   }

   return Rcpp::List::create(
     Rcpp::Named("weights") = w.diagonal(),
     Rcpp::Named("z") = z,
     Rcpp::Named("coefficients") = beta,
     Rcpp::Named("maxit") = it
   );
 }

//' Asymptotic variance of robust logistic regression estimator with Tukey's weights
//'
//' @param x a n x p matrix of design
//' @param beta a p-vector of parameter (starting values)
//' @param c tuning parameter for Tukey's weights (default value is 4.685061)
//' @export
// [[Rcpp::export]]
Rcpp::List roblogisticMqleVar(
     Eigen::MatrixXd& x,
     Eigen::VectorXd& start,
     double c = 4.685061
 ){
   unsigned int n = x.rows();
   unsigned int p = x.cols();
   unsigned int p1 = p+1;

   // data storage
   Eigen::VectorXd mu(n),eta(n);
   Eigen::MatrixXd x1(n,p+1);
   Eigen::MatrixXd w = Eigen::MatrixXd::Zero(n,n);
   Eigen::MatrixXd d = Eigen::MatrixXd::Zero(n,n);
   Eigen::MatrixXd d1 = Eigen::MatrixXd::Zero(n,n);
   Eigen::MatrixXd w2 = Eigen::MatrixXd::Zero(n,n);
   double varmu, varmup, v12, a, t1;
   Eigen::MatrixXd asym_var_rob(p1,p1), asym_var_rob1(p1,p1), asym_var_mle(p1,p1);

   // pre-computation
   x1.rightCols(p) = x;
   x1.col(0) = Eigen::VectorXd::Constant(n,1.0);

   // Compute w, z
   eta = x1 * start;
   mu = (eta).unaryExpr(Sigmoid());
   for(unsigned int j=0;j<n;++j){
     // Fisher info (MLE)
     varmu = V(mu(j));
     w2(j,j) = varmu;
     
     // Computation 
     varmup = V1(mu(j));
     v12 = std::sqrt(varmu);
     w(j,j) = varmu * Edpsi(mu(j),c) + v12 * varmup * Edpsir(mu(j),c) / 0.2e1 + v12 * varmu * dEpsi(mu(j),c);
     a = Epsi(mu(j),c);
     t1 = varmu * Epsi2(mu(j),c);
     d(j,j) = t1 - varmu * a * a;
     d1(j,j) = t1;
   }

   Eigen::MatrixXd H = (x1.transpose() * w * x1).inverse();
   asym_var_rob = H * x1.transpose() * d * x1 * H;
   asym_var_rob1 = H * x1.transpose() * d1 * x1 * H;
   asym_var_mle = (x1.transpose() * w2 * x1).inverse();
   double efficiency = asym_var_mle.trace() / asym_var_rob.trace();
   double efficiency1 = asym_var_mle.trace() / asym_var_rob1.trace();

   return Rcpp::List::create(
     Rcpp::Named("V_rob") = asym_var_rob,
     Rcpp::Named("V_rob1") = asym_var_rob1,
     Rcpp::Named("V_mle") = asym_var_mle,
     Rcpp::Named("efficiency") = efficiency,
     Rcpp::Named("efficiency1") = efficiency1
   );
 }

//' Robust logistic regression initial estimator (inconsistent) with Tukey's weights
//'
//' @param y a vector of responses
//' @param x a n x p matrix of design
//' @param beta a p-vector of parameter (starting values)
//' @param c tuning parameter for Tukey's weight (default value is 4.685061)
//' @param maxit max number of iteration for IRWLS
//' @param eps max number of iteration for IRWLS
//' @param tol tolerance for stopping criterion
//' @param verbose print info
//' @export
// [[Rcpp::export]]
Rcpp::List roblogisticWmle1(
     Eigen::ArrayXd& y,
     Eigen::MatrixXd& x,
     Eigen::VectorXd& start,
     double c = 4.685061,
     unsigned int maxit=200,
     double eps = 1e-2,
     double tol=1e-7,
     bool verbose=false
){
   unsigned int n = y.size();
   unsigned int p = x.cols();
   unsigned int p1 = p+1;
   unsigned int conv = 1;
   
   // data storage
   Eigen::VectorXd mu(n),eta(n),z(n),y_tilde(n);
   Eigen::MatrixXd x1(n,p+1);
   Eigen::MatrixXd w = Eigen::MatrixXd::Zero(n,n);
   double varmu, xp, s;
   Eigen::MatrixXd A(p1,p1);
   Eigen::VectorXd b(p1),Dbeta(p1),beta(p1);
   
   // pre-computation
   x1.rightCols(p) = x;
   x1.col(0) = Eigen::VectorXd::Constant(n,1.0);
   for(unsigned int i=0;i<n;++i){
     y_tilde(i) = eps;
     if(y(i) == 1.0){y_tilde(i) = 1.0 - eps;}
   }
   
   // Fisher scoring
   unsigned int it = 0;
   while(it < maxit){
     // Compute w, z
     eta = x1 * start;
     mu = (eta).unaryExpr(Sigmoid());
     for(unsigned int j=0;j<n;++j){
       varmu = V(mu(j));
       xp = x1.row(j).norm(); // l2-norm
       s = std::abs(y_tilde(j) - mu(j)) * xp; 
       w(j,j) = varmu * Ew(mu(j), c, xp) - varmu * Edwymu(mu(j), c, xp) * xp;
       z(j) = (y_tilde(j) - mu(j)) * wc(s, c); 
     }
     
     // Fisher scoring
     A = x1.transpose() * w * x1;
     b = x1.transpose() * z;
     Dbeta = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
     // Dbeta = A.colPivHouseholderQr().solve(b);
     beta = start + Dbeta;
     
     // Check convergence
     double relE = std::sqrt(Dbeta.squaredNorm() / std::max(1e-20, start.squaredNorm()));
     ++it;
     
     if(verbose) {
       Rcpp::Rcout << "Relative error " << relE << " at iteration " << it << std::endl;
     }
     start = beta;
     
     if(relE <= tol){
       conv = 0;
       break;
     }
   }
   
   return Rcpp::List::create(
     Rcpp::Named("weights") = w.diagonal(),
     Rcpp::Named("z") = z,
     Rcpp::Named("coefficients") = beta,
     Rcpp::Named("maxit") = it,
     Rcpp::Named("conv") = conv
   );
 }

// probably a better way to handle this
Eigen::VectorXd roblogisticWmle1b(
    Eigen::ArrayXd& y,
    Eigen::MatrixXd& x,
    Eigen::VectorXd& start,
    double c = 4.685061,
    unsigned int maxit=200,
    double tol=1e-7
){
  unsigned int n = y.size();
  unsigned int p = x.cols();
  unsigned int p1 = p+1;
  
  // data storage
  Eigen::VectorXd mu(n),eta(n),z(n);
  Eigen::MatrixXd x1(n,p+1);
  Eigen::MatrixXd w = Eigen::MatrixXd::Zero(n,n);
  double varmu, xp, s;
  Eigen::MatrixXd A(p1,p1);
  Eigen::VectorXd b(p1),Dbeta(p1),beta(p1);
  
  // pre-computation
  x1.rightCols(p) = x;
  x1.col(0) = Eigen::VectorXd::Constant(n,1.0);
  
  // Fisher scoring
  unsigned int it = 0;
  while(it < maxit){
    // Compute w, z
    eta = x1 * start;
    mu = (eta).unaryExpr(Sigmoid());
    for(unsigned int j=0;j<n;++j){
      varmu = V(mu(j));
      xp = x1.row(j).norm(); // l2-norm
      s = std::abs(y(j) - mu(j)) * xp; 
      w(j,j) = varmu * Ew(mu(j), c, xp) - varmu * Edwymu(mu(j), c, xp) * xp;
      z(j) = (y(j) - mu(j)) * wc(s, c); 
    }
    
    // Fisher scoring
    A = x1.transpose() * w * x1;
    b = x1.transpose() * z;
    Dbeta = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
    beta = start + Dbeta;
    
    // Check convergence
    double relE = std::sqrt(Dbeta.squaredNorm() / std::max(1e-20, start.squaredNorm()));
    ++it;
    
    if(relE <= tol){
      start = beta;
      break;
    }
    else {
      start = beta;
    }
  }
  
  return beta;
}

//' Robust logistic regression estimator with Tukey's weights
//'
//' @param y a vector of responses
//' @param x a n x p matrix of design
//' @param beta a p-vector of parameter (starting values)
//' @param c tuning parameter for Tukey's weight (default value is 4.685061)
//' @param maxit max number of iteration for IRWLS
//' @param tol tolerance for stopping criterion
//' @param verbose print info
//' @export
// [[Rcpp::export]]
Rcpp::List roblogisticWmle(
     Eigen::ArrayXd& y,
     Eigen::MatrixXd& x,
     Eigen::VectorXd& start,
     double c = 4.685061,
     unsigned int maxit=200,
     double tol=1e-7,
     bool verbose=false
 ){
   unsigned int n = y.size();
   unsigned int p = x.cols();
   unsigned int p1 = p+1;
   
   // data storage
   Eigen::VectorXd mu(n),eta(n),z(n);
   Eigen::MatrixXd x1(n,p+1);
   Eigen::MatrixXd w = Eigen::MatrixXd::Zero(n,n);
   double varmu, xp, s;
   Eigen::MatrixXd A(p1,p1);
   Eigen::VectorXd b(p1),Dbeta(p1),beta(p1);
   
   // pre-computation
   x1.rightCols(p) = x;
   x1.col(0) = Eigen::VectorXd::Constant(n,1.0);
   
   // Fisher scoring
   unsigned int it = 0;
   while(it < maxit){
     // Compute w, z
     eta = x1 * start;
     mu = (eta).unaryExpr(Sigmoid());
     for(unsigned int j=0;j<n;++j){
       varmu = V(mu(j));
       xp = x1.row(j).norm(); // l2-norm
       s = std::abs(y(j) - mu(j)) * xp; 
       w(j,j) = varmu * Ew(mu(j), c, xp) - varmu * Edwymu(mu(j), c, xp) * xp - varmu * dEwymu(mu(j), c, xp);
       z(j) = (y(j) - mu(j)) * wc(s, c) - Ewymu(mu(j), c, xp); 
     }
     
     // Fisher scoring
     A = x1.transpose() * w * x1;
     b = x1.transpose() * z;
     Dbeta = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
     // Dbeta = A.colPivHouseholderQr().solve(b);
     beta = start + Dbeta;
     
     // Check convergence
     double relE = std::sqrt(Dbeta.squaredNorm() / std::max(1e-20, start.squaredNorm()));
     ++it;
     
     if(verbose) {
       Rcpp::Rcout << "Relative error " << relE << " at iteration " << it << std::endl;
     }
     
     if(relE <= tol){
       start = beta;
       break;
     }
     else {
       start = beta;
     }
   }
   
   return Rcpp::List::create(
     Rcpp::Named("weights") = w.diagonal(),
     Rcpp::Named("z") = z,
     Rcpp::Named("coefficients") = beta,
     Rcpp::Named("maxit") = it
   );
 }


//' Asymptotic variance of robust logistic regression estimator with Tukey's weights
//'
//' @param x a n x p matrix of design
//' @param beta a p-vector of parameter (starting values)
//' @param c tuning parameter for Tukey's weights (default value is 4.685061)
//' @export
// [[Rcpp::export]]
Rcpp::List roblogisticWmleVar(
     Eigen::MatrixXd& x,
     Eigen::VectorXd& start,
     double c = 4.685061
 ){
   unsigned int n = x.rows();
   unsigned int p = x.cols();
   unsigned int p1 = p+1;
   
   // data storage
   Eigen::VectorXd mu(n),eta(n);
   Eigen::MatrixXd x1(n,p+1);
   Eigen::MatrixXd w = Eigen::MatrixXd::Zero(n,n);
   Eigen::MatrixXd w1 = Eigen::MatrixXd::Zero(n,n);
   Eigen::MatrixXd w2 = Eigen::MatrixXd::Zero(n,n);
   Eigen::MatrixXd d = Eigen::MatrixXd::Zero(n,n);
   Eigen::MatrixXd d1 = Eigen::MatrixXd::Zero(n,n);
   double varmu, a, xp;
   Eigen::MatrixXd asym_var_rob(p1,p1), asym_var_rob1(p1,p1), asym_var_mle(p1,p1);
   
   // pre-computation
   x1.rightCols(p) = x;
   x1.col(0) = Eigen::VectorXd::Constant(n,1.0);
   
   // Compute w, z
   eta = x1 * start;
   mu = (eta).unaryExpr(Sigmoid());
   for(unsigned int j=0;j<n;++j){
     // Fisher info (MLE)
     varmu = V(mu(j));
     w2(j,j) = varmu;
     
     // Computation 
     xp = x1.row(j).norm(); // l2-norm
     w1(j,j) = varmu * Edwymu(mu(j),c,xp) * xp - varmu * Ew(mu(j),c,xp);
     w(j,j) = w1(j,j) - varmu * dEwymu(mu(j),c,xp);
     a = Ewymu(mu(j),c,xp);
     d1(j,j) = Ewymu2(mu(j),c,xp);
     d(j,j) = d1(j,j) - a * a;
   }
   
   Eigen::MatrixXd H = (x1.transpose() * w * x1).inverse();
   Eigen::MatrixXd H1 = (x1.transpose() * w1 * x1).inverse();
   asym_var_rob = H * x1.transpose() * d * x1 * H;
   asym_var_rob1 = H1 * x1.transpose() * d1 * x1 * H1;
   asym_var_mle = (x1.transpose() * w2 * x1).inverse();
   double efficiency = asym_var_mle.trace() / asym_var_rob.trace();
   double efficiency1 = asym_var_mle.trace() / asym_var_rob1.trace();
   
   return Rcpp::List::create(
     Rcpp::Named("V_rob") = asym_var_rob,
     Rcpp::Named("V_rob1") = asym_var_rob1,
     Rcpp::Named("V_mle") = asym_var_mle,
     Rcpp::Named("efficiency") = efficiency,
     Rcpp::Named("efficiency1") = efficiency1
   );
 }

//' Simulation of logistic regression
//'
//' @param x a n x p matrix of design
//' @param beta a p-vector of parameter (starting values)
//' @param seed for random number generator
//' @export
// [[Rcpp::export]]
Eigen::ArrayXd r_logistic(
    Eigen::VectorXd& beta,
    Eigen::MatrixXd& x,
    unsigned int seed
){
  // data storage
  unsigned int p = x.cols();
  unsigned int n = x.rows();
  Eigen::ArrayXd sig(n);
  Eigen::VectorXd eta(n);
  Eigen::MatrixXd x1(n,p+1);
  Eigen::ArrayXd y(n);
  double u;
  
  // pre-computation
  x1.rightCols(p) = x;
  x1.col(0) = Eigen::VectorXd::Constant(n,1.0);
  eta = x1 * beta;
  sig = eta.unaryExpr(Sigmoid());
  std::mt19937_64 engine(seed);  // Mersenne twister random number engine
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  
  for(unsigned int i(0);i<n;++i){
    u = unif(engine);
    y(i) = (sig(i) < u) ? 0.0 : 1.0;
  }
  return y;
}

//' Iterative bootstrap for robust logistic regression with inconsistent initial estimator with Tukey's weights
//'
//' @param x a n x p matrix of design
//' @param start an inconsistent estimator (also used as starting values)
//' @param c tuning parameter for Tukey's weight (default value is 4.685061)
//' @param H number of estimators for Monte Carlo approximation
//' @param maxit max number of iteration for IRWLS
//' @param tol tolerance for stopping criterion
//' @param verbose print info
//' @param seed for random number generator
//' @export
// [[Rcpp::export]]
Rcpp::List IBroblogisticWmle1(
     Eigen::MatrixXd& x,
     Eigen::VectorXd& start,
     double c = 4.685061,
     unsigned int H = 200,
     unsigned int maxit=200,
     double tol=1e-7,
     bool verbose=false,
     unsigned int seed=321
 ){
   unsigned int p = x.cols();
   unsigned int p1 = p+1;

   // data storage
   Eigen::MatrixXd boot(p1,H);
   Eigen::VectorXd Dbeta(p1),beta(p1),beta0(p1),pi(p1);

   // Iterative bootstrap
   beta0 = start;
   unsigned int it = 0;
   while(it < maxit){
     // Monte Carlo approximation
     // #pragma omp parallel for num_threads(ncores)
     for(unsigned int i=0; i<H; ++i){
         unsigned int se = seed + i;
         Eigen::ArrayXd y = r_logistic(beta0,x,se);
         boot.col(i) = roblogisticWmle1b(y,x,beta0,c);
     }
     pi = boot.rowwise().mean();
     Dbeta = start - pi;
     beta = beta0 + Dbeta;

     // Check convergence
     double relE = std::sqrt(Dbeta.squaredNorm() / std::max(1e-20, start.squaredNorm()));
     ++it;

     if(verbose) {
       Rcpp::Rcout << "Relative error " << relE << " at iteration " << it << std::endl;
     }
     beta0 = beta;

     if(relE <= tol){
       break;
     }
    }

   return Rcpp::List::create(
     Rcpp::Named("coefficients") = beta0,
     Rcpp::Named("maxit") = it,
     Rcpp::Named("diff") = Dbeta
   );
}

//' Stochastic approximation for robust logistic regression with inconsistent initial estimator with Tukey's weights
//'
//' @param x a n x p matrix of design
//' @param start an inconsistent estimator (also used as starting values)
//' @param c tuning parameter for Tukey's weight (default value is 4.685061)
//' @param maxit max number of iteration for IRWLS
//' @param tol tolerance for stopping criterion
//' @param verbose print info
//' @param seed for random number generator
//' @param k constant for stochastic approximation (1 by default)
//' @export
// [[Rcpp::export]]
Rcpp::List StocApproblogisticWmle1(
     Eigen::MatrixXd& x,
     Eigen::VectorXd& start,
     double c = 4.685061,
     unsigned int maxit=10000,
     double tol=1e-7,
     bool verbose=false,
     unsigned int seed=321,
     double k = 1
 ){
   unsigned int p = x.cols();
   unsigned int p1 = p+1;
   
   // data storage
   Eigen::VectorXd Dbeta(p1),Dinit(p1),beta(p1),beta0(p1),boot(p1);
   double alpha;
   
   // Iterative bootstrap
   beta0 = start;
   unsigned int it = 0;
   while(it < maxit){
     // Monte Carlo approximation
     unsigned int se = seed + it;
     Eigen::ArrayXd y = r_logistic(beta0,x,se);
     boot = roblogisticWmle1b(y,x,beta0,c);
     Dinit = start - boot;
     alpha = k / (it + 1);
     beta = beta0 + alpha * Dinit;
     
     // Check convergence
     Dbeta = beta - beta0;
     double relE = std::sqrt(Dbeta.squaredNorm() / std::max(1e-20, start.squaredNorm()));
     ++it;
     
     if(verbose) {
       Rcpp::Rcout << "Relative error " << relE << " at iteration " << it << std::endl;
     }
     beta0 = beta;
     
     if(relE <= tol){
       break;
     }
   }
   
   return Rcpp::List::create(
     Rcpp::Named("coefficients") = beta0,
     Rcpp::Named("maxit") = it,
     Rcpp::Named("diff") = Dbeta
   );
 }


// --------------
// WMLE (using lbfgs and not Fisher scoring)
// --------------
class wmle_logistic: public Numer::MFuncGrad
{
private:
  const Eigen::ArrayXd y;
  const Eigen::MatrixXd x;
  const double c;
  const unsigned int n = y.size();
  const unsigned int p = x.cols();

public:
  wmle_logistic(const Eigen::ArrayXd& y_,const Eigen::MatrixXd& x_, const double& c_) :
  y(y_), x(x_), c(c_) {}
  double f_grad(Numer::Constvec& beta, Numer::Refvec grad);
};

double wmle_logistic::f_grad(
    Numer::Constvec& beta,
    Numer::Refvec grad
){
  // data storage
  Eigen::ArrayXd sig(n);
  Eigen::VectorXd eta(n), mu(n), z(n);
  Eigen::MatrixXd x1(n,p+1);
  Eigen::VectorXd psi(p+1);
  Eigen::MatrixXd Dpsi(p+1,p+1);
  Eigen::MatrixXd w = Eigen::MatrixXd::Zero(n,n);
  double varmu, xp, s;
  
  // pre-computation
  x1.rightCols(p) = x;
  x1.col(0) = Eigen::VectorXd::Constant(n,1.0);
  eta = x1 * beta;
  mu = eta.unaryExpr(Sigmoid());

  // Compute of and grad
  for(unsigned int j=0;j<n;++j){
    double k = 1.0;
    if(y(j) == 1) k = -1.0;
    varmu = V(mu(j));
    xp = x1.row(j).norm(); // l2-norm
    s = std::abs(y(j) - mu(j)) * xp; 
    w(j,j) = varmu * wc(s, c) - varmu * wc1(s, c) * (y(j) - mu(j)) * xp * k;
    z(j) = (y(j) - mu(j)) * wc(s, c); 
  }
  
  psi = x1.transpose() * z;
  Dpsi = -x1.transpose() * w * x1;
  
  // objective function
  const double f = psi.squaredNorm() / 0.2e1 / n;
  
  // gradient
  grad = Dpsi.transpose() * psi/ n;
  
  return f;
}

//' MLE for logistic regression with misclassified responses
//'
//' @param y a vector of responses
//' @param x a n x p matrix of design
//' @param c tuning parameter for Tukey's weight (default value is 4.685061)
//' @export
// [[Rcpp::export]]
Rcpp::List logistic_wmle(
     Eigen::ArrayXd& y,
     Eigen::MatrixXd& x,
     double c
 ){
   // Regress
   unsigned int p = x.cols() + 1;
   double fopt;
   Eigen::VectorXd beta(p);
   beta.setZero();
   double prop = y.mean();
   beta(0) = std::log(prop) - std::log(1.0 - prop);
   wmle_logistic f(y,x,c);
   int res = Numer::optim_lbfgs(f,beta,fopt);
   
   return Rcpp::List::create(
     Rcpp::Named("coefficients") = beta,
     Rcpp::Named("fopt") = fopt,
     Rcpp::Named("status") = res
   );
}

Rcpp::List check(
    Eigen::ArrayXd& y,
    Eigen::MatrixXd& x,
    Eigen::VectorXd& beta,
    double c
){
  double fopt;
  Eigen::VectorXd grad(beta.size());
  wmle_logistic f(y,x,c);
  fopt = f.f_grad(beta,grad);
  
  return Rcpp::List::create(
    Rcpp::Named("of") = fopt,
    Rcpp::Named("grad") = grad
  );
}

// //' Iterative bootstrap for robust logistic regression with inconsistent initial estimator with Tukey's weights
// //'
// //' @param x a n x p matrix of design
// //' @param start an inconsistent estimator (also used as starting values)
// //' @param c tuning parameter for Tukey's weight (default value is 4.685061)
// //' @param H number of estimators for Monte Carlo approximation
// //' @param maxit max number of iteration for IRWLS
// //' @param tol tolerance for stopping criterion
// //' @param verbose print info
// //' @param seed for random number generator
// //' @export
// // [[Rcpp::export]]
// Rcpp::List logistic_wmle_ib(
//      Eigen::MatrixXd& x,
//      Eigen::VectorXd& start,
//      double c = 4.685061,
//      unsigned int H = 200,
//      unsigned int maxit=200,
//      double tol=1e-7,
//      bool verbose=false,
//      unsigned int seed=321
//  ){
//    unsigned int p = x.cols();
//    unsigned int p1 = p+1;
//    
//    // data storage
//    Eigen::MatrixXd boot(p1,H);
//    Eigen::VectorXd Dbeta(p1),beta(p1),beta0(p1),pi(p1);
//    
//    // Iterative bootstrap
//    beta0 = start;
//    unsigned int it = 0;
//    while(it < maxit){
//      // Monte Carlo approximation
//      // #pragma omp parallel for num_threads(ncores)
//      for(unsigned int i=0; i<H; ++i){
//        unsigned int se = seed + i;
//        Eigen::ArrayXd y = r_logistic(beta0,x,se);
//        boot.col(i) = logistic_wmle(y,x,c);
//      }
//      pi = boot.rowwise().mean();
//      Dbeta = start - pi;
//      beta = beta0 + Dbeta;
//      
//      // Check convergence
//      double relE = std::sqrt(Dbeta.squaredNorm() / std::max(1e-20, start.squaredNorm()));
//      ++it;
//      
//      if(verbose) {
//        Rcpp::Rcout << "Relative error " << relE << " at iteration " << it << std::endl;
//      }
//      beta0 = beta;
//      
//      if(relE <= tol){
//        break;
//      }
//    }
//    
//    return Rcpp::List::create(
//      Rcpp::Named("coefficients") = beta0,
//      Rcpp::Named("maxit") = it,
//      Rcpp::Named("diff") = Dbeta
//    );
//  }


// //' Stochastic approximation for robust logistic regression with inconsistent initial estimator with Tukey's weights
// //'
// //' @param x a n x p matrix of design
// //' @param start an inconsistent estimator (also used as starting values)
// //' @param c tuning parameter for Tukey's weight (default value is 4.685061)
// //' @param H number of estimators for Monte Carlo approximation
// //' @param maxit max number of iteration for IRWLS
// //' @param tol tolerance for stopping criterion
// //' @param verbose print info
// //' @param seed for random number generator
// //' @export
// // [[Rcpp::export]]
// Rcpp::List logistic_wmle_stocapp(
//      Eigen::MatrixXd& x,
//      Eigen::VectorXd& start,
//      double c = 4.685061,
//      unsigned int maxit=10000,
//      double tol=1e-5,
//      bool verbose=false,
//      unsigned int seed=321
//  ){
//    unsigned int p = x.cols();
//    unsigned int p1 = p+1;
//    
//    // data storage
//    Eigen::VectorXd Dinit(p1),Dbeta(p1),beta(p1),beta0(p1),boot(p1);
//    
//    // Stochastic approximation
//    beta0 = start;
//    unsigned int it = 0;
//    while(it < maxit){
//      unsigned int se = seed + it;
//      Eigen::ArrayXd y = r_logistic(beta0,x,se);
//      boot = logistic_wmle(y,x,c);
//      Dinit = start - boot;
//      double alpha = 1.0 / (it + 1.0);
//      beta = beta0 + alpha * Dinit;
//      
//      // Check convergence
//      Dbeta = beta - beta0;
//      double relE = std::sqrt(Dbeta.squaredNorm() / std::max(1e-20, start.squaredNorm()));
//      ++it;
//      
//      if(verbose) {
//        Rcpp::Rcout << "Relative error " << relE << " at iteration " << it << std::endl;
//      }
//      beta0 = beta;
//      
//      if(relE <= tol){
//        break;
//      }
//    }
//    
//    return Rcpp::List::create(
//      Rcpp::Named("coefficients") = beta0,
//      Rcpp::Named("maxit") = it,
//      Rcpp::Named("diff") = Dbeta
//    );
//  }

// OLD ATTEMPTS:

// //' Asymptotic variance of robust logistic regression estimator with Tukey's weights
// //'
// //' @param x a n x p matrix of design
// //' @param beta a p-vector of parameter (starting values)
// //' @param c tuning parameter for Tukey's weights (default value is 4.685061)
// //' @export
// // [[Rcpp::export]]
// Rcpp::List roblogisticMqleVar(
//      Eigen::MatrixXd& x,
//      Eigen::VectorXd& start,
//      double c = 4.685061
//  ){
//    unsigned int n = x.rows();
//    unsigned int p = x.cols();
//    unsigned int p1 = p+1;
//    
//    // data storage
//    Eigen::VectorXd mu(n),eta(n);
//    Eigen::MatrixXd x1(n,p+1);
//    Eigen::MatrixXd w = Eigen::MatrixXd::Zero(n,n);
//    Eigen::MatrixXd b = Eigen::MatrixXd::Zero(n,n);
//    Eigen::MatrixXd w2 = Eigen::MatrixXd::Zero(n,n);
//    double t1, t2, t3, t4, alpha, varmu, varmup, c2, c4, c6, c8, v2, v3, mu3, mu_3, mu5, mu_5, mu7, mu_7, mu9, mu_9;
//    double f1, f2, pr;
//    Eigen::MatrixXd asym_var_rob(p1,p1), asym_var_mle(p1,p1);
//    
//    // pre-computation
//    x1.rightCols(p) = x;
//    x1.col(0) = Eigen::VectorXd::Constant(n,1.0);
//    
//    // Compute w, z
//    eta = x1 * start;
//    mu = (eta).unaryExpr(Sigmoid());
//    for(unsigned int j=0;j<n;++j){
//      // Fisher info (MLE)
//      varmu = V(mu(j));
//      w2(j,j) = varmu;
//      // Computation of E[alpha] (expectation)
//      t1 = std::sqrt(varmu);
//      // r = y(j) / t1 - mu(j) / t1; // Pearson's residual
//      // bool is_zero = std::abs(r) > c;
//      t2 = mu(j) / t1;
//      bool is_below0 = t2 <= c;
//      bool is_below1 = 0.1e1 / t2 <= c;
//      bool is_zero = !is_below0 || !is_below1;
//      if(!is_zero) {
//        if(is_below0) {
//          pr = 0.1e1 - mu(j);
//        } else {
//          pr = mu(j);
//        }
//        c2 = c * c;
//        c4 = c2 * c2;
//        varmup = V1(mu(j));
//        t2 = varmup * varmup;
//        f1 = varmup / varmu / t1 / c4 - 0.2e1 / c4 - 0.2e1 * varmup / c2;
//        // f1 = 0.0;
//        alpha = 0.1e1 - 1.5e1 / c4 + (0.5e1 / c4 + t2 / c2 + t2 / c4) / varmu - 0.6e1 * varmu / c2  - t2 / varmu / varmu / 0.2e1 / c4;
//        w(j,j) = alpha * varmu * pr + f1;
//        // Computation of E[psi^2] (expectation)
//        f2 = 0.4e1 * varmu / c2 - 0.2e1 / c4 - t2 / c4 / varmu;
//        // f2 = 0.0;
//        v2 = varmu * varmu;
//        v3 = v2 * varmu;
//        c6 = c4 * c2;
//        c8 = c4 * c4;
//        t1 = mu(j);
//        t2 = t1 * t1;
//        t3 = 0.1e1 - t1;
//        t4 = t3 * t3;
//        mu3 = t1 * t2;
//        mu_3 = t3 * t4;
//        mu5 = mu3 * t2;
//        mu_5 = mu_3 * t4;
//        mu7 = mu5 * t2;
//        mu_7 = mu5 * t4;
//        mu9 = mu7 * t2;
//        mu_9 = mu_7 * t4;
//        b(j,j) = pr * (varmu + (mu9 + mu_9) / c8 / v3 - 0.4e1 * (mu7 + mu_7) / c6 / v2 + 0.6e1 * (mu5 + mu_5) / c4 / varmu - 0.4e1 * (mu3 + mu_3) / c2) + f2;
//        }
//    }
//    
//    Eigen::MatrixXd H = (x1.transpose() * w * x1).inverse();
//    asym_var_rob = H * x1.transpose() * b * x1 * H;
//    asym_var_mle = (x1.transpose() * w2 * x1).inverse();
//    double efficiency = asym_var_mle.trace() / asym_var_rob.trace();
//    
//    return Rcpp::List::create(
//      Rcpp::Named("V_rob") = asym_var_rob,
//      Rcpp::Named("V_mle") = asym_var_mle,
//      Rcpp::Named("efficiency") = efficiency
//    );
//  }