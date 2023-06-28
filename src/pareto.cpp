// Robust Pareto regression (VGLM) with Tukey's weight (inconsistent)
#include "common.h"

// ------------------
// Implementation
// ------------------
struct Exp {
  Exp(){}
  const double operator()(const double& x) const {return std::exp(x);}
};

// Standard case (vs robust)

//' Pareto regression (MLE)
//'
//' @param y a vector of responses
//' @param x a n x p matrix of design
//' @param start is a p-vector of parameters (starting values)
//' @param maxit max number of iteration for IRWLS
//' @param tol tolerance for stopping criterion
//' @param verbose print info
//' @export
// [[Rcpp::export]]
Rcpp::List paretoMle(
     Eigen::ArrayXd& y,
     Eigen::MatrixXd& x,
     Eigen::VectorXd& start,
     unsigned int maxit=200,
     double tol=1e-7,
     bool verbose=false
 ){
   unsigned int n = y.size();
   unsigned int p = x.cols();
   unsigned int p1 = p+1;
   unsigned int conv = 1;
   
   // data storage
   Eigen::VectorXd mu(n),eta(n),z(n);
   Eigen::MatrixXd x1(n,p+1);
   double y0;
   Eigen::MatrixXd A(p1,p1);
   Eigen::VectorXd b(p1),Dbeta(p1),beta(p1);
   
   // pre-computation
   x1.rightCols(p) = x;
   x1.col(0) = Eigen::VectorXd::Constant(n,1.0);
   y0 = y.minCoeff();
   A = x1.transpose() * x1;
   
   // Fisher scoring
   unsigned int it = 0;
   while(it < maxit){
     // Compute w, z
     eta = x1 * start;
     mu = (eta).unaryExpr(Exp());
     for(unsigned int j=0;j<n;++j){
       // if(y0 / y(j) <= 0.0) { // handle case when y is too large
       //   z(j) = 0.0;
       // } else {
         z(j) = 0.1e1 + mu(j) * std::log(y0 / y(j)); 
       // }
     }
     
     // Fisher scoring
     b = x1.transpose() * z;
     Dbeta = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
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
     
     if(!std::isfinite(relE)) {
       if(verbose) {
         Rcpp::Rcout << "Algorithm stopped because of non-finite difference!" << std::endl;
       }
       conv = 1;
       break;
     }
   }
   
   return Rcpp::List::create(
     Rcpp::Named("coefficients") = beta,
     Rcpp::Named("scale") = y0,
     Rcpp::Named("maxit") = it,
     Rcpp::Named("conv") = conv
   );
}

// Robust case (TODO: Fisher scoring)

//' Robust Pareto regression initial estimator (inconsistent) with Tukey's weights
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
Rcpp::List paretoWmle1(
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
   unsigned int conv = 1;
   
   // data storage
   Eigen::VectorXd mu(n),eta(n),z0(n),w0(n);
   Eigen::MatrixXd x1(n,p+1);
   // Eigen::MatrixXd w = Eigen::MatrixXd::Zero(n,n);
   double xp, s, y0, t1, t11, wT, t2;
   Eigen::MatrixXd A(p1,p1);
   Eigen::VectorXd b(p1),Dbeta(p1),beta(p1);
   
   // pre-computation
   x1.rightCols(p) = x;
   x1.col(0) = Eigen::VectorXd::Constant(n,1.0);
   y0 = y.minCoeff();

   // Newton-Raphson method
   unsigned int it = 0;
   while(it < maxit){
     // Compute w, z
     eta = x1 * start;
     mu = (eta).unaryExpr(Exp());
     std::vector<int> ind;
     for(unsigned int j=0;j<n;++j){
       t2 = 1;
       xp = x1.row(j).norm(); // l2-norm but any norm would be possible
       t1 = mu(j) * std::log(y0 / y(j));
       if(t1 < -1){t2 = -1;}
       t11 = 0.1e1 + t1;
       s = std::abs(t11) * xp; 
       wT = wc(s, c);
       z0(j) = t11 * wT; 
       w0(j) = wc1(s, c) * t11 * t1 * xp * t2 + wT * t1;
       if(w0(j)<0){ind.push_back(j);}
     }
     
     // If more than 50% of the weights are okay, continue
     unsigned int m = ind.size();
     
     double ratio = (double)m / (double)n;
     if(verbose) {
       Rcpp::Rcout << "There is " << 100 - 100 * ratio << " % of negative weights at iteration " << it << std::endl;
     }
     
     if(m < std::round(0.5 * n)){
       conv = 1;
       break;
     }
     
     
     // Fisher scoring
     Eigen::MatrixXd x2(m,p+1);
     Eigen::MatrixXd w = Eigen::MatrixXd::Zero(m,m);
     Eigen::VectorXd z(m);
     for(unsigned int i=0;i<m;++i) {
       x2.row(i) = x1.row(ind[i]);
       z(i) = z0(ind[i]);
       w(i,i) = w0(ind[i]);
     }

     A = -x2.transpose() * w * x2;
     b = x2.transpose() * z;
     Dbeta = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
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
     
     if(!std::isfinite(relE)) {
       if(verbose) {
         Rcpp::Rcout << "Algorithm stopped because of non-finite difference!" << std::endl;
       }
       conv = 1;
       break;
     }
   }
   
   return Rcpp::List::create(
     Rcpp::Named("coefficients") = beta,
     Rcpp::Named("scale") = y0,
     Rcpp::Named("maxit") = it,
     Rcpp::Named("conv") = conv
   );
 }

//' Simulation of Pareto regression
//'
//' @param x a n x p matrix of design
//' @param beta a p-vector of parameter (starting values)
//' @param k minimum
//' @param seed for random number generator
//' @export
// [[Rcpp::export]]
Eigen::ArrayXd r_pareto(
     Eigen::VectorXd& beta,
     double k,
     Eigen::MatrixXd& x,
     unsigned int seed
 ){
   // data storage
   unsigned int p = x.cols();
   unsigned int n = x.rows();
   Eigen::ArrayXd mu(n);
   Eigen::VectorXd eta(n);
   Eigen::MatrixXd x1(n,p+1);
   Eigen::ArrayXd y(n);
   double u;
   
   // pre-computation
   x1.rightCols(p) = x;
   x1.col(0) = Eigen::VectorXd::Constant(n,1.0);
   eta = x1 * beta;
   mu = eta.unaryExpr(Exp());
   std::mt19937_64 engine(seed);  // Mersenne twister random number engine
   std::uniform_real_distribution<double> unif(0.0, 1.0);
   
   for(unsigned int i(0);i<n;++i){
     u = unif(engine);
     y(i) = k / std::pow(u, 1.0 / mu(i));
   }
   return y;
 }



//' Robust Pareto regression initial estimator (inconsistent) with Tukey's weights
//'
//' @param start a p-vector of parameter (starting values)
//' @param y a vector of responses
//' @param x a n x p matrix of design
//' @param c tuning parameter for Tukey's weight (default value is 4.685061)
//' @export
// [[Rcpp::export]]
double paretoWmle_of(
     Eigen::VectorXd& start,
     Eigen::ArrayXd& y,
     Eigen::MatrixXd& x,
     double c = 4.685061
 ){
   unsigned int n = y.size();
   unsigned int p = x.cols();
   unsigned int p1 = p+1;
   
   // data storage
   Eigen::VectorXd mu(n),eta(n),z(n);
   Eigen::MatrixXd x1(n,p1);
   double xp, s, y0, t1, t11, wT;
   Eigen::MatrixXd A(p1,p1);
   Eigen::VectorXd b(p1),Dbeta(p1),beta(p1);
   
   // pre-computation
   x1.rightCols(p) = x;
   x1.col(0) = Eigen::VectorXd::Constant(n,1.0);
   y0 = y.minCoeff();
   
   // Compute w, z
   eta = x1 * start;
   mu = (eta).unaryExpr(Exp());
   for(unsigned int j=0;j<n;++j){
     xp = x1.row(j).norm(); // l2-norm but any norm would be possible
     t1 = mu(j) * std::log(y0 / y(j));
     t11 = 0.1e1 + t1;
     s = std::abs(t11) * xp; 
     wT = wc(s, c);
     z(j) = t11 * wT; 
   }
   
   const double of = (x1.transpose() * z).squaredNorm();
   
   return of;
 }

//' Robust Pareto regression initial estimator (inconsistent) with Huber's weights
//'
//' @param y a vector of responses
//' @param x a n x p matrix of design
//' @param beta a p-vector of parameter (starting values)
//' @param c tuning parameter for Huber's weight (default value is 1.345)
//' @param maxit max number of iteration for IRWLS
//' @param tol tolerance for stopping criterion
//' @param verbose print info
//' @export
// [[Rcpp::export]]
Rcpp::List paretoWmle1H(
     Eigen::ArrayXd& y,
     Eigen::MatrixXd& x,
     Eigen::VectorXd& start,
     double c = 1.345,
     unsigned int maxit=200,
     double tol=1e-7,
     bool verbose=false
 ){
   unsigned int n = y.size();
   unsigned int p = x.cols();
   unsigned int p1 = p+1;
   unsigned int conv = 1;
   
   // data storage
   Eigen::VectorXd mu(n),eta(n),z(n);
   Eigen::MatrixXd x1(n,p+1);
   Eigen::MatrixXd w = Eigen::MatrixXd::Zero(n,n);
   double xp, s, y0, t1, t11, wT, t2;
   Eigen::MatrixXd A(p1,p1);
   Eigen::VectorXd b(p1),Dbeta(p1),beta(p1);
   
   // pre-computation
   x1.rightCols(p) = x;
   x1.col(0) = Eigen::VectorXd::Constant(n,1.0);
   y0 = y.minCoeff();
   
   // Newton-Raphson method
   unsigned int it = 0;
   while(it < maxit){
     // Compute w, z
     eta = x1 * start;
     mu = (eta).unaryExpr(Exp());
     for(unsigned int j=0;j<n;++j){
       t2 = 1;
       xp = x1.row(j).norm(); // l2-norm but any norm would be possible
       t1 = mu(j) * std::log(y0 / y(j));
       if(t1 < -1){t2 = -1;}
       if(t1 == -1){t2 = 0.0;}
       t11 = 0.1e1 + t1;
       s = std::abs(t11) * xp; 
       wT = wH(s, c);
       z(j) = t11 * wT; 
       w(j,j) = wH1(s, c) * t11 * t1 * xp * t2 + wT * t1;
     }
     
     
     // Fisher scoring
     A = -x1.transpose() * w * x1;
     b = x1.transpose() * z;
     Dbeta = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
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
     
     if(!std::isfinite(relE)) {
       conv = 1;
       break;
       if(verbose) {
         Rcpp::Rcout << "Algorithm stopped because of non-finite difference!" << std::endl;
       }
     }
   }
   
   return Rcpp::List::create(
     Rcpp::Named("coefficients") = beta,
     Rcpp::Named("scale") = y0,
     Rcpp::Named("maxit") = it,
     Rcpp::Named("conv") = conv
   );
 }