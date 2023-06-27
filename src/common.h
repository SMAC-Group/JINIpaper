#ifndef _COMMON
#define _COMMON
// ------------------
// Headers
// ------------------
// [[Rcpp::depends(RcppEigen,RcppNumerical,BH)]]

// Enable C++14 via the Rcpp plugin
// [[Rcpp::plugins("cpp14")]]

// Libraries
#include <RcppEigen.h>
#include <random>
#include <math.h>
#include <numeric>
#include <algorithm>
#include <vector>
#include <RcppNumerical.h>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/distributions/beta.hpp>

// using namespace Eigen;
// using namespace std;

// Add a flag to enable OpenMP at compile time
// [[Rcpp::plugins(openmp)]]

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

// Change error handling for Boost functions
// Define a specific policy:
typedef boost::math::policies::policy<
  boost::math::policies::digits10<5>,
  boost::math::policies::overflow_error<boost::math::policies::ignore_error>
> my_policy;

// Sigmoid for logistic
struct Sigmoid {
  Sigmoid(){}
  const double operator()(const double& x) const {return 0.1e1 / (0.1e1 + std::exp(-x));}
};

struct Indicator {
  Indicator(){}
  const double operator()(const double &t) const {return (t <= 0.5) ? 0.0 : 1.0;}
};

// Function for robust estimation
// Tukey's weight:
inline double wc(double x, double c){return (std::abs(x) <= c) ? 0.1e1 - 0.2e1 * x*x/c/c + x*x*x*x/c/c/c/c : 0.0;}
inline double wc1(double x, double c){return (std::abs(x) <= c) ? - 0.4e1 * x/c/c + 0.4e1 * x*x*x/c/c/c/c : 0.0;}
// Huber's weight:
inline double wH(double x, double c){return (std::abs(x) <= c) ? 0.1e1 : c / std::abs(x); }
// from https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> inline constexpr
  int signum(T x, std::false_type is_signed) {
    return T(0) < x;
}
template <typename T> inline constexpr
  int signum(T x, std::true_type is_signed) {
    return (T(0) < x) - (x < T(0));
}
template <typename T> inline constexpr
  int signum(T x) {
    return signum(x, std::is_signed<T>());
}
inline double wH1(double x, double c){return (std::abs(x) <= c) ? 0.0 : c * signum(x); }

#endif
