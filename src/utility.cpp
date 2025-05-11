#include <RcppArmadillo.h>
#include "utility.hpp"

/* UTILITY FUNCTIONS
 * util_entropic_plan : compute the entropic-regularized plan using Sinkhorn algorithm
 * 
 * 
 */

using namespace arma;
using namespace Rcpp;
using namespace std;


// util_entropic_plan ==========================================================
static double logsumexp(const arma::vec& x) {
  double xmax = x.max();
  return xmax + std::log(arma::sum(arma::exp(x - xmax)));
}
arma::mat util_entropic_plan(const arma::mat& C, const arma::vec& a, const arma::vec& b,
                             double lambda, int maxiter, double abstol){
  int M = C.n_rows;
  int N = C.n_cols;
  
  arma::vec f(M, arma::fill::zeros);
  arma::vec g(N, arma::fill::zeros);
  arma::vec logKsum(std::max(M, N));
  
  arma::vec marg_u(M,fill::zeros);
  arma::vec marg_v(N,fill::zeros);
  
  for (int it = 0; it < maxiter; it++) {
    // Update f
    for (int m = 0; m < M; m++) {
      arma::vec row = (-C.row(m).t() + g) / lambda;
      logKsum(m) = logsumexp(row);
    }
    f = lambda * (arma::log(a) - logKsum.subvec(0, M - 1));
    
    // Update g
    for (int n = 0; n < N; n++) {
      arma::vec col = (-C.col(n) + f) / lambda;
      logKsum(n) = logsumexp(col);
    }
    g = lambda * (arma::log(b) - logKsum.subvec(0, N - 1));
    
    // Marginal check
    marg_u.zeros();
    marg_v.zeros();
    for (int m = 0; m < M; m++) {
      for (int n = 0; n < N; n++) {
        double val = std::exp((f(m) + g(n) - C(m, n)) / lambda);
        marg_u(m) += val;
        marg_v(n) += val;
      }
    }
    
    double err = arma::norm(marg_u - a, 1) + arma::norm(marg_v - b, 1);
    if (err < abstol) {
      break;
    }
  }
  
  // Final plan
  arma::mat plan(M, N, arma::fill::zeros);
  for (int m = 0; m < M; m++) {
    for (int n = 0; n < N; n++) {
      plan(m, n) = std::exp((f(m) + g(n) - C(m, n)) / lambda);
    }
  }
  
  return plan;
}