// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <limits>
#include <cmath>
#include "utility.h"

/* UTILITY FUNCTIONS
 * util_plan_entropic : compute the entropic-regularized plan using Sinkhorn algorithm
 * util_plan_emd_R    : compute the exact EMD plan using lpSolve - deprecated
 * util_mvrnorm       : generate multivariate normal random variables with given mean and covariance
 */

using namespace arma;
using namespace Rcpp;
using namespace std;


// util_plan_entropic ==========================================================
static double logsumexp(const arma::vec& x) {
  double xmax = x.max();
  return xmax + std::log(arma::sum(arma::exp(x - xmax)));
}

arma::mat util_plan_entropic(const arma::vec& a, const arma::vec& b, const arma::mat& C,
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

// // util_plan_emd_R ==========================================================
// arma::mat util_plan_emd_R(const arma::vec& a, const arma::vec& b, const arma::mat& C){
//   Rcpp::Environment pkg = Rcpp::Environment::namespace_env("T4transport");
//   Rcpp::Function aux_emd = pkg["aux_emd"];
//   //Function aux_emd("aux_emd");
//   SEXP result = aux_emd(wrap(a), wrap(b), wrap(C));
//   return(Rcpp::as<arma::mat>(result));
// }

// util_mvrnorm  ============================================================
// [[Rcpp::export]]
arma::mat util_mvrnorm(const arma::vec& par_mean, 
                       const arma::mat& par_cov,
                       int num_samples){
  arma::mat X = arma::trans(arma::mvnrnd(par_mean, par_cov, num_samples));
  return(X);
}

// util_pairwise_sqdist =====================================================
// [[Rcpp::export]]
arma::mat util_pairwise_sqdist(const arma::mat& X, const arma::mat& Y){
  int N = X.n_rows;
  int M = Y.n_rows;
  arma::mat D(N,M,fill::zeros);
  
  for (int n=0; n<N; n++){
    for (int m=0; m<M; m++){
      D(n,m) = arma::accu(arma::square(X.row(n)-Y.row(m)));
    }
  }
  
  return(D);
}

// compute the pairwise distance
// [[Rcpp::export]]
arma::mat util_pairwise_dist(const arma::mat& X){
  int N = X.n_rows;
  arma::mat D(N,N,fill::zeros);
  
  for (int n=0; n<(N-1); n++){
    for (int m=(n+1); m<N; m++){
      D(n,m) = std::sqrt(arma::accu(arma::square(X.row(n)-X.row(m))));
      D(m,n) = D(n,m);
    }
  }
  
  return(D);
}