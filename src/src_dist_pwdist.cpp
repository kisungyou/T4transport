#include "RcppArmadillo.h"
#include <cmath>
#include "utility.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
Rcpp::List cpp_pwdist(const arma::mat& X, const arma::vec& p, 
                      const arma::mat& Y, const arma::vec& q, 
                      int maxiter, 
                      double abstol){
  // parameters
  int N = X.n_rows;
  int M = Y.n_rows;
  int d = Y.n_cols;
  
  // initial correspondence
  arma::mat cross_sqdist = util_pairwise_sqdist(X, Y);
  arma::mat iter_Gamma = util_plan_emd_R(p, q, cross_sqdist);
  arma::mat iter_P(d,d,fill::zeros);
  arma::mat iter_YP(M,d,fill::zeros);
  double old_cost = 100000000000.0;
  double new_cost = 0.0;
  double inc_cost = 0.0;
  
  // iteration 
  for (int it=0; it<maxiter; it++){
    // step 1. procrustes update for a current plan
    arma::mat U;
    arma::vec s;
    arma::mat V;
    arma::svd(U,s,V,arma::trans(Y)*iter_Gamma.t()*X);
    iter_P.reset();
    iter_P = U*V.t();
    
    // step 2. ot update
    iter_YP = Y*iter_P;
    for (int n=0; n<N; n++){
      for (int m=0; m<M; m++){
        cross_sqdist(n,m) = arma::accu(arma::square(X.row(n)-iter_YP.row(m)));
      }
    }
    iter_Gamma = util_plan_emd_R(p, q, cross_sqdist);
    
    // cost
    new_cost = arma::dot(cross_sqdist, iter_Gamma);
    if ((it > 2)&&(new_cost > old_cost)){
      break;
    }
    inc_cost = std::abs(new_cost-old_cost);
    old_cost = new_cost;
    if (inc_cost < abstol){
      break;
    }
  }
  
  // return the matrices
  return(Rcpp::List::create(
      Rcpp::Named("est_plan") = iter_Gamma,
      Rcpp::Named("est_P") = iter_P,
      Rcpp::Named("est_cost") = old_cost));
}