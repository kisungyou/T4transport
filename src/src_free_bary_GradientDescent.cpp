#include "RcppArmadillo.h"
#include "elementary.h"
#include "utility.h"
#include "utility_EMD.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

arma::mat aux_freebary_sqdist(const arma::mat& source,
                              const arma::mat& target){
  // get the params
  int M = source.n_rows;
  int N = target.n_rows;
  
  // squared distances
  arma::mat Dsq(M,N,fill::zeros);
  for (int m=0; m<M; m++){
    for (int n=0; n<N; n++){
      Dsq(m,n) = std::pow(arma::norm(source.row(m) - target.row(n), 2), 2.0);
    }
  }
  
  return(Dsq);
}

// [[Rcpp::export]]
Rcpp::List cpp_free_bary_gradient(const arma::field<arma::mat>& measures, // support points
                                 const arma::field<arma::vec>& marginals, // marginal measures
                                 const arma::vec& weights, // weights for each measure 
                                 int num_support, // number of support points
                                 int maxiter, // maximum number of iterations
                                 double abstol // absolute tolerance for convergence
){
  // PARAMETERS
  int N = measures.n_elem;    // number of measures
  int P = measures(0).n_cols; // problem dimensionality
  
  // INITIALIZE
  // free-support target measure by gaussian sampling
  arma::rowvec par_mean(P, fill::zeros);
  arma::mat par_cov(P,P,fill::zeros);
  for (int n=0; n<N; n++){
    par_mean += (1.0/static_cast<double>(N))*arma::mean(measures(n), 0);
    par_cov  += (1.0/static_cast<double>(N))*arma::cov(measures(n)); 
  }
  //par_cov = arma::diagmat(par_cov)*2.0; // zero-out off-diagonals and make it dispersed
  par_cov = arma::diagmat(par_cov)*2.0 + 1e-8 * arma::eye(P,P); // jitter added
  arma::mat X_old = arma::trans(arma::mvnrnd(par_mean.t(), par_cov, num_support));
  arma::mat X_new(num_support, P, fill::zeros);
  
  // uniform target weights
  arma::vec target_weights = arma::ones(num_support)/static_cast<double>(num_support);
  
  // temporary variables
  arma::field<arma::mat> now_dist2(N);
  arma::field<arma::mat> now_plans(N);
  
  // compute all pairwise distances
  for (int n=0; n<N; n++){
    now_dist2(n).reset();
    now_dist2(n) = aux_freebary_sqdist(X_old, measures(n));
  }
  // compute all plans
  for (int n=0; n<N; n++){
    now_plans(n).reset();
    //now_plans(n) = util_plan_emd_R(target_weights, marginals(n), now_dist2(n));
    now_plans(n) = util_plan_emd_C(target_weights, marginals(n), now_dist2(n));
  }
  
  // cost
  double old_cost = 0.0;
  double new_cost = 0.0;
  for (int n=0; n<N; n++){
    old_cost += weights(n) * arma::accu(now_plans(n) % now_dist2(n));
  }
  
  // iteration recorder
  arma::vec record_cost(maxiter+1, fill::zeros);
  record_cost(0) = old_cost;

  // MAIN COMPUTATION
  for (int it=0; it<maxiter; it++){
    // update each entry, starting with cleaning
    double scalar = 0.0;
    X_new.fill(0.0);
    for (int i=0; i<num_support; i++){
      // iterate over all plans
      for (int n=0; n<N; n++){
        int m_n = measures(n).n_rows; // number of support points in measure n
        for (int j=0; j<m_n; j++){
          scalar = weights(n)*(1.0/target_weights(i))*(now_plans(n)(i,j));
          X_new.row(i) += scalar*(measures(n).row(j));
        }
      }
    }
    
    // compute all pairwise distances
    for (int n=0; n<N; n++){
      now_dist2(n).reset();
      now_dist2(n) = aux_freebary_sqdist(X_new, measures(n));
    }
    // compute all plans
    for (int n=0; n<N; n++){
      now_plans(n).reset();
      //now_plans(n) = util_plan_emd_R(target_weights, marginals(n), now_dist2(n));
      now_plans(n) = util_plan_emd_C(target_weights, marginals(n), now_dist2(n));
    }
    
    // compute the updated cost
    new_cost = 0.0;
    for (int n=0; n<N; n++){
      new_cost += weights(n) * arma::accu(now_plans(n) % now_dist2(n));
    }
    
    // updating : if cost increased? stop it as an invalid update.
    if (new_cost > old_cost * (1.0 + 1e-12)){
      break;
    }
    // if (new_cost >= old_cost){
    //   break;
    // } 
    
    // updating : if decreased? replace and stop conditionally
    auto rel_change = [](double a, double b){
      double denom = std::max(1.0, std::abs(a));
      return std::abs(a - b)/denom;
    };
    if (rel_change(old_cost, new_cost) < abstol){
      X_old = X_new;
      record_cost(it+1) = new_cost;
      break;
    }
    
    // if ((std::abs(old_cost - new_cost)/std::abs(old_cost)) < abstol){
    //   X_old = X_new;
    //   record_cost(it+1) = new_cost;
    //   break;
    // }  
    
    // updating : for the rest of the case, just update
    record_cost(it+1) = new_cost;
    X_old = X_new;
    old_cost = new_cost;
  }
  
  // return
  return(Rcpp::List::create(
    Rcpp::Named("support") = X_old,
    Rcpp::Named("cost_history") = record_cost
  ));
}