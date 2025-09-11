#include "RcppArmadillo.h"
#include "elementary.h"
#include "utility.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

arma::mat aux_freemedian_sqdist(const arma::mat& source,
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
arma::mat cpp_single_barycenter(const arma::field<arma::mat>& measures, // support points
                                const arma::field<arma::vec>& marginals, // marginal measures
                                const arma::vec& weights, // weights for each measure 
                                arma::mat& init_support //initial support points)
){
  // TEMPORARY PARAM
  bool printer = false;
  int maxiter = 4;      // for quick evaluation
  double abstol = 1e-5; // for quick evaluation
  
  // PARAMETERS
  int N = measures.n_elem;    // number of measures
  int P = measures(0).n_cols; // problem dimensionality
  double NN = static_cast<double>(N);
  int num_support = init_support.n_rows;
  
  // INITIALIZATION
  arma::mat X_old = init_support;
  arma::mat X_new(num_support, P, fill::zeros);
  
  // uniform target weights
  arma::vec target_weights = arma::ones(num_support)/static_cast<double>(num_support);
  
  // temporary variables
  arma::field<arma::mat> now_dist2(N);
  arma::field<arma::mat> now_plans(N);
  
  // compute all pairwise distances
  for (int n=0; n<N; n++){
    now_dist2(n).reset();
    now_dist2(n) = aux_freemedian_sqdist(X_old, measures(n));
  }
  // compute all plans
  for (int n=0; n<N; n++){
    now_plans(n).reset();
    now_plans(n) = util_plan_emd_R(target_weights, 
              marginals(n), now_dist2(n));
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
      now_dist2(n) = aux_freemedian_sqdist(X_new, measures(n));
    }
    // compute all plans
    for (int n=0; n<N; n++){
      now_plans(n).reset();
      now_plans(n) = util_plan_emd_R(target_weights, 
                marginals(n), now_dist2(n));
    }
    
    // compute the updated cost
    new_cost = 0.0;
    for (int n=0; n<N; n++){
      new_cost += weights(n) * arma::accu(now_plans(n) % now_dist2(n));
    }
    
    // updating : if cost increased? stop it as an invalid update.
    if (new_cost >= old_cost){
      break;
    } 
    
    // updating : if decreased? replace and stop conditionally
    if ((std::abs(old_cost - new_cost)/std::abs(old_cost)) < abstol){
      X_old = X_new;
      record_cost(it+1) = new_cost;
      break;
    }  
    
    // updating : for the rest of the case, just update
    record_cost(it+1) = new_cost;
    X_old = X_new;
    old_cost = new_cost;
    
    if (printer){
      Rcpp::Rcout << "iteration " << it+1 << " complete with cost=" << old_cost << std::endl;
    }
  }
  
  // return
  return(X_old);
}

// [[Rcpp::export]]
Rcpp::List cpp_free_bary_gradient_init(const arma::field<arma::mat>& measures, // support points
                                       const arma::field<arma::vec>& marginals, // marginal measures
                                       const arma::vec& weights, // weights for each measure 
                                       int maxiter, // maximum number of iterations
                                       double abstol, // absolute tolerance for convergence
                                       arma::mat& init_support //initial support points
){
  // TEMPORARY PARAM
  bool printer = false;
  
  // PARAMETERS
  int N = measures.n_elem;    // number of measures
  int P = measures(0).n_cols; // problem dimensionality
  double NN = static_cast<double>(N);
  int num_support = init_support.n_rows;
  
  // INITIALIZATION
  arma::mat X_old = init_support;
  arma::mat X_new(num_support, P, fill::zeros);
  
  // uniform target weights
  arma::vec target_weights = arma::ones(num_support)/static_cast<double>(num_support);
  
  // temporary variables
  arma::field<arma::mat> now_dist2(N);
  arma::field<arma::mat> now_plans(N);
  
  // compute all pairwise distances
  for (int n=0; n<N; n++){
    now_dist2(n).reset();
    now_dist2(n) = aux_freemedian_sqdist(X_old, measures(n));
  }
  // compute all plans
  for (int n=0; n<N; n++){
    now_plans(n).reset();
    now_plans(n) = util_plan_emd_R(target_weights, 
              marginals(n), now_dist2(n));
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
      now_dist2(n) = aux_freemedian_sqdist(X_new, measures(n));
    }
    // compute all plans
    for (int n=0; n<N; n++){
      now_plans(n).reset();
      now_plans(n) = util_plan_emd_R(target_weights, 
                marginals(n), now_dist2(n));
    }
    
    // compute the updated cost
    new_cost = 0.0;
    for (int n=0; n<N; n++){
      new_cost += weights(n) * arma::accu(now_plans(n) % now_dist2(n));
    }
    
    // updating : if cost increased? stop it as an invalid update.
    if (new_cost >= old_cost){
      break;
    } 
    
    // updating : if decreased? replace and stop conditionally
    if ((std::abs(old_cost - new_cost)/std::abs(old_cost)) < abstol){
      X_old = X_new;
      record_cost(it+1) = new_cost;
      break;
    }  
    
    // updating : for the rest of the case, just update
    record_cost(it+1) = new_cost;
    X_old = X_new;
    old_cost = new_cost;
    
    if (printer){
      Rcpp::Rcout << "iteration " << it+1 << " complete with cost=" << old_cost << std::endl;
    }
  }
  
  // return
  return(Rcpp::List::create(
      Rcpp::Named("support") = X_old,
      Rcpp::Named("cost_history") = record_cost
  ));
}