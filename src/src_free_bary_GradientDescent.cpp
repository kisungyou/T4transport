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
  int M = source.n_rows;
  int N = target.n_rows;
  
  arma::mat Dsq(M, N, fill::zeros);
  for (int m=0; m<M; m++){
    for (int n=0; n<N; n++){
      Dsq(m,n) = std::pow(arma::norm(source.row(m) - target.row(n), 2), 2.0);
    }
  }
  return Dsq;
}

inline double aux_freebary_relchange(double a, double b){
  double denom = std::max(1.0, std::abs(a));
  return std::abs(a - b)/denom;
}

// full barycentric-projection destination (alpha = 1 target)
arma::mat aux_freebary_fullstep(const arma::field<arma::mat>& measures,
                                const arma::field<arma::mat>& plans,
                                const arma::vec& weights,
                                const arma::vec& target_weights){
  int N = measures.n_elem;
  int M = target_weights.n_elem;
  int P = measures(0).n_cols;
  
  arma::vec inv_target_weights = 1.0 / target_weights;
  arma::mat X_full(M, P, fill::zeros);
  
  for (int n=0; n<N; n++){
    // plans(n): M x m_n
    // measures(n): m_n x P
    arma::mat tmp = plans(n) * measures(n);   // M x P
    tmp.each_col() %= inv_target_weights;     // row i divided by v_i
    X_full += weights(n) * tmp;
  }
  return X_full;
}

// shared core
Rcpp::List free_bary_gradient_core(const arma::field<arma::mat>& measures,
                                   const arma::field<arma::vec>& marginals,
                                   const arma::vec& weights,
                                   const arma::mat& X_init,
                                   int maxiter,
                                   double abstol,
                                   double alpha){
  if ((alpha <= 0.0) || (alpha > 1.0)){
    Rcpp::stop("free_bary_gradient_core: 'alpha' should lie in (0,1].");
  }
  
  int N = measures.n_elem;
  int M = X_init.n_rows;
  int P = X_init.n_cols;
  
  arma::mat X_old = X_init;
  arma::mat X_full(M, P, fill::zeros);
  arma::mat X_prop(M, P, fill::zeros);
  
  arma::vec target_weights = arma::ones(M)/static_cast<double>(M);
  
  arma::field<arma::mat> now_dist2(N);
  arma::field<arma::mat> now_plans(N);
  
  arma::field<arma::mat> prop_dist2(N);
  arma::field<arma::mat> prop_plans(N);
  
  double old_cost = 0.0;
  for (int n=0; n<N; n++){
    now_dist2(n) = aux_freebary_sqdist(X_old, measures(n));
    now_plans(n) = util_plan_emd_C(target_weights, marginals(n), now_dist2(n));
    old_cost += weights(n) * arma::accu(now_plans(n) % now_dist2(n));
  }
  
  arma::vec record_cost(maxiter+1, fill::zeros);
  record_cost(0) = old_cost;
  int last_filled = 0;
  
  for (int it=0; it<maxiter; it++){
    // alpha = 1 target
    X_full = aux_freebary_fullstep(measures, now_plans, weights, target_weights);
    
    // damped proposal
    X_prop = (1.0 - alpha)*X_old + alpha*X_full;
    
    // evaluate proposal
    double new_cost = 0.0;
    for (int n=0; n<N; n++){
      prop_dist2(n) = aux_freebary_sqdist(X_prop, measures(n));
      prop_plans(n) = util_plan_emd_C(target_weights, marginals(n), prop_dist2(n));
      new_cost += weights(n) * arma::accu(prop_plans(n) % prop_dist2(n));
    }
    
    // reject cost-increasing move
    if (new_cost > old_cost * (1.0 + 1e-12)){
      break;
    }
    
    // accept
    X_old = X_prop;
    now_dist2 = prop_dist2;
    now_plans = prop_plans;
    
    record_cost(it+1) = new_cost;
    last_filled = it+1;
    
    if (aux_freebary_relchange(old_cost, new_cost) < abstol){
      old_cost = new_cost;
      break;
    }
    
    old_cost = new_cost;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("support")      = X_old,
    Rcpp::Named("cost_history") = record_cost.subvec(0, last_filled),
    Rcpp::Named("alpha")        = alpha,
    Rcpp::Named("niter")        = last_filled
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_free_bary_gradient_damped(const arma::field<arma::mat>& measures,
                                         const arma::field<arma::vec>& marginals,
                                         const arma::vec& weights,
                                         int num_support,
                                         int maxiter,
                                         double abstol,
                                         double alpha){
  if (num_support < 2){
    Rcpp::stop("cpp_free_bary_gradient_damped: 'num_support' should be at least 2.");
  }
  
  int N = measures.n_elem;
  int P = measures(0).n_cols;
  
  // same random Gaussian initialization as your current code
  arma::rowvec par_mean(P, fill::zeros);
  arma::mat par_cov(P, P, fill::zeros);
  for (int n=0; n<N; n++){
    par_mean += (1.0/static_cast<double>(N)) * arma::mean(measures(n), 0);
    par_cov  += (1.0/static_cast<double>(N)) * arma::cov(measures(n));
  }
  par_cov = arma::diagmat(par_cov)*2.0 + 1e-8 * arma::eye(P, P);
  
  arma::mat X_init = arma::trans(arma::mvnrnd(par_mean.t(), par_cov, num_support));
  
  return free_bary_gradient_core(
    measures, marginals, weights, X_init, maxiter, abstol, alpha
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_free_bary_gradient_damped_init(const arma::field<arma::mat>& measures,
                                              const arma::field<arma::vec>& marginals,
                                              const arma::vec& weights,
                                              int maxiter,
                                              double abstol,
                                              double alpha,
                                              const arma::mat& init_support){
  int P = measures(0).n_cols;
  
  if (init_support.n_cols != static_cast<arma::uword>(P)){
    Rcpp::stop("cpp_free_bary_gradient_damped_init: 'init_support' has wrong number of columns.");
  }
  if (init_support.n_rows < 2){
    Rcpp::stop("cpp_free_bary_gradient_damped_init: need at least 2 support points.");
  }
  
  return free_bary_gradient_core(
    measures, marginals, weights, init_support, maxiter, abstol, alpha
  );
}