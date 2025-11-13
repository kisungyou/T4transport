#include "RcppArmadillo.h"
#include "elementary.h"
#include "utility.h"
#include "utility_EMD.h"

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
      now_dist2(n) = aux_freemedian_sqdist(X_new, measures(n));
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
      now_dist2(n) = aux_freemedian_sqdist(X_new, measures(n));
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

// [[Rcpp::export]]
Rcpp::List cpp_free_median_PF(const arma::field<arma::mat>& measures,   // X_n: (m_n x P)
                              const arma::field<arma::vec>& marginals,  // b_n: (m_n)
                              const arma::vec& weights,                 // pi_n (size N, sum to 1)
                              int maxiter,
                              double abstol,
                              arma::mat& init_support)
{
  // OPTINAL CONSTANTS
  const double irls_eps = 1e-12; // distance clip
  const double damping = 1.0;    // damping factor
  const bool verbose = false;     // print the progress
  
  // BASIC SIZES
  const int N = measures.n_elem;             // number of input measures
  const int P = init_support.n_cols;         // dimension
  const int m = init_support.n_rows;         // number of median support points
  
  // TARGET WEIGHTS v (uniform free-support by default)
  arma::vec v = arma::ones(m) / static_cast<double>(m);
  
  // STATE
  arma::mat Z_old = init_support;            // (m x P)
  arma::mat Z_new(m, P, fill::zeros);
  
  // scratch containers per measure
  arma::field<arma::mat> dist2(N);           // cost matrices C_n (m x m_n)
  arma::field<arma::mat> plans(N);           // optimal plans Γ_n (m x m_n)
  arma::vec distances(N, fill::zeros);       // d_n = W2(ν, μ_n)
  arma::vec irls_w(N, fill::zeros);          // w_n = π_n / max(d_n, δ)
  arma::vec tilde_w(N, fill::zeros);         // normalized IRLS weights
  
  // INITIAL: compute costs, plans, distances, objective
  for (int n = 0; n < N; ++n) {
    dist2(n) = aux_freemedian_sqdist(Z_old, measures(n));                 // m x m_n
    //plans(n) = util_plan_emd_R(v, marginals(n), dist2(n));                // m x m_n
    plans(n) = util_plan_emd_C(v, marginals(n), dist2(n));                // m x m_n
    distances(n) = std::sqrt(arma::accu(plans(n) % dist2(n)));            // W2 via plan & cost
  }
  double old_cost = arma::dot(weights, distances);                        // Φ(ν)
  
  arma::vec cost_hist(maxiter + 1, fill::zeros);
  cost_hist(0) = old_cost;
  
  // MAIN LOOP
  double wsum = 0.0;

  for (int it = 0; it < maxiter; ++it) {
    // (1) IRLS WEIGHTS: w_n = π_n / max(d_n, δ), normalize to sum 1
    wsum = 0.0;
    for (int n = 0; n < N; ++n) {
      irls_w(n) = weights(n)/std::max(distances(n), irls_eps);
      wsum += irls_w(n);
    }
    tilde_w = irls_w / wsum;   // sum_n tilde_w(n) = 1
    
    // (Optional "capture" rule): if any d_n < δ, give it all the weight.
    // This prevents divide-by-zero and follows Weiszfeld capture.
    for (int n = 0; n < N; ++n) {
      if (distances(n) <= irls_eps) {
        tilde_w.zeros();
        tilde_w(n) = 1.0;
        break;
      }
    }
    
    // (2) WEISZFELD UPDATE VIA BARYCENTRIC PROJECTIONS
    //    Z_tmp(i,:) = sum_n tilde_w(n) * (1/v_i) * sum_j Γ_{ij}^{(n)} x_{n,j}
    Z_new.zeros();
    for (int i = 0; i < m; ++i) {
      // accumulate the weighted barycentric projections across n
      arma::rowvec Zi_new(P, fill::zeros);
      
      for (int n = 0; n < N; ++n) {
        const arma::mat& Xn = measures(n);    // (m_n x P)
        const arma::mat& Gn = plans(n);       // (m x m_n)
        
        // barycentric projection T_n(z_i): (1/v_i) * sum_j Γ_{ij} x_{n,j}
        arma::rowvec Tn_i = (Gn.row(i) * Xn) / v(i);  // (1 x P)
        
        Zi_new += tilde_w(n) * Tn_i;
      }
      
      // optional damping: z^{+} = (1-α) z + α * Zi_new
      Z_new.row(i) = (1.0 - damping) * Z_old.row(i) + damping * Zi_new;
    }
    
    // (3) RE-SOLVE PLANS FROM UPDATED SUPPORT, RECOMPUTE DISTANCES & COST
    for (int n = 0; n < N; ++n) {
      dist2(n)    = aux_freemedian_sqdist(Z_new, measures(n));            // m x m_n
      //plans(n)    = util_plan_emd_R(v, marginals(n), dist2(n));           // m x m_n
      plans(n)    = util_plan_emd_C(v, marginals(n), dist2(n));           // m x m_n
      distances(n)= std::sqrt(arma::accu(plans(n) % dist2(n)));
    }
    double new_cost = arma::dot(weights, distances);
    
    // (4) ACCEPT/REJECT + STOPPING
    if (new_cost >= old_cost) {
      // no descent => reject and stop
      if (verbose) Rcpp::Rcout << "iteration " << it+1 << " rejected (no descent). "
                               << "cost=" << new_cost << " >= " << old_cost << "\n";
      break;
    }
    
    cost_hist(it + 1) = new_cost;
    
    // relative improvement
    double rel = std::abs(old_cost - new_cost) / std::max(old_cost, 1e-16);
    Z_old = Z_new;
    old_cost = new_cost;
    
    if (verbose) {
      Rcpp::Rcout << "iteration " << it+1
                  << " complete; cost=" << old_cost
                  << ", rel.decr=" << rel << "\n";
    }
    
    if (rel < abstol) break;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("support") = Z_old,
    Rcpp::Named("cost_history") = cost_hist
  );
}