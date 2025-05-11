#include "RcppArmadillo.h"
#include "elementary.h"
#include <cmath>

/*
 * ROUTINES FOR IMAGE AND HISTOGRAMS
 * (1) routine_bary14C 
 * (2) routine_bary15B
 */

// (1) routine_bary14C  ========================================================
// [[Rcpp::export]]
arma::vec routine_bary14C(arma::mat& dxy, arma::field<arma::vec>& marginals, arma::vec weights,
                          double p, double lambda, int maxiter, double abstol, bool printer, arma::vec initvec, int nthread){
  // PARAMETERS
  int K = weights.n_elem;
  int N = dxy.n_rows;
  
  // DATA TRANSFORMATION
  arma::mat cxy = arma::pow(dxy, p);
  arma::mat exK = arma::exp(-lambda*cxy);
  if (arma::all(arma::vectorise(exK) < 1e-15)){
    Rcpp::stop("* Routine for bary14C : regularization parameter 'lambda' is too small. Please use a larger number.");
  }
  
  // INITIALIZATION
  arma::vec a_hat_old = initvec;
  arma::vec a_hat_new(N,fill::zeros); a_hat_new.fill((1.0/static_cast<double>(N)));
  arma::vec a_til(N,fill::zeros); a_til.fill((1.0/static_cast<double>(N)));
  arma::vec aiter(N,fill::zeros);
  arma::vec agrad(N,fill::zeros);
  
  arma::mat agradmat(N,K,fill::zeros);
  
  double t0 = 2.0;
  double beta = 0.0;
  double a_hat_inc = 0.0;
  
  // MAIN ITERATION
  for (int it=0; it<maxiter; it++){
    beta  = (static_cast<double>(it)+2.0)/2.0;
    aiter = (1.0-(1.0/beta))*a_hat_old + (1.0/beta)*a_til;
    
    agrad.fill(0.0);
    for (int k=0; k<K; k++){
      agrad += weights(k)*cpp_subgrad_weight(aiter, marginals(k), cxy, lambda);
    }

    a_til = a_til%arma::exp(-t0*beta*agrad);
    a_til = a_til/arma::accu(a_til);
    a_hat_new = (1.0 - (1.0/beta))*a_hat_old + (1.0/beta)*a_til;
    a_hat_new = a_hat_new/arma::accu(a_hat_new);
    if (a_hat_new.has_nan()||arma::any(a_hat_new < 0)){
      return(a_hat_old);
    }
    a_hat_inc = arma::norm(a_hat_new - a_hat_old,2);
    a_hat_old = a_hat_new;
    if ((a_hat_inc < abstol)&&(it>0)){
      break;
    }
    if (printer){
      Rcpp::Rcout << "* Routine for bary14C : iteration " << it+1 << "/" << maxiter << " complete; abs. increment =" << a_hat_inc << std::endl;
    }
  }
  
  // RETURN
  return(a_hat_old);
}


// (2) routine_bary15B =========================================================
// [[Rcpp::export]]
arma::vec routine_bary15B(arma::mat& dxy, arma::field<arma::vec>& marginals, arma::vec weights,
                              double p, double lambda, int maxiter, double abstol, bool printer, arma::vec initvec, int nthread){
  // PARAMETERS
  int K = weights.n_elem;
  int N = dxy.n_rows;
  
  // DATA TRANSFORMATION
  arma::mat G = arma::exp(-arma::pow(dxy, p)/lambda);
  if (arma::all(arma::vectorise(G) < 1e-15)){
    Rcpp::stop("* Routine for bary15B : regularization parameter 'lambda' is too small. Please use a larger number.");
  }
  
  // INITIALIZATION
  arma::vec p_old = initvec;
  arma::vec p_new(N,fill::zeros);
  double    p_inc = 0.0;
  arma::mat vec_v = arma::ones<arma::mat>(N,K)/static_cast<double>(N);
  arma::mat vec_u = arma::ones<arma::mat>(N,K)/static_cast<double>(N);
  arma::cube plan_old(N,N,K,fill::ones);
  arma::cube plan_new(N,N,K,fill::zeros);
  arma::vec  plan_dist(K,fill::zeros);
  
  // MAIN ITERATION
  for (int it=0; it<maxiter; it++){
    // 1. update u
    for (int k=0; k<K; k++){
      vec_u.col(k) = p_old/(G*vec_v.col(k));
    }
    // 3. update v
    for (int k=0; k<K; k++){
      vec_v.col(k) = marginals(k)/(G.t()*vec_u.col(k));
    }
    // 2. update p & plan
    p_new.fill(1.0);
    for (int k=0; k<K; k++){
      p_new = p_new%arma::pow(vec_u.col(k)%(G*vec_v.col(k)), weights(k));
      plan_new.slice(k) = arma::diagmat(vec_u.col(k))*G*arma::diagmat(vec_v.col(k));
      plan_dist(k)      = arma::norm(plan_old.slice(k)-plan_new.slice(k),"fro");
    }
    p_new /= arma::accu(p_new);
    // 4. updater
    p_inc = arma::norm(p_old-p_new,2);
    p_old = p_new;
    plan_old = plan_new;
    if ((p_inc < abstol)&&(it>0)&&(plan_dist.max()<abstol)){
      break;
    }
    if (printer){
      Rcpp::Rcout << "* Routine for bary15B : iteration " << it+1 << "/" << maxiter << " complete; abs. increment =" << plan_dist.max() << std::endl;
    }
  }
  
  // RETURN
  return(p_old);
}
