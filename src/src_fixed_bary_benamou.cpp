#include "RcppArmadillo.h"
#include "elementary.h"
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;



// cpp_barybregman15  : Benamou et al. (2015) ==================================
// [[Rcpp::export]]
arma::vec cpp_barybregman15(arma::field<arma::mat>& listdXY, arma::field<arma::vec>& marginals, arma::vec weights,
                            double p, double lambda, int maxiter, double abstol, bool printer, arma::vec initvec){
  // PARAMETERS
  int K = weights.n_elem;
  int N = listdXY(0).n_rows; // total number of atoms

  // DATA TRANSFORMATION
  arma::field<arma::mat> list_G(K);
  for (int k=0; k<K; k++){
    list_G(k) = arma::exp(-arma::pow(listdXY(k), p)/lambda);
    if (list_G(k).min() < 1e-300){
      Rcpp::stop("* barybregman15 : regularization parameter 'lambda' is too small. Please use a larger number.");
    }
  }
  
  // INITIALIZATION
  arma::vec p_old = initvec;
  arma::vec p_new(N,fill::zeros);
  double    p_inc = 0.0;
  
  arma::field<arma::vec> vec_v(K);
  arma::field<arma::vec> vec_u(K);
  for (int k=0; k<K; k++){
    vec_v(k) = arma::ones<arma::vec>(marginals(k).n_elem)/(static_cast<double>(marginals(k).n_elem));
  }
  arma::field<arma::mat> plans_old(K);
  arma::field<arma::mat> plans_new(K);
  for (int k=0; k<K; k++){
    plans_old(k) = arma::ones<arma::mat>(N,marginals(k).n_elem);
  }
  arma::vec  plans_dist(K,fill::zeros);
  
  // MAIN ITERATION
  for (int it=0; it<maxiter; it++){
    // 1. update u
    for (int k=0; k<K; k++){
      vec_u(k) = p_old/(list_G(k)*vec_v(k));
    }
    // 2. update v
    for (int k=0; k<K; k++){
      vec_v(k) = marginals(k)/(arma::trans(list_G(k))*vec_u(k));
    }
    // 3. update p and plan
    p_new.fill(1.0);
    for (int k=0; k<K; k++){
      // stability update
      arma::vec prod = vec_u(k) % (list_G(k) * vec_v(k));
      prod = arma::clamp(prod, 1e-14, arma::datum::inf); // optional safeguard
      p_new = p_new % arma::pow(prod, weights(k));       // accumulate
      //p_new = p_new%arma::pow(vec_u(k)%(list_G(k)*vec_v(k)), weights(k));
      
      plans_new(k)  = arma::diagmat(vec_u(k))*(list_G(k))*arma::diagmat(vec_v(k));
      plans_dist(k) = arma::norm(plans_old(k)-plans_new(k),"fro");
    }
    
    p_new = arma::clamp(p_new, 0.0, arma::datum::inf);
    p_new /= arma::accu(p_new);
    if (p_new.has_nan()||arma::any(p_new < 0)||p_new.has_inf()){
      return(p_old);
    }
    
    // 4. updater
    p_inc = arma::norm(p_old-p_new,2);
    p_old = p_new;
    plans_old = plans_new;
    if ((p_inc < abstol)&&(it>0)&&(plans_dist.max()<abstol)){
      break;
    }
    if (printer){
      Rcpp::Rcout << "* barybregman15 : iteration " << it+1 << "/" << maxiter << " complete; absolute increment=" << plans_dist.max() << std::endl;
    }
  }
  
  // RETURN
  return(p_old);
  
  
  // // PARAMETERS
  // int K = weights.n_elem;
  // int N = listdXY(0).n_rows;
  // double NN = static_cast<double>(N);
  // 
  // // DATA TRANSFORMATION
  // arma::field<arma::mat> list_G(K);
  // for (int k = 0; k < K; k++) {
  //   arma::mat Ck = arma::pow(listdXY(k), p);
  //   arma::mat exponent = -Ck / lambda;
  //   exponent = arma::clamp(exponent, -700.0, 0.0); // prevent underflow
  //   list_G(k) = arma::exp(exponent);
  //   if (list_G(k).min() < 1e-300) {
  //     Rcpp::stop("* barybregman15 : lambda too small or distances too large; Gibbs kernel underflowed.");
  //   }
  // }
  // 
  // // INITIALIZATION
  // arma::vec p_old = initvec;
  // arma::vec p_new(N, fill::zeros);
  // double p_inc = 0.0;
  // 
  // arma::field<arma::vec> vec_v(K);
  // arma::field<arma::vec> vec_u(K);
  // for (int k = 0; k < K; k++) {
  //   vec_v(k) = arma::ones<arma::vec>(marginals(k).n_elem) / static_cast<double>(marginals(k).n_elem);
  // }
  // 
  // arma::field<arma::mat> plans_old(K);
  // arma::field<arma::mat> plans_new(K);
  // for (int k = 0; k < K; k++) {
  //   plans_old(k) = arma::ones<arma::mat>(N, marginals(k).n_elem);
  // }
  // 
  // arma::vec plans_dist(K, fill::zeros);
  // 
  // // MAIN ITERATION
  // for (int it = 0; it < maxiter; it++) {
  //   // 1. update u
  //   for (int k = 0; k < K; k++) {
  //     vec_u(k) = p_old / (list_G(k) * vec_v(k) + 1e-12);
  //   }
  //   
  //   // 2. update v
  //   for (int k = 0; k < K; k++) {
  //     vec_v(k) = marginals(k) / (arma::trans(list_G(k)) * vec_u(k) + 1e-12);
  //   }
  //   
  //   // 3. update p and plan
  //   p_new.fill(1.0);
  //   for (int k = 0; k < K; k++) {
  //     arma::vec prod = vec_u(k) % (list_G(k) * vec_v(k));
  //     prod = arma::clamp(prod, 1e-14, arma::datum::inf);  // stabilize
  //     p_new = p_new % arma::pow(prod, weights(k));        // multiplicative update
  //     
  //     plans_new(k) = arma::diagmat(vec_u(k)) * list_G(k) * arma::diagmat(vec_v(k));
  //     plans_dist(k) = arma::norm(plans_old(k) - plans_new(k), "fro");
  //   }
  //   
  //   p_new = arma::clamp(p_new, 1e-14, arma::datum::inf);
  //   p_new /= arma::accu(p_new);
  //   
  //   if (p_new.has_nan() || arma::any(p_new < 0) || p_new.has_inf()) {
  //     Rcpp::Rcout << "* barybregman15 : numerical error — returning last valid iterate.\n";
  //     return p_old;
  //   }
  //   
  //   // 4. updater
  //   p_inc = arma::norm(p_old - p_new, 2);
  //   p_old = p_new;
  //   plans_old = plans_new;
  //   
  //   if ((p_inc < abstol) && (it > 0) && (plans_dist.max() < abstol)) {
  //     break;
  //   }
  //   
  //   if (printer) {
  //     Rcpp::Rcout << "* barybregman15 : iteration " << it + 1 << "/" << maxiter
  //                 << " complete; max plan diff = " << plans_dist.max() << std::endl;
  //   }
  // }
  // 
  // // RETURN
  // return p_old;
}