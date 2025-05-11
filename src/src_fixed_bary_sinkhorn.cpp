#include "RcppArmadillo.h"
#include "elementary.h"
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// method by Cuturi and Doucet (2014) ==========================================
// NOTE       : codes are written in the paper's notation (1/lambda)*h(P)
// [[Rcpp::export]]
arma::vec cpp_fixed_sinkhorn14(arma::field<arma::mat>& listdXY, arma::field<arma::vec>& marginals, arma::vec weights,
                               double p, double lambda, int maxiter, double abstol, bool printer, arma::vec initvec){
  // PARAMETERS
  int K = weights.n_elem;
  int N = listdXY(0).n_rows; // total number of atoms
  double NN = static_cast<double>(N);
  
  // DATA TRANSFORMATION
  arma::field<arma::mat> listcXY(K);
  for (int k=0; k<K; k++){
    listcXY(k) = arma::pow(listdXY(k), p);
    arma::mat kernel = arma::exp(-lambda * listcXY(k));
    if (kernel.min() < 1e-300) {
      Rcpp::stop("* barysinkhorn14 : lambda is too small, causing underflow in exp(-C/lambda).");
    }
  }
  
  // INITIALIZATION
  arma::vec a_hat_old = initvec;
  arma::vec a_hat_new(N,fill::zeros); a_hat_new.fill(1.0/NN);
  arma::vec a_til(N,fill::zeros); a_til.fill(1.0/NN);
  arma::vec aiter(N,fill::zeros);
  arma::vec agrad(N,fill::zeros);
  
  double t0 = 2.0;
  double beta = 0.0;
  double a_hat_inc = 0.0;
  
  // MAIN ITERATION
  for (int it=0; it<maxiter; it++){
    beta  = (static_cast<double>(it)+2.0)/2.0;
    aiter = (1.0-(1.0/beta))*a_hat_old + (1.0/beta)*a_til;
    agrad.fill(0.0);
    for (int k=0; k<K; k++){
      agrad += weights(k)*cpp_subgrad_weight(aiter, marginals(k), listcXY(k), lambda);
    }
    //a_til = a_til%arma::exp(-t0*beta*agrad);
    arma::vec delta = arma::clamp(-t0 * beta * agrad, -100.0, 100.0);
    a_til = a_til % arma::exp(delta);
    a_til /= arma::accu(a_til);  // ensure normalization
    
    a_hat_new = (1.0 - (1.0/beta))*a_hat_old + (1.0/beta)*a_til;
    a_hat_new = arma::clamp(a_hat_new, 1e-14, arma::datum::inf);
    a_hat_new /= arma::accu(a_hat_new); // ensure normalization
    if (a_hat_new.has_nan()||arma::any(a_hat_new < 0)||a_hat_new.has_inf()){
      return(a_hat_old);
    }
    a_hat_inc = arma::norm(a_hat_new - a_hat_old,2);
    a_hat_old = a_hat_new;
    if ((a_hat_inc < abstol)&&(it>0)){
      break;
    }
    if (printer){
      Rcpp::Rcout << "* barysinkhorn14 : iteration " << it+1 << "/" << maxiter << " complete; absolute increment=" << a_hat_inc << std::endl;
    }
  }
  
  // RETURN
  return(a_hat_old);
}