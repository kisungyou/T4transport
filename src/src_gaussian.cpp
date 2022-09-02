#include "RcppArmadillo.h"
#include <cmath>

using namespace Rcpp;
using namespace arma;
using namespace std;


// C++ ROUTINES FOR GAUSSIAN DISTRIBUTIONS
// gauss_wassdist   : SPD matrices and their Wasserstein distance of order 2
// gauss_weiszfeld  : weiszfeld algorithm
// gauss_spdbary16A : wasserstein barycenter for SPDs by Alvarez-Esteban
// gauss_spdmed22Y  : my proposal for Wasserstein Median


// gauss_wassdist --------------------------------------------------------------
double gauss_wassdist(arma::mat x, arma::mat y){
  arma::mat xhalf = arma::sqrtmat_sympd(x);
  double d2 = arma::trace(x) + arma::trace(y) - 2.0*arma::trace(arma::sqrtmat_sympd(xhalf*y*xhalf));
  return(std::sqrt(d2));
}

// gauss_weiszfeld -------------------------------------------------------------
// [[Rcpp::export]]
arma::rowvec gauss_weiszfeld(arma::mat &X, arma::vec &weights, double abstol, int maxiter){
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  
  // prepare
  arma::rowvec xold(P,fill::zeros);
  for (int n=0; n<N; n++){
    xold += X.row(n)*weights(n);
  }
  arma::rowvec xtmp(P,fill::zeros);
  arma::rowvec xnew(P,fill::zeros);
  arma::vec dists(N,fill::zeros);
  double xtmp2 = 0.0;
  double xinc  = 0.0;
  double norm2 = 0.0;
  
  // iteration
  double epsnum = 100*arma::datum::eps;
  for (int it=0;it<maxiter;it++){
    // step 1. compute distance
    for (int n=0;n<N;n++){
      norm2 = arma::norm(X.row(n)-xold, 2);
      if (norm2 < epsnum){
        dists(n) = norm2 + epsnum;
      } else {
        dists(n) = norm2;
      }
    }
    // step 2. compute numerator and denominator
    xtmp.fill(0.0);
    xtmp2 = 0.0;
    for (int n=0;n<N;n++){
      xtmp  += weights(n)*X.row(n)/dists(n);
      xtmp2 += weights(n)/dists(n);
    }
    xnew = xtmp/xtmp2;
    
    // step 3. updating information
    xinc = arma::norm(xold-xnew,2);
    xold = xnew;
    if (xinc < abstol){
      break;
    }
  }
  
  // return
  return(xold);
}

// gauss_spdbary16A ------------------------------------------------------------
// [[Rcpp::export]]
arma::mat gauss_spdbary16A(arma::cube array3d, arma::vec weight, double abstol, int maxiter){
  // PREPARE
  int p = array3d.n_rows;
  int N = array3d.n_slices;
  
  double S_inc = 10000.0;
  arma::mat S_old = arma::mean(array3d, 2);
  int S_old_rank = arma::rank(S_old);
  if (S_old_rank < p){
    S_old.fill(0.0);
    for (int n=0; n<N; n++){
      S_old = arma::logmat_sympd(array3d.slice(n))/static_cast<double>(N);
    }
    S_old = arma::expmat_sym(S_old);
  }
  
  arma::mat S_tmp(p,p,fill::zeros);
  arma::mat S_new(p,p,fill::zeros);
  
  arma::mat S0_half(p,p,fill::zeros);
  arma::mat S0_hinv(p,p,fill::zeros);
  
  arma::vec S_weight = weight/arma::accu(weight);
  
  // MAIN
  int itcount = 0;
  for (int it=0; it<maxiter; it++){
    // MAIN-1. initialize 'tmp' and compute auxiliary ones
    itcount += 1;
    S_tmp.fill(0.0);
    S0_half = arma::sqrtmat_sympd(S_old);
    S0_hinv = arma::inv_sympd(S0_half);
    
    // MAIN-2. compute the intermediate summation
    for (int n=0; n<N; n++){
      S_tmp += S_weight(n)*arma::sqrtmat_sympd(S0_half*array3d.slice(n)*S0_half);
    }
    
    // MAIN-3. finally attain an iterate
    S_new = S0_hinv*S_tmp*S_tmp*S0_hinv;
    
    // MAIN-3. update
    S_inc = arma::norm(S_old-S_new,"fro");
    S_old = S_new;
    if (S_inc < abstol){
      break;
    }
  }
  
  // RETURN
  return(S_old);
}


// gauss_spdmed22Y -------------------------------------------------------------
// [[Rcpp::export]]
arma::mat gauss_spdmed22Y(arma::cube &array3d, arma::vec &weight, double abstol, int maxiter){
  // parameters
  int p = array3d.n_rows;
  int N = array3d.n_slices;
  
  // initialize
  arma::mat cov_old = gauss_spdbary16A(array3d, weight, 1e-10, 100);
  arma::mat cov_new(p,p,fill::zeros);
  
  // quantities
  double nearness = 100.0*arma::datum::eps;
  double cov_inc  = 10.0;
  double tmp_dist = 0.0;
  arma::vec tmp_weight(N,fill::zeros);
  
  // iterate
  for (int it=0; it<maxiter; it++){
    // it-1. update the relative weights
    tmp_weight.fill(0.0);
    for (int n=0; n<N; n++){
      tmp_dist = gauss_wassdist(array3d.slice(n), cov_old);
      if (tmp_dist < nearness){
        return(array3d.slice(n));
      }
      tmp_weight(n) = weight(n)/tmp_dist;
    }
    
    // it-2. normalize the weight
    tmp_weight /= arma::accu(tmp_weight);
    
    // it-3. update the SPD matrix
    cov_new = gauss_spdbary16A(array3d, tmp_weight, 1e-10, 100);
    
    // it-4. error & update
    cov_inc = arma::norm(cov_new-cov_old, "fro");
    cov_old = cov_new;
    cov_new.fill(0.0);
    if (cov_inc < abstol){
      break;
    }
  }
  
  // return
  return(cov_old);
}