#include "elementary.h"
#include <cmath>

// (01) compute_pdist2 : compute pairwise distance matrix ======================
// [[Rcpp::export]]
arma::mat compute_pdist2(arma::mat& X, arma::mat& Y){
  int M = X.n_rows;
  int N = Y.n_rows;

  arma::mat output(M,N,fill::zeros);
  for (int m=0; m<M; m++){
    for (int n=0; n<N; n++){
      output(m,n) = arma::norm(X.row(m)-Y.row(n), 2);
    }
  }
  return(output);
}
// (02) cpp_subgrad_weight =====================================================
arma::vec cpp_subgrad_weight(arma::vec a, arma::vec b, arma::mat M, double lambda){
  int n = a.n_elem; double nn = static_cast<double>(n);
  
  arma::vec a_clamped = arma::clamp(a, 1e-12, 1.0);  // avoids division by 0
  arma::mat K = arma::exp(-lambda*M);
  arma::mat Ktil = arma::diagmat(1.0/a_clamped)*K;
  
  arma::vec uold(n,fill::zeros); uold.fill(1.0/nn);
  arma::vec unew(n,fill::zeros); unew.fill(1.0/nn);
  double uinc = 10000.0;
  
  int par_iter   = 100;  // arbitrarily set by KY
  double par_tol = 1e-4; 
  
  for (int it=0; it<par_iter; it++){
    uold = 1.0/(Ktil*(b/(K.t()*uold)));
    if (it%10 == 0){
      unew = 1.0/(Ktil*(b/(K.t()*uold)));
      uinc = arma::norm(uold-unew,2);
      uold = unew;
      if (uinc < par_tol){
        break;
      }
    }
  }
  arma::vec vold = b/(K.t()*uold);
  arma::vec alpha = ((1/lambda)*(arma::log(uold))) - ((arma::accu(arma::log(uold)))/(lambda*nn));
  return(alpha);
}
arma::mat cpp_subgrad_plan(arma::vec a, arma::vec b, arma::mat M, double lambda){
  int n = a.n_elem; double nn = static_cast<double>(n);
  arma::mat K = arma::exp(-lambda*M);
  arma::mat Ktil = arma::diagmat(1.0/a)*K;
  
  arma::vec uold(n,fill::zeros); uold.fill(1.0/nn);
  arma::vec unew(n,fill::zeros); unew.fill(1.0/nn);
  double uinc = 10000.0;
  
  int par_iter   = 100;  // arbitrarily set by KY
  double par_tol = 1e-4; 
  
  for (int it=0; it<par_iter; it++){
    uold = 1.0/(Ktil*(b/(K.t()*uold)));
    if (it%10 == 0){
      unew = 1.0/(Ktil*(b/(K.t()*uold)));
      uinc = arma::norm(uold-unew,2);
      uold = unew;
      if (uinc < par_tol){
        break;
      }
    }
  }
  arma::vec vold = b/(K.t()*uold);
  arma::vec alpha = ((1/lambda)*(arma::log(uold))) - ((arma::accu(arma::log(uold)))/(lambda*nn));
  arma::mat theplan = arma::diagmat(uold)*K*arma::diagmat(vold);
  return(theplan);
}
arma::field<arma::mat> cpp_subgrad_both(arma::vec a, arma::vec b, arma::mat M, double lambda){
  int n = a.n_elem; double nn = static_cast<double>(n);
  arma::mat K = arma::exp(-lambda*M);
  arma::mat Ktil = arma::diagmat(1.0/a)*K;
  
  arma::vec uold(n,fill::zeros); uold.fill(1.0/nn);
  arma::vec unew(n,fill::zeros); unew.fill(1.0/nn);
  double uinc = 10000.0;
  
  int par_iter   = 100;  // arbitrarily set by KY
  double par_tol = 1e-4; 
  
  for (int it=0; it<par_iter; it++){
    uold = 1.0/(Ktil*(b/(K.t()*uold)));
    if (it%10 == 0){
      unew = 1.0/(Ktil*(b/(K.t()*uold)));
      uinc = arma::norm(uold-unew,2);
      uold = unew;
      if (uinc < par_tol){
        break;
      }
    }
  }
  arma::vec vold = b/(K.t()*uold);
  arma::vec alpha = ((1/lambda)*(arma::log(uold))) - ((arma::accu(arma::log(uold)))/(lambda*nn));
  arma::field<arma::mat> output(2);
  output(0) = arma::reshape(alpha, n, 1);
  output(1) = arma::diagmat(uold)*K*arma::diagmat(vold);
  return(output);
}


arma::mat cpp_sinkhorn_getmap(arma::mat c, arma::mat p, arma::mat q, double lambda, int maxiter, double abstol){
  // parameters
  int m = p.n_elem;
  int n = q.n_elem;
  
  // prepare
  // Gibbs kernel
  arma::mat K  = arma::exp(-c/lambda);
  arma::mat Kt = arma::trans(K); 
  
  // updaters
  arma::vec a_old(m,fill::ones);
  arma::vec a_new(m,fill::zeros);
  arma::vec b_old(n,fill::zeros);
  arma::vec b_new(n,fill::zeros);
  
  // stopping criterion
  double a_inc = 100.0;
  double b_inc = 100.0;
  
  // iteration
  for (int it=0; it<maxiter; it++){
    // update b & a once
    b_new = q/(Kt*a_old);
    a_new = p/(K*b_new);
    
    // stopping
    a_inc = arma::norm(a_old-a_new,2);
    b_inc = arma::norm(b_old-b_new,2);
    a_old = a_new;
    b_old = b_new;
    
    // stop if small
    if ((a_inc < abstol)&&(b_inc < abstol)){
      break;
    }
  }
  
  // compute the estimated coupling
  arma::mat output = arma::diagmat(a_old)*K*arma::diagmat(b_old);
  return(output);
}
 
// (04) cpp_mvrnorm ============================================================
// [[Rcpp::export]]
arma::mat cpp_mvrnorm(int n, const arma::vec& mu, const arma::mat& cov){
  arma::mat X = arma::mvnrnd(mu, cov, n);
  return(arma::trans(X));
}
