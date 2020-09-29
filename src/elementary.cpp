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
