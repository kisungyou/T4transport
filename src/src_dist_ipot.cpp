#include "RcppArmadillo.h"
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

double cpp_cost(arma::mat C, arma::mat P){
  arma::mat transposed = arma::trans(C)*P;
  return(arma::as_scalar(arma::accu(transposed.diag())));
}
// [[Rcpp::export]]
Rcpp::List cpp_ipot20(const arma::vec a, const arma::vec b, const arma::mat& dab, 
                      double lambda, double p, int maxiter, double abstol, int L){
  // Initialize
  int M = dab.n_rows;
  int N = dab.n_cols;
  double MM = static_cast<double>(M);
  double NN = static_cast<double>(N);
  
  arma::mat costm = arma::pow(dab, p);
  arma::mat G     = arma::exp(-costm/lambda);
  if (G.min() < 1e-20){
    Rcpp::stop("* ipot : regularization parameter 'lambda' is too small. Please use a larger number.");
  }
  
  arma::vec vec_a(M,fill::zeros); vec_a.fill(1.0/MM);
  arma::vec vec_b(N,fill::zeros); vec_b.fill(1.0/NN);
  
  arma::mat plan_old(M,N,fill::ones);
  arma::mat plan_new(M,N,fill::zeros);
  double    plan_inc = 0.0;
  arma::mat Q(M,N,fill::zeros);
  
  int it = 0;
  for (it=0; it<maxiter; it++){
    Q = G%plan_old;
    for (int l=0; l<L; l++){
      //vec_a = a/(Q*vec_b);
      //vec_b = b/(Q.t()*vec_a);
      
      // added support for overflow prevention
      vec_a = a / (Q * vec_b + 1e-12);
      vec_b = b / (Q.t() * vec_a + 1e-12);
    }
    if (vec_a.has_nan() || vec_a.has_inf()){
      Rcpp::stop("* ipot : regularization parameter 'lambda' is too small. Please use a larger number.");  
    }
    if (vec_b.has_nan() || vec_b.has_inf()){
      Rcpp::stop("* ipot : regularization parameter 'lambda' is too small. Please use a larger number.");  
    }
    
    plan_new = arma::diagmat(vec_a)*Q*arma::diagmat(vec_b);
    plan_inc = arma::norm(plan_old-plan_new,"fro");
    plan_old = plan_new;
    
    if ((it > 1) && (plan_inc / std::sqrt(static_cast<double>(M * N)) < abstol)) {
      break;
    }
  }
  
  // compute the total cost
  double cost = arma::accu(plan_old % costm);
  double dist = std::pow(cost, 1.0/p);
  
  // return
  return(Rcpp::List::create(
      Rcpp::Named("distance") = dist,
      Rcpp::Named("plan") = plan_old));
}