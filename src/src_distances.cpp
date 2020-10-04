#include "RcppArmadillo.h"
#include <cmath>

using namespace Rcpp;
using namespace arma;
using namespace std;

// extra function
double cpp_cost(arma::mat C, arma::mat P){
  arma::mat transposed = arma::trans(C)*P;
  return(arma::as_scalar(arma::accu(transposed.diag())));
}

/* CPP CODE FOR DISTANCE COMPUTATIONS
 * (1) cpp_sinkhorn13 : basic method by Cuturi's (2013)
 */

// (1) cpp_sinkhorn13 : basic method by Cuturi's (2013) ========================
//     code is taken mostly from author's website https://marcocuturi.net/SI.html
// [[Rcpp::export]]
Rcpp::List cpp_sinkhorn13(arma::vec a, arma::vec b, arma::mat dab, double lambda, double p, int maxiter, double abstol){
  // Initialize
  int M = dab.n_rows;
  int N = dab.n_cols;
  double MM = static_cast<double>(M);
  double NN = static_cast<double>(N);
  
  arma::mat costm(M,N,fill::zeros);
  for (int m=0; m<M; m++){
    for (int n=0; n<N; n++){
      costm(m,n) = std::pow(dab(m,n), p);
    }
  }
  
  arma::mat G = arma::exp(-costm/lambda);
  
  if (arma::all(arma::vectorise(G) < 1e-20)){
    Rcpp::stop("* sinkhorn : regularization parameter 'lambda' is too small. Please use a larger number.");
  }
  
  arma::vec vold(N,fill::zeros); vold.fill(1.0/NN);
  arma::vec vnew(N,fill::zeros); vnew.fill(1.0/NN);
  arma::vec uold(M,fill::zeros); uold.fill(1.0/MM);
  arma::vec unew(M,fill::zeros);
  
  arma::vec updater_inc(2,fill::zeros);
  arma::mat plan_old = arma::diagmat(uold)*G*arma::diagmat(vold);
  arma::mat plan_new(M,N,fill::zeros);

  unsigned int it = 0;
  for (it=0; it<maxiter; it++){
    unew = a/(G*vold);
    vnew = b/(G.t()*unew);
    if (unew.has_nan() || unew.has_inf()){
      Rcpp::stop("* sinkhorn : regularization parameter 'lambda' is too small. Please use a larger number.");  
    }
    if (vnew.has_nan() || vnew.has_inf()){
      Rcpp::stop("* sinkhorn : regularization parameter 'lambda' is too small. Please use a larger number.");  
    }
    
    plan_new = arma::diagmat(unew)*G*arma::diagmat(vnew);
    if (arma::approx_equal(plan_old, plan_new, "absdiff", abstol)){
      break;
    } 
    
    plan_old = plan_new;
    uold     = unew;
    vold     = vnew;
  }
  
  arma::mat myplan = arma::diagmat(uold)*G*arma::diagmat(vold);
  double    mydist = std::pow(cpp_cost(myplan, costm), 1.0/p);
  
  Rcpp::List output;
  output["distance"] = mydist;
  output["iteration"] = it;
  output["plan"] = myplan;
  return(output);
}

// (2) cpp_ipot20     : inexact proximal by Xie  (2020) ========================
// [[Rcpp::export]]
Rcpp::List cpp_ipot20(arma::vec a, arma::vec b, arma::mat dab, double lambda, double p, int maxiter, double abstol, int L){
  // Initialize
  int M = dab.n_rows;
  int N = dab.n_cols;
  double MM = static_cast<double>(M);
  double NN = static_cast<double>(N);
  
  arma::mat costm = arma::pow(dab, p);
  arma::mat G     = arma::exp(-costm/lambda);
  if (arma::all(arma::vectorise(G) < 1e-20)){
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
      vec_a = a/(Q*vec_b);
      vec_b = b/(Q.t()*vec_a);
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
    
    if ((it>1)&&(plan_inc<abstol)){
      break;
    }
  }
  
  double mydist = std::pow(cpp_cost(plan_old, costm), 1.0/p);
  
  Rcpp::List output;
  output["distance"] = mydist;
  output["iteration"] = it;
  output["plan"] = plan_old;
  return(output);
}