#include "RcppArmadillo.h"
#include "elementary.h"
#include <cmath>

//#ifdef _OPENMP
//#include <omp.h>
//#endif

using namespace Rcpp;
using namespace arma;
using namespace std;

/* BARYCENTER OF IMAGES
 * (1) image_barysinkhorn14 : modification of "cpp_barysinkhorn14"
 */

// (1) image_barysinkhorn14 ====================================================
// [[Rcpp::export]]
arma::vec image_barysinkhorn14(arma::mat& dxy, arma::field<arma::vec>& marginals, arma::vec weights,
                               double p, double lambda, int maxiter, double abstol, bool printer, arma::vec initvec, int nthread){
  // PARAMETERS
  int K = weights.n_elem;
  int N = dxy.n_rows;
  
  // DATA TRANSFORMATION
  arma::mat cxy = arma::pow(dxy, p);
  arma::mat exK = arma::exp(-lambda*cxy);
  if (arma::all(arma::vectorise(exK) < 1e-15)){
    Rcpp::stop("* image C14 : regularization parameter 'lambda' is too small. Please use a larger number.");
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
  int k;
  
  // MAIN ITERATION
  for (int it=0; it<maxiter; it++){
    beta  = (static_cast<double>(it)+2.0)/2.0;
    aiter = (1.0-(1.0/beta))*a_hat_old + (1.0/beta)*a_til;
    
    agrad.fill(0.0);
    for (int k=0; k<K; k++){
      agrad += weights(k)*cpp_subgrad_weight(aiter, marginals(k), cxy, lambda);
    }

//#ifdef _OPENMP
//#pragma omp parallel num_threads(nthread) default(none) shared(K, aiter, marginals, cxy, lambda) //private(k)
//{
//#pragma omp for
//  for (k=0; k<K; k++){
//    agradmat.col(k) = weights(k)*cpp_subgrad_weight(aiter, marginals(k), cxy, lambda);
//  }
//#pragma omp single
//  agrad = arma::sum(agradmat, 1);
//} 
//#else
//    agrad.fill(0.0);
//    for (k=0; k<K; k++){
//      agrad += weights(k)*cpp_subgrad_weight(aiter, marginals(k), cxy, lambda);
//    }
//#endif
    
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
      Rcpp::Rcout << "* image14C : iteration " << it+1 << "/" << maxiter << " complete; absolute increment=" << a_hat_inc << std::endl;
    }
  }
  
  // RETURN
  return(a_hat_old);
}
