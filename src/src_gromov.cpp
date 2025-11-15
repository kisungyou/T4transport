#include "RcppArmadillo.h"
#include <cmath>
#include "utility.h"
#include "utility_EMD.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

/*
 * aux_unreg_cost : GW^2 value for the unregularized GW problem
 * aux_unreg_grad : compute the gradient given the current iterate
 * 
 * cpp_gwdist : main computation for the distance
 *  - cpp_gwdist_mm
 *  - cpp_gwdist_pg
 *  - cpp_gwdist_fw 
 *  - cpp_gwdist_**_plan : just return all the plans
 * cpp_gwbary : main computation for the barycenter
 */

double aux_unreg_cost(const arma::mat& Dx, 
                      const arma::mat& Dy, 
                      const arma::mat& P){
  int nx = Dx.n_rows;
  int ny = Dy.n_rows;
  
  double cost = 0.0;
  double diff = 0.0;
  for (int i1=0; i1<nx; i1++){
    for (int j1=0; j1<ny; j1++){
      for (int i2=0; i2<nx; i2++){
        for (int j2=0; j2<ny; j2++){
          diff = Dx(i1,i2) - Dy(j1,j2);
          cost += diff * diff * P(i1,j1) * P(i2,j2);
        }
      }
    }
  }
  return(cost);
}

arma::mat aux_unreg_grad(const arma::mat& Dx,
                         const arma::mat& Dy,
                         const arma::mat& P, 
                         const arma::mat& a, 
                         const arma::mat& b){
  // n = |X|, m = |Y|
  const arma::uword n = Dx.n_rows;
  const arma::uword m = Dy.n_rows;

  // Ones vectors
  arma::vec ones_x = arma::ones<arma::vec>(n);  // 1_n
  arma::vec ones_y = arma::ones<arma::vec>(m);  // 1_m
  
  // Ensure a, b are column vectors
  arma::vec avec = arma::conv_to<arma::vec>::from(a);  // n x 1
  arma::vec bvec = arma::conv_to<arma::vec>::from(b);  // m x 1
  
  // Elementwise squares of distance matrices: Dx^(2), Dy^(2)
  arma::mat Dx2 = Dx % Dx;  // n x n
  arma::mat Dy2 = Dy % Dy;  // m x m
  
  // First term: Dx^(2) * a * 1_m^T  -> n x m
  arma::mat term1 = Dx2 * avec * ones_y.t();
  
  // Second term: 1_n * (Dy^(2) * b)^T -> n x m
  arma::vec Dy2_b = Dy2 * bvec;       // m x 1
  arma::mat term2 = ones_x * Dy2_b.t(); // n x m
  
  // Cross term: Dx * P * Dy^T -> n x m
  arma::mat cross = Dx * P * Dy.t();
  
  // Gradient (up to a constant factor; this is the usual form)
  arma::mat grad = term1 + term2 - 2.0 * cross;
  
  return grad;
}

// GW DISTANCE ==============================================
Rcpp::List cpp_gwdist_mm(const arma::mat& Dx, 
                         const arma::mat& Dy, 
                         const arma::vec& p, 
                         const arma::vec& q, 
                         int maxiter, 
                         double abstol){
  // parameters
  int nx = Dx.n_rows;
  int ny = Dy.n_rows;
  
  // initialization
  arma::mat old_P = p*q.t();
  arma::mat new_P(nx, ny, fill::zeros);
  
  double old_cost = aux_unreg_cost(Dx, Dy, old_P);
  double new_cost = 321.0;
  double inc_cost = 123.0;
  double inc_P = 100.0;
  
  arma::mat now_grad(nx,ny,fill::zeros);
  arma::mat now_plan(nx,ny,fill::zeros);
  
  // iteration
  for (int it=0; it<maxiter; it++){
    // compute the gradient
    now_grad.reset();
    now_grad = aux_unreg_grad(Dx, Dy, old_P, p, q);
    now_grad -= now_grad.min();
    now_grad += 1e-12;
    
    // solve the linear OT subproblem 
    now_plan.reset();
    now_plan = util_plan_emd_C(p, q, now_grad);
    
    // mm : just use the EMD solution directly
    new_P = now_plan;
    new_cost = aux_unreg_cost(Dx, Dy, new_P);
    
    // updater
    if ((it > 0)&&(new_cost > old_cost)){
      break;
    }
    inc_cost = std::abs(new_cost - old_cost);
    inc_P = arma::norm(new_P - old_P, "fro");
    
    old_P = new_P;
    old_cost = new_cost;
    if ((it > 0)&&((inc_cost < abstol)||(inc_P < abstol))){
      break;
    }
  }
  
  // return
  return(Rcpp::List::create(
      Rcpp::Named("est_plan") = old_P,
      Rcpp::Named("est_cost") = old_cost));
}

arma::mat cpp_gwdist_mm_plan(const arma::mat& Dx, 
                         const arma::mat& Dy, 
                         const arma::vec& p, 
                         const arma::vec& q, 
                         int maxiter, 
                         double abstol){
  // parameters
  int nx = Dx.n_rows;
  int ny = Dy.n_rows;
  
  // initialization
  arma::mat old_P = p*q.t();
  arma::mat new_P(nx, ny, fill::zeros);
  
  double old_cost = aux_unreg_cost(Dx, Dy, old_P);
  double new_cost = 321.0;
  double inc_cost = 123.0;
  double inc_P = 100.0;
  
  arma::mat now_grad(nx,ny,fill::zeros);
  arma::mat now_plan(nx,ny,fill::zeros);
  
  // iteration
  for (int it=0; it<maxiter; it++){
    // compute the gradient
    now_grad.reset();
    now_grad = aux_unreg_grad(Dx, Dy, old_P, p, q);
    now_grad -= now_grad.min();
    now_grad += 1e-12;
    
    // solve the linear OT subproblem 
    now_plan.reset();
    now_plan = util_plan_emd_C(p, q, now_grad);
    
    // mm : just use the EMD solution directly
    new_P = now_plan;
    new_cost = aux_unreg_cost(Dx, Dy, new_P);
    
    // updater
    if ((it > 0)&&(new_cost > old_cost)){
      break;
    }
    inc_cost = std::abs(new_cost - old_cost);
    inc_P = arma::norm(new_P - old_P, "fro");
    
    old_P = new_P;
    old_cost = new_cost;
    if ((it > 0)&&((inc_cost < abstol)||(inc_P < abstol))){
      break;
    }
  }
  
  // return
  return(old_P);
}

Rcpp::List cpp_gwdist_fw(const arma::mat& Dx, 
                         const arma::mat& Dy, 
                         const arma::vec& p, 
                         const arma::vec& q, 
                         int maxiter, 
                         double abstol){
  int nx = Dx.n_rows;
  int ny = Dy.n_rows;
  
  // initialization: outer product p q^T
  arma::mat old_P = p * q.t();
  arma::mat new_P(nx, ny, fill::zeros);
  
  double old_cost = aux_unreg_cost(Dx, Dy, old_P);
  double new_cost = 321.0;
  double inc_cost = 123.0;
  double inc_P    = 100.0;
  
  arma::mat now_grad(nx, ny, fill::zeros);
  arma::mat ot_cost(nx, ny, fill::zeros);
  arma::mat now_plan(nx, ny, fill::zeros);
  arma::mat direction(nx, ny, fill::zeros);
  
  for (int it = 0; it < maxiter; it++){
    // 1. gradient at current plan
    now_grad.reset();
    now_grad = aux_unreg_grad(Dx, Dy, old_P, p, q);
    
    // 2. linear OT subproblem: Q^k = argmin_{Q} <grad, Q>
    ot_cost = now_grad;
    double cmin = ot_cost.min();
    ot_cost -= cmin;       // shift to be nonnegative
    ot_cost += 1e-12;
    
    now_plan.reset();
    now_plan = util_plan_emd_C(p, q, ot_cost);
    
    // 3. FW update: P^{k+1} = (1 - gamma_k) P^k + gamma_k Q^k
    double gamma = 2.0 / static_cast<double>(it + 2); // gamma_0 = 1
    direction = now_plan - old_P;
    new_P = old_P + gamma * direction;
    
    // 4. evaluate cost
    new_cost = aux_unreg_cost(Dx, Dy, new_P);
    
    // 5. stopping checks (same style as MM)
    if ((it > 0) && (new_cost > old_cost)){
      // cost increased; revert and stop
      break;
    }
    inc_cost = std::abs(new_cost - old_cost);
    inc_P    = arma::norm(new_P - old_P, "fro");
    
    old_P   = new_P;
    old_cost = new_cost;
    
    if ((it > 0) && ((inc_cost < abstol) || (inc_P < abstol))){
      break;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("est_plan") = old_P,
    Rcpp::Named("est_cost") = old_cost
  );
}

arma::mat cpp_gwdist_fw_plan(const arma::mat& Dx, 
                         const arma::mat& Dy, 
                         const arma::vec& p, 
                         const arma::vec& q, 
                         int maxiter, 
                         double abstol){
  int nx = Dx.n_rows;
  int ny = Dy.n_rows;
  
  // initialization: outer product p q^T
  arma::mat old_P = p * q.t();
  arma::mat new_P(nx, ny, fill::zeros);
  
  double old_cost = aux_unreg_cost(Dx, Dy, old_P);
  double new_cost = 321.0;
  double inc_cost = 123.0;
  double inc_P    = 100.0;
  
  arma::mat now_grad(nx, ny, fill::zeros);
  arma::mat ot_cost(nx, ny, fill::zeros);
  arma::mat now_plan(nx, ny, fill::zeros);
  arma::mat direction(nx, ny, fill::zeros);
  
  for (int it = 0; it < maxiter; it++){
    // 1. gradient at current plan
    now_grad.reset();
    now_grad = aux_unreg_grad(Dx, Dy, old_P, p, q);
    
    // 2. linear OT subproblem: Q^k = argmin_{Q} <grad, Q>
    ot_cost = now_grad;
    double cmin = ot_cost.min();
    ot_cost -= cmin;       // shift to be nonnegative
    ot_cost += 1e-12;
    
    now_plan.reset();
    now_plan = util_plan_emd_C(p, q, ot_cost);
    
    // 3. FW update: P^{k+1} = (1 - gamma_k) P^k + gamma_k Q^k
    double gamma = 2.0 / static_cast<double>(it + 2); // gamma_0 = 1
    direction = now_plan - old_P;
    new_P = old_P + gamma * direction;
    
    // 4. evaluate cost
    new_cost = aux_unreg_cost(Dx, Dy, new_P);
    
    // 5. stopping checks (same style as MM)
    if ((it > 0) && (new_cost > old_cost)){
      // cost increased; revert and stop
      break;
    }
    inc_cost = std::abs(new_cost - old_cost);
    inc_P    = arma::norm(new_P - old_P, "fro");
    
    old_P   = new_P;
    old_cost = new_cost;
    
    if ((it > 0) && ((inc_cost < abstol) || (inc_P < abstol))){
      break;
    }
  }
  
  return(old_P);
}

Rcpp::List cpp_gwdist_pg(const arma::mat& Dx, 
                         const arma::mat& Dy, 
                         const arma::vec& p, 
                         const arma::vec& q, 
                         int maxiter, 
                         double abstol){
  int nx = Dx.n_rows;
  int ny = Dy.n_rows;
  
  // initialization
  arma::mat old_P = p * q.t();
  arma::mat new_P(nx, ny, fill::zeros);
  
  double old_cost = aux_unreg_cost(Dx, Dy, old_P);
  double new_cost = 321.0;
  double inc_cost = 123.0;
  double inc_P    = 100.0;
  
  arma::mat now_grad(nx, ny, fill::zeros);
  arma::mat R(nx, ny, fill::zeros);
  arma::mat proj_cost(nx, ny, fill::zeros);
  arma::mat now_plan(nx, ny, fill::zeros);
  
  for (int it = 0; it < maxiter; it++){
    // 1. gradient at current plan
    now_grad.reset();
    now_grad = aux_unreg_grad(Dx, Dy, old_P, p, q);
    
    // 2. step size schedule (simple diminishing)
    double stepsize = 1.0 / std::sqrt(static_cast<double>(it) + 1.0);
    
    // 3. unconstrained gradient step
    R = old_P - stepsize * now_grad;
    
    // 4. approximate projection via OT:
    //    P^{k+1} = argmax_Q <R, Q> = argmin_Q <-R, Q>
    proj_cost = -R;
    double cmin = proj_cost.min();
    proj_cost -= cmin;
    proj_cost += 1e-12;
    
    now_plan.reset();
    now_plan = util_plan_emd_C(p, q, proj_cost);
    
    new_P   = now_plan;
    new_cost = aux_unreg_cost(Dx, Dy, new_P);
    
    // 5. stopping checks (same style as MM)
    if ((it > 0) && (new_cost > old_cost)){
      // cost increased; revert and stop
      break;
    }
    inc_cost = std::abs(new_cost - old_cost);
    inc_P    = arma::norm(new_P - old_P, "fro");
    
    old_P   = new_P;
    old_cost = new_cost;
    
    if ((it > 0) && ((inc_cost < abstol) || (inc_P < abstol))){
      break;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("est_plan") = old_P,
    Rcpp::Named("est_cost") = old_cost
  );
}
arma::mat cpp_gwdist_pg_plan(const arma::mat& Dx, 
                         const arma::mat& Dy, 
                         const arma::vec& p, 
                         const arma::vec& q, 
                         int maxiter, 
                         double abstol){
  int nx = Dx.n_rows;
  int ny = Dy.n_rows;
  
  // initialization
  arma::mat old_P = p * q.t();
  arma::mat new_P(nx, ny, fill::zeros);
  
  double old_cost = aux_unreg_cost(Dx, Dy, old_P);
  double new_cost = 321.0;
  double inc_cost = 123.0;
  double inc_P    = 100.0;
  
  arma::mat now_grad(nx, ny, fill::zeros);
  arma::mat R(nx, ny, fill::zeros);
  arma::mat proj_cost(nx, ny, fill::zeros);
  arma::mat now_plan(nx, ny, fill::zeros);
  
  for (int it = 0; it < maxiter; it++){
    // 1. gradient at current plan
    now_grad.reset();
    now_grad = aux_unreg_grad(Dx, Dy, old_P, p, q);
    
    // 2. step size schedule (simple diminishing)
    double stepsize = 1.0 / std::sqrt(static_cast<double>(it) + 1.0);
    
    // 3. unconstrained gradient step
    R = old_P - stepsize * now_grad;
    
    // 4. approximate projection via OT:
    //    P^{k+1} = argmax_Q <R, Q> = argmin_Q <-R, Q>
    proj_cost = -R;
    double cmin = proj_cost.min();
    proj_cost -= cmin;
    proj_cost += 1e-12;
    
    now_plan.reset();
    now_plan = util_plan_emd_C(p, q, proj_cost);
    
    new_P   = now_plan;
    new_cost = aux_unreg_cost(Dx, Dy, new_P);
    
    // 5. stopping checks (same style as MM)
    if ((it > 0) && (new_cost > old_cost)){
      // cost increased; revert and stop
      break;
    }
    inc_cost = std::abs(new_cost - old_cost);
    inc_P    = arma::norm(new_P - old_P, "fro");
    
    old_P   = new_P;
    old_cost = new_cost;
    
    if ((it > 0) && ((inc_cost < abstol) || (inc_P < abstol))){
      break;
    }
  }
  
  return(old_P);
}



// [[Rcpp::export]]
Rcpp::List cpp_gwdist(const arma::mat& Dx, 
                      const arma::mat& Dy, 
                      const arma::vec& p, 
                      const arma::vec& q, 
                      int maxiter, 
                      double abstol, 
                      std::string method){
  if (method=="mm"){
    return(cpp_gwdist_mm(Dx, Dy, p, q, maxiter, abstol));
  } else if (method=="fw"){
    return(cpp_gwdist_fw(Dx, Dy, p, q, maxiter, abstol));
  } else if (method=="pg"){
    return(cpp_gwdist_pg(Dx, Dy, p, q, maxiter, abstol));
  }{
    Rcpp::stop("Unknown method for Gromov-Wasserstein distance computation: " + method);
  }
}

arma::mat cpp_gwdist_plan(const arma::mat& Dx, 
                      const arma::mat& Dy, 
                      const arma::vec& p, 
                      const arma::vec& q, 
                      int maxiter, 
                      double abstol, 
                      std::string method){
  if (method=="mm"){
    return(cpp_gwdist_mm_plan(Dx, Dy, p, q, maxiter, abstol));
  } else if (method=="fw"){
    return(cpp_gwdist_fw_plan(Dx, Dy, p, q, maxiter, abstol));
  } else if (method=="pg"){
    return(cpp_gwdist_pg_plan(Dx, Dy, p, q, maxiter, abstol));
  }{
    Rcpp::stop("Unknown method for Gromov-Wasserstein distance computation: " + method);
  }
}

// [[Rcpp::export]]
arma::mat cpp_gwbary(const arma::field<arma::mat>& Ds,
                     const arma::field<arma::vec>& marginals, // a^(k)
                     const arma::vec& weights,                // λ_k
                     const arma::mat& init_D,                 // initial barycenter distance
                     int maxiter,
                     double abstol,
                     std::string method){
  // PARAMETERS
  int N = Ds.n_elem;                 // number of input measures
  int num_support = init_D.n_rows;   // barycenter support size m
  
  arma::vec weight_bary = arma::ones<arma::vec>(num_support)
    / static_cast<double>(num_support); // b (uniform)
  
  // sanity checks
  if (weights.n_elem != static_cast<arma::uword>(N)){
    Rcpp::stop("Length of 'weights' must match number of measures.");
  }
  for (int n = 0; n < N; n++){
    const arma::mat& Dk = Ds(n);
    if (Dk.n_rows != Dk.n_cols){
      Rcpp::stop("Ds(" + std::to_string(n) + ") is not square.");
    }
    const arma::vec& ak = marginals(n);
    if (ak.n_elem != Dk.n_rows){
      Rcpp::stop("marginals(" + std::to_string(n) + ") length does not match Ds("
                   + std::to_string(n) + ").");
    }
  }
  if (init_D.n_rows != init_D.n_cols){
    Rcpp::stop("init_D must be a square matrix.");
  }
  
  // precompute denominator matrix B = b b^T (constant over k and iterations)
  arma::mat B = weight_bary * weight_bary.t(); // m x m
  
  // INITIALIZE
  arma::field<arma::mat> couplings(N);       // P^{(k)} matrices, each m x n_k
  arma::mat old_D = init_D;                  // current barycenter distance
  arma::mat new_D(num_support, num_support); // updated barycenter distance
  
  // OUTER ITERATION
  for (int it = 0; it < maxiter; it++){
    // ------------------------------------
    // Step 1. update couplings P^{(k)}
    // ------------------------------------
    for (int n = 0; n < N; n++){
      const arma::mat& Dk = Ds(n);         // n_k x n_k
      const arma::vec& ak = marginals(n);  // n_k x 1
      
      couplings(n) = cpp_gwdist_plan(old_D, Dk,
                weight_bary, ak,
                maxiter, abstol, method);
    }
    
    // ------------------------------------
    // Step 2. update barycenter distance D_Y
    // D_Y(i,i') = sum_k λ_k [P_k D_k P_k^T]_{i,i'} / (b_i b_{i'})
    // ------------------------------------
    arma::mat num(num_support, num_support, arma::fill::zeros);
    
    for (int n = 0; n < N; n++){
      const arma::mat& Dk = Ds(n);        // n_k x n_k
      const arma::mat& Pk = couplings(n); // m   x n_k
      double wk = weights(n);             // λ_k
      
      arma::mat Mk = Pk * Dk * Pk.t();    // m x m
      num += wk * Mk;
    }
    
    // elementwise update of new_D using denominator B = b b^T
    new_D.zeros();
    for (int i = 0; i < num_support; i++){
      for (int j = 0; j < num_support; j++){
        double Bij = B(i,j);  // = b_i b_j
        if (Bij > 0.0){
          new_D(i,j) = num(i,j) / Bij;
        } else {
          // shouldn't happen if b has positive entries,
          // but keep fallback
          new_D(i,j) = old_D(i,j);
        }
      }
    }
    
    // enforce symmetry and zero diagonal
    new_D = 0.5 * (new_D + new_D.t());
    new_D.diag().zeros();
    
    // ------------------------------------
    // Step 3. check outer convergence
    // ------------------------------------
    double diff = arma::norm(new_D - old_D, "fro");
    old_D = new_D;
    
    if ((it > 0) && (diff < abstol)){
      break;
    }
  }
  
  return old_D;
}