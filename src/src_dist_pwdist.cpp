#include "RcppArmadillo.h"
#include <cmath>
#include "utility.h"
#include "utility_EMD.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

bool all_finite(const arma::mat& A){
  return A.is_finite();
}
void clip_nonneg_inplace(arma::mat& A){
  A.elem(arma::find_nonfinite(A)).zeros();
  A.elem(arma::find(A < 0.0)).zeros();
}
double safe_cost(const arma::mat& C, const arma::mat& G){
  arma::mat Gc = G;
  clip_nonneg_inplace(Gc);
  // C is squared distances, should be >= 0; clip just in case
  arma::mat Cc = C;
  Cc.elem(arma::find_nonfinite(Cc)).zeros();
  Cc.elem(arma::find(Cc < 0.0)).zeros();
  double val = arma::accu(Cc % Gc);
  return (val > 0.0) ? val : 0.0;
}
void normalize_masses(arma::vec& w){
  w.elem(arma::find_nonfinite(w)).zeros();
  w.elem(arma::find(w < 0.0)).zeros();
  double s = arma::accu(w);
  if (s > 0.0) w /= s;
}

// Center to mean and scale so max row-norm is 1
arma::mat normalize_cloud(const arma::mat& X) {
  arma::mat Y = X;
  arma::rowvec c = arma::mean(Y, 0);
  Y.each_row() -= c;
  arma::vec r = arma::sqrt(arma::sum(arma::square(Y), 1));
  double maxr = r.max();
  if (maxr > 0.0) Y /= maxr;
  return Y;
}

// Pairwise squared Euclidean distances (dense)
arma::mat pairwise_sqdist(const arma::mat& X) {
  arma::vec nrm = arma::sum(arma::square(X), 1);
  arma::mat D = arma::repmat(nrm, 1, X.n_rows);
  D += arma::repmat(nrm.t(), X.n_rows, 1);
  D -= 2.0 * (X * X.t());
  D.elem(arma::find(D < 0)).zeros();
  return D;
}

// Median of strictly upper triangular entries (off-diagonal)
double median_offdiag(const arma::mat& D) {
  const arma::uword n = D.n_rows;
  if (n < 2) return 0.0;
  arma::vec v; v.set_size(n * (n - 1) / 2);
  arma::uword p = 0;
  for (arma::uword i = 0; i < n; ++i)
    for (arma::uword j = i + 1; j < n; ++j)
      v(p++) = D(i, j);
  return arma::median(v);
}

// Fully connected RBF weights with median-bandwidth; zero diagonal
arma::mat rbf_weights(const arma::mat& X) {
  arma::mat D2 = pairwise_sqdist(X);
  arma::mat D  = arma::sqrt(D2);
  double sigma = median_offdiag(D);
  if (!(sigma > 0.0) || !std::isfinite(sigma)) sigma = 1.0;
  const double inv2s2 = 1.0 / (2.0 * sigma * sigma);
  
  arma::mat W = arma::exp(-D2 * inv2s2);
  W.diag().zeros(); // no self-loops
  // Enforce symmetry (numerical safety)
  W = 0.5 * (W + W.t());
  return W;
}

// Unnormalized Laplacian L = D - W (dense)
arma::mat laplacian_from_weights(const arma::mat& W) {
  arma::vec deg = arma::sum(W, 1);      // row sums
  arma::mat L = arma::diagmat(deg) - W; // dense Laplacian
  return L;
}

// Returns 2nd smallest eigenvector of L, normalized to [0,1]
arma::vec fiedler_rbf_vector(const arma::mat& X) {
  const arma::uword n = X.n_rows;
  if (n < 2) return arma::zeros<arma::vec>(n);
  
  arma::mat W = rbf_weights(X);
  arma::mat L = laplacian_from_weights(W);
  
  arma::vec eval;
  arma::mat evec;
  // Dense symmetric eigen-decomposition (sorted ascending)
  bool ok = arma::eig_sym(eval, evec, L, "dc");
  if (!ok || eval.n_elem < 2) return arma::zeros<arma::vec>(n);
  
  arma::vec vf = evec.col(1); // 2nd smallest eigenvector
  double vmin = vf.min(), vmax = vf.max();
  if (!(vmax > vmin)) return arma::zeros<arma::vec>(n);
  return (vf - vmin) / (vmax - vmin);
}


struct OT1DResult {
  double cost;
  arma::mat Gamma; // dense n x m plan (mostly zeros but dense by request)
};


// 1D OT between fx (n) and fy (m) with cost |x-y| (uniform masses)
OT1DResult ot1d_uniform_dense(const arma::vec& fx, const arma::vec& fy) {
  const arma::uword n = fx.n_elem, m = fy.n_elem;
  const double ax = 1.0 / static_cast<double>(n);
  const double by = 1.0 / static_cast<double>(m);
  
  arma::uvec ix = arma::sort_index(fx); // ascending
  arma::uvec iy = arma::sort_index(fy);
  
  arma::mat G(n, m, arma::fill::zeros);
  
  arma::uword i = 0, j = 0;
  double rx = ax, ry = by;
  double cost = 0.0;
  
  while (i < n && j < m) {
    double flow = std::min(rx, ry);
    arma::uword ii = ix(i), jj = iy(j);
    G(ii, jj) += flow;
    cost += flow * std::abs(fx(ii) - fy(jj));
    rx -= flow; ry -= flow;
    if (rx <= 1e-15) { ++i; rx = ax; }
    if (ry <= 1e-15) { ++j; ry = by; }
  }
  OT1DResult res; res.cost = cost; res.Gamma = std::move(G);
  return res;
}

// Procrustes: P = U V^T where SVD(Y^T Γ^T X)
arma::mat procrustes_from_plan_dense(const arma::mat& X,
                                     const arma::mat& Y,
                                     const arma::mat& Gamma) {
  arma::mat M = Y.t() * Gamma.t() * X; // d x d
  arma::mat U, V;
  arma::vec s;
  arma::svd(U, s, V, M);
  return U * V.t();
}

// Returns orthogonal P (d x d) for PW alignment
arma::mat fiedler_rbf_initializer(const arma::mat& X_raw,
                                  const arma::mat& Y_raw) {
  arma::mat X = normalize_cloud(X_raw);
  arma::mat Y = normalize_cloud(Y_raw);
  
  // Fiedler coordinates (normalized to [0,1])
  arma::vec fx = fiedler_rbf_vector(X);
  arma::vec fy = fiedler_rbf_vector(Y);
  
  // Try both signs for fy (eigenvectors are sign-ambiguous)
  OT1DResult r_pos = ot1d_uniform_dense(fx,  fy);
  OT1DResult r_neg = ot1d_uniform_dense(fx, -fy);
  const OT1DResult& best = (r_pos.cost <= r_neg.cost) ? r_pos : r_neg;
  
  // Orthogonal initializer
  arma::mat P = procrustes_from_plan_dense(X, Y, best.Gamma);
  
  // Safety: if SVD failed (rare), fall back to identity
  if (P.n_rows != X.n_cols || P.n_cols != X.n_cols) {
    return arma::eye<arma::mat>(X.n_cols, X.n_cols);
  }
  return P;
}


// [[Rcpp::export]]
Rcpp::List cpp_pwdist(const arma::mat& X, const arma::vec& p, 
                      const arma::mat& Y, const arma::vec& q, 
                      int maxiter, 
                      double abstol){
  // parameters
  int N = X.n_rows;
  int M = Y.n_rows;
  int d = Y.n_cols;
  
  // initial correspondence
  arma::mat init_P = fiedler_rbf_initializer(X, Y);
  arma::mat cross_sqdist = util_pairwise_sqdist(X, Y*init_P);
  cross_sqdist.elem(arma::find_nonfinite(cross_sqdist)).zeros();
  cross_sqdist.elem(arma::find(cross_sqdist < 0.0)).zeros(); // paranoia guard
  //arma::mat iter_Gamma = util_plan_emd_R(p, q, cross_sqdist);
  arma::mat iter_Gamma = util_plan_emd_C(p, q, cross_sqdist);
  clip_nonneg_inplace(iter_Gamma);
  arma::mat iter_P(d,d,fill::zeros);
  arma::mat iter_YP(M,d,fill::zeros);
  double old_cost = safe_cost(cross_sqdist, iter_Gamma);
  double new_cost = 0.0;
  double inc_cost = 0.0;
  
  // iteration 
  for (int it=0; it<maxiter; it++){
    // step 1. procrustes update for a current plan
    arma::mat U;
    arma::vec s;
    arma::mat V;
    arma::svd(U,s,V,arma::trans(Y)*iter_Gamma.t()*X);
    iter_P.reset();
    iter_P = U*V.t();
    
    // step 2. ot update
    iter_YP = Y*iter_P;
    for (int n=0; n<N; n++){
      for (int m=0; m<M; m++){
        cross_sqdist(n,m) = arma::accu(arma::square(X.row(n)-iter_YP.row(m)));
      }
    }
    cross_sqdist.elem(arma::find_nonfinite(cross_sqdist)).zeros();
    cross_sqdist.elem(arma::find(cross_sqdist < 0.0)).zeros(); // paranoia guard
    //iter_Gamma = util_plan_emd_R(p, q, cross_sqdist);
    iter_Gamma = util_plan_emd_C(p, q, cross_sqdist);
    clip_nonneg_inplace(iter_Gamma);
    
    // cost
    new_cost = safe_cost(cross_sqdist, iter_Gamma);
    if ((it > 0)&&(new_cost > old_cost)){
      break;
    }
    inc_cost = std::abs(new_cost-old_cost);
    old_cost = new_cost;
    if ((it > 0)&&(inc_cost < abstol)){
      break;
    }
  }
  
  // return the matrices
  return(Rcpp::List::create(
      Rcpp::Named("est_plan") = iter_Gamma,
      Rcpp::Named("est_P") = iter_P,
      Rcpp::Named("est_cost") = old_cost));
}