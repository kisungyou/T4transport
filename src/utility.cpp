// [[Rcpp::depends(RcppArmadillo, BH)]]

#include <RcppArmadillo.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/successive_shortest_path_nonnegative_weights.hpp>
#include <boost/graph/find_flow_cost.hpp>
#include <limits>
#include <cmath>
#include "utility.h"

/* UTILITY FUNCTIONS
 * util_plan_entropic : compute the entropic-regularized plan using Sinkhorn algorithm
 * util_plan_emd_BH   : compute the exact EMD plan using Boost
 * util_plan_emd_R    : compute the exact EMD plan using lpSolve
 * 
 * util_mvrnorm       : generate multivariate normal random variables with given mean and covariance
 */



using namespace arma;
using namespace Rcpp;
using namespace std;


// util_plan_entropic ==========================================================
static double logsumexp(const arma::vec& x) {
  double xmax = x.max();
  return xmax + std::log(arma::sum(arma::exp(x - xmax)));
}

arma::mat util_plan_entropic(const arma::vec& a, const arma::vec& b, const arma::mat& C,
                             double lambda, int maxiter, double abstol){
  int M = C.n_rows;
  int N = C.n_cols;
  
  arma::vec f(M, arma::fill::zeros);
  arma::vec g(N, arma::fill::zeros);
  arma::vec logKsum(std::max(M, N));
  
  arma::vec marg_u(M,fill::zeros);
  arma::vec marg_v(N,fill::zeros);
  
  for (int it = 0; it < maxiter; it++) {
    // Update f
    for (int m = 0; m < M; m++) {
      arma::vec row = (-C.row(m).t() + g) / lambda;
      logKsum(m) = logsumexp(row);
    }
    f = lambda * (arma::log(a) - logKsum.subvec(0, M - 1));
    
    // Update g
    for (int n = 0; n < N; n++) {
      arma::vec col = (-C.col(n) + f) / lambda;
      logKsum(n) = logsumexp(col);
    }
    g = lambda * (arma::log(b) - logKsum.subvec(0, N - 1));
    
    // Marginal check
    marg_u.zeros();
    marg_v.zeros();
    for (int m = 0; m < M; m++) {
      for (int n = 0; n < N; n++) {
        double val = std::exp((f(m) + g(n) - C(m, n)) / lambda);
        marg_u(m) += val;
        marg_v(n) += val;
      }
    }
    
    double err = arma::norm(marg_u - a, 1) + arma::norm(marg_v - b, 1);
    if (err < abstol) {
      break;
    }
  }
  
  // Final plan
  arma::mat plan(M, N, arma::fill::zeros);
  for (int m = 0; m < M; m++) {
    for (int n = 0; n < N; n++) {
      plan(m, n) = std::exp((f(m) + g(n) - C(m, n)) / lambda);
    }
  }
  
  return plan;
}

// util_plan_emd_BH =========================================================
// Type definitions for Boost graph
typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS>::vertex_descriptor Vertex;
typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS>::edge_descriptor Edge;

typedef boost::adjacency_list<
  boost::vecS, boost::vecS, boost::directedS,
  boost::no_property,
  boost::property<boost::edge_capacity_t, long,
                  boost::property<boost::edge_residual_capacity_t, long,
                                  boost::property<boost::edge_reverse_t, Edge,
                                                  boost::property<boost::edge_weight_t, long> > > >
> Graph;

typedef boost::property_map<Graph, boost::edge_capacity_t>::type CapacityMap;
typedef boost::property_map<Graph, boost::edge_weight_t>::type WeightMap;
typedef boost::property_map<Graph, boost::edge_reverse_t>::type ReverseMap;
typedef boost::property_map<Graph, boost::edge_residual_capacity_t>::type ResidualMap;

struct FlowEdge {
  int from, to;
  Edge edge;
};

// Add edge helper
void addEdge(int from, int to, long capacity, long cost,
             Graph &g, CapacityMap &cap, WeightMap &weight,
             ReverseMap &rev, std::vector<FlowEdge> *transportEdges = nullptr) {
  Edge e, rev_e;
  bool success;
  boost::tie(e, success) = boost::add_edge(from, to, g);
  boost::tie(rev_e, success) = boost::add_edge(to, from, g);
  cap[e] = capacity;
  weight[e] = cost;
  cap[rev_e] = 0;
  weight[rev_e] = -cost;
  rev[e] = rev_e;
  rev[rev_e] = e;
  
  if (transportEdges && from < to && from >= 0 && to >= 0)
    transportEdges->push_back({from, to, e});
}

arma::mat util_plan_emd_BH(const arma::vec& a, const arma::vec& b, const arma::mat& C) {
  int n = a.n_elem, m = b.n_elem;
  int N = n + m + 2, src = n + m, sink = n + m + 1;
  long scale = 1e6;
  
  // === Runtime Checks ===
  if (C.n_rows != (size_t)n || C.n_cols != (size_t)m) {
    stop("Dimension mismatch: cost matrix C must be n x m where length(a) = n and length(b) = m.");
  }
  if (!a.is_finite() || !b.is_finite() || !C.is_finite()) {
    stop("Non-finite values detected in input: 'a', 'b', or 'C'.");
  }
  if (arma::any(a < 0) || arma::any(b < 0) || arma::any(vectorise(C) < 0)) {
    stop("Negative values are not allowed in 'a', 'b', or 'C'.");
  }
  double sa = arma::accu(a), sb = arma::accu(b);
  if (std::abs(sa - sb) > 1e-8) {
    stop("Total mass mismatch: sum(a) != sum(b).");
  }
  if (arma::max(vectorise(C)) * scale > static_cast<double>(std::numeric_limits<long>::max())) {
    stop("Scaled cost matrix has values too large for 64-bit integer. Try reducing the scale factor.");
  }
  
  // === Graph Construction ===
  Graph g(N);
  CapacityMap cap = boost::get(boost::edge_capacity, g);
  WeightMap weight = boost::get(boost::edge_weight, g);
  ReverseMap rev = boost::get(boost::edge_reverse, g);
  ResidualMap rescap = boost::get(boost::edge_residual_capacity, g);
  std::vector<FlowEdge> flowEdges;
  
  for (int i = 0; i < n; ++i)
    addEdge(src, i, static_cast<long>(a(i) * scale), 0, g, cap, weight, rev);
  
  for (int j = 0; j < m; ++j)
    addEdge(n + j, sink, static_cast<long>(b(j) * scale), 0, g, cap, weight, rev);
  
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j)
      addEdge(i, n + j, std::numeric_limits<long>::max(),
              static_cast<long>(C(i, j) * scale),
              g, cap, weight, rev, &flowEdges);
  
  // === Solve min-cost flow ===
  boost::successive_shortest_path_nonnegative_weights(g, src, sink);
  
  // === Extract transport plan ===
  arma::mat P(n, m, arma::fill::zeros);
  for (const auto& fe : flowEdges) {
    long f = cap[fe.edge] - rescap[fe.edge];
    if (f > 0)
      P(fe.from, fe.to - n) = static_cast<double>(f) / scale;
  }
  
  return P;
}
// util_plan_emd_R ==========================================================
arma::mat util_plan_emd_R(const arma::vec& a, const arma::vec& b, const arma::mat& C){
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("T4transport");
  Rcpp::Function aux_emd = pkg["aux_emd"];
  //Function aux_emd("aux_emd");
  SEXP result = aux_emd(wrap(a), wrap(b), wrap(C));
  return(Rcpp::as<arma::mat>(result));
}

// util_mvrnorm  ============================================================
// [[Rcpp::export]]
arma::mat util_mvrnorm(const arma::vec& par_mean, 
                       const arma::mat& par_cov,
                       int num_samples){
  arma::mat X = arma::trans(arma::mvnrnd(par_mean, par_cov, num_samples));
  return(X);
}

// util_pairwise_sqdist =====================================================
arma::mat util_pairwise_sqdist(const arma::mat& X, const arma::mat& Y){
  int N = X.n_rows;
  int M = Y.n_rows;
  arma::mat D(N,M,fill::zeros);
  
  for (int n=0; n<N; n++){
    for (int m=0; m<M; m++){
      D(n,m) = arma::accu(arma::square(X.row(n)-Y.row(m)));
    }
  }
  
  return(D);
}