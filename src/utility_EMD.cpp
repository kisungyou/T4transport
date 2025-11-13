// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "network_simplex_simple.h"
#include "full_bipartitegraph.h"
#include "utility_EMD.h"

using namespace Rcpp;
using namespace arma;
using namespace lemon;

// We only need Arc; don't pull in DIGRAPH_TYPEDEFS to avoid unused typedef warnings.
typedef FullBipartiteDigraph Digraph;
typedef Digraph::Arc Arc;

// -----------------------------------------------------------------------------
// Exact EMD using network_simplex (balanced OT)
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
arma::mat util_plan_emd_C(const arma::vec &a,
                          const arma::vec &b,
                          const arma::mat &C) {
  const int n1 = static_cast<int>(a.n_elem);
  const int n2 = static_cast<int>(b.n_elem);
  
  // Basic dimension checks
  if (static_cast<int>(C.n_rows) != n1 || static_cast<int>(C.n_cols) != n2) {
    Rcpp::stop("util_plan_emd_C: dimension mismatch between a, b, and C.");
  }
  
  const double sum_a = arma::accu(a);
  const double sum_b = arma::accu(b);
  if (std::abs(sum_a - sum_b) > 1e-12) {
    Rcpp::stop("util_plan_emd_C: a and b must have the same total mass.");
  }
  
  // ---------------------------------------------------------------------------
  // Compress supports: keep only indices where weights > 0
  // ---------------------------------------------------------------------------
  std::vector<int>    indI;
  std::vector<int>    indJ;
  std::vector<double> w1;  // supplies (positive)
  std::vector<double> w2;  // demands (negative, as expected by supplyMap)
  
  indI.reserve(n1);
  indJ.reserve(n2);
  w1.reserve(n1);
  w2.reserve(n2);
  
  for (int i = 0; i < n1; ++i) {
    double val = a(i);
    if (val > 0.0) {
      indI.push_back(i);
      w1.push_back(val);
    }
  }
  for (int j = 0; j < n2; ++j) {
    double val = b(j);
    if (val > 0.0) {
      indJ.push_back(j);
      w2.push_back(-val);  // demands as negative supply
    }
  }
  
  const int n = static_cast<int>(indI.size());
  const int m = static_cast<int>(indJ.size());
  
  if (n == 0 || m == 0) {
    // Trivial case: at least one measure is all zeros
    return arma::mat(n1, n2, arma::fill::zeros);
  }
  
  // ---------------------------------------------------------------------------
  // Build bipartite graph and configure network simplex
  // ---------------------------------------------------------------------------
  Digraph di(n, m);  // n left (supply) nodes, m right (demand) nodes
  
  // maxIter: you can tune this if desired
  const std::uint64_t maxIter = 100000000;
  
  NetworkSimplexSimple<Digraph, double, double> net(
      di,
      true,                 // arc_mixing (as in original code)
      n + m,                // number of nodes
      n * m,                // number of arcs
      maxIter);
  
  // Set supplies and demands
  net.supplyMap(w1.data(), n, w2.data(), m);
  
  // Set costs on each arc: arcs are ordered lexicographically over (i,j)
  std::size_t arcId = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      const double cij = C(indI[i], indJ[j]);
      if (cij < 0.0) {
        Rcpp::stop("util_plan_emd_C: negative costs are not allowed.");
      }
      net.setCost(di.arcFromId(static_cast<int>(arcId)), cij);
      ++arcId;
    }
  }
  
  // ---------------------------------------------------------------------------
  // Solve
  // ---------------------------------------------------------------------------
  NetworkSimplexSimple<Digraph,double,double>::ProblemType status = net.run();
  if (status == NetworkSimplexSimple<Digraph,double,double>::UNBOUNDED) {
    Rcpp::stop("util_plan_emd_C: network simplex unbounded (check cost matrix C).");
  }
  
  // Original vendor code does not do this:
  // // Accept both OPTIMAL and MAX_ITER_REACHED, like POT’s EMD_wrap
  // if (ret != (int)net.OPTIMAL && ret != (int)net.MAX_ITER_REACHED) {
  //   Rcpp::stop("util_plan_emd_C: network simplex failed (status %d).", ret);
  // }
  // Here is a little control for infeasible and unbounded solutions:
  // if (ret == NetworkSimplexSimple<Digraph,double,double>::INFEASIBLE) {
  //   Rcpp::stop("util_plan_emd_C: network simplex infeasible.");
  // }
  // if (ret == NetworkSimplexSimple<Digraph,double,double>::UNBOUNDED) {
  //   Rcpp::stop("util_plan_emd_C: network simplex unbounded.");
  // }
  
  // ---------------------------------------------------------------------------
  // Extract optimal flow and assemble full transport plan
  // ---------------------------------------------------------------------------
  arma::mat G(n1, n2, arma::fill::zeros);
  
  arcId = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      Arc a_edge = di.arcFromId(static_cast<int>(arcId));
      ++arcId;
      
      double flow = net.flow(a_edge);
      if (flow <= 0.0) {
        continue;
      }
      
      const int I = indI[i];  // original row index
      const int J = indJ[j];  // original col index
      G(I, J) = flow;
    }
  }
  
  return G;
}
