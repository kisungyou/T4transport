#ifndef T4TRANSPORT_UTILITY_EMD_H
#define T4TRANSPORT_UTILITY_EMD_H

#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>

arma::mat util_plan_emd_C(const arma::vec &a,
                          const arma::vec &b,
                          const arma::mat &C);
Rcpp::List util_dual_emd_C(const arma::vec &a,
                           const arma::vec &b,
                           const arma::mat &C,
                           const bool return_plan);

#endif // T4TRANSPORT_UTILITY_EMD_H
