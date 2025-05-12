#ifndef T4TRANSPORT_UTILITY_H
#define T4TRANSPORT_UTILITY_H

#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>

arma::mat util_entropic_plan(const arma::mat& C, const arma::vec& a, const arma::vec& b,
                             double lambda, int maxiter, double abstol);

#endif // T4TRANSPORT_UTILITY_H