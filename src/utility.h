#ifndef T4TRANSPORT_UTILITY_H
#define T4TRANSPORT_UTILITY_H

#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>

arma::mat util_plan_entropic(const arma::mat& C, const arma::vec& a, const arma::vec& b,
                             double lambda, int maxiter, double abstol);
arma::mat util_plan_emd_BH(const arma::mat& C, const arma::vec& a, const arma::vec& b);
arma::mat util_plan_emd_R(const arma::mat& C, const arma::vec& a, const arma::vec& b);

#endif // T4TRANSPORT_UTILITY_H