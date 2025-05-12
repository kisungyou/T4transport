#ifndef T4TRANSPORT_UTILITY_H
#define T4TRANSPORT_UTILITY_H

#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>

arma::mat util_plan_entropic(const arma::vec& a, const arma::vec& b, const arma::mat& C,
                             double lambda, int maxiter, double abstol);
arma::mat util_plan_emd_BH(const arma::vec& a, const arma::vec& b, const arma::mat& C);
arma::mat util_plan_emd_R(const arma::vec& a, const arma::vec& b, const arma::mat& C);

#endif // T4TRANSPORT_UTILITY_H