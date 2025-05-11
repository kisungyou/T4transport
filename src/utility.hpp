#pragma once

#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>

arma::mat util_entropic_plan(const arma::mat& C, const arma::vec& a, const arma::vec& b,
                             double lambda, int maxiter, double abstol);