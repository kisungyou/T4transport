# Wasserstein Distance via Inexact Proximal Point Method

The Inexact Proximal Point Method (IPOT) offers a computationally
efficient approach to approximating the Wasserstein distance between two
empirical measures by iteratively solving a series of regularized
optimal transport problems. This method replaces the entropic
regularization used in Sinkhorn's algorithm with a proximal formulation
that avoids the explicit use of entropy, thereby mitigating numerical
instabilities.

Let \\C := \\X_m - Y_n\\^p\\ be the cost matrix, where \\X_m\\ and
\\Y_n\\ are the support points of two discrete distributions \\\mu\\ and
\\\nu\\, respectively. The IPOT algorithm solves a sequence of
optimization problems: \$\$ \Gamma^{(t+1)} = \arg\min\_{\Gamma \in
\Pi(\mu, \nu)} \langle \Gamma, C \rangle + \lambda D(\Gamma \\
\Gamma^{(t)}), \$\$ where \\\lambda \> 0\\ is the proximal
regularization parameter and \\D(\cdot \\ \cdot)\\ is the
Kullback–Leibler divergence. Each subproblem is solved approximately
using a fixed number of inner iterations, making the method inexact.

Unlike entropic methods, IPOT does not require \\\lambda \rightarrow 0\\
for convergence to the unregularized Wasserstein solution. It is
therefore more robust to numerical precision issues, especially for
small regularization parameters, and provides a closer approximation to
the true optimal transport cost with fewer artifacts.

## Usage

``` r
ipot(X, Y, p = 2, wx = NULL, wy = NULL, lambda = 1, ...)

ipotD(D, p = 2, wx = NULL, wy = NULL, lambda = 1, ...)
```

## Arguments

- X:

  an \\(M\times P)\\ matrix of row observations.

- Y:

  an \\(N\times P)\\ matrix of row observations.

- p:

  an exponent for the order of the distance (default: 2).

- wx:

  a length-\\M\\ marginal density that sums to \\1\\. If `NULL`
  (default), uniform weight is set.

- wy:

  a length-\\N\\ marginal density that sums to \\1\\. If `NULL`
  (default), uniform weight is set.

- lambda:

  a regularization parameter (default: 0.1).

- ...:

  extra parameters including

  maxiter

  :   maximum number of iterations (default: 496).

  abstol

  :   stopping criterion for iterations (default: 1e-10).

  L

  :   small number of inner loop iterations (default: 1).

- D:

  an \\(M\times N)\\ distance matrix \\d(x_m, y_n)\\ between two sets of
  observations.

## Value

a named list containing

- distance:

  \\\mathcal{W}\_p\\ distance value

- plan:

  an \\(M\times N)\\ nonnegative matrix for the optimal transport plan.

## References

Xie Y, Wang X, Wang R, Zha H (2020-07-22/2020-07-25). “A Fast Proximal
Point Method for Computing Exact Wasserstein Distance.” In Adams RP,
Gogate V (eds.), *Proceedings of the 35th Uncertainty in Artificial
Intelligence Conference*, volume 115 of *Proceedings of Machine Learning
Research*, 433–453.

## Examples

``` r
# \donttest{
#-------------------------------------------------------------------
#  Wasserstein Distance between Samples from Two Bivariate Normal
#
# * class 1 : samples from Gaussian with mean=(-1, -1)
# * class 2 : samples from Gaussian with mean=(+1, +1)
#-------------------------------------------------------------------
## SMALL EXAMPLE
set.seed(100)
m = 20
n = 30
X = matrix(rnorm(m*2, mean=-1),ncol=2) # m obs. for X
Y = matrix(rnorm(n*2, mean=+1),ncol=2) # n obs. for Y

## COMPARE WITH WASSERSTEIN 
outw = wasserstein(X, Y)
ipt1 = ipot(X, Y, lambda=1)
ipt2 = ipot(X, Y, lambda=10)

## VISUALIZE : SHOW THE PLAN AND DISTANCE
pmw = paste0("Exact plan\n dist=",round(outw$distance,2))
pm1 = paste0("IPOT (lambda=1)\n dist=",round(ipt1$distance,2))
pm2 = paste0("IPOT (lambda=10)\n dist=",round(ipt2$distance,2))

opar <- par(no.readonly=TRUE)
par(mfrow=c(1,3), pty="s")
image(outw$plan, axes=FALSE, main=pmw)
image(ipt1$plan, axes=FALSE, main=pm1)
image(ipt2$plan, axes=FALSE, main=pm2)

par(opar)
# }
```
