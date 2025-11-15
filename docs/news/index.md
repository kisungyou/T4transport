# Changelog

## T4transport 0.1.6

- Added
  [`img2measure()`](https://www.kisungyou.com/T4transport/reference/img2measure.md)
  to extract a discrete measure from a gray-scale image matrix.

## T4transport 0.1.5

CRAN release: 2025-11-13

- Replaced `lpSolve` package’s EMD solver with Nicolas Bonneel’s
  highly-optimized [C
  library](https://github.com/nbonneel/network_simplex/).
- Added functions for Procrustes-Wasserstein geometry;
  [`pwdist()`](https://www.kisungyou.com/T4transport/reference/pwdist.md)
  and
  [`pwbary()`](https://www.kisungyou.com/T4transport/reference/pwbary.md).

## T4transport 0.1.4

- `rbarygd()` for a Riemannian gradient descent to compute the
  free-support barycenter is added.

## T4transport 0.1.3

CRAN release: 2025-05-29

- Changed the structure of the package.
- Applied the `log-sum-exp` trick for numerical stability.
- [`rbary23L()`](https://www.kisungyou.com/T4transport/reference/rbary23L.md)
  added for simple free-support Wasserstein barycenter computation.
- Two ad hoc median routines for images and histograms are removed.

## T4transport 0.1.2

CRAN release: 2023-04-11

- [`gaussmed1d()`](https://www.kisungyou.com/T4transport/reference/gaussmed1d.md)
  and
  [`gaussmedpd()`](https://www.kisungyou.com/T4transport/reference/gaussmedpd.md)
  now fully respects the product manifold perspective.
- [`gwdist()`](https://www.kisungyou.com/T4transport/reference/gwdist.md)
  added for Sliced-Wasserstein distance computation.

## T4transport 0.1.1

CRAN release: 2022-09-03

- Support for special types of data. See the reference for more details.

## T4transport 0.1.0

CRAN release: 2020-10-09

- Added a `NEWS.md` file to track changes to the package.
- Initial release.
