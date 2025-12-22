# T4transport 0.1.7

* Removed dependence on CVXR package.

# T4transport 0.1.6

* API references restructured.
* Added a number of functions for histograms, including distance, barycenter, interpolation, and median.
* Added `img2measure()` to extract a discrete measure from a gray-scale image matrix.
* New family of routines for Gromov-Wasserstein geometry, `gwdist()` and `gwbary()`, added.

# T4transport 0.1.5

* Replaced `lpSolve` package's EMD solver with Nicolas Bonneel's highly-optimized [C library](https://github.com/nbonneel/network_simplex/).
* Added functions for Procrustes-Wasserstein geometry: `pwdist()` and `pwbary()`.

# T4transport 0.1.4

* `rbarygd()` for a Riemannian gradient descent to compute the free-support barycenter is added.

# T4transport 0.1.3

* Changed the structure of the package. 
* Applied the `log-sum-exp` trick for numerical stability to the entropic regularization solver.
* `rbary23L()` added for simple free-support Wasserstein barycenter computation.
* Two ad hoc median routines for images and histograms are removed.

# T4transport 0.1.2

* `gaussmed1d()` and `gaussmedpd()` now fully respects the product manifold perspective.
* `gwdist()` added for Sliced-Wasserstein distance computation.

# T4transport 0.1.1

* Support for special types of data. See the reference for more details.

# T4transport 0.1.0

* Added a `NEWS.md` file to track changes to the package.
* Initial release.
