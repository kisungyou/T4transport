# Package index

## \[1\] Core OT for Empirical Measures

### Distances

- [`ipot()`](https://www.kisungyou.com/T4transport/reference/ipot.md)
  [`ipotD()`](https://www.kisungyou.com/T4transport/reference/ipot.md) :
  Wasserstein Distance via Inexact Proximal Point Method
- [`pwdist()`](https://www.kisungyou.com/T4transport/reference/pwdist.md)
  : Procrustes-Wasserstein Distance
- [`sinkhorn()`](https://www.kisungyou.com/T4transport/reference/sinkhorn.md)
  [`sinkhornD()`](https://www.kisungyou.com/T4transport/reference/sinkhorn.md)
  : Wasserstein Distance via Entropic Regularization and Sinkhorn
  Algorithm
- [`swdist()`](https://www.kisungyou.com/T4transport/reference/swdist.md)
  : Sliced Wasserstein Distance
- [`wassboot()`](https://www.kisungyou.com/T4transport/reference/wassboot.md)
  : Wasserstein Distance Estimation with Boostrapping
- [`wasserstein()`](https://www.kisungyou.com/T4transport/reference/wasserstein.md)
  [`wassersteinD()`](https://www.kisungyou.com/T4transport/reference/wasserstein.md)
  : Wasserstein Distance via Linear Programming

### Free-Support Centroids

- [`pwbary()`](https://www.kisungyou.com/T4transport/reference/pwbary.md)
  : Procrustes-Wasserstein Barycenter
- [`rbary23L()`](https://www.kisungyou.com/T4transport/reference/rbary23L.md)
  : Free-Support Barycenter by von Lindheim (2023)
- [`rbaryGD()`](https://www.kisungyou.com/T4transport/reference/rbaryGD.md)
  : Free-Support Barycenter by Riemannian Gradient Descent
- [`rmedIRLS()`](https://www.kisungyou.com/T4transport/reference/rmedIRLS.md)
  : Free-Support Median by IRLS
- [`rmedWB()`](https://www.kisungyou.com/T4transport/reference/rmedWB.md)
  : Free-Support Median by Weiszfeld Update with Barycentric Projection

### Fixed-Support Centroids

- [`fbary14C()`](https://www.kisungyou.com/T4transport/reference/fbary14C.md)
  [`fbary14Cdist()`](https://www.kisungyou.com/T4transport/reference/fbary14C.md)
  : Fixed-Support Barycenter by Cuturi & Doucet (2014)
- [`fbary15B()`](https://www.kisungyou.com/T4transport/reference/fbary15B.md)
  [`fbary15Bdist()`](https://www.kisungyou.com/T4transport/reference/fbary15B.md)
  : Fixed-Support Barycenter by Benamou et al. (2015)

## \[2\] Modality-Specific

### Cumulative Distribution Functions

- [`ecdfbary()`](https://www.kisungyou.com/T4transport/reference/ecdfbary.md)
  : Barycenter of Empirical CDFs
- [`ecdfmed()`](https://www.kisungyou.com/T4transport/reference/ecdfmed.md)
  : Wasserstein Median of Empirical CDFs

### Gaussian Measures

- [`gaussbary1d()`](https://www.kisungyou.com/T4transport/reference/gaussbary1d.md)
  : Barycenter of Gaussian Distributions in \\\mathbb{R}\\
- [`gaussbarypd()`](https://www.kisungyou.com/T4transport/reference/gaussbarypd.md)
  : Barycenter of Gaussian Distributions in \\\mathbb{R}^p\\
- [`gaussmed1d()`](https://www.kisungyou.com/T4transport/reference/gaussmed1d.md)
  : Wasserstein Median of Gaussian Distributions in \\\mathbb{R}\\
- [`gaussmedpd()`](https://www.kisungyou.com/T4transport/reference/gaussmedpd.md)
  : Wasserstein Median of Gaussian Distributions in \\\mathbb{R}^p\\

### Histograms

- [`histbary()`](https://www.kisungyou.com/T4transport/reference/histbary.md)
  : Barycenter of Histograms
- [`histbary14C()`](https://www.kisungyou.com/T4transport/reference/histbary14C.md)
  : Barycenter of Histograms by Cuturi and Doucet (2014)
- [`histbary15B()`](https://www.kisungyou.com/T4transport/reference/histbary15B.md)
  : Barycenter of Histograms by Benamou et al. (2015)
- [`histdist()`](https://www.kisungyou.com/T4transport/reference/histdist.md)
  : Distance between Histograms
- [`histinterp()`](https://www.kisungyou.com/T4transport/reference/histinterp.md)
  : Interpolation between Histograms
- [`histmed()`](https://www.kisungyou.com/T4transport/reference/histmed.md)
  : Wasserstein Median of Histograms

### Images

- [`imagebary()`](https://www.kisungyou.com/T4transport/reference/imagebary.md)
  : Barycenter of Images
- [`imagebary14C()`](https://www.kisungyou.com/T4transport/reference/imagebary14C.md)
  : Barycenter of Images according to Cuturi & Doucet (2014)
- [`imagebary15B()`](https://www.kisungyou.com/T4transport/reference/imagebary15B.md)
  : Barycenter of Images according to Benamou et al. (2015)
- [`imagedist()`](https://www.kisungyou.com/T4transport/reference/imagedist.md)
  : Wasserstein Distance between Two Images
- [`imageinterp()`](https://www.kisungyou.com/T4transport/reference/imageinterp.md)
  : Interpolation between Images
- [`imagemed()`](https://www.kisungyou.com/T4transport/reference/imagemed.md)
  : Wasserstein Median of Images

## \[3\] Gromov–Wasserstein

- [`gwbary()`](https://www.kisungyou.com/T4transport/reference/gwbary.md)
  : Gromov-Wasserstein Barycenter
- [`gwdist()`](https://www.kisungyou.com/T4transport/reference/gwdist.md)
  : Gromov-Wasserstein Distance

## \[4\] Data Sets

- [`digit3`](https://www.kisungyou.com/T4transport/reference/digit3.md)
  : MNIST Images of Digit 3
- [`digits`](https://www.kisungyou.com/T4transport/reference/digits.md)
  : MNIST Images of All Digits

## \[5\] Utilities

- [`fiedler()`](https://www.kisungyou.com/T4transport/reference/fiedler.md)
  : Compute the fiedler vector of a point cloud
- [`gaussvis2d()`](https://www.kisungyou.com/T4transport/reference/gaussvis2d.md)
  : Sampling from a Bivariate Gaussian Distribution for Visualization
- [`img2measure()`](https://www.kisungyou.com/T4transport/reference/img2measure.md)
  : Extract a discrete measure from a gray-scale image matrix
