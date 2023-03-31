
<!-- README.md is generated from README.Rmd. Please edit that file -->

# msireg

<!-- badges: start -->
<!-- badges: end -->

The goal of *msireg* is to co-register mass spectrometry images (MSI)
with microscopic optical images using `SimpleIK`.

## Installation

You can install the development version of msireg from
[GitHub](https://github.com/) with:

``` r
if ( !require(devtools) ) { # if not installed already 
    install.packages("devtools") 
}

devtools::install_github("sslakkimsetty/msireg")
```

## Example usage

# Details

*msireg* co-registers MS images with optical images (H&E stained, AF)
using the family of registration methods from the `SimpleITK` package.
The package offers processing and summarization of MS imaging data and
optical images. This package uses spatial shrunken centroids method from
the `Cardinal` package to filter out unimportant features.

# References

\[1\] Bemis, K. D., Harry, A., Eberlin, L. S., Ferreira, C., van de Ven,
S. M., Mallick, P., Stolowitz, M., and Vitek, O.”, “Cardinal: an R
package for statistical analysis of mass spectrometry-based imaging
experiments.

\[2\] L.J.P. van der Maaten and G.E. Hinton. Visualizing
High-Dimensional Data Using t-SNE. Journal of Machine Learning Research
9(Nov):2579-2605, 2008.

\[3\] Gregoire Pau, Florian Fuchs, Oleg Sklyar, Michael Boutros, and
Wolfgang Huber (2010):“,”EBImage - an R package for image processing
with applications to cellular phenotypes.”, “Bioinformatics

\[4\] Beare, R., Lowekamp, B., & Yaniv, Z. (2018). Image Segmentation,
Registration and Characterization in R with SimpleITK. Journal of
Statistical Software, 86(8), 1–35.
<https://doi.org/10.18637/jss.v086.i08>

\[5\] Reference for data (in data/)
