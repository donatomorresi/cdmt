
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Change Detection by Multispectral Trends (CDMT)

<!-- badges: start -->

[![R-CMD-check](https://github.com/donato-morresi/cdmt/workflows/R-CMD-check/badge.svg)](https://github.com/donato-morresi/cdmt/actions)
<!-- badges: end -->

## About

The `cdmt` *R* package provides an implementation of CDMT, which is a
Landsat time series-based algorithm for mapping inter-annual changes in
linear trends at the pixel level. CDMT was designed to detect abrupt and
gradual spectral changes associated with forest disturbance dynamics. A
detailed description of the algorithm will be provided in an upcoming
paper. Landsat time series can be univariate, i.e. include a single
spectral band/index, or multivariate, i.e. include multiple spectral
bands/indices. A modified version of the High Dimensional Trend
Segmentation (HiTS) procedure proposed by Maeng (2019) is at the core of
CDMT. `cdmt` relies on the `terra` package for managing raster data and
parallelising computations.

## Installation

You can install the development version of `cdmt` with:

``` r
# install.packages("devtools")
devtools::install_github("donatomorresi/cdmt")
```

### References

Maeng, H. (2019). Adaptive multiscale approaches to regression and trend
segmentation. Ph.D. thesis, London School of Economics and Political
Science.
