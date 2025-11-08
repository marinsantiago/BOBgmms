# BOBgmms <img src="man/figures/BOBgmms.png" alt="Description of image" width="130" height="150" align="right"> 

<!-- badges: start -->

[![R-CMD-check](https://github.com/marinsantiago/BOBgmms/workflows/R-CMD-check/badge.svg)](https://github.com/marinsantiago/BOBgmms/workflows/R-CMD-check/badge.svg)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

</br>

## Overview

This repository contains the R package `BOBgmms` (developer's version), which 
implements the [Bayesian optimized bootstrap for approximate posterior sampling in Gaussian mixture models](https://doi.org/10.1007/s11222-025-10763-y) (Marin et al., 2026).

## Installation

You can install the developer's version via `devtools` as:

``` r
# install.packages("devtools")
devtools::install_github("marinsantiago/BOBgmms")
```

On the other hand, if you wish to install the package from the `BOBgmms.zip` file in the supplementary materials to Marin et al. (2026):

  1. Decompress the zip file `BOBgmms.zip`. The folder `BOBgmms` should result.
  
  2. In R, set your working directory to the folder `BOBgmms`.
  
  3. Run the following R code:
  
``` r
# install.packages("devtools")
devtools::build()
devtools::install()
```

## <a name="system"></a> System Requirements

Parallelization over multiple CPU workers is conducted via *forking* 
(rather than *sockets*), so it is only available on POSIX systems (e.g., macOS, 
Linux, Unix, BSD), not on Windows. On non-POSIX platforms (such as Windows), 
the package functions are still operational but the number of CPU workers will 
be automatically set to one. For further details, see 
[Package '`parallel`'](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf). 

One can verify the OS type by running the following R code:

``` r
.Platform$OS.type
```

## Examples

Detailed guidelines for using the package functions are referred to their help pages in R. Additional examples are available at [https://github.com/marinsantiago/BOBgmms-examples](https://github.com/marinsantiago/BOBgmms-examples).

## <a name="cite"></a> Citation

If you use any part of this package in your work, please consider citing our *Statistics and Computing* paper:

```
@article{marin_bob,
  title   = {BOB: Bayesian optimized bootstrap for approximate posterior sampling in Gaussian mixture models},
  author  = {Santiago Marin and Bronwyn Loong and Anton H. Westveld},
  journal = {Statistics and Computing},
  volume  = {36},
  pages   = {14},
  year    = {2026},
  doi     = {10.1007/s11222-025-10763-y}
}
```

## Disclaimer

The software is provided "as is", without warranty of any kind, express or implied,
including but not limited to the warranties of merchantability, fitness for a particular
purpose and noninfringement. In no event shall the authors or copyright holders be liable
for any claim, damages, or other liability, whether in an action of contract, 
tort or otherwise, arising from, out of, or in connection with the software or the use
or other dealings in the software.

## <a name="refs"></a> References

Marin, S., Loong, B., and Westveld, A. H. (2026), "BOB: Bayesian optimized bootstrap for approximate posterior sampling in Gaussian mixture models." *Statistics and Computing*, **36**, 14. [doi:10.1080/10618600.2025.2572327](https://doi.org/10.1007/s11222-025-10763-y)
