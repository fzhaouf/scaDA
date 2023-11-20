
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scaDA

<!-- badges: start -->
<!-- badges: end -->

Single-cell ATAC-seq sequencing data (scATAC-seq) has been a widely
adopted technology to investigate chromatin accessibility on the
single-cell level. Analyzing scATAC-seq can provide valuable insights
into identifying cell populations and revealing the epigenetic
heterogeneity across cell populations in different biological contexts.
One important aspect of scATAC-seq data analysis is performing
differential chromatin accessibility (DA) analysis, which will help
identify cell populations and reveal epigenetic heterogeneity. While
numerous differential expression methods have been proposed for
single-cell RNA sequencing data, DA methods for scATAC-seq data are
underdeveloped and remain a major challenge due to the high sparsity and
high dimensionality of the data. To fill the gap, we introduce a novel
and robust zero-inflated negative binomial framework named scaDA for DA
analysis. The model links the prevalence, mean and dispersion parameters
to covariates such as cell populations, treatment conditions, and batch
effect. The statistical inference is based on the EM algorithm and the
dispersion parameter is shrunk using an empirical Bayes method to
stabilize the parameter estimation by leveraging information from other
accessible chromatin regions in the genome. Consequently, we performed
both simulation studies and real data applications, which demonstrate
the superiority of scaDA compared to existing approaches.

## Installation

You can install the development version of scaDA from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fzhaouf/scaDA")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(scaDA)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
