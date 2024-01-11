
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scaDA

Single-cell ATAC-seq sequencing data (scATAC-seq) has been widely used
to investigate chromatin accessibility on the single-cell level. One
important application of scATAC-seq data analysis is performing
differential chromatin accessibility (DA) analysis. The data
characteristics of scATAC-seq such as excessive zeros and large
variability of chromatin accessibility across cells impose a unique
challenge for DA analysis. Existing statistical methods focus on
detecting the difference in mean of the chromatin accessible regions and
treat the dispersion and prevalence as nuisances. Motivated by real data
exploration, where dispersion and prevalence demonstrate distribution
differences among cell types, we introduce a novel composite statistical
test named “scaDA”, which is based on the zero-inflated negative
binomial regression model, for differential distribution analysis of
scATAC-seq by jointly testing the abundance, prevalence, and dispersion
simultaneously. scaDA further adopts an empirical Bayes shrinkage
technique and iterative estimation procedure to refine the estimates for
all three parameters.

# Installation

You can install scaDA from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fzhaouf/scaDA")
```

or from CRAN:

``` r
install.packages('scaDA')
```

# Examples

One example datasets are included in the package to illustrate the
application of scaDA. To test for DA peaks between two specific groups
of cells, please specify the group.1 and group.2 parameters in
“estParams” function. If the group.2 parameter is omitted or set to
NULL, the pipeline will test for DA peaks between the group specified by
group.1 parameter and all other cell types in the dataset. The first
example using Human Brain 3K, demonstrates how scaDA perform DA test
between interested cell type and all other cells. The second example
using the same dataset showcases the method’s capability to determine DA
regions between two specific groups in the data. The same setting can be
used to conduct DA test for same cell type between two distinct
conditions (e.g. normal/disease).

## DA test between specific cell type and all other types in Human Brain 3K

``` r
library(scaDA)
data("HumanBrain", package = "scaDA") # load human brain dataset
```

The cell type composition is shown as following and we choose
“microglia” as the interested cell type for demonstration purpose.

| cell type       | cell size | Proportion |
|-----------------|-----------|------------|
| granule neuron  | 636       | 22%        |
| oligodendrocyte | 552       | 20%        |
| cOPC            | 365       | 13%        |
| bergmann glia   | 344       | 12%        |
| ependymal       | 204       | 7%         |
| purkinje cell   | 144       | 5%         |
| astrocytes      | 125       | 4%         |
| microglia       | 104       | 4%         |

### construct scaDAdataset object

To construct scaDAdataset object using “scaDAdatasetFromMatrix”
function, it requires a peak-by-cell read counts matrix from the
scATAC-seq experiment along with cell labels information.

``` r
counts = HumanBrain@assays$ATAC@counts
coldata = HumanBrain@meta.data$celltype
scaDA.obj <- scaDAdatasetFromMatrix(count = as.matrix(counts), colData = data.frame(coldata))
```

### DA analysis pipeline

Specify the group.1 parameter in “estParams” function to interested cell
type and omit the group.2 parameter or set it to NULL.

``` r
scaDA.obj <- estParams(scaDA.obj, group.1 = "microglia")
#> start initial parameter estiamte
scaDA.obj <- shrinkDisp(scaDA.obj)
#> start shrink dispersion parameter
scaDA.obj <- optParams(scaDA.obj)
#> start optimize parameter estimates
# report results in a dataframe
results1 = scaDA.obj@result
print(results1[c(1:10),])
#>       tstats         pval          FDR      log2fc
#> 1  22.838016 4.364731e-05 0.0004275972 -0.73181994
#> 2  11.306194 1.018032e-02 0.0148401123  0.75685565
#> 3  22.917815 4.200760e-05 0.0004275972 -0.72107317
#> 4  12.055718 7.194712e-03 0.0113302551  0.16586482
#> 5   5.710075 1.266000e-01 0.1328436311 -0.08226962
#> 6  16.238320 1.013277e-03 0.0028146580  1.07895272
#> 7  11.369542 9.886599e-03 0.0145060396 -0.54729025
#> 8   9.398709 2.443370e-02 0.0303147655  0.43849966
#> 9   7.989589 4.622736e-02 0.0524714596 -0.60699970
#> 10 20.045826 1.660706e-04 0.0008364590 -0.93735478
```

## DA test between two specific cell types in Human Brain 3K

Specify both group.1 and group.2 parameters in “estParams” function to
interested cell types for DA test. microglia and astrocytes are used for
demonstration. The same setting can be used for case-control DA test.

``` r

counts = HumanBrain@assays$ATAC@counts
coldata = HumanBrain@meta.data$celltype
scaDA.obj <- scaDAdatasetFromMatrix(count = as.matrix(counts), colData = data.frame(coldata))
scaDA.obj <- estParams(scaDA.obj, group.1 = "microglia", group.2 = "astrocytes")
#> start initial parameter estiamte
scaDA.obj <- shrinkDisp(scaDA.obj)
#> start shrink dispersion parameter
scaDA.obj <- optParams(scaDA.obj)
#> start optimize parameter estimates
# report results in a dataframe
results2 = scaDA.obj@result
print(results2[c(1:10),])
#>       tstats         pval          FDR      log2fc
#> 1  15.829929 1.228754e-03 0.0065797653 -0.01083656
#> 2  17.223335 6.357910e-04 0.0046408101  1.25667661
#> 3  18.707597 3.142201e-04 0.0031422008 -0.53539134
#> 4  25.152561 1.434696e-05 0.0004628050  1.14661976
#> 5  27.696727 4.205273e-06 0.0001966811  0.67369349
#> 6   6.542731 8.799246e-02 0.0964829641  0.71616728
#> 7   9.619585 2.209262e-02 0.0328556281 -0.31919845
#> 8   8.992854 2.938605e-02 0.0399267006 -0.32318888
#> 9  15.167677 1.678822e-03 0.0073956904  0.01849312
#> 10 11.490644 9.348168e-03 0.0180914164 -0.71459130
```

## ver-1.0.1 update: add optParamsParallel

There is a parallel version based on R “parallel” library that
significantly speeds up the final optimization process.

``` r

counts = HumanBrain@assays$ATAC@counts
coldata = HumanBrain@meta.data$celltype
scaDA.obj <- scaDAdatasetFromMatrix(count = as.matrix(counts), colData = data.frame(coldata))
scaDA.obj <- estParams(scaDA.obj, group.1 = "microglia", group.2 = "astrocytes")
#> start initial parameter estiamte
scaDA.obj <- shrinkDisp(scaDA.obj)
#> start shrink dispersion parameter
scaDA.obj <- optParamsParallel(scaDA.obj)
#> start optimize parameter estimates
# report results in a dataframe
results3 = scaDA.obj@result
print(results2[c(1:10),])
#>       tstats         pval          FDR      log2fc
#> 1  15.829929 1.228754e-03 0.0065797653 -0.01083656
#> 2  17.223335 6.357910e-04 0.0046408101  1.25667661
#> 3  18.707597 3.142201e-04 0.0031422008 -0.53539134
#> 4  25.152561 1.434696e-05 0.0004628050  1.14661976
#> 5  27.696727 4.205273e-06 0.0001966811  0.67369349
#> 6   6.542731 8.799246e-02 0.0964829641  0.71616728
#> 7   9.619585 2.209262e-02 0.0328556281 -0.31919845
#> 8   8.992854 2.938605e-02 0.0399267006 -0.32318888
#> 9  15.167677 1.678822e-03 0.0073956904  0.01849312
#> 10 11.490644 9.348168e-03 0.0180914164 -0.71459130
```
