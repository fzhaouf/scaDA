
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scaDA

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

## Examples

Two example datasets are included in the package to illustrate the
application of scaDA. The first dataset, Human Brain, exemplify how
scaDA identifies DA regions across various cell type from normal
samples. The AD dataset showcases the method’s capability to determine
DA regions within the same cell type under distinct states or
conditions.

### Human Brain

``` r
library(scaDA)
data("Hbrain", package = "scaDA") # load human brain dataset
```

The cell type composition is shown as following and we choose “Inh L1-4
LAMP5 LCP2” as interested cell type for demonstration purpose.

| cell type               | cell size | prop  |
|-------------------------|-----------|-------|
| Inh L4-6 SST B3GAT2     | 533       | 0.188 |
| Exc L5-6 RORB TTC12     | 502       | 0.177 |
| Exc L5-6 FEZF2 ABO      | 418       | 0.147 |
| Exc L5-6 FEZF2 EFTUD1P1 | 350       | 0.123 |
| Exc L3-5 RORB ESR1      | 204       | 0.072 |
| Exc L4-5 RORB FOLH1B    | 144       | 0.051 |
| Inh L1-4 LAMP5 LCP2     | 121       | 0.042 |
| Exc L5-6 THEMIS C1QL3   | 105       | 0.037 |
| Inh L2-3 VIP CASC6      | 105       | 0.037 |
| Exc L4-6 RORB SEMA3E    | 100       | 0.035 |

#### construct scaDAdataset object

``` r
counts = Hbrain@assays$ATAC@counts
coldata = Hbrain@meta.data$celltype
scaDA.obj <- scaDAdatasetFromMatrix(count = as.matrix(counts), colData = data.frame(coldata))
```

#### DA analysis pipeline

``` r
scaDA.obj <- estParams(scaDA.obj,celltype2 = "Inh L1-4 LAMP5 LCP2")
#> start initial parameter estiamte
scaDA.obj <- shrinkDisp(scaDA.obj)
#> start shrink dispersion parameter
scaDA.obj <- optParams(scaDA.obj)
#> start optimize parameter estimates
# report results in a dataframe
result = scaDA.obj@params$result
print(result[c(1:10),])
#>       tstats         pval          FDR     log2fc
#> 1  22.215805 5.882247e-05 0.0003001146  0.7336385
#> 2  23.821714 2.721499e-05 0.0001755806  0.7841482
#> 3  20.962088 1.072021e-04 0.0004305305  0.2124737
#> 4   8.949328 2.997224e-02 0.0316831263  0.2123390
#> 5  12.593176 5.604318e-03 0.0070967017 -0.4007064
#> 6  20.519291 1.324684e-04 0.0004987510 -0.9389558
#> 7  15.590417 1.375691e-03 0.0022928178 -1.4311649
#> 8  18.136637 4.122105e-04 0.0009699070 -0.5684360
#> 9  13.825139 3.153093e-03 0.0043671653 -0.2941101
#> 10  5.427206 1.430579e-01 0.1434884022  0.2183336
```

### Alzheimer’s disease dataset

The included example dataset is for determining DA regions within the
same cell type under distinct states or conditions. This dataset
contains 3 batches of Normal/AD samples and for each sample sequencing
data is available for 6 cell types: ASC, EX, IHN, MG, ODC, OPC. For
demonstration purpose, we use MG cell type from “Sample-100” of controls
and “Sample-43” of ADs.

| Batch |            | Control   |
|-------|------------|-----------|
| B1    | sample-96  | sample-43 |
| B1    | sample-100 | sample-45 |
| B2    | sample-82  | sample-46 |
| B2    | sample-66  | sample-40 |
| B2    | sample-101 | sample-27 |
| B2    | \\         | sample-50 |
| B2    | \\         | sample-22 |
| B3    | sample-90  | sample-37 |
| B3    | sample-52  | sample-33 |
| B3    | sample-58  | sample-47 |
| B3    | \\         | sample-17 |
| B3    | \\         | sample-19 |

``` r
data("ADdata", package = "scaDA") # load AD dataset
table(ADdata@meta.data$Cell.Type) # cell types and their size in AD/Normal
#> 
#>     ASC      EX     INH      MG     ODC     OPC PER.END 
#>     670     971     382     619    4156     286      67

selec.celltype="MG" # select MG as cell of interest
AD.samples =ADdata@meta.data$Diagnosis=='AD'
control.samples=ADdata@meta.data$Diagnosis=='Control'
celltype.loc=ADdata@meta.data$Cell.Type==selec.celltype
ADcell.loc= which(AD.samples & celltype.loc)
contrlcell.loc=which(control.samples & celltype.loc)
print(paste0('AD sample size for ',selec.celltype,': ',length(ADcell.loc)))
#> [1] "AD sample size for MG: 290"
print(paste0('Control sample size for ',selec.celltype,': ',length(contrlcell.loc)))
#> [1] "Control sample size for MG: 329"

counts = ADdata@assays$ATAC@counts 
coldata = ADdata@meta.data$Cell.Type
scaDA.obj <- scaDAdatasetFromMatrix(count = as.matrix(counts), colData = data.frame(coldata))
scaDA.obj <- estParams2(scaDA.obj,contrl.loc=contrlcell.loc,case.loc=ADcell.loc)
#> start initial parameter estiamte
scaDA.obj <- shrinkDisp(scaDA.obj)
#> start shrink dispersion parameter
scaDA.obj <- optParams(scaDA.obj)
#> start optimize parameter estimates
# report results in a dataframe
result = scaDA.obj@params$result
print(result[c(1:10),])
#>      tstats         pval          FDR     log2fc
#> 1  25.61022 1.150848e-05 1.490736e-05 0.13568317
#> 2  53.21647 1.648524e-11 1.485157e-10 1.08666235
#> 3  40.95068 6.698611e-09 1.824398e-08 0.99799114
#> 4  44.30547 1.299712e-09 4.592623e-09 0.08885513
#> 5  46.28596 4.930511e-10 2.012453e-09 0.41448943
#> 6  45.73556 6.455348e-10 2.561646e-09 0.78501542
#> 7  39.53505 1.336901e-08 3.252800e-08 0.72939522
#> 8  17.58483 5.356512e-04 6.093759e-04 0.41432035
#> 9  53.19752 1.663933e-11 1.485654e-10 0.53591747
#> 10 23.90239 2.617998e-05 3.264336e-05 1.04432405
```
