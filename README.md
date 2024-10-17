IReNA: integrated regulatory network analysis of single-cell
transcriptomes
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![](https://img.shields.io/badge/r-version4.04-green.svg)](https://www.r-project.org)
[![](https://img.shields.io/badge/Seurat-version4.01-red.svg)](https://satijalab.org/seurat/articles/get_started.html)
[![](https://img.shields.io/badge/monocle-version2.18-blue.svg)](http://cole-trapnell-lab.github.io/monocle-release)
[![](https://img.shields.io/badge/publication-iscience-purple.svg)](https://www.cell.com/iscience/pdf/S2589-0042(22)01631-5.pdf)

IReNA (Integrated Regulatory Network Analysis) is an R package to
perform regulatory network analysis. IReNA contains two methods to
reconstruct gene regulatory networks. The first is using single-cell RNA
sequencing (scRNA-seq) data alone. The second is integrating scRNA-seq
data and chromatin accessibility profiles from Assay for Transposase
Accessible Chromatin using sequencing (scATAC-seq or bulk ATAC-seq).
IReNA performs modular regulatory network to reveal key transcription
factors and significant regulatory relationships among modules,
providing biological insights on regulatory mechanisms.

## Workflow

<img src="docs/Readme%20figure/Workflow1.png"
style="width:80.0%;height:80.0%" />

## Installation

IReNA needs R version 4.0 or higher, and
[Bioconductor](http://bioconductor.org/) version 3.12.

First, install Bioconductor, open R platform and run:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(version = "3.12")
```

Next, install several Bioconductor dependencies:

``` r
BiocManager::install(c('Rsamtools', 'ChIPseeker', 'monocle',
                       'RcisTarget', 'RCy3', 'clusterProfiler'))
```

Then, install IReNA from GitHub:

``` r
install.packages("devtools")
devtools::install_github("jiang-junyao/IReNA")
```

Finally, check whether IReNA was installed correctly, restart R session
and run:

``` r
library(IReNA)
```

## Quick start

The primary novelty of IReNA lies in its ability to decode regulatory
relationships among modules using a hypergeometric test. Consequently,
if you possess pre-built networks and gene groups, you can directly
execute IReNA with the following code.

``` r
library(IReNA)
### load test data
load(system.file("extdata", "qucik_start_test.rda", package = "IReNA"))
### Structure of input data. Please note that colnames of your table should be 
### same as the following test data
print(head(grn_test))
#>                    TF TFGroup          Target TargetGroup Correlation
#> 3370  ENSG00000114315       2 ENSG00000013441           1   0.7201089
#> 3404  ENSG00000118263       2 ENSG00000013441           1   0.7146384
#> 4323  ENSG00000108055       3 ENSG00000013441           1  -0.7023170
#> 10725 ENSG00000025156       1 ENSG00000025156           1   1.0000000
#> 18768 ENSG00000064703       1 ENSG00000064703           1   1.0000000
#> 21449 ENSG00000066422       1 ENSG00000066422           1   1.0000000
print(head(group_test))
#>                 KmeansGroup
#> ENSG00000011007           1
#> ENSG00000013441           1
#> ENSG00000015479           1
#> ENSG00000023516           1
#> ENSG00000025156           1
#> ENSG00000028839           1
### IReNA analysis
IReNA_result <- network_analysis(grn_test,group_test)
#> [1] "Total TFs: 227"
#> [1] "Enriched TFs: 76"
#> [1] "Significant regulations: 5"
### Here is enriched TFs that regulate other module
IReNA_result$TF_module_regulation[1:3,]
#>                     TF TFGroup LogFDR TargetGroup RegulationType
#> out2   ENSG00000124766       1    Inf      Group4       Negative
#> out2.1 ENSG00000126003       1    Inf      Group3       Negative
#> out2.2 ENSG00000196757       1    Inf      Group3       Negative
### Here is the network of enriched TFs
IReNA_result$TF_network[1:3,]
#>                     TF TFGroup TFMinNlogfdr TFMinGroup SigActModules
#> 262739 ENSG00000124766       1          Inf         N4            NA
#> 270782 ENSG00000126003       1          Inf         N3            NA
#> 742638 ENSG00000196757       1          Inf         N3            NA
#>        SigRepModules          Target TargetGroup Correlation Regulation
#> 262739             4 ENSG00000124766           1           1   Positive
#> 270782             3 ENSG00000126003           1           1   Positive
#> 742638             3 ENSG00000196757           1           1   Positive
### Here is the simplified netowrk
IReNA_result$intramodular_network[1:3,]
#>                    TFGroup TargetGroup Regulation        Correlation
#> Regulation12Pnum.4       2           2   Positive  0.808763726245598
#> Regulation12Nnum.5       2           3   Negative -0.719932554995428
#> Regulation21Nnum.3       3           2   Negative -0.719932554995428
#>                    NumberRegulation       Pvalue  NlogFdr
#> Regulation12Pnum.4   180;184;66;184 1.755047e-51 50.23283
#> Regulation12Nnum.5      16;16;16;16 0.000000e+00      Inf
#> Regulation21Nnum.3      16;16;16;16 0.000000e+00      Inf
```

## Full tutorials

- [Regulatory network analysis through only scRNA-seq
  data](https://jiang-junyao.github.io/IReNA/only-scRNA)

- [Regulatory network analysis through intergrating scRNA-seq data and
  scATAC-seq data](https://jiang-junyao.github.io/IReNA/scATAC+scRNA)

- [Regulatory network analysis through intergrating scRNA-seq data and
  bulk ATAC-seq
  data](https://jiang-junyao.github.io/IReNA/bulk-ATAC+scRNA)

## External links

An example for [using IReNA to identify transcription factors critical
for retinal
regeneration](https://github.com/jiewwwang/Single-cell-retinal-regeneration)

## Citation

Official publication: [IReNA: integrated regulatory network analysis of
single-cell
transcriptomes](https://www.cell.com/iscience/pdf/S2589-0042(22)01631-5.pdf)

## Help and Suggestion

If you have any question, comment or suggestion, please use github issue
tracker to report issues of IReNA or contact <jyjiang@link.cuhk.edu.hk>.
I will answer you timely, and please remind me again if you have not
received response more than three days.
