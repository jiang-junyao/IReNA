IReNA: integrated regulatory network analysis of single-cell
transcriptomes
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![](https://img.shields.io/badge/r-version4.04-green.svg)](https://www.r-project.org)
[![](https://img.shields.io/badge/Seurat-version4.01-red.svg)](https://satijalab.org/seurat/articles/get_started.html)
[![](https://img.shields.io/badge/monocle-version2.18-blue.svg)](http://cole-trapnell-lab.github.io/monocle-release)
[![](https://img.shields.io/badge/Preprint-biorxiv-purple.svg)](https://doi.org/10.1101/2021.11.22.469628)

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
### IReNA analysis
IReNA_result <- network_analysis(grn_test,group_test,2,2,ModuleFDR = 0.05)
#> [1] "Total TFs: 227"
#> [1] "Enriched TFs: 107"
#> [1] "Significant regulations: 6"
### Here is enriched TFs that regulate other module
IReNA_result$TF_module_regulation[1:5,]
#>                     TF TFGroup           LogFDR TargetGroup RegulationType
#> out2   ENSG00000070061       1 3.14614893110143      Group1       Positive
#> out2.1 ENSG00000071564       1 2.02272256522902      Group4       Negative
#> out2.2 ENSG00000124766       1              Inf      Group4       Negative
#> out2.3 ENSG00000126003       1              Inf      Group1       Positive
#> out2.4 ENSG00000126003       1              Inf      Group3       Negative
### Here is the network of enriched TFs
IReNA_result$TF_network[1:5,]
#>                     TF TFGroup TFMinNlogfdr TFMinGroup SigActModules
#> 29492  ENSG00000070061       1     3.146149         P1             1
#> 31393  ENSG00000182944       3    22.352923         P3             3
#> 34854  ENSG00000071564       1     2.022723         N4            NA
#> 35366  ENSG00000087510       2          Inf         N3             2
#> 262739 ENSG00000124766       1          Inf         N4            NA
#>        SigRepModules          Target TargetGroup Correlation Regulation
#> 29492             NA ENSG00000070061           1   1.0000000   Positive
#> 31393              2 ENSG00000070061           1  -0.7386442   Negative
#> 34854              4 ENSG00000071564           1   1.0000000   Positive
#> 35366              3 ENSG00000071564           1   0.7346326   Positive
#> 262739             4 ENSG00000124766           1   1.0000000   Positive
### Here is the simplified netowrk
IReNA_result$intramodular_network[1:5,]
#>                    TFGroup TargetGroup Regulation        Correlation
#> Regulation12Pnum         1           1   Positive                  1
#> Regulation12Pnum.4       2           2   Positive  0.804181395829719
#> Regulation12Nnum.5       2           3   Negative -0.719932554995428
#> Regulation21Nnum.3       3           2   Negative -0.719932554995428
#> Regulation12Pnum.7       3           3   Positive  0.854643538073071
#>                    NumberRegulation       Pvalue   NlogFdr
#> Regulation12Pnum        9;19;388;19 2.297460e-10  9.270775
#> Regulation12Pnum.4  290;306;101;306 6.269939e-56 54.658669
#> Regulation12Nnum.5      16;17;17;16 0.000000e+00       Inf
#> Regulation21Nnum.3      16;16;18;17 0.000000e+00       Inf
#> Regulation12Pnum.7     43;44;363;44 4.465846e-60 58.681089
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
