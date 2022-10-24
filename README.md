IReNA: integrated regulatory network analysis of single-cell
transcriptomes
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![](https://img.shields.io/badge/r-version4.04-green.svg)](https://www.r-project.org)
[![](https://img.shields.io/badge/Seurat-version4.01-red.svg)](https://satijalab.org/seurat/articles/get_started.html)
[![](https://img.shields.io/badge/monocle-version2.18-blue.svg)](http://cole-trapnell-lab.github.io/monocle-release)
[![](https://img.shields.io/badge/Published-iscience-purple.svg)](https://doi.org/10.1101/2021.11.22.469628)

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

<img src="docs/Readme%20figure/Workflow1.png" style="width:80.0%;height:80.0%" />

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

## Examples and tutorials

-   [Regulatory network analysis through only scRNA-seq
    data](https://jiang-junyao.github.io/IReNA/only-scRNA)

-   [Regulatory network analysis through intergrating scRNA-seq data and
    scATAC-seq data](https://jiang-junyao.github.io/IReNA/scATAC+scRNA)

-   [Regulatory network analysis through intergrating scRNA-seq data and
    bulk ATAC-seq
    data](https://jiang-junyao.github.io/IReNA/bulk-ATAC+scRNA)

## External links

An example for [using IReNA to identify transcription factors critical
for retinal
regeneration](https://github.com/jiewwwang/Single-cell-retinal-regeneration)

## Citation

Preprint: [IReNA: integrated regulatory network analysis of single-cell
transcriptomes](https://doi.org/10.1101/2021.11.22.469628)

[Hoang T, Wang J, Boyd P, et al. Gene regulatory networks controlling
vertebrate retinal regeneration. Science 2020;
370(6519):eabb8598](https://www.science.org/doi/10.1126/science.abb8598)

## Help and Suggestion

If you have any question, comment or suggestion, please use github issue
tracker to report issues of IReNA or contact <jiangjunyao789@163.com>. I
will answer you timely, and please remind me again if you have not
received response more than three days.
