IReNA: integrated regulatory network analysis of single-cell
transcriptomes
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

IReNA (Integrated Regulatory Network Analysis) is an R package to
reconstruct regulatory networks through integrating scRNA-seq and
ATAC-seq data. Compared with other regulatory network analysis method
(SCENIC), IReNA provides modularized regulatory network analysis to
discover the biological significance of transcription factors and the
regulatory role of each module throughout the process.

IReNA is still under testing, so there may be some changes for this
package in the coming weeks.

## Installation

IReNA needs R version 4.0 or higher,
[Bioconductor](http://bioconductor.org/) version 3.12.

First install Bioconductor, open R and run:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(version = "3.12")
```

Next, install a few Bioconductor dependencies that aren’t automatically
installed:

``` r
BiocManager::install(c('Rsamtools', 'ChIPseeker', 'monocle',
                       'RcisTarget', 'RCy3', 'clusterProfiler'))
```

Then, install IReNA from GitHub:

``` r
install.packages("devtools")
devtools::install_github("jiang-junyao/IReNA")
```

Finally, check whether IReNA was installed correctly, start a new R
session and run:

``` r
library(IReNA)
```

## Examples and tutorials

If you use **both ATAC-seq data and scRNA-seq or bulk RNA-seq data** to
reconstruct regulatory network, please run [ATAC-seq raw data
preprocessing](https://jiang-junyao.github.io/IReNA/ATAC-seq-preprocessing)
first to get bam file, peaks file and footprints file, then run codes of
**part 1**, **part 3** and **part 4** in [IReNA
tutorial](https://jiang-junyao.github.io/IReNA/tutorial). If you **only
use scRNA-seq or bulk RNA-seq data** to reconstruct regulatory network,
just directly run codes of **part 1**, **part 2** and **part 4** in
[IReNA tutorial](https://jiang-junyao.github.io/IReNA/tutorial).

-   [ATAC-seq raw data
    preprocessing](https://jiang-junyao.github.io/IReNA/ATAC-seq-preprocessing)

-   [IReNA tutorial](https://jiang-junyao.github.io/IReNA/tutorial)

## Workflow

![workflow](docs/Readme%20figure/Workflow.png)

## How to cite this package

If you use IReNA package, please cite the following Science
paper: <https://science.sciencemag.org/content/370/6519/eabb8598>.

## Help and Suggestion

If you have any question, comment or suggestion, please use github issue
tracker to report coding related issues of IReNA or contact
<jiangjunyao789@163.com>. I will answer you timely, and please remind me
again if you have not received response more than three days.
