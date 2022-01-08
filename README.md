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
Accessible Chromatin using sequencing (ATAC-seq). IReNA performs modular
regulatory network to reveal key transcription factors and significant
regulatory relationships among modules, providing biological insights on
regulatory mechanisms.

## Installation

Please see the Installation part of [IReNA
tutorial](https://jiang-junyao.github.io/IReNA/tutorial#installation)

## Examples and tutorials

First of all, run [scRNA-seq preprocessing
pipeline](https://jiang-junyao.github.io/IReNA/scRNA-seq-preprocessing)
to get seurat object with pseudotime, or use your own methods to
calculate the pseudotime of each cell, and then add it to the metadata
of seurat object.

If **only scRNA-seq (or bulk RNA-seq data)** are used to reconstruct
regulatory network, just directly run the steps of **part 1**, **part
2** and **part 4** in [IReNA
tutorial](https://jiang-junyao.github.io/IReNA/tutorial).

If **both ATAC-seq data and scRNA-seq (or bulk RNA-seq data)** are used
to reconstruct regulatory network, please run [ATAC-seq raw data
preprocessing](https://jiang-junyao.github.io/IReNA/ATAC-seq-preprocessing)
first to get bam file, peaks file and footprint file, then run the steps
of **part 1**, **part 3** and **part 4** in [IReNA
tutorial](https://jiang-junyao.github.io/IReNA/tutorial).

-   [scRNA-seq raw data
    preprocessing](https://jiang-junyao.github.io/IReNA/scRNA-seq-preprocessing)

-   [ATAC-seq raw data
    preprocessing](https://jiang-junyao.github.io/IReNA/ATAC-seq-preprocessing)

-   [IReNA tutorial](https://jiang-junyao.github.io/IReNA/tutorial)

## Workflow

<img src="docs/Readme%20figure/Workflow.png" style="width:30.0%;height:30.0%" />

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
