# GapClust
## Detecting rare cells from expression profiles by single cell RNA-seq.

## Introduction
GapClust takes advantage of the gap between minor cluster and neighbouring abundant cluster to let rare cells within minor cluster stand out through delicately designed statistics. Meanwhile, GapClust does not
struggle to search for rare cell informative genes like most of the competitors, but learns the cluster size as well as rare cells using simple arithmetic calculation.

## Installation
R users  can easily install GapClust by running following code in R console.
```
devtools::install_github("fabotao/GapClust")
```
## R packages prerequisites
```
Seurat
rflann
irlba
e1071
```
## Publications

## Copyright
This software package is distributed under GNU GPL v3.

This work is free to use for academic and research purposes. Please contact maintainer for commercial use of this work.
