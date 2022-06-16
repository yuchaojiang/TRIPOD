# TRIPOD
Detecting Transcriptional Regulatory Relationships in Single-Cell RNA and Chromatin Accessibility Multiomic Data

## Author
Yuriko Harigaya, Nancy R. Zhang, Yuchao Jiang

## Maintainer
Yuriko Harigaya <harigaya@email.unc.edu>

## Description
TRIPOD is a statistical framework for detecting three-way regulatory relationships between a cis-regulatory region, a transcription factor, and a target gene, which we call a "trio,"
using single-cell ATAC/RNA multiomic data. The main functionality of this package is as follows.

* Infers trio regulatory relationships using robust nonparametric models
* Builds RNA prediction models from ATAC-seq data
* Identifies cell-types in which a trio regulatory relationship is active by estimating the influence of a data point in linear regression

## Installation
```r
install.packages("devtools")
devtools::install_github("yharigaya/TRIPOD/package")
```

## Vignettes

The vignettes for the PBMC data provide a brief overview of the data analysis
using the TRIPOD package, whereas those for the mouse embryonic brain data
contain detailed descriptions of the procedure.

* [Preprocessing - PBMC](http://htmlpreview.github.io/?https://github.com/yharigaya/TRIPOD/blob/main/vignettes/preprocessing_pbmc.html)
* [TRIPOD - PBMC](http://htmlpreview.github.io/?https://github.com/yharigaya/TRIPOD/blob/main/vignettes/TRIPOD_pbmc.html)
* [Preprocessing - mouse embryonic brain](http://htmlpreview.github.io/?https://github.com/yharigaya/TRIPOD/blob/main/vignettes/preprocessing_e18.html)
* [TRIPOD - mouse embryonic brain](http://htmlpreview.github.io/?https://github.com/yharigaya/TRIPOD/blob/main/vignettes/TRIPOD_e18.html)

##  Manuscript code

* [10X Genomics PBMC](https://github.com/yharigaya/TRIPOD/tree/main/scripts/10x_pbmc/)
* [10X Genomics mouse embryonic brain](https://github.com/yharigaya/TRIPOD/tree/main/scripts/10x_e18/)
* [SHARE-seq mouse skin](https://github.com/yharigaya/TRIPOD/tree/main/scripts/share_seq_skin/)
* [SNARE-seq adult mouse brain](https://github.com/yharigaya/TRIPOD/tree/main/scripts/snare_seq_mbrain/)


##  Common questions

* [What if my TF/motif of interest is not included in the default annotation?](https://github.com/yharigaya/TRIPOD/blob/main/instruction/custom.md)
* [Annotated motifs in addition to the JASPAR database?](https://github.com/yharigaya/TRIPOD/blob/main/instruction/database.md)
* [How to handle multiple TFs with the same motif?](https://github.com/yharigaya/TRIPOD/blob/main/instruction/multiple.md)

## Reference
Yuriko Harigaya, Zhaojun Zhang, Hongpan Zhang, Chongzhi Zang, Nancy R Zhang, Yuchao Jiang. Nonparametric Interrogation of Transcriptional Regulation in Single-Cell RNA and Chromatin Accessibility Multiomic Data. ***bioRxiv***, 2021. ([link](https://www.biorxiv.org/content/10.1101/2021.09.22.461437v1))
