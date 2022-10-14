# TRIPOD
Nonparametric interrogation of transcriptional regulation in single-cell RNA and chromatin accessibility multiomic data

## Author
Yuchao Jiang, Yuriko Harigaya, Nancy R. Zhang

## Maintainer
Yuchao Jiang <yuchaoj@email.unc.edu>

## Description
TRIPOD is a statistical framework for detecting three-way regulatory relationships between a cis-regulatory region, a transcription factor, and a target gene, which we call a "trio,"
using single-cell ATAC/RNA multiomic data. The main functionality of this package is as follows.

* Infers trio regulatory relationships using robust nonparametric models
* Builds RNA prediction models from ATAC-seq data
* Identifies influential cell types in which a trio regulatory relationship is active

## Installation
```r
install.packages("devtools")
devtools::install_github("yuchaojiang/TRIPOD/package")
```

## Vignettes

The vignettes for the PBMC data provide a brief overview of the data analysis
using the TRIPOD package, whereas those for the mouse embryonic brain data
contain detailed descriptions of the procedure.

* [Preprocessing - PBMC](http://htmlpreview.github.io/?https://github.com/yuchaojiang/TRIPOD/blob/main/vignettes/preprocessing_pbmc.html)
* [TRIPOD - PBMC](http://htmlpreview.github.io/?https://github.com/yuchaojiang/TRIPOD/blob/main/vignettes/TRIPOD_pbmc.html)
* [Preprocessing - mouse embryonic brain](http://htmlpreview.github.io/?https://github.com/yuchaojiang/TRIPOD/blob/main/vignettes/preprocessing_e18.html)
* [TRIPOD - mouse embryonic brain](http://htmlpreview.github.io/?https://github.com/yuchaojiang/TRIPOD/blob/main/vignettes/TRIPOD_e18.html)

## Issues & bugs

Please do not email but post [here](https://github.com/yuchaojiang/TRIPOD/issues) on GitHub for fastest reply/help.

##  Common questions

* [What if my TF/motif of interest is not included in the default annotation?](https://github.com/yuchaojiang/TRIPOD/blob/main/instruction/custom.md)
* [Annotated motifs in addition to the JASPAR database?](https://github.com/yuchaojiang/TRIPOD/blob/main/instruction/database.md)
* [How to handle multiple TFs with the same motif?](https://github.com/yuchaojiang/TRIPOD/blob/main/instruction/multiple.md)

##  Manuscript code

* [10X Genomics PBMC](https://github.com/yuchaojiang/TRIPOD/tree/main/scripts/10x_pbmc/)
* [10X Genomics mouse embryonic brain](https://github.com/yuchaojiang/TRIPOD/tree/main/scripts/10x_e18/)
* [SHARE-seq mouse skin](https://github.com/yuchaojiang/TRIPOD/tree/main/scripts/share_seq_skin/)
* [SNARE-seq adult mouse brain](https://github.com/yuchaojiang/TRIPOD/tree/main/scripts/snare_seq_mbrain/)

## Reference
Yuchao Jiang, Yuriko Harigaya, Zhaojun Zhang, Hongpan Zhang, Chongzhi Zang, Nancy R Zhang. Nonparametric single-cell multiomic characterization of three-way relationships between transcription factors, target genes, and chromatin accessibilities. ***Cell Systems***, 13 (9), 737-751, 2022. ([link](https://www.sciencedirect.com/science/article/abs/pii/S2405471222003489))
