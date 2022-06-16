## What if my TF/motif of interest is not included in the default annotation?
  
Currently, `TRIPOD` can be used with two comprehensive TF binding
motif databases,
[JASPAR](https://jaspar.genereg.net/)
and [HOCOMOCO](https://hocomoco11.autosome.ru/downloads_v11).

If your TF/motif of interest is not included in the databases,
you can manually create a position frequency matrix (PFM) and add it to
the motif data. You can also replace the existing motif information
for a given TF with a custom PFM.

In the following example script, the PFM for the mouse Olig2 TF
is replaced with a manually-created PFM corresponding to the `CAGCTG` motif
in a list of PFMs from the JASPAR2020 database. 
Note that 633 TF binding motifs in human (taxon ID: 9696) are used
because motifs are highly conserved 
and the list of mouse TF motifs from the database contains only 107
entries.
Also, note that `<path/to/seurat>` needs to be replaced with a path to
a file containing a Seurat object processed from multiome data.
An example file, which has been processed from the 10X Genomics mouse
embryonic data, can be downloaded
[here](https://www.dropbox.com/s/flzriz8m8uiaopn/e18.qc.rds?dl=0).

```r
library(tidyverse)
library(Seurat)
library(Signac)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)

# read in the Seurat object
object <- <path/to/seurat> %>% readRDS()

# get a motif object from the JASPAR database
pfm.set <- getMatrixSet(
    x = JASPAR2020,
    opts = list(species = 9606, all_versions = FALSE)
)

# construct a PFM of the Olig2 binding motif
pfm.olig2 <- PFMatrix(
    ID = "PFM.Olig2",
    name = "Olig2",
    strand = "+",
    bg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
    tags = list(
        species = "10090",
        tax_group = "vertebrates"),
        profileMatrix = matrix(
        c(0L, 1L, 0L, 0L, 0L, 0L, # A
        1L, 0L, 0L, 1L, 0L, 0L, # C
        0L, 0L, 1L, 0L, 0L, 1L, # G
        0L, 0L, 0L, 0L, 1L, 0L), # T
        byrow = TRUE, nrow = 4,
        dimnames = list(c("A", "C", "G", "T")))
)

# remove the original PFM of the Olig2 motif from the PFM set
pfm.set <- pfm.set[-which(sapply(pfm.set, function(x) x@name == "OLIG2"))]

# add the custom PFM
pfm.set$PFM.Olig2 <- pfm.olig2

# creates a motif matrix
motif.matrix <- CreateMotifMatrix(
    features = granges(object),
    pwm = pfm.set,
    genome = "mm10",
    use.counts = FALSE
)

# add the PFM set to a Seurat object
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pfm.set)
object <- SetAssayData(object, assay = "ATAC", slot = "motifs", 
    new.data = motif.object)
```

See also section 3.6 in [the preprocessing vignettes for the mouse embryonic
brain data](http://htmlpreview.github.io/?https://github.com/yharigaya/TRIPOD/blob/main/vignettes/preprocessing_e18.html).
