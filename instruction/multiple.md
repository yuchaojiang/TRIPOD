## How to handle multiple TFs with the same motif?
  
Our original workflow, which is described in the vignettes, excludes
cases where a given motif corresponds to multiple TFs as well as
those where a given TF corresponds to multiple TFs. 
To allow for such cases, one can use an alternative set of functions
(`getObjectsForModelFit2()`, `filterSeuratObject2()`, 
`getXYMatrices2()`, and `fitModel2()`) as shown in the
following example code.
Note that `<path/to/seurat>` needs to be replaced with a path to
a file containing a Seurat object processed from multiome data.
An example file, which has been processed from the 10X Genomics mouse
embryonic data, can be downloaded
[here](https://www.dropbox.com/s/4afi9rp4t5d5km0/e18.chromvar.rds?dl=0).

```r
library(tidyverse)
library(Seurat)
library(Signac)
library(GenomicRanges)
library(BiocParallel)
library(TRIPOD)
library(nbpMatching)

# read in the Seurat object
object <- <path/to/seurat> %>% readRDS()

# get objects required for model fitting
tripod <- getObjectsForModelFit2(
    object = object,
    chr = paste0("chr", 1:19), # mouse
    convert = TRUE
)

object <- filterSeuratObject2(
    object = object, tripod.object = tripod
)

object <- processSeuratObject(
    object = object,
    dim.rna = 1:30,
    dim.atac = 2:50,
    verbose = FALSE
)

res <- 15
object <- getClusters(
    object = object,
    graph.name = "SCT_snn",
    resolution = res,
    verbose = FALSE
)

metacell <- getMetacellMatrices(
    object = object,
    cluster.name = "seurat_clusters",
    min.num = 20
)

object <- removeSmallMetacells(object = object, min.num = 20)

# fit models for representative genes
genes <- c("Neurog2", "Eomes", "Olig2", "Sox10")
ext.upstream <- 2e5
xymats.list <- bplapply(
    genes,
    getXYMatrices2,
    ext.upstream = ext.upstream,
    tripod = tripod,
    metacell = metacell
)
names(xymats.list) <- genes

xymats.m.list <- bplapply(
    xymats.list,
    fitModel2,
    model.name = "marginal"
)
names(xymats.m.list) <- genes

xymats.c.list <- bplapply(
    xymats.list,
    fitModel2,
    model.name = "conditional"
)
names(xymats.c.list) <- genes

xymats.i.list <- bplapply(
    xymats.list,
    fitModel2,
    model.name = "interaction"
)
names(xymats.i.list) <- genes

xymats.tX.list <- bplapply(
    xymats.list,
    fitModel2,
    model.name = "TRIPOD",
    match.by = "Xp"
)
names(xymats.tX.list) <- genes

xymats.tY.list <- bplapply(
    xymats.list,
    fitModel2,
    model.name = "TRIPOD",
    match.by = "Yt"
)
names(xymats.tY.list) <- genes
```


