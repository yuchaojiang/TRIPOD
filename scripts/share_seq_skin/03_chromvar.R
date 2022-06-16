
library(Seurat)
library(Signac)

skin=readRDS('skin.frag.rds')

## run chromVAR
DefaultAssay(skin) <- "ATAC"
levels(skin)
#  [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14"
# [16] "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27"
celltype.names <- levels(skin)

library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)

DefaultAssay(skin) <- "ATAC"
# https://satijalab.org/signac/articles/motif_vignette.html
# get a list of motif position frequency matrices from the jaspar database
# https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=info&id=10090
# 10090 is the species ID for mouse
# 9606 is the species ID for human
# https://satijalab.org/signac/articles/motif_vignette.html uses 9606 for mouse
# and this link gives the justification: https://github.com/timoast/signac/issues/58 

pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(skin), 
                                  pwm = pwm_set, 
                                  genome = 'mm10', 
                                  use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
skin <- SetAssayData(skin, assay = 'ATAC', slot = 'motifs', new.data = motif.object)
skin <- RunChromVAR( # note that this step can take 30-60 minutes
  object = skin,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

saveRDS(skin, file='skin.chromvar.rds')


