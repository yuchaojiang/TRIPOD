
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79) # For mouse skin
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)
library(presto)

library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)

library(glmnet)
library(pheatmap)
library(fields)
library(qvalue)
library(gplots)
library(patchwork)
library(olsrr)
library(pdftools)
library(ape)
library(dendextend)
library(gplots)
library(nbpMatching)
library(fields)
library(RColorBrewer)
library(gridBase)
library(grid)

load('ATAC_RNA_WNN_branch_integrated_high_var.rda')

# Remove all previously loaded functions
rm(list=lsf.str())
# Source updated functions
sapply(list.files(path='./R/',pattern='.R', full.names = TRUE),source,.GlobalEnv)

skin@assays$ATAC@motifs # This is the motif object
DefaultAssay(skin)='ATAC'
granges(skin)
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.positions <- matchMotifs(
  pwms = pwm_set,
  subject = granges(skin),
  out = 'positions',
  genome = 'mm10'
)

motif.matrix <- CreateMotifMatrix(features = granges(skin), pwm = pwm_set, genome = 'mm10', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set, positions=motif.positions)
motif.object 
skin <- SetAssayData(skin, assay = 'ATAC', slot = 'motifs', new.data = motif.object)

# Need to re-define the directory of the fragment file
tmp <- Fragments(skin)[[1]]
tmp<-UpdatePath(tmp, new.path="../../Multiome_Shared/data/share_seq/skin/skin.atac.fragments.tsv.gz")
skin <- SetAssayData(skin, slot = "fragments", new.data = tmp)
rm(tmp); rm(motif.matrix); rm(motif.object); rm(pwm_set)

# Footprinting 
skin=Footprint(object = skin, genome=BSgenome.Mmusculus.UCSC.mm10,
               assay='ATAC', motif.name=toupper(motifxTF[,2])) # These are all the TFs, which takes quite long

save.image(file='ATAC_RNA_WNN_footprint.rda')
