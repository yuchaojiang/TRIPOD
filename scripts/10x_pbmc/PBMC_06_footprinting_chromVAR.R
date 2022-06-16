
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)

library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

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

# Remove all previously loaded functions
rm(list=lsf.str())
# Source updated functions
sapply(list.files(path='./R/',pattern='.R', full.names = TRUE),source,.GlobalEnv)

load('ATAC_RNA_WNN.rda') 

VlnPlot(pbmc, features = 'rna_CCR7', slot = "counts", log = TRUE, group.by = 'citeseq.celltype') + ggtitle('TF Expression')+ylab('Expression')
# This is chromVAR deviation scores
VlnPlot(pbmc, features = motifxTF[motifxTF[,2]=='CEBPB',1], slot = "data", assay='chromvar', group.by = 'citeseq.celltype') + ggtitle('chromVAR motif accessibility deviation')+ylab('Deviation')

pbmc@assays$ATAC@motifs # This is the motif object
DefaultAssay(pbmc)='ATAC'
granges(pbmc)
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.positions <- matchMotifs(
  pwms = pwm_set,
  subject = granges(pbmc),
  out = 'positions',
  genome = 'hg38'
)

motif.matrix <- CreateMotifMatrix(features = granges(pbmc), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set, positions=motif.positions)
motif.object # A Motif object containing 633 motifs in 105913 regions
pbmc <- SetAssayData(pbmc, assay = 'ATAC', slot = 'motifs', new.data = motif.object)

# Need to re-define the directory of the fragment file
tmp <- Fragments(pbmc)[[1]]
tmp<-UpdatePath(tmp, new.path="/proj/yuchaojlab/yuchaoj/multiomics/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz")
pbmc <- SetAssayData(pbmc, slot = "fragments", new.data = tmp)
pbmc=Footprint(object = pbmc, genome=BSgenome.Hsapiens.UCSC.hg38,
          assay='ATAC', motif.name=motifxTF[,2]) # These are all the TFs, which takes quite long
rm(tmp); rm(motif.matrix); rm(motif.object); rm(pwm_set)

PlotFootprint(pbmc, features = c("EBF1"), group.by='citeseq.celltype') |
  PlotFootprint(pbmc, features = c("TCF4"), group.by='citeseq.celltype') 

save.image(file='ATAC_RNA_WNN_footprint.rda')



