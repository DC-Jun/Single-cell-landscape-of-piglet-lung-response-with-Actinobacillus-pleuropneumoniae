
## Cyto TRACE
library(reticulate)
use_condaenv("cytoTRACE")
library(CytoTRACE)
library(SeuratData)

ifnb <- LoadData("ACCumap")
cds <- readRDS("./monocle_cells.rds")
sce <- ACCumap[,ACCumap$celltype %in% c("epithelial","fibroblast_like")]
mat<- as.matrix(sce@assays$RNA@counts)
results <- CytoTRACE(mat = mat)
phe <- sce$celltype
phe = as.character(phe)
names(phe) <- rownames(sce@meta.data)
plotCytoGenes(results, numOfGenes = 10)

mat<- as.matrix(sce@assays$RNA@counts)
results <- CytoTRACE(mat = mat)
phe <- sce$celltype
phe = as.character(phe)
names(phe) <- rownames(sce@meta.data)
plotCytoGenes(results, numOfGenes = 10)
cds$cytotrace=results[["CytoTRACE"]]
cols <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
pdf("cytotrace.pdf",6,6)
plot_cell_trajectory(cds,color_by="cytotrace", size=1,show_backbone=TRUE)+
  scale_colour_gradientn(colours = cols)
dev.off()