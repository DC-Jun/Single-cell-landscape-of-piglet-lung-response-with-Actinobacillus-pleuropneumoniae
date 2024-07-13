

##AUCell score

library(ggraph)
library(ggplot2)
library(Seurat)
library(AUCell)
library(SeuratData)

ACCumap <- LoadData("ACCumap")
cells_rankings <- AUCell_buildRankings(ACCumap@assays$RNA@data) 
Hallmarker <- read.gmt("c5.go.bp.v7.5.1.symbols.gmt") 
geneSets <- lapply(unique(Hallmarker$term), function(x){print(x);Hallmarker$gene[Hallmarker$term == x]})
names(geneSets) <- unique(Hallmarker$term)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

##set gene set of interest here for plotting
geneSet <- "GOBP_INFLAMMATORY_RESPONSE"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
ACCumap$AUC  <- aucs

neu8 <- subset(ACCumap,new.cluster=="8")
VlnPlot(neu8,features = "AUC",pt.size = 0,group.by = "orig.ident",sort = 'decreasing',cols = c( "#455cc4","#ff2d51"))+stat_compare_means(method = "t.test")
  ggsave('AUcell_GOBP INFLAMMATORY RESPONSE_Neu.pdf',w=3.96,h=3.55) 
  
mono3<- subset(ACCumap,new.cluster=="3")  
VlnPlot(mono3,features = "AUC",pt.size = 0,group.by = "orig.ident",sort = 'decreasing',cols = c( "#455cc4","#ff2d51"))+stat_compare_means(method = "t.test")
  ggsave('AUcell_GOBP INFLAMMATORY RESPONSE_mono.pdf',w=3.96,h=3.55) 

AM4 <- subset(ACCumap,new.cluster=="4")  
VlnPlot(AM4,features = "AUC",pt.size = 0,group.by = "orig.ident",sort = 'increasing',cols = c(  "#455cc4","#ff2d51"))+stat_compare_means(method = "t.test")
  ggsave('AUcell_GOBP INFLAMMATORY RESPONSE_AM6.pdf',w=3.96,h=3.55) 

IM6  <- subset(ACCumap,new.cluster=="6")  
VlnPlot(IM6,features = "AUC",pt.size = 0,group.by = "orig.ident",sort = 'increasing',cols = c( "#455cc4","#ff2d51"))+stat_compare_means(method = "t.test")
  ggsave('AUcell_GOBP INFLAMMATORY RESPONSE_IM6 .pdf',w=3.96,h=3.55) 

pDC12  <- subset(ACCumap,new.cluster=="12")  
VlnPlot(pDC12,features = "AUC",pt.size = 0,group.by = "orig.ident",sort = 'decreasing',cols = c( "#455cc4","#ff2d51"))+stat_compare_means(method = "t.test")
  ggsave('AUcell_GOBP INFLAMMATORY RESPONSE_pDC.pdf',w=3.96,h=3.55) 

