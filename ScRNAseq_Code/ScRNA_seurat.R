
########## Seurat analysis  
options(stringsAsFactors = F)
library(dplyr)
library(tibble)
library(cowplot)
library(Seurat)
library(ggplot2)
library(tidyverse)

merge_All <- readRDS('seurat_process_files.rds')

###  Load datasets and create Seurat objects 
gf=gzfile('/APP/featuresV.tsv.gz','rt')
data<-read.table(gf)
ensemblQ=data$V2[!grepl('ENSSSCG', data$V2)]
counts_matrix = GetAssayData(merge_All, slot="counts")
merge_All_sub <- subset(merge_All, features = ensemblQ)
counts_matrixQ = GetAssayData(merge_All_sub, slot="counts")
table(grepl("^ENSSS",rownames(merge_All_sub)))
saveRDS(merge_All_sub, file = "/seurat_enss.rds")
merge_All<- merge_All_sub

rownames(merge_All)[grepl('^ND[1-6]', rownames(merge_All))]
rownames(merge_All)[grepl('^COX[1-3]$', rownames(merge_All))]
rownames(merge_All)[grepl('^ATP[68]$', rownames(merge_All))]
rownames(merge_All)[grepl('CYTB', rownames(merge_All))]

# Calculate mitochondrial results separately and sum
mt <- cbind(
    PercentageFeatureSet(merge_All, pattern = "^ND[1-6]"),
    cox <- PercentageFeatureSet(merge_All, pattern = "^COX[1-3]$"),
    atp <- PercentageFeatureSet(merge_All, pattern = "^ATP[68]"),
    cytb <- PercentageFeatureSet(merge_All, pattern = "CYTB")
  ) 
merge_All[["percent.nd"]] <- apply(mt, 1, sum)
head(merge_All@meta.data)

#Zone the selected cells
theme_set(theme_classic())
ggplot(merge_All@meta.data, aes(nFeature_RNA, percent.nd)) +
  geom_point(size = 0.8, colour = "black") +
  geom_hline(yintercept = 20, linetype = "dotdash", colour = "red",size = 1.2)+ #上界
  geom_hline(yintercept = 0, linetype = "dotdash", colour = "red",size = 1.2)+ #下界
  geom_vline(xintercept = 4000, linetype = "dotdash", colour = "red",size = 1.2)+ #右界
  geom_vline(xintercept = 300, linetype = "dotdash", colour = "red",size = 1.2) #左界

ggplot(merge_All@meta.data, aes(nCount_RNA, percent.nd)) +
  geom_point(size = 0.8, colour = "black") +
  geom_hline(yintercept = 5, linetype = "dotdash", colour = "red",size = 1.2)+ #上界
  geom_hline(yintercept = 0, linetype = "dotdash", colour = "red",size = 1.2)+ #下界
  geom_vline(xintercept = 50000, linetype = "dotdash", colour = "red",size = 1.2)+ #右界
  geom_vline(xintercept = 2000, linetype = "dotdash", colour = "red",size = 1.2) #左界

### filter low quality cells
merge_All <- subset(merge_All, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 10 )
merge_All <- NormalizeData(merge_All, normalization.method = "LogNormalize", scale.factor = 100000)
all.genes <- rownames(merge_All)
merge_All <- ScaleData(merge_All, features = all.genes)#如失败则改成features=NULL
merge_All <- FindVariableFeatures(merge_All, selection.method = "vst", nfeatures = 3000)
merge_All <- RunPCA(merge_All, features = VariableFeatures(object = merge_All))

ElbowPlot(ACC, ndims = 15, reduction = "pca")
merge_All <- FindNeighbors(merge_All, dims = 1:15)
merge_All <- FindClusters(merge_All, resolution = 0.2)

merge_All <- RunUMAP(merge_All, min.dist =0.8, local.connectivity = 10L,  repulsion.strength = 10, negative.sample.rate = 5, 
                   n.epochs = 200, n.neighbors = 60L, dims = c(1:15))

DimPlot(object = merge_All, 
                 reduction = 'umap',
               cols = c( "#006E82", "#AA0A3C", "#8214A0", "#00A0FA", "#FFA500","#ADFF2F",
                                "#CC79A7", "#20B2AA", "#AB82FF", "#FF34B3","#0072B2",  "#00CD00",  '#FF0033', "#CD5C5C","#008B8B",
                                             "#FA7850", "#14D2DC", "#FA78FA","#FA5078" ),
                 label = F, label.box = F,cols.highlight = "#DE2D26", pt.size = 0.1) +NoLegend()

###plot
 DimPlot(object = merge_All, 
                          reduction = 'umap',
                           cols = c( "#ff2d51","#455cc4"),
                           label = F, label.box = F,cols.highlight = "#DE2D26", pt.size = 0.1,group.by = 'orig.ident')+NoLegend()
ggsave('umap_orig.ident.pdf',w=6.2,h=5.7) 
saveRDS(merge_All,"merge_All.rds")

#cluster biomarkers
merge_All_markers <- FindAllMarkers(merge_All, only.pos = TRUE, min.pct = 0.6, logfc.threshold = 0.25)
write.csv(x = merge_All_markers, file = "merge_all_markers.csv", quote = FALSE)

##########
##Boxplot

library(ggplot2)
library(cowplot)
library(tidyverse)
library(ggsci)
library(ggpubr)
library(cowplot)

consistentcolors = colors <- c("#ff2d51","#455cc4")
df_all <- data.frame( TissueType = ACCumap_NEW$TissueType, orig.ident = ACCumap_NEW$orig.ident,
                      cell_all = ACCumap_NEW$rename)
df_summed <- df_all %>% group_by(orig.ident, TissueType, cell_all) %>% tally()
df_summed <-  df_summed %>% group_by(TissueType) %>% mutate(freq = n/sum(n))
head(df_summed)
df_summed$orig.ident <- factor(df_summed$orig.ident,levels=c("CON","APP"))

ggboxplot(df_summed, x = "cell_all", y = "freq", color = "orig.ident", add = "jitter") + 
   ylim(0, 0.85) + stat_compare_means(aes(group  = freq), label = "p", method = "t.test") + 
   theme(axis.text.x = element_text(angle = 90, hjust = 1))+ scale_colour_manual(values = consistentcolors)
ggsave(paste0("boxplot_frequency.pdf"), width =7.6, height = 5.5)


