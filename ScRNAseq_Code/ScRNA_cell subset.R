
###T cell subsets analyzed
library(ggplot2)
library(Seurat)
library(SeuratData)

ACCumap <- LoadData("ACCumap")
t3 <- subset(ACCumap, new.cluster%in% c('2','11','16'))
t3 <- NormalizeData(t3, normalization.method = "LogNormalize", scale.factor = 1e4) 
t3 <- FindVariableFeatures(t3, selection.method = 'vst', nfeatures = 2000)
all.genes <- rownames(t3)
t3 <- ScaleData(t3, features = all.genes)
t3 <- RunPCA(t3, features = VariableFeatures(object = t3))
t3 <- JackStraw(t3, dims = 15, num.replicate = 100, maxit = 10000)
t3 <- ScoreJackStraw(t3, dims = 1:15)
JackStrawPlot(t3, dims = 1:15)
ElbowPlot(t3, ndims = 20, reduction = "pca")
t3 <- FindNeighbors(t3, dims = 1:10)
t3 <- FindClusters(t3, resolution = 0.3) 
t3umap <- RunUMAP(t3,min.dist =0.8, local.connectivity = 10L,  repulsion.strength = 10, negative.sample.rate = 5, 
n.epochs = 200, n.neighbors = 60L, dims = 1:10) 

DimPlot(object = t3umap, 
        reduction = 'umap',
        cols = c('#CC3399','#993366','#009933', '#0066CC','#003399','#9966CC','#99CC33','#003399',
        "#FAE6BE", "#0072B2", "#A0FA82", "#F0F032", "#0AB45A"),
        label = T, label.box = F,cols.highlight = "#DE2D26", pt.size = 0.5)
ggsave('t3_subset_umap.pdf',w=5.9,h=4.5)
        
t3_cell.markers <- FindAllMarkers(t3umap, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
write.csv(x = t3_0.3cell.markers , file = "t3_cell.markers.csv", quote = FALSE)
top5 <- t3_cell.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(t3, features = top5$gene, slot = 'scale.data',size = 3)
ggsave('t3_heatmap_markers.pdf',w=8,h=7)
saveRDS(t3, file = "/t3umap.rds")


#Stacked violin plot
my36colors <- c( '#C5DEBA', '#58A4C3', '#E4C755',
         '#AA9A59','#E39A35', '#625D9E','#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6')

modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
       p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
               xlab("") + ylab(feature) + ggtitle("") +
               theme(legend.position = "none",
               axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.x = element_blank(),
               axis.ticks.y = element_line(),
               axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
               plot.margin = plot.margin )
       return(p)}

## main function
StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
       plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
            plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
            theme(axis.text.x=element_text(), axis.ticks.x = element_line())
       p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
       return(p)}
       
StackedVlnPlot(t3umap, c('CD3E','CD8A','CD4','GNLY', 'NKG7','CXCR3','CXCR6','XCL1','EOMES','GATA3','TRDC','TOP2A','PCLAF','CCR7','IL7R','LEF1','GZMB','GZMM','MARCO'), pt.size=0, cols=my36colors)
ggsave('T cell subset violinPlot.pdf',w=6,h=20)


