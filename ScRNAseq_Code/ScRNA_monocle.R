
## Monocle analysis

options(stringsAsFactors = F)
library(tibble)
library(ggplot2)
library(Seurat)
library(cowplot)
library(dplyr)
library(monocle)

subset.merge <- readRDS("merge_All.rds")
## Construct the cds object
exp <- subset.merge@assays$RNA@counts
gene <- data.frame(gene_short_name = rownames(exp), row.names=rownames(exp))
sample <- subset.merge@meta.data 

CDS <- newCellDataSet(
    as.matrix(exp), 
    phenoData = new('AnnotatedDataFrame', data = sample), 
    featureData = new('AnnotatedDataFrame', data = gene)
  )
pData(CDS) %>% head()

CDS <- estimateSizeFactors(CDS)
CDS <- estimateDispersions(CDS) 

## choose gene
disp_table <- dispersionTable(CDS)
ordering_genes <- subset(disp_table,mean_expression>=0.01&dispersion_empirical>=1*dispersion_fit)$gene_id
CDS <- setOrderingFilter(CDS, ordering_genes)
plot_ordering_genes(CDS)
  

CDS <- reduceDimension(CDS, max_components = 2, method = 'DDRTree') 
CDS <- orderCells(CDS, reverse = T)
saveRDS(CDS,"monocle_cells.rds")

 # plot  
  plot_cell_trajectory(CDS, color_by = "State", cell_size =0.5)+scale_color_manual(values = c("#FA7850",  "#0AB45A", "#0072B2",  
    "#FA78FA"))   
 ggsave('mono.cells_state.pdf', w=3.68,h=3.06)    
 
  plot_cell_trajectory(CDS, color_by = "orig.ident", cell_size = 0.5) +
    facet_wrap(~orig.ident, nrow = 1) +
    theme(text = element_text(size=25))+scale_color_manual(values = c("#ff2d51", "#455cc4"))   
  ggsave('mono.cells_group.pdf', w=5.51,h=3.56)
  
  plot_cell_trajectory(CDS, color_by = "rename", cell_size = 0.5) +scale_color_manual(values = c("#FFA500","#CC79A7","#00A0FA"))
   ggsave('mono.cells_rename.pdf', w=3.68,h=3.06)
   
    ##show_branch_points = FALSE
  plot_cell_trajectory(CDS, color_by = "TissueType", cell_size = 0.5) +
    facet_wrap(~TissueType, nrow = 2) +
    theme(text = element_text(size=25))+scale_color_manual(values = c("#FF9400", "#4790F5",  "#47BD3A", "#F553BA",  "#2B7E94","#952B6E"))
  ggsave('mono.cells_tissuetype.pdf', w=8.27,h=5.17)                               
  
  plot_cell_trajectory(CDS, color_by = "Pseudotime",cell_size = 0.5,show_branch_points = FALSE) + 
    theme(text = element_text(size = 25))
  ggsave('mono.cells_pseudotime.pdf', w=5.08,h=4.5)
  
 
##BEAM(Branched expression analysis modeling)

disp_table <- dispersionTable(CDS)
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
CDS_sub <- CDS[disp.genes,]
plot_cell_trajectory(CDS_sub, color_by = "State")
beam_res <- BEAM(CDS_sub, branch_point = 1, cores = 8)

beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
write.csv(beam_res,"beam_res_siggene.csv",row.names = F)

CDS_sub_beam <- CDS_sub[row.names(subset(beam_res, qval < 1e-4)),]
p5 <- plot_genes_branched_heatmap(CDS_sub_beam,  branch_point = 1, num_clusters = 3, show_rownames = T,
      branch_colors = c("#FA7850",  "#0AB45A", "#0072B2", "#FA78FA"))  
ggsave('BEAM analysis.pdf', w=7.61,h=8.81)    
    