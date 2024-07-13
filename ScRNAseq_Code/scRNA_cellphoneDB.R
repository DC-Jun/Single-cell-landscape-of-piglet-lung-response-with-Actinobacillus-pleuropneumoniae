

####Cellphonedb_analysis
##APP
write.table(as.matrix(APP@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(APP@meta.data), APP@meta.data[,'new.rename', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown"
write.table(meta_data, 'cellphonedb_APPmeta.txt', sep='\t', quote=F, row.names=F)

###Control
write.table(as.matrix(CON@assays$RNA@data), 'cellphonedb_CONcount.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(CON@meta.data), CON@meta.data[,'new.rename', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown"
write.table(meta_data, 'cellphonedb_CONmeta.txt', sep='\t', quote=F, row.names=F)

###APP cellphoneDB
cellphonedb method statistical_analysis  cellphonedb_APP1021meta.txt  cellphonedb_APP1021count.txt      --counts-data=gene_name  
cellphonedb plot dot_plot 
cellphonedb plot heatmap_plot cellphonedb_APP1021meta.txt   
tree out/


##########Control cellphoneDB
cellphonedb method statistical_analysis  cellphonedb_CON1021meta.txt  cellphonedb_CON1021count.txt      --counts-data=gene_name 
cellphonedb plot dot_plot 
cellphonedb plot heatmap_plot cellphonedb_CON1021meta.txt   
tree out/