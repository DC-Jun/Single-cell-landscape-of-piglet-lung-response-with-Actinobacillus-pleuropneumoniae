
##differentially expressed genes (DEG)
##For example neutrophil cell

library(clusterProfiler)
library(ggplot2)
library(Seurat)
library(SeuratData)
ACCumap <- LoadData("ACCumap")

neu <- subset(ACCumap, new.cluster=="8")
diff_neu <- FindMarkers(neu, 
                              min.pct = 0.25,logfc.threshold = 0.25,group.by = "orig.ident",ident.1 ="APP",ident.2 ="CON")
write.csv(diff_neu ,file = "/neu_DEG.csv")

## GO analysis
diff_neu<- read.csv("/neu_DEG.csv")
ac <- subset(diff_neu,avg_logFC>0);
gene <- bitr(ac$X, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Ss.eg.db")
geneid <- gene$ENTREZID
go <- enrichGO(geneid, OrgDb = org.Ss.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
go <- clusterProfiler::simplify(go, cutoff=0.7, by="p.adjust", select_fun=min) 

barplot(go,showCategory=30,drop=T,order= T,label_format = 100) +scale_fill_continuous(low = "#03899C", high = "#FFD373", name = "p.adjust")
+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave('neu_GO.pdf',w=9.6,h=7.3)


