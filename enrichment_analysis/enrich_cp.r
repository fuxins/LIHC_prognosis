
enrich_path <- "/home/fux/fux/github/LIHC_prognosis/enrichment_analysis"


symbol <- readr::read_rds(file.path(enrich_path,"symbol.rds"))


source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")
biocLite("topGO")
biocLite("org.Hs.eg.db")
biocLite("Rgraphviz")

library(clusterProfiler)

library(org.Hs.eg.db)
library(topGO)
library(Rgraphviz)


stringr::str_split(string = symbol$`1`,
                   pattern = "\\|",
                   simplify = T)[,1]->surv.rna_entr_03


surv.rna_entr <- bitr(expr$symbol,fromType = "SYMBOL",toType = "ENTREZID",
     OrgDb="org.Hs.eg.db")

head(surv.rna_entr)

genelist <- as.character(surv.rna_entr$ENTREZID)

#GO----------------------------------------------------------------------------------
ego <- enrichGO( gene = as.character(surv.rna_entr$ENTREZID),
                OrgDb="org.Hs.eg.db",
                keyType = "ENTREZID",
                pvalueCutoff = 0.05,
                readable=TRUE)


dotplot(ego, showCategory=30)
enrichMap(ego, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
cnetplot(ego1, foldChange=surv.rna_entr$ENTREZID)
plotGOgraph(ego)
#cnetplot(ego, foldChange=surv.rna_entr$ENTREZID)

ego1 <- setReadable(ego, OrgDb = org.Hs.eg.db)

write.csv(as.data.frame(ego),file.path(enrich_path,"G-enrich.csv"),row.names =F)

#KEGG--------------------------------------------------------------------------------------
ekk <- enrichKEGG(genelist, organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", 
                 qvalueCutoff=0.5)


dotplot(ekk, showCategory=30)
enrichMap(ekk, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
write.csv(as.data.frame(ekk),file.path(enrich_path,"KEGG-enrich.csv"),row.names =F)

head(ekk)
str(ekk)

#GSEA--------------------------------------------------------------------------------
gsecc <- gseGO(genelist, ont="CC", OrgDb=org.Hs.eg.db, verbose=F)


