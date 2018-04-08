
enrich_path <- "/home/fux/fux/github/LIHC_prognosis/enrichment_analysis"


symbol <- readr::read_rds(file.path(enrich_path,"symbol.rds"))


source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")

biocLite("org.Hs.eg.db")

library(clusterProfiler)

library(org.Hs.eg.db)




stringr::str_split(string = symbol$`1`,
                   pattern = "\\|",
                   simplify = T)[,1]->surv.rna_entr_03


surv.rna_entr <- bitr(surv.rna_entr_03,fromType = "SYMBOL",toType = "ENTREZID",
     OrgDb="org.Hs.eg.db")

head(surv.rna_entr)

genelist <- as.character(surv.rna_entr$ENTREZID)

#GO----------------------------------------------------------------------------------
ego <- enrichGO( gene = as.character(surv.rna_entr$ENTREZID),
                OrgDb="org.Hs.eg.db",
                keyType = "ENTREZID",
                pvalueCutoff = 1,
                readable=TRUE)


dotplot(ego, showCategory=30)
enrichMap(ego, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
#cnetplot(ego, foldChange=surv.rna_entr$ENTREZID)

ego1 <- setReadable(ego, OrgDb = org.Hs.eg.db)

write.csv(as.data.frame(ego),file.path(enrich_path,"G-enrich.csv"),row.names =F)

#KEGG--------------------------------------------------------------------------------------
ekk <- enrichKEGG(genelist, organism="hsa", pvalueCutoff=1, pAdjustMethod="BH", 
                 qvalueCutoff=0.5)

dotplot(ekk, showCategory=30)

write.csv(as.data.frame(ekk),file.path(enrich_path,"KEGG-enrich.csv"),row.names =F)

head(ekk)
str(ekk)

#GSEA--------------------------------------------------------------------------------
gsecc <- gseGO(genelist, ont="CC", OrgDb=org.Hs.eg.db, verbose=F)


