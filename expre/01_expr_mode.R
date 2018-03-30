
list_path <- "/data/fux/github/LIHC_prognosis/survival_analysis"
expr_path <- "/data/fux/github/LIHC_prognosis/expre"

symbol <- readr::read_tsv(file.path(list_path,"surv_gene_01"))
mirlist <- readr::read_tsv(file.path(list_path,"mirlist.txt"))
  mirlist[30,1] <- colnames(mirlist)
  colnames(mirlist) <- "mir_id"
  
mir_matrix <- readr::read_tsv(file.path(list_path,"surv_matrix_mirna.tsv"))
rna_matrix <- readr::read_tsv(file.path(list_path,"surv_matrix_rna.tsv"))
clinical <- readr::read_rds(file.path(list_path,"03_clinical.rds"))

mir_matrix%>%
  select(Patient_id,mirlist$mir_id)->mir_expr

rna_matrix%>%
  select(Patient_id,symbol$'1')->rna_expr

#-------------------------------------------------------------------------------
#做热图分析表达量
#mirna
mir_expr%>%
  select(-`hsa-miR-151a-3p`,-`hsa-miR-210-3p`,-`hsa-miR-194-5p`,-`hsa-miR-21-5p`,-`hsa-miR-9-5p`)->mir_expr_01



mir_expr_01 <- as.matrix(mir_expr_01)
row.names(mir_expr_01) <- mir_expr_01[,1]
mir_expr_02 <- mir_expr_01[,-1]

class(mir_expr_02) <- "numeric"

library(gplots)
heatmap.2(mir_expr_02)

typeof(mir_expr_02)

#rna
rna_expr_01 <- as.matrix(rna_expr)
row.names(rna_expr_01) <- rna_expr_01[,1]
rna_expr_02 <- rna_expr_01[,-1]

class(rna_expr_02) <- "numeric"
rna_expr_02[,20:40]->rna_expr_03

typeof(rna_expr_02)
heatmap.2(rna_expr_03)
#-----------------------------------------------------------------------------
#用箱线图看mir表达量的分布
library(ggplot2)

reshape2::melt(mir_expr)->mir_box
mir_box%>%
  filter(value>10000)%>%
  count(variable)->high_mir
mir_box%>%
  filter(variable != "hsa-miR-21-5p")->mir_box_03

mir_box%>%
  filter(variable %in% high_mir$variable)->mir_box_01
mir_box%>%
  anti_join(mir_box_01)->mir_box_02

ggplot(mir_box_02,aes(x=variable,y=value))+
  geom_boxplot()+
  #geom_hline(aes(yintercept=10000),col="red")+
  xlab("miRNA")+
  ylab("expr")

#hsa-miR-210-3p在一个样本中表达量异常高
#mir_box中不包含21-5p,共29条mirna
mir_box%>%
  filter(variable =="hsa-miR-210-3p",value>10000)->mir_210_3p
mir_box%>%
  filter(Patient_id == mir_210_3p$Patient_id)%>%
  filter(value>5000)
  ggplot(.,aes(x=variable,y=value))+
  geom_bar(stat = "identity")+
  xlab("miRNA")+
  ylab("expr")+
  labs(title = "Patient:TCGA-CC-5263")
  
#---------------------------------------------------------------------------  
