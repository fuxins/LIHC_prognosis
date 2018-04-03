row_rna_path <- "/home/fux/fux/lihc/LIHC/gdac.broadinstitute.org_LIHC.mRNAseq_Preprocess.Level_3.2016012800.0.0"
row_mir_path <- "/home/fux/fux/lihc/LIHC/gdac.broadinstitute.org_LIHC.miRseq_Mature_Preprocess.Level_3.2016012800.0.0"
out_path <- "/home/fux/fux/github/LIHC_prognosis/survival_analysis"

row_rna <- readr::read_tsv(file.path(row_rna_path,"LIHC.uncv2.mRNAseq_raw_counts.txt"))
row_mir <- readr::read_tsv(file.path(row_mir_path,"LIHC.miRseq_mature_RPM.txt"))


filter_valid_data <- function(obs)
{
#  obs_0 <- ifelse(is.na(obs),0,obs)#缺失值太多了
  zero_rate <- sum(obs_0==0,na.rm = T)/length(obs_0)
  score <- ifelse(zero_rate>0.5,0,1)
  return(score)
}

val_rna <- mutate(row_rna,score=0)

for(i in 1:nrow(row_rna)){
  val_rna[i,]$score <- filter_valid_data(row_rna[i,])
}

val_rna%>%
  filter(score==1)%>%
  select(-score)->val_rna_01

val_mir <- mutate(row_mir,score=0)

for (i in 1:nrow(row_mir)){
  val_mir[i,]$score <- filter_valid_data(row_mir[i,])
}


val_mir%>%
  filter(score==1)%>%
  select(-score)->val_mir_01

readr::write_rds(val_mir_01,file.path(out_path,"row_miRNA.rds"))
readr::write_rds(val_rna_01,file.path(out_path,"row_rna.rds"))
