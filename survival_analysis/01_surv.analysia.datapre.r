rm(list=ls())
#LIB------------------------------------------------------
library(magrittr)
library(ggplot2)
library(dplyr)

#PATH----------------------------------------------------------

data_path <- "/home/fux/fux/data/LIHC_prognosis"
out_path <- "/home/fux/fux/github/LIHC_prognosis/survival_analysis"


#INFILE------------------------------------------------------
clinical <- readr::read_tsv(file.path(data_path,"clinical_information_for_survival"))%>%
  select("gene_id","Months","dead/alive")
colnames(clinical) <- c("Patient_id","Time","Status")
mir <- readr::read_tsv(file.path(data_path,"zong_desease_mirna_normalized_rawcount.idchange.txt"))
rna <- readr::read_tsv(file.path(data_path,"zong_desease_rna_normalized_rawcount.xls"))

#Make survival metrix---------------------------------------------------- 

# for miRNA

fun_barcode <- function(.b){
  stringr::str_sub(
    string = .b,
    start = 1,
    end = 12
  )
} #get short barcode from long barcode

fun_tn_type <- function(.b){
  type <- .b %>% 
    stringr::str_split(pattern = "\\.", simplify = T) %>% 
    .[, 4] %>% 
    stringr::str_sub(1, 2)
} # get tumor and normal info

mir%>%
   filter(!is.na(mir$gene_id),mir$gene_id!="1.62975")->mir_01

   as.data.frame(t(data.frame(mir_01,row.names = 1)))->mir_02
   
mir_02%>%
   mutate(barcode=rownames(mir_02))%>%
   mutate(type=fun_tn_type(.$barcode))%>%
   filter(type=="01")%>%
   mutate(Patient_id=fun_barcode(.$barcode))%>%
   select(-barcode,-type)->mir_03

mir_03$Patient_id <- stringr::str_replace_all(mir_03$Patient_id,"\\.","-")

clinical%>%
  inner_join(mir_03)->surv_matrix

readr::write_tsv(surv_matrix,file.path(out_path,"surv_matrix_mirna.tsv"))


#for rna
as.data.frame(t(data.frame(rna,row.names = 1)))->rna_01
rna_01%>%
  mutate(Patient_id=row.names(rna_01))->rna_02

rna_02%>%
  mutate(type=fun_tn_type(Patient_id))%>%
  filter(type =="01")%>%
  select(-type)->rna_03

rna_03$Patient_id <- fun_barcode(rna_03$Patient_id)
rna_03$Patient_id <- stringr::str_replace_all(rna_03$Patient_id,"\\.","-")

clinical%>%
  inner_join(rna_03)->rna_matrix

readr::write_tsv(rna_matrix,file.path(out_path,"surv_matrix_rna.tsv"))
