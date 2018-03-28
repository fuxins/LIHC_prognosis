
list_path <- "/data/fux/github/LIHC_prognosis/survival_analysis"
expr_path <- "/data/fux/github/LIHC_prognosis/expre"

symbol <- readr::read_tsv(file.path(list_path,"surv_gene_01"))
mirlist <- readr::read_tsv(file.path(list_path,"mirlist.txt"))
  mirlist[30,1] <- colnames(mirlist)
  colnames(mirlist) <- "mir_id"
  
mir_matrix <- readr::read_tsv(file.path(list_path,"surv_matrix_mirna.tsv"))
rna_matrix <- readr::read_tsv(file.path(list_path,"surv_matrix_rna.tsv"))

mir_matrix%>%
  select(Patient_id,mirlist$mir_id)->mir_expr

rna_matrix%>%
  select(Patient_id,symbol$'1')->rna_expr


