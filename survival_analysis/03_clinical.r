
library(dplyr)

clinical_path <- "/home/fux/fux/lihc/LIHC/gdac.broadinstitute.org_LIHC.Clinical_Pick_Tier1.Level_4.2016012800.0.0"
mir_matrix_path <- "/data/fux/github/LIHC_prognosis/survival_analysis"

clinical <- readr::read_tsv(file.path(clinical_path,"LIHC.clin.merged.picked.txt"))
mir_matrix <- readr::read_tsv(file.path(mir_matrix_path,"surv_matrix_mirna.tsv"))

mir_matrix$Patient_id
upper_barcode <- toupper(colnames(clinical))
upper_barcode[1] <- colnames(clinical)[1]

colnames(clinical) <- upper_barcode

readr::write_rds(clinical,file.path(mir_matrix_path,"03_clinical.rds"))

                 