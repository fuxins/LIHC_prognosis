library(dplyr)

data_path <- "/home/fux/fux/github/LIHC_prognosis/survival_analysis"
diff_path <- "/home/fux/fux/github/LIHC_prognosis/differential_expression"

diff.mir <- readr::read_rds(file.path(data_path,"row_miRNA.rds"))
diff.rna <- readr::read_rds(file.path(data_path,"row_rna.rds"))

diff.mir_01 <- as.data.frame(t(data.frame(diff.mir,row.names = 1)))
diff.rna_01 <- as.data.frame(t(data.frame(row_rna,row.names = 1)))

calculate_fC_P <- function(df){
  
  sample <- df%>%
    mutate(sample = stringr::str_sub(
      string = row.names(df),
      start = 1,
      end = 12
    ),
    type = stringr::str_split(string = row.names(df),
                              pattern = "\\.",
                              simplify = T)[,4]%>%
      stringr::str_sub(1,2))%>%
   
    filter(type %in% c("01","11"))%>%
    group_by(sample)%>%
    filter(n()==2,length(unique(type))==2)%>%
    ungroup()
  
}

diff.mir_02 <- calculate_fC_P(diff.mir_01)



#for rna-------------------------------------------------------------------------------------------------

diff.rna_02 <- calculate_fC_P(diff.rna_01)#Something wrong with it,no idea

diff.rna_02 <- diff.rna_01%>%
  mutate(sample = stringr::str_sub(
    string = row.names(diff.rna_01),
    start = 1,
    end = 12
  ),
  type = stringr::str_split(string = row.names(diff.rna_01),
                            pattern = "\\.",
                            simplify = T)[,4]%>%
    stringr::str_sub(1,2))

diff.rna_02%>%
  filter(type %in% c("01","11"))%>%
  group_by(sample)%>%
  filter(n()==2,length(unique(type))==2)%>%
  ungroup()->diff.rna_03

#-------------------------------------------------------------------------------------------
readr::write_rds(diff.mir_02,file.path(diff_path,"miRNA_pair49.rds"))
readr::write_rds(diff.rna_03,file.path(diff_path,"RNA_pair50.rds"))
#-------------------------------------------------------------------------------------------

diff.mir_03 <- reshape2::melt(diff.mir_02,id=c("sample","type"))
diff.rna_04 <- reshape2::melt(diff.rna_03,id=c("sample","type"))

#miRNA
diff.mir_03 %>%
  tidyr::drop_na() %>%
  group_by(sample,variable)%>%
  filter(n()==2,length(unique(type)==2))%>%
  ungroup()->diff.mir_05



diff.mir_05%>%
  group_by(variable)%>%
  do(broom::tidy(
    
    tryCatch(
      t.test(value ~ type, data = .),
      error = function(e) {1}
    )
    
    )) %>% 
  dplyr::ungroup() %>% 
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(variable, p.value, fdr)%>%
  filter(p.value<0.05)-> df_pvalue

#RNA
diff.rna_04%>%
  tidyr::drop_na() %>%
  group_by(sample,variable)%>%
  filter(n()==2,length(unique(type)==2))%>%
  ungroup()%>%
  group_by(variable)%>%
  do(broom::tidy(
    
    tryCatch(
      t.test(value ~ type, data = .),
      error = function(e) {1}
    )
    
  )) %>% 
  dplyr::ungroup() %>% 
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(variable, p.value, fdr)%>%
  filter(p.value<0.05)-> df_pvalue_rna
  
df_pvalue_rna%>%
  filter(p.value<0.001)->df_pvalue_rna_01
#-----------------------------------------------------------------------------------------
readr::write_rds(df_pvalue,file.path(diff_path,"miRNA_267.rds"))
readr::write_rds(df_pvalue_rna_01,file.path(diff_path,"RNA_9001.rds"))
#t-test use formular---------------------------------------------------------------------------
calculate_fC_P_01 <- function(df){
  df%>%
    filter(is.na(value))->sample_na
  df%>%
    filter(sample %in% sample_na$sample)->sample_na_01
  df%>%
    anti_join(sample_na_01)->df_01
  return(df_01)
}

diff.mir_04 <- diff.mir_03%>%
  filter(sample=="0")

for (var in unique(diff.mir_03$variable)) {
  diff.mir_03%>%
    filter(variable == var)->df
  diff.mir_04%>%
    full_join(calculate_fC_P_01(df))->diff.mir_04
  warnings('off')
}



#------------------------------------------------------------------------------------



