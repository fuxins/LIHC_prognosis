library(dplyr)


data_path <- "/home/fux/fux/github/LIHC_prognosis/survival_analysis"
diff_path <- "/home/fux/fux/github/LIHC_prognosis/differential_expression"

diff.mir <- readr::read_rds(file.path(data_path,"row_miRNA.rds"))

diff.mir_01 <- as.data.frame(t(data.frame(diff.mir,row.names = 1)))


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
#-------------------------------------------------------------------------------------------
readr::write_rds(diff.mir_02,file.path(diff_path,"miRNA_pair49.rds"))
#-------------------------------------------------------------------------------------------

diff.mir_03 <- reshape2::melt(diff.mir_02,id=c("sample","type")) 



#t-test
diff.mir_03%>%
  group_by(sample)%>%
  tidyr::drop_na(value)%>%
  do(broom::tidy(t.test(value ~ type,data=.))) %>%
  ungroup()%>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(variable,p.value, fdr) -> df_pvalue


#------------------------------------------------------------------------------------

