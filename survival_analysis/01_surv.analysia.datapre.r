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


#fun_expr_survival_merge <- function(filter_expr, clinical){


fun_draw_survival <- function(symbol, p.value, cancer_types, expr_clinical_ready){
  gene <- symbol
  p_val <- signif(-log10(p.value), digits = 3)
  fig_name <- paste(cancer_types, gene, p_val, "pdf", sep = ".")
  # print(fig_name)
  .d <- 
    expr_clinical_ready %>% 
    dplyr::filter(symbol == gene)
  
  .d_diff <- survival::survdiff(survival::Surv(time, status) ~ group, data = .d)
  
  kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.d$group))) - 1)
  
  if(kmp > 0.05) {return(NA)} else{
    fit_x <- survival::survfit(survival::Surv(time, status) ~ group, data = .d , na.action = na.exclude)
    survminer::ggsurvplot(fit_x, data = .d, pval=T, pval.method = T,
                          title = paste(paste(cancer_types, gene, sep = "-"), "Coxph =", signif(p.value, 2)),
                          xlab = "Survival in days",
                          ylab = 'Probability of survival')
    ggsave(filename = fig_name, device = "pdf", path = file.path(survival_path, "boxplot"), width = 6, height = 6)
  }
}
fun_clinical_test <- function(expr_clinical_ready, cancer_types){
  if(nrow(expr_clinical_ready) < 1){return(tibble::tibble())}
  print(cancer_types)
  expr_clinical_ready %>% 
    dplyr::group_by(symbol) %>% 
    dplyr::do(
      broom::tidy(
        tryCatch(
          survival::coxph(survival::Surv(time, status) ~ expr, data = ., na.action = na.exclude),
          error = function(e){1},
          warning = function(e){1})
      )
    ) %>%
    dplyr::ungroup() %>% 
    dplyr::filter(p.value < 0.05) %>% 
    dplyr::select(symbol, estimate, p.value) %>% 
    dplyr::mutate(status = ifelse(estimate > 0, "H", "L"))-> d
  
  d %>% 
    dplyr::select(symbol, p.value) %>% 
    purrr::pwalk(fun_draw_survival, cancer_types = cancer_types, expr_clinical_ready = expr_clinical_ready) 
  
  return(d)
}

# expr_clinical %>%
#   dplyr::filter(cancer_types == "KIRC") %>% 
#   dplyr::mutate(merged_clean = purrr::map2(filter_expr, clinical, fun_expr_survival_merge)) %>%
#   dplyr::select(-filter_expr, -clinical) %>%
#   dplyr::mutate(diff_pval = purrr::map2(merged_clean, cancer_types, fun_clinical_test)) %>%
#   dplyr::select(-merged_clean) %>%
#   tidyr::unnest(diff_pval) -> expr_clinical_sig_pval

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
expr_clinical %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_library("ggplot2") %>%
  multidplyr::cluster_assign_value("fun_barcode", fun_barcode)  %>%
  multidplyr::cluster_assign_value("fun_tn_type", fun_tn_type) %>% 
  multidplyr::cluster_assign_value("fun_expr_survival_merge", fun_expr_survival_merge) %>% 
  multidplyr::cluster_assign_value("fun_clinical_test", fun_clinical_test) %>% 
  multidplyr::cluster_assign_value("fun_draw_survival", fun_draw_survival) %>% 
  multidplyr::cluster_assign_value("survival_path", survival_path) %>% 
  dplyr::mutate(merged_clean = purrr::map2(filter_expr, clinical, fun_expr_survival_merge)) %>%
  dplyr::select(-filter_expr, -clinical) %>%
  dplyr::mutate(diff_pval = purrr::map2(merged_clean, cancer_types, fun_clinical_test)) %>%
  dplyr::select(-merged_clean) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) %>% 
  tidyr::unnest(diff_pval) -> expr_clinical_sig_pval
on.exit(parallel::stopCluster(cluster))

#---------------------------------------------------------------------------------------------

fun_rank_cancer <- function(pattern){
  pattern %>% 
    dplyr::summarise_if(.predicate = is.numeric, dplyr::funs(sum(., na.rm = T))) %>%
    tidyr::gather(key = cancer_types, value = rank) %>%
    dplyr::arrange(dplyr::desc(rank))
} #get cancer rank
fun_rank_gene <- function(pattern){
  pattern %>% 
    dplyr::rowwise() %>%
    dplyr::do(
      symbol = .$symbol,
      rank =  unlist(.[-1], use.names = F) %>% sum(na.rm = T)
    ) %>%
    dplyr::ungroup() %>%
    tidyr::unnest() %>%
    dplyr::arrange(rank)
} # get gene rank

expr_clinical_sig_pval %>% 
  dplyr::select(cancer_types, symbol) %>% 
  dplyr::mutate(n = 1) %>% 
  tidyr::spread(key = cancer_types, value = n) -> pattern

cancer_rank <- pattern %>% fun_rank_cancer()
gene_rank <- 
  pattern %>% 
  fun_rank_gene() %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::filter( rank >= 5) %>% 
  dplyr::mutate(color = plyr::revalue(status, replace = c('a' = "#e41a1c", "l" = "#377eb8", "i" = "#4daf4a", "p" = "#984ea3"))) %>% 
  dplyr::arrange(color, rank)

expr_clinical_sig_pval %>% 
  ggplot(aes(x = cancer_types, y = symbol, color = status)) +
  geom_point(aes(size = -log10(p.value))) +
  scale_x_discrete(limit = cancer_rank$cancer_types) +
  scale_y_discrete(limit = gene_rank$symbol) +
  scale_size_continuous(name = "P-value") +
  theme(
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    ),
    
    axis.title = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text.y = element_text(color = gene_rank$color),
    
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black")
  ) +
  ggthemes::scale_color_gdocs(name = "Surivival Worse")-> p
ggsave(
  filename = "fig_03_d_survival_sig_genes_serminar.pdf",
  plot = p,
  device = "pdf",
  width = 14,
  height = 9,
  path = survival_path
)



save.image(file = file.path(survival_path, ".rda_03_d_survival_gene_expr.rda"))
load(file = file.path(survival_path, ".rda_03_d_survival_gene_expr.rda"))
