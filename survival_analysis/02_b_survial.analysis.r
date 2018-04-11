rm(list = ls())
update.packages("multidplyr")
install.packages("devtools")
library(devtools)
devtools::install_github("hadley/multidplyr")
library(multidplyr)

data_path <- "/home/fux/fux/github/LIHC_prognosis/survival_analysis"


row_mir <- readr::read_rds(file.path(data_path,"row_miRNA.rds"))
row_rna <- readr::read_rds(file.path(data_path,"row_rna.rds"))
clinical <- readr::read_tsv(file.path(data_path,"clinical_information_for_survival"))
#exclude gene/mir number(expression==0)>0.5
colnames(row_mir)[1] <- "symbol"
colnames(row_rna)[1] <- "symbol"


expr_clinical_ready%>%
  count(symbol)

fun_barcode <- function(.b){
  stringr::str_sub(
    string = .b,
    start = 1,
    end = 12
  )
} #get short barcode from long barcode

fun_tn_type <- function(.b){
  type <- .b %>% 
    stringr::str_split(pattern = "-", simplify = T) %>% 
    .[, 4] %>% 
    stringr::str_sub(1, 2)
} # get tumor and normal info

fun_symbol_brif <- function(.b){
  stringr::str_split(string = .b,
                     pattern = "\\|",
                     simplify = T)[,1]->symbol
  return(symbol)
}

row_mir$symbol <- fun_symbol_brif(row_mir$symbol)
row_rna$symbol <- fun_symbol_brif(row_rna$symbol)
expr_clinical_sig_pval_rna$symbol <- fun_symbol_brif(expr_clinical_sig_pval_rna$symbol)


fun_expr_survival_merge <- function(filter_expr,clinical){
  # merge clinical and expr data
  filter_expr %>% 
  #  dplyr::select(-entrez_id) %>% 
    tidyr::gather(key = barcode, value = expr, -symbol) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type == "01") %>% 
    dplyr::mutate(barcode = fun_barcode(barcode)) %>% 
    dplyr::select(symbol, barcode, expr)  %>% 
    dplyr::inner_join(clinical, by = c("barcode" = "gene_id")) %>% 
    dplyr::select(symbol, barcode, expr, time = days, status = vital_status) %>% 
    dplyr::filter(!is.na(time), time > 0, !is.na(status)) %>% 
    dplyr::mutate(status = plyr::revalue(status, replace = c("Alive" = 0, "Dead" = 1))) %>%
    dplyr::mutate(status = as.numeric(status)) %>% 
    dplyr::mutate(expr = log2(expr + 1)) %>% 
    tidyr::drop_na(expr) %>% 
    dplyr::group_by(symbol) %>% 
    dplyr::mutate(group = as.factor(ifelse(expr <= median(expr),"Low", "High"))) %>% 
    dplyr::ungroup() -> expr_clinical_ready
} 


#--------------------------------------------------------------------------------------------------
  
fun_draw_survival <- function(symbol,p.value,expr_clinical_ready){
  gene <- symbol
  p_val <- signif(-log10(p.value), digits = 3)
  fig_name <- paste(gene, "pdf", sep = ".")
  # print(fig_name)
  .d <- 
    expr_clinical_ready %>% 
    dplyr::filter(symbol == gene)
  
  .d_diff <- survival::survdiff(survival::Surv(time, status) ~ group, data = .d)
  
  kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.d$group))) - 1)
  
  if(kmp > 0.05) {return(NA)} else{
    fit_x <- survival::survfit(survival::Surv(time, status) ~ group, data = .d , na.action = na.exclude)
    survminer::ggsurvplot(fit_x, data = .d, pval=T, pval.method = T,
                          title = paste(paste(gene, sep = "-"), "Coxph =", signif(p.value, 2)),
                          xlab = "Survival in days",
                          ylab = 'Probability of survival')
    ggsave(filename = fig_name, device = "pdf", path = file.path(data_path, "boxplot_rna"), width = 6, height = 6)
  }
}


fun_clinical_test <- function(expr_clinical_ready){
  if(nrow(expr_clinical_ready) < 1){return(tibble::tibble())}
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
    purrr::pwalk(fun_draw_survival, expr_clinical_ready = expr_clinical_ready) 
  
  return(d)
}


#mir--------------------------------------------------------------------------------------------
row_mir%>%
  fun_expr_survival_merge(clinical)%>%
  fun_clinical_test()->expr_clinical_sig_pval

readr::write_csv(expr_clinical_sig_pval,file.path(data_path,"surv_coxph_mir.csv"))

#rna----------------------------------------------------------------
row_rna%>%
  fun_expr_survival_merge(clinical)%>%
  fun_clinical_test()->expr_clinical_sig_pval_rna
readr::write_csv(expr_clinical_sig_pval_rna,file.path(data_path,"surv_coxph_rna.csv"))

#cl <- parallel::detectCores()
#cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
#row_mir %>%
#  multidplyr::partition(cluster = cluster) %>%
#  multidplyr::cluster_library("magrittr") %>%
#  multidplyr::cluster_library("ggplot2") %>%
#  multidplyr::cluster_assign_value("fun_barcode", fun_barcode)  %>%
#  multidplyr::cluster_assign_value("fun_tn_type", fun_tn_type) %>% 
#  multidplyr::cluster_assign_value("fun_expr_survival_merge", fun_expr_survival_merge) %>% 
#  multidplyr::cluster_assign_value("fun_clinical_test", fun_clinical_test) %>% 
#  multidplyr::cluster_assign_value("fun_draw_survival", fun_draw_survival) %>% 
#  multidplyr::cluster_assign_value("surv.mir_path", surv.mir_path) %>%
#  multidplyr::cluster_assign_value("row_mir",row_mir)%>%
#  multidplyr::cluster_assign_value("clinical",clinical)%>%
 # dplyr::mutate(merged_clean = purrr::map2(row_mir, clinical, fun_expr_survival_merge)) %>%
 # fun_expr_survival_merge(clinical)%>%
  #dplyr::select(-row_mir, -clinical) %>%
  #dplyr::mutate(diff_pval = purrr::map2(merged_clean,fun_clinical_test)) %>%
  #dplyr::select(-merged_clean) %>%
#  fun_clinical_test()%>%
#  dplyr::collect() %>%
#  dplyr::as_tibble() %>%
#  dplyr::ungroup() %>%
#  dplyr::select(-PARTITION_ID) %>% 
#  tidyr::unnest(diff_pval) -> expr_clinical_sig_pval
#on.exit(parallel::stopCluster(cluster))

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
  ggplot(aes(x = "LIHC",y = symbol, color = status)) +
  geom_point(aes(size = -log10(p.value))) +
#  scale_x_discrete(limit = cancer_rank$cancer_types) +
 # scale_y_discrete(limit = gene_rank$symbol) +
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
 #   axis.text.y = element_text(color = gene_rank$color),
    
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black")
  ) +
  ggthemes::scale_color_gdocs(name = "Surivival Worse")-> p

#ggsave(
 # filename = "fig_03_d_survival_sig_genes_serminar.pdf",
#  plot = p,
#  device = "pdf",
#  width = 14,
#  height = 9,
#  path = survival_path
#)

#genelist-------------------------------------
  