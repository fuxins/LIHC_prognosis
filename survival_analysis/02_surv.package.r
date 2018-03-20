library(survival)
library(ggplot2)

#Path-----------------------------
surv_path <- "/home/fux/fux/github/LIHC_prognosis/survival_analysis"
out_path <- "/home/fux/fux/github/LIHC_prognosis/survival_analysis/jpeg"

#INFILE----------------------------------------------
mir_matrix <- readr::read_tsv(file.path(surv_path,"surv_matrix_mirna.tsv"))


#-------------------------------------------------------
fun_sa.matrix <- function(i,sa_matrix){
  expr_mean <- mean(as.numeric(unlist(sa_matrix[,i])),na.rm = T)
  sa_matrix%>%
    select(Patient_id,Time,Status,i)%>%
    mutate(expr=ifelse(surv_matrix[,i]<expr_mean,"L","H"))->sa_matrix_01
  my.surv <- Surv(sa_matrix_01$Time,sa_matrix_01$Status)
  my.fit <- survfit(my.surv~sa_matrix_01$expr)
  my.diff <- survdiff(my.surv~sa_matrix_01$expr,rho = 1)
  p.value <- 1 - pchisq(my.diff$chisq,length(my.diff$n) - 1)
  if(p.value<=0.05){
   plot(my.fit,main = colnames(mir_matrix[,i]),col=c("red","blue"))
  }
    
}

setwd(out_path)

for (i in 4:length(mir_matrix)) {
  
  if(!is.na(sum(as.numeric(unlist(mir_matrix[,i]))))){
    jpeg(filename = colnames(mir_matrix[,i]))  
    fun_sa.matrix(i,sa_matrix = mir_matrix)
    dev.off()
  }
  
}


#TEST---------------------------------------------

expr_mean <- mean(unlist(mir_matrix[,i]))
mir_matrix%>%
  select(Patient_id,Time,Status,i)%>%
  mutate(expr=ifelse(mir_matrix[,i]<expr_mean,"L","H"))->mir_matrix_01
my.surv <- Surv(mir_matrix_01$Time,mir_matrix_01$Status)
my.fit <- survfit(my.surv~mir_matrix_01$expr)
my.diff <- survdiff(my.surv~mir_matrix_01$expr,rho = 1)
p.value <- 1 - pchisq(my.diff$chisq,length(my.diff$n) - 1)
if(p.value<=0.05){
  plot(my.fit,main = colnames(mir_matrix[,i]))
}

mir_matrix_01
expr_mean
mean(unlist(mir_matrix[,i]),na.rm = T)


      