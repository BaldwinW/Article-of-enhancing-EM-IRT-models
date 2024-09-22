###############################################
############## 2PL EM-IRT Model ###############
###############################################

## Date: 20230518

## Load Packages
Packages <- c("truncnorm", "mirt", "MASS", "mirtCAT", "stringr")

for (BAO in Packages) {
  if (BAO %in% rownames(installed.packages()) == FALSE){
    install.packages(BAO)
  }
  if (BAO %in% (.packages()) == FALSE){
    library(BAO, character.only = TRUE)
  }
}


## Main Function
RG_IRT_Model_EMIRT <- function(
  csv_table,
  Path_Save_Result){
  
  mydata <- read.csv(file = csv_table)
  
  responses <- mydata[1:2000, 1:50]
  theta_true <- mydata[1:2000, 51]
  a_true <- mydata[2001, 1:50]
  rownames(a_true) <- NULL
  b_true <- mydata[2002, 1:50]
  rownames(b_true) <- NULL
  
  tpl <- mirt(responses, 
              model=1, 
              itemtype='2PL', 
              SE=TRUE)
  
  itempars_2PL <- coef(tpl, 
                       printSE=TRUE, 
                       IRTpars=TRUE,
                       as.data.frame=TRUE)
  
  ## Extracting theta estimates
  theta_2PL <- fscores(tpl,
                       method="EAP",
                       full.scores=TRUE,
                       full.scores.SE=FALSE)
  
  ## Save Item parameters
  lis <- seq(1, by=4, length.out=50)
  a_est <- itempars_2PL[lis, 1]
  b_est <- itempars_2PL[lis+1, 1]

  theta_values <- cbind(theta_true, theta_2PL)
  names(theta_values) <- c("theta_true", "theta_est")
  write.csv(theta_values,
            file = paste0(Path_Save_Result, "/EMIRT/Theta/Theta_", substr(csv_table, str_length(csv_table)-49, str_length(csv_table)-4), ".csv"), 
            row.names = F)
  
  a_values <- data.frame(t(data.frame(rbind(a_true, a_est))))
  names(a_values) <- c("a_true", "a_est")
  write.csv(a_values, 
            file = paste0(Path_Save_Result, "/EMIRT/a/a_", substr(csv_table, str_length(csv_table)-49, str_length(csv_table)-4), ".csv"), 
            row.names = F)
  
  b_values <- data.frame(t(data.frame(rbind(b_true, b_est))))
  names(b_values) <- c("b_true", "b_est")
  write.csv(b_values, 
            file = paste0(Path_Save_Result, "/EMIRT/b/b_", substr(csv_table, str_length(csv_table)-49, str_length(csv_table)-4), ".csv"), 
            row.names = F)
}



