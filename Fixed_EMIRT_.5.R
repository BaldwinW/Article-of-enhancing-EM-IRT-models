###############################################
######   Fixed item parameters Model     ######
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
RG_IRT_Model_Fixed2PL <- function(
  csv_table,
  Path_Save_Result){
  
  mydata <- read.csv(file = csv_table)
  
  responses <- mydata[1:2000, 1:50]
  theta_true <- mydata[1:2000, 51]
  a_true <- mydata[2001, 1:50]
  rownames(a_true) <- NULL
  b_true <- mydata[2002, 1:50]
  rownames(b_true) <- NULL
  
  N_Sim <- 2000
  perc_unmo_peop <- 0.5 # 0.1; 0.25; or 0.50 
  num_unmo <- N_Sim*perc_unmo_peop #the number of selected unmotivated people
  
  ## Two steps approach No.1
  ## Extracting the fixed item parameter estimates based on data without guessing behavior
  data_fix <- responses[(num_unmo+1):2000, 1:50]
  
  tplfix <- mirt(data_fix,
                 model=1, 
                 itemtype='2PL',
                 SE=TRUE)
  
  itempars_fix <- coef(tplfix,
                       printSE=TRUE,
                       IRTpars=TRUE,
                       as.data.frame=TRUE)
  
  ## Save Item parameters
  lis <- seq(1, by=4, length.out=50)
  a_est_f <- as.numeric(itempars_fix[lis,1])
  b_est_f <- as.numeric(itempars_fix[lis+1,1])
  d_est_f <- -a_est_f*b_est_f
  
  item_param_est_f <- data.frame(a1 = a_est_f, d = d_est_f)
  
  ## Generate IRT model based on fixed item parameters using package mirtCAT
  tplfix_1 <- generate.mirt_object(item_param_est_f, 
                                   itemtype = '2PL')
  
  ## Get the new theta estimates
  theta_2PL_1 <- fscores(tplfix_1, 
                         method ="EAP",
                         response.pattern = responses, 
                         full.scores.SE=FALSE)
  
  theta_values <- data.frame(theta_true, theta_2PL_1[, 1])
  names(theta_values) <- c("theta_true", "theta_est")
  write.csv(theta_values, 
            file = paste0(Path_Save_Result, "/Fixed/Theta/Theta_", substr(csv_table, 51, str_length(csv_table)-4), ".csv"), 
            row.names = F)
  
  
  a_values <- data.frame(t(data.frame(rbind(a_true, a_est_f))))
  names(a_values) <- c("a_true", "a_est")
  write.csv(a_values, 
            file = paste0(Path_Save_Result, "/Fixed/a/a_", substr(csv_table, 51, str_length(csv_table)-4), ".csv"), 
            row.names = F)
  
  
  b_values <- data.frame(t(data.frame(rbind(b_true, b_est_f))))  
  names(b_values) <- c("b_true", "b_est")
  write.csv(b_values, 
            file = paste0(Path_Save_Result, "/Fixed/b/b_", substr(csv_table, 51, str_length(csv_table)-4), ".csv"), 
            row.names = F)
}
