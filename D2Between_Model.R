###############################################
#########   D2 Between 2PL Model     ##########
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
RG_IRT_Model_MBet2PL <- function(
  csv_table,
  Path_Save_Result){
  
  mydata <- read.csv(file = csv_table)
  
  responses <- mydata[1:2000, 1:50]
  theta_true <- mydata[1:2000, 51]
  a_true <- mydata[2001, 1:50]
  rownames(a_true) <- NULL
  b_true <- mydata[2002, 1:50]
  rownames(b_true) <- NULL
  
  ## Rapid-guessing indicator matrix
  data_replace_r <- ifelse(is.na(responses),1,0)
  
  ## Adding rapid-guessing indicator matrix to the response data
  data_replace_m <- cbind(responses, data_replace_r)
  
  ## Rename item 1-100
  item_id_m <- paste('item', 1:100, sep = '_')
  colnames(data_replace_m) <- item_id_m
  
  ## Find the location of the minimum theta of simulee
  min_theta <- which.min(theta_true)
  
  ## Replace the response data in item 51-100 to 1 and item 1-50 to NA 
  data_replace_m[min_theta, 51:100] = rep(1, times=50)
  data_replace_m[min_theta, 1:50] = NA
  
  
  ## Fit two dimensional 2PL IRT between item model
  model_m <- mirt.model ('
          F1 = 1-50
          F2 = 51-100
          COV = F1*F2
          CONSTRAIN = (51-100, a1)')
  
  mtwopl1 <- mirt(data_replace_m, model_m, '2PL')
  
  est_itempars_mirt1 <- coef(mtwopl1, 
                             allpars = TRUE,
                             as.data.frame=TRUE)
  
  MB1 <- MDIFF(mtwopl1)[1:50]
  b_values <- data.frame(t(data.frame(rbind(b_true, MB1))))
  names(b_values) <- c("b_true", "b_est")
  write.csv(b_values, 
            file = paste0(Path_Save_Result, "/D2Between/b/b_", substr(csv_table, 51, str_length(csv_table)-4), ".csv"), 
            row.names = F)
  
  MA1 <- MDISC(mtwopl1)[1:50]
  a_values <- data.frame(t(data.frame(rbind(a_true, MA1))))
  names(a_values) <- c("a_true", "a_est")
  write.csv(a_values, 
            file = paste0(Path_Save_Result, "/D2Between/a/a_", substr(csv_table, 51, str_length(csv_table)-4), ".csv"), 
            row.names = F)
  
  ## Estimation of true ability and propensity of rapid guessing
  est_theta_mirt1 <- fscores(mtwopl1)
  
  theta_values <- data.frame(theta_true, est_theta_mirt1)
  write.csv(theta_values, 
            file = paste0(Path_Save_Result, "/D2Between/Theta/Theta_", substr(csv_table, 51, str_length(csv_table)-4), ".csv"), 
            row.names = F)

}
