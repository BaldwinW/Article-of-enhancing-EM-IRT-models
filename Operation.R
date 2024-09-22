###############################################
########### Running all Functions #############
###############################################

## Date: 20230518
## Authors: Bowen, Hailey


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

## Settings
getwd()
setwd("/Users/bowenwang/Documents/UF/Rapid_Guessing_MIRT-main")

## Generate Data
#source("Data_Generate.R")
#RG_IRT_Data_Generation(
  #perc_unmo_peop_list = c(.1, .25, .5),
  #rela_list = c(-1, -.5, -.2, 0), 
  #percent_guess_items_list = c(.1, .3, .5, .7),
  #guessing_pattern_list = c("a","b","c","d"),
  #reps = 2,
  #N_Sim = 2000,
  #N_items = 50,
  #Path_Save_Data = "/Users/bowenwang/Documents/UF/res_mac")

## Modeling
## Location of the List of Data(.csv) Files
Path_Save_Result <- "/Users/bowenwang/Documents/UF/Rapid_Guessing_MIRT-main/Test_Model"

csv_path <- "/Users/bowenwang/Documents/UF/Test_Data_0518"
#csv_path <- "/Users/bowenwang/Documents/UF/Test_Data_2"
csv_files <- list.files(path = csv_path, pattern = "Data_punmo_0\\.5_rela_-0.5_gitem_0.5.*\\.csv$", full.names = T)

## Model 1 
source("EMIRT_Model.R")
for (file in csv_files) {

  cat("Starting at:")
  print(Sys.time())
  start_time <- Sys.time()
  
  RG_IRT_Model_EMIRT(file,
                     Path_Save_Result)
  
  end_time <- Sys.time()
  cat("Ending at:")
  print(Sys.time())
  cat("Time for Running:", end_time-start_time)
}

## Model 2 
source("Fixed_EMIRT.R")
for (file in csv_files) {
  
  cat("Starting at:")
  print(Sys.time())
  start_time <- Sys.time()
  
  RG_IRT_Model_Fixed2PL(file,
                     Path_Save_Result)
  
  end_time <- Sys.time()
  cat("Ending at:")
  print(Sys.time())
  cat("Time for Running:", end_time-start_time)
}


## Model 3 
source("D2Between_Model.R")
for (file in csv_files) {
  
  cat("Starting at:")
  print(Sys.time())
  start_time <- Sys.time()
  
  RG_IRT_Model_MBet2PL(file, 
                       Path_Save_Result)
  
  end_time <- Sys.time()
  cat("Ending at:")
  print(Sys.time())
  cat("Time for Running:", end_time-start_time)
}

## Model 4
source("D2Within_Model.R")
for (file in csv_files) {
  
  cat("Starting at:")
  print(Sys.time())
  start_time <- Sys.time()
  
  RG_IRT_Model_MWti2PL(file, 
                       Path_Save_Result)
  
  end_time <- Sys.time()
  cat("Ending at:")
  print(Sys.time())
  cat("Time for Running:", end_time-start_time)
}






