###############################################
#### Generate Data for Rapid_Guessing #########
###############################################

## Date: 20230518

## Load Packages
Packages <- c("truncnorm", "MASS")

for (BAO in Packages) {
  if (BAO %in% rownames(installed.packages()) == FALSE){
    install.packages(BAO)
  }
  if (BAO %in% (.packages()) == FALSE){
    library(BAO, character.only = TRUE)
  }
}


## Self-defined 2PL Function
twopl <- function(item_param, theta, i){
  (exp(item_param[i,3]*
         (theta-item_param[i,2])))/(1+(exp(item_param[i,3]*(theta-item_param[i,2]))))
}


## Data generation Function (Main Function)
RG_IRT_Data_Generation <- function(
  perc_unmo_peop_list,
  rela_list, 
  percent_guess_items_list,
  guessing_pattern_list,
  reps,
  N_Sim,
  N_items,
  Path_Save_Data){ 
  
  cat("Starting at:")
  print(Sys.time())
  start_time <- Sys.time()
  
  ## percentage of unmotivated people (people who made rapid guessing)
  for (perc_unmo_peop in perc_unmo_peop_list){
    cat("Percentage of unmotivated people:", perc_unmo_peop, "\n")
    
    #rela iis the mean of theta for unmotivated people, which are smaller than 0, 
    ##suggesting people with low true ability intend to make rapid guessing
    ###in terms of MIRT, rela indicates the correlation between true ability and propensity of rapid guessing
    for (rela in rela_list){  
      cat("Relation:", rela, "\n")
      
      #percentage of guessing items
      for (percent_guess_items in percent_guess_items_list){
        
        #Three guessing patterns
        for (guessing_pattern in guessing_pattern_list){
          
          #number of replications for each condition
          for (r in 1:reps) {
            
            #Generate item an people parameters for EM-IRT
            beta <- rtruncnorm(n=50, a=-1.25, b=1.71, mean=0.28, sd=0.61)  #creating item parameters based on 
            
            #alpha <- rtruncnorm(n=50, a=0.35, b=1.79, mean=0.97, sd=0.29)
            alpha <- runif(n=50, min=0.35, max=1.79)
            
            item_param <- data.frame("item"=c(1:50),
                                     "b"=beta,
                                     "a"=alpha)
            
            num_guess_items <- N_items*percent_guess_items #the number of the guessing items
            
            num_unmo <- N_Sim*perc_unmo_peop #the number of selected unmotivated people
            
            new <- c()
            new2 <- c()
            
            theta_unmo <- rnorm(num_unmo,rela,1)   #theta of unmotivated people
            
            theta_mo <- rnorm(N_Sim-(num_unmo),0,1) #theta of motivated people
            
            theta <- c(theta_unmo,theta_mo)#the selected unmotivated people + the rest motivated people
            
            new <- c(new,theta)
            new2 <- c(new2,theta_unmo)#theta of unmotivated people
            
            # Generating item responses 
            data <- c()
            
            for(j in 1:N_items) {
              data <- cbind(data,
                            ifelse(twopl(item_param,theta,j) >= runif(length(theta)),1,0))
            }
            
            #the probabilities of correct responding to an item for each unmotivated people 
            
            p_selected <- c()
            
            for(j in 1:N_items){
              p_correct<-(twopl(item_param,new2,j))
              
              p_selected <- c(p_selected,p_correct)
              
            }
            
            #the probability of correct responding to 50 items for each unmotivated people
            p_unmo <- matrix(data=p_selected, nrow=num_unmo, ncol=50)
            
            #create a sequence (1, 51, 101, 151,...)
            #person1 responding 50 items; person2 responding 50 items...
            g <- seq(1, by=50,length.out=num_unmo) 
            
            #Guessing pattern a
            if(guessing_pattern == "a"){
              
              p_new <- c()
              for (m in g){
                n <- m+49
                result1 <- order(t(p_unmo)[m:n],decreasing=F)[1:num_guess_items]##choose n most difficult items 
                rp <- replace(t(p_unmo)[m:n], result1, 0.25)                    ##replace the probability to .25
                p_new <- c(p_new, rp)
              }
              
              
            }  
            
            #Guessing pattern b
            if(guessing_pattern == "b"){
              
              p_new <- c()
              for (m in g){
                n <- m+49
                l <- 50-num_guess_items+1
                rp <- replace(t(p_unmo)[m:n], (l:50), 0.25)#replace the probability of responding last n items to .25
                p_new <- c(p_new, rp)
              }
              
              
            }
            
            #Guessing pattern c
            
            if(guessing_pattern == "c"){        #split 50 items into 5 bins (10,10,10,10,10)
              p_new <-c()
              for (m in g){ #set up sequences 
                g1 <- m+9 #g=(1,51,101,151...)
                d <- m+10
                d1 <- m+19
                e <- m+20
                e1 <- m+29
                w <- m+30
                w1 <- m+39
                s <- m+40
                s1 <- m+49
                if(percent_guess_items == .1){    #5 items out of 50 items are from guessing (0,0,1,2,2)
                  random1 <- sample(1:10,1)
                  random2 <- sample(1:10,2)
                  rp1 <- t(p_unmo)[m:g1]  #keep 1st bin same
                  rp2 <- t(p_unmo)[d:d1] #keep 2nd bin same
                  rp3 <- replace(t(p_unmo)[e:e1], random1, 0.25) #randomly replace 1 item in 3rd bin to .25
                  rp4 <- replace(t(p_unmo)[w:w1], random2, 0.25) #randomly replace 2 items in 4th bin to .25
                  rp5 <- replace(t(p_unmo)[s:s1], random2, 0.25) #randomly replace 2 items in 5th bin to .25
                  p_new <- c(p_new, rp1, rp2, rp3, rp4, rp5)
                }
                
                
                if(percent_guess_items == .3){   #15 items out of 50 items are from guessing (1,2,3,4,5)
                  random1 <- sample(1:10,1)
                  random2 <- sample(1:10,2)
                  random3 <- sample(1:10,3)
                  random4 <- sample(1:10,4)
                  random5 <- sample(1:10,5)
                  rp1 <- replace(t(p_unmo)[m:g1], random1, 0.25)#randomly replace 1 item in 1st bin to .25
                  rp2 <- replace(t(p_unmo)[d:d1], random2, 0.25)#randomly replace 2 item in 2nd bin to .25
                  rp3 <- replace(t(p_unmo)[e:e1], random3, 0.25)#randomly replace 3 item in 3rd bin to .25
                  rp4 <- replace(t(p_unmo)[w:w1], random4, 0.25)#randomly replace 4 item in 4th bin to .25
                  rp5 <- replace(t(p_unmo)[s:s1], random5, 0.25)#randomly replace 5 item in 5th bin to .25
                  p_new <- c(p_new, rp1, rp2, rp3, rp4, rp5)
                }
                
                if(percent_guess_items == .5){   #25 items out of 50 items are from guessing (3,4,5,6,7)
                  random1 <- sample(1:10,3)
                  random2 <- sample(1:10,4)
                  random3 <- sample(1:10,5)
                  random4 <- sample(1:10,6)
                  random5 <- sample(1:10,7)
                  rp1 <- replace(t(p_unmo)[m:g1], random1, 0.25)#randomly replace 3 item in 1st bin to .25
                  rp2 <- replace(t(p_unmo)[d:d1], random2, 0.25)#randomly replace 4 item in 2nd bin to .25
                  rp3 <- replace(t(p_unmo)[e:e1], random3, 0.25)#randomly replace 5 item in 3rd bin to .25
                  rp4 <- replace(t(p_unmo)[w:w1], random4, 0.25)#randomly replace 6 item in 4th bin to .25
                  rp5 <- replace(t(p_unmo)[s:s1], random5, 0.25)#randomly replace 7 item in 5th bin to .25
                  p_new <- c(p_new, rp1, rp2, rp3, rp4, rp5)
                }
                
                if(percent_guess_items == .7){     #35 items out of 50 items are from guessing (5,6,7,8,9)
                  random1 <- sample(1:10,5)
                  random2 <- sample(1:10,6)
                  random3 <- sample(1:10,7)
                  random4 <- sample(1:10,8)
                  random5 <- sample(1:10,9)
                  rp1 <- replace(t(p_unmo)[m:g1], random1, 0.25)#randomly replace 5 item in 1st bin to .25
                  rp2 <- replace(t(p_unmo)[d:d1], random2, 0.25)#randomly replace 6 item in 2nd bin to .25
                  rp3 <- replace(t(p_unmo)[e:e1], random3, 0.25)#randomly replace 7 item in 3rd bin to .25
                  rp4 <- replace(t(p_unmo)[w:w1], random4, 0.25)#randomly replace 8 item in 4th bin to .25
                  rp5 <- replace(t(p_unmo)[s:s1], random5, 0.25)#randomly replace 9 item in 5th bin to .25
                  p_new <- c(p_new, rp1, rp2, rp3, rp4, rp5)
                }
                
                
              } ## closing sequences loop
            } ## closing guessing pattern c 
            
            ## Guessing pattern d
            if(guessing_pattern == "d"){
              
              p_new <- c()
              for (m in g){
                n <- m+49
                random <- sample(1:50,num_guess_items)#randomly select n items from guessing
                rp <- replace(t(p_unmo)[m:n], random, 0.25)#replace the probability of responding those items to 0.25
                p_new <- c(p_new, rp) 
              }
            } ##close "d" 
            
            
            ## Flip data file so each row is a simulee and each column is an item
            p_new <- t(matrix(data = p_new, nrow = 50, ncol = num_unmo))
            
            ## Getting the location of guessing responses 
            guessing_place <- which(p_new == 0.25, arr.ind = T)
            
            ## Replacing guessing data to NA and getting the final data
            data_replace <- replace(data, guessing_place, NA)
            
            # print(data_replace)
            # print(is.data.frame(data_replace))
            # print(is(list(data_replace)))
            
            ## Giving items names for analysis
            item_id <- paste('item', 1:50, sep = '_')
            
            colnames(data_replace) <- item_id
            
            data_replace <- data.frame(data_replace)
            data_replace$theta_true <- theta
            data_replace[2001, 1:50] <- alpha
            data_replace[2002, 1:50] <- beta
            
            filename <- paste0(Path_Save_Data, "/Data_punmo_", perc_unmo_peop, 
                   "_rela_", rela, 
                   "_gitem_", percent_guess_items, 
                   "_gpattern_", guessing_pattern, "_rep_", r, ".csv")
            
            
            write.csv(data_replace, file = filename, row.names = F)
            
          }  ## closing replications 
        } #closing guessing pattern loop
      }  #closing percent_guess_item loop
    }  #closing rela loop
  } #closing percent_unmo_people loop
  
  end_time <- Sys.time()
  cat("Ending at:")
  print(Sys.time())
  cat("Time for Running:", end_time-start_time)
  
} ## closing the main function

