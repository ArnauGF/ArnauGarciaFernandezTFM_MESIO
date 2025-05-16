library(JMbayes2)
library(MASS)
library(lattice)
library(dplyr)
library(tidyr)

pow2 <- function(x){
  x^2
}

# we load the lists of data sets
load("G:/La meva unitat/TFM/MFPCA_simuls/R/datasets_MFPCAvsSL_29apr2025.RData")


repl <- 100
list_DF_infCens <- list()
list_DF_test_infCens <- list()


for(count in 1:repl){
  n <- 175 # number of subjects
  n_test <- 175
  K <- 10 # number of measurements per subject
  
  # we set DF ans DF_test to the corresponding element of the list
  DF <- list_DF[[count]]
  DF_test <- list_DF_test[[count]]
  
  ##saving true times
  trueTimes <- DF[!duplicated(DF$id),]$Time
  
  ##saving true times
  trueTimes_test <- DF_test[!duplicated(DF_test$id),]$Time
  
  ####Let us simulate informative censoring:
  
  param1 <- -4.2
  param2 <- 0.01
  probs_vec1 <- numeric(n*K)
  probs_vec2 <- numeric(n*K)
  probs_vec3 <- numeric(n*K)
  probs_vec4 <- numeric(n*K)
  
  kk <- 1
  for(i in 1:n){
    for(j in 1:K){
      if((kk-1)%%K==0){
        probs_vec1[kk] <- 0
        probs_vec2[kk] <- 0
        probs_vec3[kk] <- 0
        probs_vec4[kk] <- 0
      } else{
        thet1 <- param1 + param2*DF[DF$id==i,]$y1[j]
        probs_vec1[kk] <- 1/(1+exp(-thet1))
        thet2 <- param1 + param2*DF[DF$id==i,]$y2[j]
        probs_vec2[kk] <- 1/(1+exp(-thet2))
        thet3 <- param1 + param2*DF[DF$id==i,]$y3[j]
        probs_vec3[kk] <- 1/(1+exp(-thet3))
        thet4 <- param1 + param2*DF[DF$id==i,]$y4[j]
        probs_vec4[kk] <- 1/(1+exp(-thet4))
      }
      kk <- kk + 1
    }
  }
  
  
  #L=1
  DF$probs_drop_y1 <- probs_vec1
  DF$dropout_y1 <- rbinom(n*K, size = 1, prob = probs_vec1)
  
  DF_test$probs_dropout_y1 <- probs_vec1
  DF_test$dropout_y1 <- rbinom(n*K, size = 1, prob = probs_vec1)
  
  #L=2
  DF$probs_drop_y2 <- probs_vec2
  DF$dropout_y2 <- rbinom(n*K, size = 1, prob = probs_vec2)
  
  DF_test$probs_dropout_y2 <- probs_vec2
  DF_test$dropout_y2 <- rbinom(n*K, size = 1, prob = probs_vec2)
  
  #L=3
  DF$probs_drop_y3 <- probs_vec3
  DF$dropout_y3 <- rbinom(n*K, size = 1, prob = probs_vec3)
  
  DF_test$probs_dropout_y3 <- probs_vec3
  DF_test$dropout_y3 <- rbinom(n*K, size = 1, prob = probs_vec3)
  
  #L=4
  DF$probs_drop_y4 <- probs_vec3
  DF$dropout_y4 <- rbinom(n*K, size = 1, prob = probs_vec3)
  
  DF_test$probs_dropout_y4 <- probs_vec3
  DF_test$dropout_y4 <- rbinom(n*K, size = 1, prob = probs_vec3)
  
  DF <- DF %>%
    group_by(id) %>%
    # Find first dropout for y1, y2, y3, y4 and store corresponding time in Ctimes
    mutate(Ctimes_y1 = ifelse(any(dropout_y1 == 1), time[which(dropout_y1 == 1)[1]], Time),
           Ctimes_y2 = ifelse(any(dropout_y2 == 1), time[which(dropout_y2 == 1)[1]], Time),
           Ctimes_y3 = ifelse(any(dropout_y3 == 1), time[which(dropout_y3 == 1)[1]], Time),
           Ctimes_y4 = ifelse(any(dropout_y4 == 1), time[which(dropout_y4 == 1)[1]], Time)) %>%
    mutate(C = pmin(Ctimes_y1, Ctimes_y2, Ctimes_y3, Ctimes_y4)) %>%
    ungroup()
  
  Ctimes <- DF[!duplicated(DF$id), ]$C
  
  DF_test <- DF_test %>%
    group_by(id) %>%
    # Find first dropout for y1, y2, y3, y4 and store corresponding time in Ctimes
    mutate(Ctimes_y1 = ifelse(any(dropout_y1 == 1), time[which(dropout_y1 == 1)[1]], Time),
           Ctimes_y2 = ifelse(any(dropout_y2 == 1), time[which(dropout_y2 == 1)[1]], Time),
           Ctimes_y3 = ifelse(any(dropout_y3 == 1), time[which(dropout_y3 == 1)[1]], Time),
           Ctimes_y4 = ifelse(any(dropout_y4 == 1), time[which(dropout_y4 == 1)[1]], Time)) %>%
    mutate(C = pmin(Ctimes_y1, Ctimes_y2, Ctimes_y3, Ctimes_y4)) %>%
    ungroup()
  
  Ctimes_test <- DF_test[!duplicated(DF_test$id), ]$C
  
  Time <- pmin(trueTimes, Ctimes)
  event <- as.numeric(trueTimes <= Ctimes) # event indicator
  # we keep the longitudinal measurements before the event times
  DF$Time <- Time[DF$id]
  DF$event <- event[DF$id]
  DF <- DF[DF$time <= DF$Time, ]
  
  ## Testing data:
  Time <- pmin(trueTimes_test, Ctimes_test)
  event <- as.numeric(trueTimes_test <= Ctimes_test) # event indicator
  # we keep the longitudinal measurements before the event times
  DF_test$Time <- Time[DF_test$id]
  DF_test$event <- event[DF_test$id]
  DF_test <- DF_test[DF_test$time <= DF_test$Time, ]
  
  
  ## We save the datasets
  list_DF_infCens <- append(list_DF_infCens, list(DF))
  list_DF_test_infCens <- append(list_DF_test_infCens, list(DF_test))
  
  if(count == repl){
    save(list_DF_infCens, list_DF_test_infCens,
         file="G:/La meva unitat/TFM/MFPCA_simuls/R/datasets_informativeCens_MFPCA_29apr2025.RData")
  }
  
  print(count)
}