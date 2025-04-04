###############################################################################
###############################################################################
## Threshold MAR 50% of dropout
###############################################################################
###############################################################################
## load script with preliminars and auxiliar functions
source("G:/TFM/SLinJointModels/simul_MFPCA_missing/preliminars_functions.R")

list_pred_mfpca1 <- list()
list_pred_mfpca1_test <- list()  

list_pred_mfpca1_y1 <- list()
list_pred_mfpca1_y2 <- list()
list_pred_mfpca1_y3 <- list()
list_pred_mfpca1_y1_test <- list()
list_pred_mfpca1_y2_test <- list()
list_pred_mfpca1_y3_test <- list()



list_preds_tmod_y1_miss <- list()
list_preds_tmod_y1_miss2 <- list()
list_preds_tmod_y1_miss_test <- list()
list_preds_tmod_y2_miss <- list()
list_preds_tmod_y2_miss2 <- list()
list_preds_tmod_y2_miss_test <- list()
list_preds_tmod_y3_miss <- list()
list_preds_tmod_y3_miss2 <- list()
list_preds_tmod_y3_miss_test <- list()
list_rhat_true_model_miss <- list()

list_DF <- list()
list_DF_miss <- list()
list_DF_test <- list()
list_DF_test_miss <- list()

missing_ratio_train <- numeric(num_datasets)
missing_ratio_test <- numeric(num_datasets)

############################################################
############################################################
## I generate the data frames: Threshold MAR 30% of dropout
##.  - participant drops at t_ij with a prob determined by
##.    a logistic model with indicator Y_i(j-1)>nu as a
##.    predictor
############################################################
############################################################
for(count in 1:num_datasets){
  n <- 200 # number of subjects
  n_test <- 200
  K <- 10 # number of measurements per subject

  DF <- data.frame(id = rep(seq_len(n), each = K),
                   time = c(replicate(n, c(0,1:(K-1) + runif(K-1, -0.5,0.5) ))),
                   visit = c(replicate(n, seq(1,K))),
                   sex = rep(gl(2, n/2, labels = c("male", "female")), each = K),
                   treatment = treat)
  DF_test <- data.frame(id = rep(seq_len(n_test), each = K),
                        time = c(replicate(n_test, c(0,1:(K-1) + runif(K-1, -0.5,0.5) ))),
                        visit = c(replicate(n, seq(1,K))),
                        sex = rep(gl(2, n_test/2, labels = c("male", "female")), each = K),
                        treatment = treat)
  
  
  # design matrices for the fixed and random effects for longitudinal submodels
  X1 <- model.matrix(~ sex * time, data = DF)
  Z1 <- model.matrix(~ time, data = DF)

  X3 <- model.matrix(~ sqrt(time), data = DF )
  Z3 <- model.matrix(~ sqrt(time), data = DF)
  
  X5 <- model.matrix(~ treatment + time, data = DF)
  Z5 <- model.matrix(~ time, data = DF)
  
  ####for test data set
  X1_test <- model.matrix(~ sex * time, data = DF_test)
  Z1_test <- model.matrix(~ time, data = DF_test)
  
  X3_test <- model.matrix(~ sqrt(time), data = DF_test )
  Z3_test <- model.matrix(~ sqrt(time), data = DF_test)
  
  X5_test <- model.matrix(~ treatment + time, data = DF_test)
  Z5_test <- model.matrix(~ time, data = DF_test)
  
  #Simulate random effects
  #################################################
  b <-  mvrnorm(n = n, mu = rep(0,6), Sigma = sigmaa%*%t(sigmaa))
  b_test <-  mvrnorm(n = n, mu = rep(0,6), Sigma = sigmaa%*%t(sigmaa))
  
  
  #Simulate first longitudinal outcome
  ###############################################
  betas1 <- c(2.2, 1.6787, 1.24, -0.05) # fixed effects coefficients
  sigma1 <- 0.125 # errors sd
  
  # random effects
  b1 <-  b[, c(1,2)]
  # linear predictor
  eta_y1 <- as.vector(X1 %*% betas1 + rowSums(Z1 * b1[DF$id, ]))
  # we simulate normal longitudinal data
  DF$y1 <- rnorm(n * K, mean = eta_y1, sd = sigma1)
  # we assume that values below 4 are not observed, and set equal to 0
  DF$ind1 <- as.numeric(DF$y1 < -4)
  DF$y1 <- pmax(DF$y1, -4)
  
  #Test data:
  # we simulate random effects
  b1_test <- b_test[, c(1,2)]
  # linear predictor
  eta_y1_test <- as.vector(X1_test %*% betas1 + rowSums(Z1_test * b1_test[DF_test$id, ]))
  # we simulate normal longitudinal data
  DF_test$y1 <- rnorm(n_test * K, mean = eta_y1_test, sd = sigma1)
  # we assume that values below 4 are not observed, and set equal to 0
  DF_test$ind1 <- as.numeric(DF_test$y1 < -4)
  DF_test$y1 <- pmax(DF_test$y1, -4)
  
  #Simulate 2d longitudinal outcome
  ###############################################
  betas3 <- c(-2.5, 0.0333) # fixed effects coefficients
  betas3 <- c(2.5, 1) # fixed effects coefficients
  sigma3 <- 0.25 # errors sd
  
  
  # we simulate random effects
  b3 <- b[, c(3,4)]
  # linear predictor
  eta_y3 <- as.vector(X3 %*% betas3 +
                        rowSums(Z3 * b3[DF$id, ]))
  # we simulate normal longitudinal data
  DF$y3 <- rnorm(n * K, mean = eta_y3, sd = sigma3)
  # we assume that values below 0 are not observed, and set equal to 0
  
  
  ########Test data
  # we simulate random effects
  b3_test <- b_test[, c(3,4)]
  # linear predictor
  eta_y3_test <- as.vector(X3_test %*% betas3 + rowSums(Z3_test * b3_test[DF_test$id, ]))
  # we simulate normal longitudinal data
  DF_test$y3 <- rnorm(n_test * K, mean = eta_y3_test, sd = sigma3)
 
  #Simulate 3rd longitudinal outcome
  ###############################################
  betas5 <- c(1, 0.155, 0.12345) # fixed effects coefficients
  betas5 <- c(70, -5.44, -3.761289) # fixed effects coefficients
  
  # we simulate random effects
  b5 <- b[, c(5,6)]
  # linear predictor
  eta_y5 <- as.vector(X5 %*% betas5 + rowSums(Z5 * b5[DF$id, ]))
  DF$y5 <- rnorm(n_test * K, mean = eta_y5, sd = sigma3)
  
  ####Test data
  # we simulate random effects
  b5_test <- b_test[, c(5,6)]
  # linear predictor
  eta_y5_test <- as.vector(X5_test %*% betas5 + rowSums(Z5_test * b5_test[DF_test$id, ]))
  DF_test$y5 <- rnorm(n_test * K, mean = eta_y5_test, sd = sigma3)

  
  ### threshold MAR dropout with 50% intensity
  nu1 <- 15.5
  nu2 <- 7.5
  nu3 <- 42
  
  DF <- DF %>%
    group_by(id) %>%                       # Group by each individual
    mutate(
      # Check dropout conditions on the previous time point
      dropout_point = which.max(lag(y1) > nu1 | lag(y3) > nu2 | lag(y5) < nu3),
      
      # Set NA to y1, y2, y3 after the dropout point (ignoring the first visit)
      y1_ind = ifelse(row_number() >= dropout_point & dropout_point > 1, 1, 0),
      y3_ind = ifelse(row_number() >= dropout_point & dropout_point > 1, 1, 0),
      y5_ind = ifelse(row_number() >= dropout_point & dropout_point > 1, 1, 0)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point)
  
  DF_test <- DF_test %>%
    group_by(id) %>%                       # Group by each individual
    mutate(
      # Check dropout conditions on the previous time point
      dropout_point = which.max(lag(y1) > nu1 | lag(y3) > nu2 | lag(y5) < nu3),
      
      # Set NA to y1, y2, y3 after the dropout point (ignoring the first visit)
      y1_ind = ifelse(row_number() >= dropout_point & dropout_point > 1, 1, 0),
      y3_ind = ifelse(row_number() >= dropout_point & dropout_point > 1, 1, 0),
      y5_ind = ifelse(row_number() >= dropout_point & dropout_point > 1, 1, 0)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point)
  
  param1 <- -2.8
  param2 <- 0.675
  probs_vec1 <- numeric(n*K)
  probs_vec2 <- numeric(n*K)
  probs_vec3 <- numeric(n*K)
  
  kk <- 1
  for(i in 1:n){
    for(j in 1:K){
      if((kk-1)%%K==0){
        probs_vec1[kk] <- 0
      } else{
        thet1 <- param1 + param2*DF[DF$id==i,]$y1_ind[j]
        probs_vec1[kk] <- 1/(1+exp(-thet1))
      }
      kk <- kk + 1
    }
  }
  
  kk <- 1
  for(i in 1:n){
    for(j in 1:K){
      if((kk-1)%%K==0){
        probs_vec2[kk] <- 0
      } else{
        thet2 <- param1 + param2*DF[DF$id==i,]$y3_ind[j]
        probs_vec2[kk] <- 1/(1+exp(-thet2))
      }
      kk <- kk + 1
    }
  }
  
  kk <- 1
  for(i in 1:n){
    for(j in 1:K){
      if((kk-1)%%K==0){
        probs_vec3[kk] <- 0
      } else{
        thet3 <- param1 + param2*DF[DF$id==i,]$y5_ind[j]
        probs_vec3[kk] <- 1/(1+exp(-thet3))
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
  DF$probs_drop_y3 <- probs_vec2
  DF$dropout_y3 <- rbinom(n*K, size = 1, prob = probs_vec2)
  
  DF_test$probs_dropout_y3 <- probs_vec2
  DF_test$dropout_y3 <- rbinom(n*K, size = 1, prob = probs_vec2)
  
  #L=3
  DF$probs_drop_y5 <- probs_vec3
  DF$dropout_y5 <- rbinom(n*K, size = 1, prob = probs_vec3)
  
  DF_test$probs_dropout_y5 <- probs_vec3
  DF_test$dropout_y5 <- rbinom(n*K, size = 1, prob = probs_vec3)
  
  ###Putting NAs after dropout
  #DF_complete <- DF
  DF_miss <- DF %>%
    group_by(id) %>%                     # Group by each individual
    mutate(
      # Check dropout conditions on the previous time point
      dropout_point = which.max(dropout_y1 == 1 | dropout_y3 == 1 | dropout_y5 == 1),
      
      # Set NA to y1, y2, y3 after the dropout point (ignoring the first visit)
      y1 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y1),
      y3 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y3),
      y5 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y5)
    ) %>%
    ungroup() %>%
    dplyr::select(id, time, visit, sex, treatment, y1, ind1, y3, y5)
  
  DF_test_miss <- DF_test %>%
    group_by(id) %>%                     # Group by each individual
    mutate(
      # Check dropout conditions on the previous time point
      dropout_point = which.max(dropout_y1 == 1 | dropout_y3 == 1 | dropout_y5 == 1),
      
      # Set NA to y1, y2, y3 after the dropout point (ignoring the first visit)
      y1 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y1),
      y3 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y3),
      y5 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y5)
    ) %>%
    ungroup() %>%
    dplyr::select(id, time, visit, sex, treatment, y1, ind1, y3, y5)
  
  
  
  
  
  missing_ratio_train[count] <- counting_missing_data(DF_miss, n, K)
  missing_ratio_test[count] <- counting_missing_data(DF_test_miss, n, K)
  
  try(list_DF <- append(list_DF, list(DF)))
  try(list_DF_miss <- append(list_DF_miss, list(DF_miss)))
  try(list_DF_test <- append(list_DF_test, list(DF_test)))
  try(list_DF_test_miss <- append(list_DF_test_miss, list(DF_test_miss)))
  
  try(if(count==100){
    setwd("D:/La meva unitat/TFM/ResultsMMvsMFPCA")
    str100 <- "dataframes_MAR_threshold_50.RData"
    save(list_DF, list_DF_miss, list_DF_test, list_DF_test_miss,
         missing_ratio_train, missing_ratio_test, file=str100)
  })
  
  print(count)
}

####################################################
####################################################
## Fitting MM and MFPCA models and doing predictions
####################################################
####################################################

for(count in 1:num_datasets){
  ######################################
  #Fitting multivariate MM models
  # - We use snat_mvmer() from rstanarm
  # - We fit the true model plus two
  #.  (or 3) more models
  ######################################
  try(true_model_miss <- stan_mvmer(
    formula = list(
      y1 ~ sex * time + (time | id),
      y3 ~ sqrt(time) + (sqrt(time) | id),
      y5 ~ treatment + time + (time | id)),
    data = list_DF_miss[[count]],
    family = list(gaussian, gaussian, gaussian),
    chains = 3, cores = 8, seed = 12345, iter = 1750))
  
  try(rhat_true_model_miss <- rhat(true_model_miss))

  #PREDICTIONS (for three outcomes separately)
  try(preds_tmod_y1_miss <- posterior_predict(true_model_miss, m=1))
  try(preds_tmod_y2_miss <- posterior_predict(true_model_miss, m=2))
  try(preds_tmod_y3_miss <- posterior_predict(true_model_miss, m=3))
  
  #since those matrices take up a lot of memory, we summary the draws computing
  #the mean by columns and we and up with a vector of 2000 elements
  try(preds_tmod_y1_miss <- colMeans(preds_tmod_y1_miss))
  try(preds_tmod_y2_miss <- colMeans(preds_tmod_y2_miss))
  try(preds_tmod_y3_miss <- colMeans(preds_tmod_y3_miss))
  
  
  #PREDICTIONS (for three outcomes separately)
  # also predictiing for the missings
  DF_miss2 <- list_DF_miss[[count]] %>% dplyr::select(id, time, visit, sex, treatment)
  
  try(preds_tmod_y1_miss2 <- posterior_predict(true_model_miss, m=1,
                                               newdata = DF_miss2))
  try(preds_tmod_y2_miss2 <- posterior_predict(true_model_miss, m=2,
                                               newdata = DF_miss2))
  try(preds_tmod_y3_miss2 <- posterior_predict(true_model_miss, m=3,
                                               newdata = DF_miss2))
  
  #since those matrices take up a lot of memory, we summary the draws computing
  #the mean by columns and we and up with a vector of 2000 elements
  try(preds_tmod_y1_miss2 <- colMeans(preds_tmod_y1_miss2))
  try(preds_tmod_y2_miss2 <- colMeans(preds_tmod_y2_miss2))
  try(preds_tmod_y3_miss2 <- colMeans(preds_tmod_y3_miss2))
  
  #transforming to a matrix 200*10
  try(preds_tmod_y1_miss2 <- matrix(preds_tmod_y1_miss2, nrow = 200, ncol = 10, 
                                    byrow = TRUE))
  try(preds_tmod_y2_miss2 <- matrix(preds_tmod_y2_miss2, nrow = 200, ncol = 10, 
                                    byrow = TRUE))
  try(preds_tmod_y3_miss2 <- matrix(preds_tmod_y3_miss2, nrow = 200, ncol = 10, 
                                    byrow = TRUE))
  
  ## predictive posterior with TEST data set
  DF_test_miss2 <- list_DF_test_miss[[count]] %>% dplyr::select(id, time, visit, sex, treatment)
  
  try(preds_tmod_y1_miss_test <- posterior_predict(true_model_miss, m=1,
                                                   newdata = DF_test_miss2))
  try(preds_tmod_y2_miss_test <- posterior_predict(true_model_miss, m=2,
                                                   newdata = DF_test_miss2))
  try(preds_tmod_y3_miss_test <- posterior_predict(true_model_miss, m=3,
                                                   newdata = DF_test_miss2))
  
  #since those matrices take up a lot of memory, we summary the draws computing
  #the mean by columns and we and up with a vector of 2000 elements
  try(preds_tmod_y1_miss_test <- colMeans(preds_tmod_y1_miss_test))
  try(preds_tmod_y2_miss_test <- colMeans(preds_tmod_y2_miss_test))
  try(preds_tmod_y3_miss_test <- colMeans(preds_tmod_y3_miss_test))
  
  #transforming to a matrix 200*10
  try(preds_tmod_y1_miss_test <- matrix(preds_tmod_y1_miss_test, nrow = 200, ncol = 10, 
                                        byrow = TRUE))
  try(preds_tmod_y2_miss_test <- matrix(preds_tmod_y2_miss_test, nrow = 200, ncol = 10, 
                                        byrow = TRUE))
  try(preds_tmod_y3_miss_test <- matrix(preds_tmod_y3_miss_test, nrow = 200, ncol = 10, 
                                        byrow = TRUE))
  
  
  ######################################
  #Fitting MFPCA models
  # - We use mfpca R package
  # - We fit three different approxs
  #.  by changing PVE used
  # - We also fit different models
  #.  using or nor weights
  ######################################
  
  #We compute grid longitudinal data:
  grid_long <- grid_longitudinal_data(list_DF_miss[[count]], n, K)
  #univariate functional data
  obs_time <- grid_long[[4]]
  f1 <- funData(obs_time, as.matrix(grid_long[[1]][,-1]))
  f2 <- funData(obs_time, as.matrix(grid_long[[2]][,-1]))
  f3 <- funData(obs_time, as.matrix(grid_long[[3]][,-1]))
  #multivariate functional data class
  m1 <- multiFunData(list(f1,f2,f3))
  
  
  # Grid longitudinal data in TEST:
  grid_long_test <- grid_longitudinal_data(list_DF_test_miss[[count]], n, K)
  #univariate functional data TEST case
  obs_time_test <- grid_long_test[[4]]
  f1_test <- funData(obs_time_test, as.matrix(grid_long_test[[1]][,-1]))
  f2_test <- funData(obs_time_test, as.matrix(grid_long_test[[2]][,-1]))
  f3_test <- funData(obs_time_test, as.matrix(grid_long_test[[3]][,-1]))
  #multivariate functional data class TEST case
  m1_test <- multiFunData(list(f1_test,f2_test,f3_test))
  
  #now we use MFPCA function from mfpca R package
  
  #We fit with NO weights and pve=0.65
  try(mfpca1 <- MFPCA(m1, M = 3, 
                      uniExpansions = list(list(type="uFPCA", pve = 0.65),
                                           list(type ="uFPCA", pve=0.65),
                                           list(type="uFPCA", pve=0.65)),
                      fit = TRUE))
  
  
  ###OBS: we have to adjust M and pve depending on % of missing data we have 
  
  # predict in training data:
  try(pred_mfpca1 <- predict(mfpca1))
  
  try(pred_mfpca1_y1 <- NULL)
  try(pred_mfpca1_y1 <- pred_mfpca1[[1]]@X)
  try(pred_mfpca1_y2 <- NULL)
  try(pred_mfpca1_y2 <- pred_mfpca1[[2]]@X)
  try(pred_mfpca1_y3 <- NULL)
  try(pred_mfpca1_y3 <- pred_mfpca1[[3]]@X)
  
  #plot(pred_mfpca1)
  
  # predict in testing data:
  # value "vector" returned by MFPCA seems to be the cms in paper (and ML code)
  # we compute scores for testing data and we use meanf and psi from training
  try(mfpca1_test <- MFPCA(m1_test, M = 3, 
                           uniExpansions = list(list(type="uFPCA", pve = 0.8),
                                                list(type ="uFPCA", pve=0.8),
                                                list(type="uFPCA", pve=0.8)),
                           fit = TRUE))
  
  try(pred_mfpca1_test <- predict(mfpca1, scores = mfpca1_test$scores))
  
  try(pred_mfpca1_y1_test <- NULL)
  try(pred_mfpca1_y1_test <- pred_mfpca1_test[[1]]@X)
  try(pred_mfpca1_y2_test <- NULL)
  try(pred_mfpca1_y2_test <- pred_mfpca1_test[[2]]@X)
  try(pred_mfpca1_y3_test <- NULL)
  try(pred_mfpca1_y3_test <- pred_mfpca1_test[[3]]@X)
  
  #########################################################
  #########################################################
  ## saving data
  #########################################################
  #########################################################
  try(list_pred_mfpca1 <- append(list_pred_mfpca1, list(pred_mfpca1)))
  try(list_pred_mfpca1_test <- append(list_pred_mfpca1_test, list(pred_mfpca1_test)))
  
  try(list_pred_mfpca1_y1 <- append(list_pred_mfpca1_y1, list(pred_mfpca1_y1)))
  try(list_pred_mfpca1_y2 <- append(list_pred_mfpca1_y2, list(pred_mfpca1_y2)))
  try(list_pred_mfpca1_y3 <- append(list_pred_mfpca1_y3, list(pred_mfpca1_y3)))
  
  try(list_pred_mfpca1_y1_test <- append(list_pred_mfpca1_y1_test, 
                                         list(pred_mfpca1_y1_test)))
  try(list_pred_mfpca1_y2_test <- append(list_pred_mfpca1_y2_test, 
                                         list(pred_mfpca1_y2_test)))
  try(list_pred_mfpca1_y3_test <- append(list_pred_mfpca1_y3_test, 
                                         list(pred_mfpca1_y3_test)))
  
  try(list_preds_tmod_y1_miss <- append(list_preds_tmod_y1_miss, list(preds_tmod_y1_miss)))
  try(list_preds_tmod_y1_miss2 <- append(list_preds_tmod_y1_miss2, list(preds_tmod_y1_miss2)))
  try(list_preds_tmod_y1_miss_test <- append(list_preds_tmod_y1_miss_test, list(preds_tmod_y1_miss_test)))
  try(list_preds_tmod_y2_miss <- append(list_preds_tmod_y2_miss, list(preds_tmod_y2_miss)))
  try(list_preds_tmod_y2_miss2 <- append(list_preds_tmod_y2_miss2, list(preds_tmod_y2_miss2)))
  try(list_preds_tmod_y2_miss_test <- append(list_preds_tmod_y2_miss_test, list(preds_tmod_y2_miss_test)))
  try(list_preds_tmod_y3_miss <- append(list_preds_tmod_y3_miss, list(preds_tmod_y3_miss)))
  try(list_preds_tmod_y3_miss2 <- append(list_preds_tmod_y3_miss2, list(preds_tmod_y3_miss2)))
  try(list_preds_tmod_y3_miss_test <- append(list_preds_tmod_y3_miss_test, list(preds_tmod_y3_miss_test)))
  try(list_rhat_true_model_miss <- append(list_rhat_true_model_miss, list(rhat_true_model_miss)))
  
  
  
  
  setwd("D:/La meva unitat/TFM/ResultsMMvsMFPCA")
  try(if(count==1){
    strr <- "results_1_mar50_thr_14feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr)
  })
  try(if(count==2){
    strr2 <- "results_2_mar50_thr_14feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr2)
  })
  try(if(count==5){
    str1 <- "results_5_mar50_thr_14feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str1)
  })
  try(if(count==10){
    str10 <- "results_10_mar50_thr_14feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str10)
  })
  try(if(count==20){
    str20 <- "results_20_mar50_thr_14feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str20)
  })
  try(if(count==30){
    str30 <- "results_30_mar50_thr_14feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str30)
  })
  try(if(count==40){
    str40 <- "results_40_mar50_thr_14feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss,file=str40)
  })
  try(if(count==50){
    str50 <- "results_50_mar50_thr_14feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str50)
  })
  try(if(count==75){
    str75 <- "results_75_mar50_thr_14feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_mo111del_miss, file=str75)
  })
  try(if(count==100){
    str100 <- "results_100_mar50_thr_14feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str100)
  })
  
  print(count)
  
  
  
  ######################################
  #END FOR
  ######################################
}


