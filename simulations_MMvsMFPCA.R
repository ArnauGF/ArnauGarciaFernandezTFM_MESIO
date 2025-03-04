#Simulation study multivariate MM vs MFPCA

###########################################################################
###########################################################################
# General settings
# - L=3 (stan_mvmer() accepts max 3 outcomes)
# - Observations times in a fixed grid every 0.5. Visits every 5 or 3 or 1 year
#   with a noise to have realistic data (each possibility is a different scenario)
#.  by the moment the noise is doing summing U(-0.5,0.5)
# - We focus on dropout as a missingness mechanism, we use 30% and 60% of dropout
# - stan_mvmer() function from rstanarm package to fit multivariate MM under
#   the Bayesian framework
# - MFPCA R package to conduct MFPCA anlaysis
###########################################################################
###########################################################################

###########################################################################
#SCENARIO I:
# - MCAR: individual drops out at time t_{ij} with a probability determined
#   by a logistic model with t_{ij} as a predictor
# - Visits every year with +-0.5 noise, 
# - 30% of dropout
###########################################################################
library(dplyr)
library(tidyr)
library(MFPCA)
library(rstanarm)
library(MASS)
library(bayesplot)
num_datasets <- 100
n <- 200 # number of subjects
n_test <- 200
K <- 10
#####D matrix for the random effects:
set.seed(12345)
#we prepare 6x6 matrix for cov-var matrix of random-effects
sigmaa <- matrix(c(runif(1,0,1.5), runif(6, -0.005, 0.005), 
                   runif(1,0,1.5), runif(6, -0.005, 0.005),
                   runif(1,0,1.5), runif(6, -0.005, 0.005),
                   runif(1,0,1.5), runif(6, -0.005, 0.005),
                   runif(1,0,1.5), runif(6, -0.005, 0.005),
                   runif(1,0,1.5))
                 ,nrow = 6)
treat <- rep(as.factor(sample(c("A", "B"), size = n, replace = TRUE))
             ,each = K)

####### Function to prepare the longitudinal data in the fixed grid format
grid_longitudinal_data <- function(DF, n, K){
  DF_mfpca <- DF
  DF_mfpca <- DF_mfpca %>%
    mutate(roundtime = round(2 * time) / 2)
  max_len <- as.integer(max(DF_mfpca$roundtime) * 2 + 1)
  obs_time <- seq(0, max_len / 2, by = 0.5)
  
  #we create an empty data frame, full of NA in the longitudinal outcomes
  #we create an empty data frame
  DF_mfpca2 <- data.frame(id = rep(seq_len(n), each = K*2+1),
                          time = c(replicate(n, obs_time)),
                          y1_2 = c(replicate(n, rep(NA, K*2+1))),
                          y3_2 = c(replicate(n, rep(NA, K*2+1))),
                          y5_2 = c(replicate(n, rep(NA, K*2+1))))
  #we put the values
  # we have to group by id and fill the NAs with values whenever
  #there is info in that time of the grid
  df_filled <- DF_mfpca2 %>%
    left_join(DF_mfpca %>% dplyr::select(id, roundtime, y1, y3, y5), by = c("id", "time" = "roundtime")) %>%
    mutate(
      y1_2 = ifelse(is.na(y1_2), y1, y1_2),
      y3_2 = ifelse(is.na(y3_2), y3, y3_2),
      y5_2 = ifelse(is.na(y5_2), y5, y5_2)
    ) %>%
    dplyr::select(-y1, -y3, -y5)
  
  # now we save in 3 different data sets with long format in order to use
  # funData()
  df_y1_wide <- df_filled %>%
    dplyr::select(id, time, y1_2) %>%
    distinct(id, time, .keep_all = TRUE) %>%   
    pivot_wider(names_from = time, 
                values_from = y1_2, 
                names_prefix = "y1_time_")
  
  df_y3_wide <- df_filled %>%
    dplyr::select(id, time, y3_2) %>%
    distinct(id, time, .keep_all = TRUE) %>%   
    pivot_wider(names_from = time, 
                values_from = y3_2, 
                names_prefix = "y3_time_")
  
  df_y5_wide <- df_filled %>%
    dplyr::select(id, time, y5_2) %>%
    distinct(id, time, .keep_all = TRUE) %>%   
    pivot_wider(names_from = time, 
                values_from = y5_2, 
                names_prefix = "y5_time_")
  return(list(df_y1_wide, df_y3_wide, df_y5_wide, obs_time))
}

simulate_dropout <- function(data) {
  # Find the first time when dropout == 1
  dropout_time <- data %>% filter(dropout == 1) %>% slice(1) %>% pull(time)
  
  # If a dropout is found, set all subsequent longitudinal variables to NA
  if (length(dropout_time) > 0) {
    data <- data %>%
      mutate(across(c(y1, y3, y5), ~ ifelse(time >= dropout_time, NA, .)))
  }
  
  return(data)
}

counting_missing_data <- function(data, n, K){
  return(sum(is.na(data$y1))/(n*K))
}


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



######################################################
######################################################
##Starting the loop: MCAR 15% of dropout ratio
######################################################
######################################################

for(count in 1:num_datasets){
  ############################
  #We generate the dataset
  ###########################
  n <- 200 # number of subjects
  n_test <- 200
  K <- 10 # number of measurements per subject
  t_max <- 10 # maximum follow-up time
  min_sep <- 0.5
  
  # we construct a data frame with the design:
  # everyone has a baseline measurement, and then measurements at random 
  # follow-up times up to t_max 
  # !! for Machine Learning purposes: two consecutive measurements will have
  #more than 0.5 of
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
  
  #X2 <- model.matrix(~ treatment * time, data = DF )
  #Z2 <- model.matrix(~ time, data = DF)
  
  #X3 <- model.matrix(~ time, data = DF )
  #Z3 <- model.matrix(~ time, data = DF)
  
  X3 <- model.matrix(~ sqrt(time), data = DF )
  Z3 <- model.matrix(~ sqrt(time), data = DF)
  
  #X4 <- model.matrix(~ sex + time, data = DF)
  #Z4 <- model.matrix(~ time, data = DF)
  
  #X4 <- model.matrix(~ sex + time , data = DF)
  #Z4 <- model.matrix(~ time, data = DF)
  
  #X5 <- model.matrix(~ time, data = DF)
  #Z5 <- model.matrix(~ time, data = DF)
  
  X5 <- model.matrix(~ treatment + time, data = DF)
  Z5 <- model.matrix(~ time, data = DF)
  
  ####for test data set
  X1_test <- model.matrix(~ sex * time, data = DF_test)
  Z1_test <- model.matrix(~ time, data = DF_test)
  
  #X2_test <- model.matrix(~ treatment * time, data = DF_test )
  #Z2_test <- model.matrix(~ time, data = DF_test)
  
  #X3_test <- model.matrix(~ time, data = DF_test )
  #Z3_test <- model.matrix(~ time, data = DF_test)
  
  X3_test <- model.matrix(~ sqrt(time), data = DF_test )
  Z3_test <- model.matrix(~ sqrt(time), data = DF_test)
  
  #X4_test <- model.matrix(~ sex + time, data = DF_test)
  #Z4_test <- model.matrix(~ time, data = DF_test)
  
  #X4_test <- model.matrix(~ sex + time, data = DF_test)
  #Z4_test <- model.matrix(~ time, data = DF_test)
  
  #X5_test <- model.matrix(~ time, data = DF_test)
  #Z5_test <- model.matrix(~ time, data = DF_test)
  
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
  
  #Simulate second longitudinal outcome
  ###############################################
  #betas2 <- c(-1.8, -0.06, 0.5, 0.06) # fixed effects coefficients
  #betas2 <- c(10, -0.02, 2.222, 0.07) # fixed effects coefficients
  #sigma2 <- 0.25 # errors sd
  
  # we simulate random effects
  #b2 <- b[, c(3,4)]
  # linear predictor
  #eta_y2 <- as.vector(X2 %*% betas2 + rowSums(Z2 * b2[DF$id, ]))
  # we simulate normal longitudinal data
  #DF$y2 <- rnorm(n * K, mean = eta_y2, sd = sigma2)
  # we assume that values below 0 are not observed, and set equal to 0
  #DF$ind2 <- as.numeric(DF$y2 < -4)
  #DF$y2 <- pmax(DF$y2, -4)
  
  ####Test data
  # we simulate random effects
  #b2_test <- b_test[, c(3,4)]
  # linear predictor
  #eta_y2_test <- as.vector(X2_test %*% betas2 + rowSums(Z2_test * b2[DF_test$id, ]))
  # we simulate normal longitudinal data
  #DF_test$y2 <- rnorm(n_test * K, mean = eta_y2_test, sd = sigma2)
  # we assume that values below 0 are not observed, and set equal to 0
  #DF_test$ind2 <- as.numeric(DF_test$y2 < -4)
  #DF_test$y2 <- pmax(DF_test$y2, -4)
  
  #Simulate third longitudinal outcome
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
  
  #Simulate forth longitudinal outcome
  ###############################################
  #betas4 <- c(0.01, 0.5, -0.31416) # fixed effects coefficients
  #betas4 <- c(80, -4.25, -3.789) # fixed effects coefficients
  
  # we simulate random effects
  #b4 <- b[, c(7,8)]
  # linear predictor
  #eta_y4 <- as.vector(X4 %*% betas4 + rowSums(Z4 * b4[DF$id, ]))
  # mean of the binomial distribution
  #mu_y4 <- plogis(eta_y4)
  # we simulate binomial longitudinal data
  #DF$y4 <- rbinom(n * K, size = 1, prob = mu_y4)
  #DF$y4 <- rnorm(n_test * K, mean = eta_y4, sd = sigma3)
  
  #####Test data
  # we simulate random effects
  #b4_test <- b_test[, c(7,8)]
  # linear predictor
  #eta_y4_test <- as.vector(X4_test %*% betas4 + rowSums(Z4_test * b4_test[DF_test$id, ]))
  # mean of the binomial distribution
  #mu_y4_test <- plogis(eta_y4_test)
  # we simulate binomial longitudinal data
  #DF_test$y4 <- rbinom(n_test * K, size = 1, prob = mu_y4_test)
  #DF_test$y4 <- rnorm(n_test * K, mean = eta_y4_test, sd = sigma3)
  
  #Simulate fifth longitudinal outcome
  ###############################################
  betas5 <- c(1, 0.155, 0.12345) # fixed effects coefficients
  betas5 <- c(70, -5.44, -3.761289) # fixed effects coefficients
  
  # we simulate random effects
  b5 <- b[, c(5,6)]
  # linear predictor
  eta_y5 <- as.vector(X5 %*% betas5 + rowSums(Z5 * b5[DF$id, ]))
  # mean of the binomial distribution
  #mu_y5 <- plogis(eta_y5)
  # we simulate binomial longitudinal data
  #DF$y5 <- rbinom(n * K, size = 1, prob = mu_y5)
  DF$y5 <- rnorm(n_test * K, mean = eta_y5, sd = sigma3)
  
  ####Test data
  # we simulate random effects
  b5_test <- b_test[, c(5,6)]
  # linear predictor
  eta_y5_test <- as.vector(X5_test %*% betas5 + rowSums(Z5_test * b5_test[DF_test$id, ]))
  # mean of the binomial distribution
  #mu_y5_test <- plogis(eta_y5_test)
  # we simulate binomial longitudinal data
  #DF_test$y5 <- rbinom(n_test * K, size = 1, prob = mu_y5_test)
  
  DF_test$y5 <- rnorm(n_test * K, mean = eta_y5_test, sd = sigma3)
  
  #If using admin censoring:
  #AdminCens <- 10
  
  # we keep the longitudinal measurements before the event times: TRAIN
  #DF <- DF[DF$time <= AdminCens, ]

  # we keep the longitudinal measurements before the event times: TEST
  #DF_test <- DF_test[DF_test$time <= AdminCens, ]
  
  #Data set with just one obs per individual
  #DF.id <- DF[!duplicated(DF$id),]
  #DF_test.id <- DF_test[!duplicated(DF_test$id),]
  
  
  ### Here, we should program the DROPOUT!!!
  param1 <- -4
  param2 <- 0.2
  probs_vec <- numeric(n*K)
  kk <- 1
  for(i in 1:n){
    for(j in 1:K){
      if((kk-1)%%K==0){
        probs_vec[kk] <- 0
      } else{
        thet <- param1 + param2*DF[DF$id==i,]$time[j]
        probs_vec[kk] <- 1/(1+exp(-thet))
      }
      kk <- kk + 1
    }
  }
  
  DF$probs_drop <- probs_vec
  DF$dropout <- rbinom(n*K, size = 1, prob = probs_vec)
  
  DF_test$probs_dropout <- probs_vec
  DF_test$dropout <- rbinom(n*K, size = 1, prob = probs_vec)
  
  ###Putting NAs after dropout
  #DF_complete <- DF
  DF_miss <- DF %>%
    group_by(id) %>%
    group_modify(~ simulate_dropout(.x)) %>%
    ungroup()
  
  DF_test_miss <- DF_test %>%
    group_by(id) %>%
    group_modify(~ simulate_dropout(.x)) %>%
    ungroup()
  
  
  try(missing_ratio_train[count] <- counting_missing_data(DF_miss, n, K))
  try(missing_ratio_test[count] <- counting_missing_data(DF_test_miss, n, K))
  
  
  ###OBS!!!!!
  ## DF_miss contiene el data frame con missing data!!
  
  ######################################
  #Fitting multivariate MM models
  # - We use snat_mvmer() from rstanarm
  # - We fit the true model plus two
  #.  (or 3) more models
  ######################################
  
  #checking time of computation using data frame with missings 
  # when n=300, k=10, we have 21mins, max rhat 1.67
  #      n=250, k=10: 17mins, max rhat 1.51
  #      n=200, k=10:  8.6 mins, max rhat 1.09
  try(true_model_miss <- stan_mvmer(
    formula = list(
      y1 ~ sex * time + (time | id),
      y3 ~ sqrt(time) + (sqrt(time) | id),
      y5 ~ treatment + time + (time | id)),
    data = DF_miss,
    family = list(gaussian, gaussian, gaussian),
    chains = 3, cores = 8, seed = 12345, iter = 1500))
  
  try(rhat_true_model_miss <- rhat(true_model_miss))
  
  #OJO!!! para que devuelva las predictions a pesar de estar pasandole NAs 
  #hay que pasarle el data frame sin los NAs
  
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
  DF_miss2 <- DF_miss %>% dplyr::select(id, time, visit, sex, treatment)
  
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
  DF_test_miss2 <- DF_test_miss %>% dplyr::select(id, time, visit, sex, treatment)
  
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
  grid_long <- grid_longitudinal_data(DF_miss, n, K)
  #univariate functional data
  obs_time <- grid_long[[4]]
  f1 <- funData(obs_time, as.matrix(grid_long[[1]][,-1]))
  f2 <- funData(obs_time, as.matrix(grid_long[[2]][,-1]))
  f3 <- funData(obs_time, as.matrix(grid_long[[3]][,-1]))
  #multivariate functional data class
  m1 <- multiFunData(list(f1,f2,f3))
  
  
  # Grid longitudinal data in TEST:
  grid_long_test <- grid_longitudinal_data(DF_test, n, K)
  #univariate functional data TEST case
  obs_time_test <- grid_long_test[[4]]
  f1_test <- funData(obs_time_test, as.matrix(grid_long_test[[1]][,-1]))
  f2_test <- funData(obs_time_test, as.matrix(grid_long_test[[2]][,-1]))
  f3_test <- funData(obs_time_test, as.matrix(grid_long_test[[3]][,-1]))
  #multivariate functional data class TEST case
  m1_test <- multiFunData(list(f1_test,f2_test,f3_test))
  
  #now we use MFPCA function from mfpca R package
  
  #We fit with NO weights and pve=0.9
  try(mfpca1 <- MFPCA(m1, M = 3, 
                 uniExpansions = list(list(type="uFPCA", pve = 0.9),
                                      list(type ="uFPCA", pve=0.9),
                                      list(type="uFPCA", pve=0.9)),
                 fit = TRUE))
  
  
  ###OBS: we have to adjust M and pve depending on % of missing data we have 
  
  # predict in training data:
  try(pred_mfpca1 <- predict(mfpca1))
  
  try(pred_mfpca1_y1 <- pred_mfpca1[[1]]@X)
  try(pred_mfpca1_y2 <- pred_mfpca1[[2]]@X)
  try(pred_mfpca1_y3 <- pred_mfpca1[[3]]@X)
  
  #plot(pred_mfpca1)
  
  # predict in testing data:
  # value "vector" returned by MFPCA seems to be the cms in paper (and ML code)
  # we compute scores for testing data and we use meanf and psi from training
  try(mfpca1_test <- MFPCA(m1_test, M = 3, 
                  uniExpansions = list(list(type="uFPCA", pve = 0.9),
                                       list(type ="uFPCA", pve=0.9),
                                       list(type="uFPCA", pve=0.9)),
                  fit = TRUE))
  
  try(pred_mfpca1_test <- predict(mfpca1, scores = mfpca1_test$scores))
  
  try(pred_mfpca1_y1_test <- pred_mfpca1_test[[1]]@X)
  try(pred_mfpca1_y2_test <- pred_mfpca1_test[[2]]@X)
  try(pred_mfpca1_y3_test <- pred_mfpca1_test[[3]]@X)
  
  #plot(pred_mfpca1_test)
  
  #AQUI YA TENGO TODA LA INFO QUE QUIERO! las predicciones de esto último se 
  #hacen para cada uno de los puntos de la grid!! También para aquellos donde
  #hay NA! Por tanto para aquellos puntos estoy haciendo una predicción!
  #Y es esa predicción la que me tengo que quedar y entrar a comparar con 
  #la predicción que me puedan devolver los MM
  
  #Cuando tenemos el predict(MFPCAfit), el objeto que pred_mfpca1[[1]]@X es
  # una matriz n X max_length +1, ie cada fila es un individuo y cada columna
  # es uno de los puntos de la fixed time grid. Y esto para el outcome 
  # longitudinal 1 
  # Lo que tenemos aquí es una predicción para cada punto de la fixed grid
  # del outcome longitudinal 1.
  
  
  
  #####################################################
  ## Comparing predictions by means of standardized
  ## RMSE
  #####################################################
  
  
  
  
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
  
  
  try(list_DF <- append(list_DF, list(DF)))
  try(list_DF_miss <- append(list_DF_miss, list(DF_miss)))
  try(list_DF_test <- append(list_DF_test, list(DF_test)))
  try(list_DF_test_miss <- append(list_DF_test_miss, list(DF_test_miss)))
  
  
  setwd("D:/La meva unitat/TFM/ResultsMMvsMFPCA")
  try(if(count==1){
    strr <- "results_1_MCAR_28jan2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, 
         list_DF, list_DF_miss, list_DF_test, list_DF_test_miss, 
         missing_ratio_train, missing_ratio_test, file=strr)
  })
  try(if(count==2){
    strr2 <- "results_2_MCAR_28jan2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss,
         list_DF, list_DF_miss, list_DF_test, list_DF_test_miss, 
         missing_ratio_train, missing_ratio_test, file=strr2)
  })
  try(if(count==5){
    str1 <- "results_5_MCAR_28jan2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, 
         list_DF, list_DF_miss, list_DF_test, list_DF_test_miss, 
         missing_ratio_train, missing_ratio_test, file=str1)
  })
  try(if(count==10){
    str10 <- "results_10_MCAR_28jan2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, 
         list_DF, list_DF_miss, list_DF_test, list_DF_test_miss, 
         missing_ratio_train, missing_ratio_test, file=str10)
  })
  try(if(count==20){
    str20 <- "results_20_MCAR_28jan2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, 
         list_DF, list_DF_miss, list_DF_test, list_DF_test_miss, 
         missing_ratio_train, missing_ratio_test, file=str20)
  })
  try(if(count==30){
    str30 <- "results_30_MCAR_28jan2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, 
         list_DF, list_DF_miss, list_DF_test, list_DF_test_miss, 
         missing_ratio_train, missing_ratio_test, file=str30)
  })
  try(if(count==40){
    str40 <- "results_40_MCAR_28jan2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, 
         list_DF, list_DF_miss, list_DF_test, list_DF_test_miss, 
         missing_ratio_train, missing_ratio_test, file=str40)
  })
  try(if(count==50){
    str50 <- "results_50_MCAR_28jan2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, 
         list_DF, list_DF_miss, list_DF_test, list_DF_test_miss, 
         missing_ratio_train, missing_ratio_test, file=str50)
  })
  try(if(count==75){
    str75 <- "results_75_MCAR_28jan2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, 
         list_DF, list_DF_miss, list_DF_test, list_DF_test_miss, 
         missing_ratio_train, missing_ratio_test, file=str75)
  })
  try(if(count==100){
    str100 <- "results_100_MCAR_28jan2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, 
         list_DF, list_DF_miss, list_DF_test, list_DF_test_miss, 
         missing_ratio_train, missing_ratio_test, file=str100)
  })
  
  
  
  ######################################
  #END FOR
  ######################################
}



###############################################################
## Let us do some checks with the longitudinal data simulated 
## via graphics
###############################################################
library(lattice)
xyplot(y1 ~ time, groups = id,
       data = DF[1:825,],
       type = "l" ,xlab="Time",ylab="Y1")

xyplot(y1 ~ time, groups = id,
       data = DF_miss[1:825,],
       type = "l" ,xlab="Time",ylab="Y1")

xyplot(y1 ~ time, groups = id,
       data = DF,
       type = "l" ,xlab="Time",ylab="Y1")

#xyplot(y2 ~ time, groups = id,
#       data = DF[1:825,],
#       type = "l" ,xlab="Time",ylab="Y2")

xyplot(y3 ~ time, groups = id,
       data = DF[1:825,],
       type = "l" ,xlab="Time",ylab="Y3")

xyplot(y3 ~ time, groups = id,
       data = DF_miss[1:825,],
       type = "l" ,xlab="Time",ylab="Y3")

xyplot(y3 ~ time, groups = id,
       data = DF,
       type = "l" ,xlab="Time",ylab="Y3")

#xyplot(y4 ~ time, groups = id,
#       data = DF[1:825,],
#       type = "l" ,xlab="Time",ylab="Y4")

xyplot(y5 ~ time, groups = id,
       data = DF[1:825,],
       type = "l" ,xlab="Time",ylab="Y5")

xyplot(y5 ~ time, groups = id,
       data = DF_miss[1:825,],
       type = "l" ,xlab="Time",ylab="Y5")

xyplot(y5 ~ time, groups = id,
       data = DF,
       type = "l" ,xlab="Time",ylab="Y5")

xyplot(y5 ~ time, groups = id,
       data = DF_miss,
       type = "l" ,xlab="Time",ylab="Y5")

## checking spaghetti plot in test data
xyplot(y1 ~ time, groups = id,
       data = DF_test[1:825,],
       type = "l" ,xlab="Time",ylab="Y1")

xyplot(y3 ~ time, groups = id,
       data = DF_test[1:825,],
       type = "l" ,xlab="Time",ylab="Y3")

xyplot(y5 ~ time, groups = id,
       data = DF_test[1:825,],
       type = "l" ,xlab="Time",ylab="Y5")


# Comparisons simulated data with MFPCA
xyplot(y1 ~ time, groups = id,
       data = DF,
       type = "l" ,xlab="Time",ylab="Y1")
plot(pred_mfpca1[[1]], add = TRUE, lty = 2)

plot(pred_mfpca1[[1]])




#####Observation generation function, to be used if we do not have enough data
# because time of obs lie within a time window less than 0.5, and then a lot
# of data is ignored
obstime_gen <- function(K, t_max, min_sep){
  obstime <- numeric(K)
  obstime <- c(0, sort(runif(K - 1, 0, t_max)))
  for(i in 2:K){
    alpha <- obstime[i] - obstime[i-1]
    if(alpha < min_sep){
      obstime[i] <- obstime[i] + 0.51 - alpha
    }
  }
  return(obstime)
}

##### Function to put in the correct format the longitudinal data to use MFPCA
get_numpy <- function(df, long = c("Y1", "Y2", "Y3"), base = c("X1", "X2"), obstime = "obstime", max_len = NULL) {
  
  # Step 1: Assign a new id ('id_new') to each subject
  df <- df %>%
    group_by(id) %>%
    mutate(id_new = as.numeric(as.factor(id)) - 1) # 'id_new' assigned from 0 to num subjects
  
  # Step 2: Round 'obstime' to the nearest 0.5 and store it in 'roundtime'
  df <- df %>%
    mutate(roundtime = round(2 * !!sym(obstime)) / 2)
  
  # Step 3: Create 'visit' column if it doesn't exist
  if (!"visit" %in% names(df)) {
    df <- df %>%
      group_by(id) %>%
      mutate(visit = row_number() - 1)
  }
  
  # Step 4: Number of unique subjects
  I <- n_distinct(df$id)
  
  # Step 5: Calculate max_len if it's NULL
  if (is.null(max_len)) {
    max_len <- as.integer(max(df$roundtime) * 2 + 1) # based on 0.5 rounding
  }
  
  # Step 6: Initialize 'x_long' (3D array) and 'x_base' (2D matrix)
  x_long <- array(NA, dim = c(I, max_len, length(long)))
  x_base <- matrix(0, nrow = I, ncol = length(base))
  
  # Step 7: Populate 'x_long' and 'x_base'
  df %>%
    group_by(id) %>%
    rowwise() %>%
    do({
      ii <- .$id_new[1] + 1 # R is 1-based indexing
      jj <- as.integer(.$roundtime[1] * 2) + 1 # Adjust for R's 1-based indexing
      
      if (!is.na(.$visit) && .$visit == 0) {
        x_base[ii, ] <- unlist(.[base])
      }
      
      x_long[ii, jj, ] <- unlist(.[long])
      return(NULL)
    })
  
  
  # Step 9: Create 'obs_time' (sequence from 0 to max_len / 2, by 0.5)
  obs_time <- seq(0, max_len / 2, by = 0.5)
  
  return(list(x_long = x_long, x_base = x_base, obs_time = obs_time))
}

#################################################################
## Calibrating parameters for dropout
#################################################################
param1 <- -2
param2 <- 0.005
probs_vec <- numeric(n*K)
kk <- 1
for(i in 1:n){
  for(j in 1:K){
    thet <- param1 + param2*DF[DF$id==i,]$time[j]
    probs_vec[kk] <- 1/(1+exp(-thet))
    kk <- kk + 1
  }
}
probs_vec[1:15]

simusimu <- rbinom(n*K, size = 1, prob = probs_vec)

DF$dropout <- simusimu



########################################################################
#####################################################################
## stuff out of the loop
########################################################################
#####################################################################



# fitting the true model
true_model <- stan_mvmer(
  formula = list(
    y1 ~ sex * time + (time | id),
    y3 ~ sqrt(time) + (sqrt(time) | id),
    y5 ~ treatment + time + (time | id)),
  data = DF,
  family = list(gaussian, gaussian, gaussian),
  chains = 3, cores = 8, seed = 12345, iter = 1000)

#Save convergence stuff
#vector with rhats
rhat_true_model <- rhat(true_model)

#we should check how all params have rhat<1.1


#PREDICTIONS (for three outcomes separately)
preds_tmod <- posterior_predict(true_model, m=1)
preds_tmod_y2 <- posterior_predict(true_model, m=2)
preds_tmod_y3 <- posterior_predict(true_model, m=3)

## we have to use posterior_predict() function to do the predictions, also
# we have to use m=l to indicate the longitudinal outcome we want to predict
## with this we have 1500 (in this case) samples of the 4500 longitudinal values
## for each longitudinal outcome. We can get the mean, and then work 
## with the point estimate of the prediction. Without mean we have the predictive
## posterior distribution of the values.

#posteriorss <- as.matrix(true_model)

#mcmc_trace(posteriorss[, "y1|(Intercept)"])

## We must save metrics to check the convergence of the model (rhat and others
## maybe even some chains)

## We do predictions with training (DF) and testing (DF_test) data

#using predict() function

mMM1 <- stan_mvmer(
  formula = list(
    y1 ~ sex + time + (time | id),
    y3 ~ treatment + time^2 + (time | id),
    y5 ~ time + (time | id)),
  data = DF,
  family = list(gaussian, gaussian, gaussian),
  chains = 3, cores = 2, seed = 12345, iter = 1000)

## We must save metrics to check the convergence of the model (rhat and others
## maybe even some chains)

## We do predictions with training (DF) and testing (DF_test) data

#using predict() function

mMM2 <- stan_mvmer(
  formula = list(
    y1 ~ time + (time | id),
    y3 ~ treatment + (time | id),
    y5 ~ sex + (time | id)),
  data = DF,
  family = list(gaussian, gaussian, gaussian),
  chains = 3, cores = 2, seed = 12345, iter = 1000)

## We must save metrics to check the convergence of the model (rhat and others
## maybe even some chains)

## We do predictions with training (DF) and testing (DF_test) data

#using predict() function

mMM3 <- stan_mvmer(
  formula = list(
    y1 ~ treatment*time + (time | id),
    y3 ~ sex + sqrt(time) + (time | id),
    y5 ~ sex + time + (time | id)),
  data = DF,
  family = list(gaussian, gaussian, gaussian),
  chains = 3, cores = 2, seed = 12345, iter = 1000)

## We must save metrics to check the convergence of the model (rhat and others
## maybe even some chains)

## We do predictions with training (DF) and testing (DF_test) data

#using predict() function

load("D:/La meva unitat/TFM/ResultsMMvsMFPCA/results_2_MCAR_23jan2025.RData")







###############################################################################
###############################################################################
##Starting the loop, MCAR 30% of dropout
###############################################################################
###############################################################################
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
## I generate the data frames
############################################################
############################################################
for(count in 1:num_datasets){
  n <- 200 # number of subjects
  n_test <- 200
  K <- 10 # number of measurements per subject
  t_max <- 10 # maximum follow-up time
  min_sep <- 0.5
  
  # we construct a data frame with the design:
  # everyone has a baseline measurement, and then measurements at random 
  # follow-up times up to t_max 
  # !! for Machine Learning purposes: two consecutive measurements will have
  #more than 0.5 of
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
  
  #X2 <- model.matrix(~ treatment * time, data = DF )
  #Z2 <- model.matrix(~ time, data = DF)
  
  #X3 <- model.matrix(~ time, data = DF )
  #Z3 <- model.matrix(~ time, data = DF)
  
  X3 <- model.matrix(~ sqrt(time), data = DF )
  Z3 <- model.matrix(~ sqrt(time), data = DF)
  
  #X4 <- model.matrix(~ sex + time, data = DF)
  #Z4 <- model.matrix(~ time, data = DF)
  
  #X4 <- model.matrix(~ sex + time , data = DF)
  #Z4 <- model.matrix(~ time, data = DF)
  
  #X5 <- model.matrix(~ time, data = DF)
  #Z5 <- model.matrix(~ time, data = DF)
  
  X5 <- model.matrix(~ treatment + time, data = DF)
  Z5 <- model.matrix(~ time, data = DF)
  
  ####for test data set
  X1_test <- model.matrix(~ sex * time, data = DF_test)
  Z1_test <- model.matrix(~ time, data = DF_test)
  
  #X2_test <- model.matrix(~ treatment * time, data = DF_test )
  #Z2_test <- model.matrix(~ time, data = DF_test)
  
  #X3_test <- model.matrix(~ time, data = DF_test )
  #Z3_test <- model.matrix(~ time, data = DF_test)
  
  X3_test <- model.matrix(~ sqrt(time), data = DF_test )
  Z3_test <- model.matrix(~ sqrt(time), data = DF_test)
  
  #X4_test <- model.matrix(~ sex + time, data = DF_test)
  #Z4_test <- model.matrix(~ time, data = DF_test)
  
  #X4_test <- model.matrix(~ sex + time, data = DF_test)
  #Z4_test <- model.matrix(~ time, data = DF_test)
  
  #X5_test <- model.matrix(~ time, data = DF_test)
  #Z5_test <- model.matrix(~ time, data = DF_test)
  
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
  
  #Simulate second longitudinal outcome
  ###############################################
  #betas2 <- c(-1.8, -0.06, 0.5, 0.06) # fixed effects coefficients
  #betas2 <- c(10, -0.02, 2.222, 0.07) # fixed effects coefficients
  #sigma2 <- 0.25 # errors sd
  
  # we simulate random effects
  #b2 <- b[, c(3,4)]
  # linear predictor
  #eta_y2 <- as.vector(X2 %*% betas2 + rowSums(Z2 * b2[DF$id, ]))
  # we simulate normal longitudinal data
  #DF$y2 <- rnorm(n * K, mean = eta_y2, sd = sigma2)
  # we assume that values below 0 are not observed, and set equal to 0
  #DF$ind2 <- as.numeric(DF$y2 < -4)
  #DF$y2 <- pmax(DF$y2, -4)
  
  ####Test data
  # we simulate random effects
  #b2_test <- b_test[, c(3,4)]
  # linear predictor
  #eta_y2_test <- as.vector(X2_test %*% betas2 + rowSums(Z2_test * b2[DF_test$id, ]))
  # we simulate normal longitudinal data
  #DF_test$y2 <- rnorm(n_test * K, mean = eta_y2_test, sd = sigma2)
  # we assume that values below 0 are not observed, and set equal to 0
  #DF_test$ind2 <- as.numeric(DF_test$y2 < -4)
  #DF_test$y2 <- pmax(DF_test$y2, -4)
  
  #Simulate third longitudinal outcome
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
  
  #Simulate forth longitudinal outcome
  ###############################################
  #betas4 <- c(0.01, 0.5, -0.31416) # fixed effects coefficients
  #betas4 <- c(80, -4.25, -3.789) # fixed effects coefficients
  
  # we simulate random effects
  #b4 <- b[, c(7,8)]
  # linear predictor
  #eta_y4 <- as.vector(X4 %*% betas4 + rowSums(Z4 * b4[DF$id, ]))
  # mean of the binomial distribution
  #mu_y4 <- plogis(eta_y4)
  # we simulate binomial longitudinal data
  #DF$y4 <- rbinom(n * K, size = 1, prob = mu_y4)
  #DF$y4 <- rnorm(n_test * K, mean = eta_y4, sd = sigma3)
  
  #####Test data
  # we simulate random effects
  #b4_test <- b_test[, c(7,8)]
  # linear predictor
  #eta_y4_test <- as.vector(X4_test %*% betas4 + rowSums(Z4_test * b4_test[DF_test$id, ]))
  # mean of the binomial distribution
  #mu_y4_test <- plogis(eta_y4_test)
  # we simulate binomial longitudinal data
  #DF_test$y4 <- rbinom(n_test * K, size = 1, prob = mu_y4_test)
  #DF_test$y4 <- rnorm(n_test * K, mean = eta_y4_test, sd = sigma3)
  
  #Simulate fifth longitudinal outcome
  ###############################################
  betas5 <- c(1, 0.155, 0.12345) # fixed effects coefficients
  betas5 <- c(70, -5.44, -3.761289) # fixed effects coefficients
  
  # we simulate random effects
  b5 <- b[, c(5,6)]
  # linear predictor
  eta_y5 <- as.vector(X5 %*% betas5 + rowSums(Z5 * b5[DF$id, ]))
  # mean of the binomial distribution
  #mu_y5 <- plogis(eta_y5)
  # we simulate binomial longitudinal data
  #DF$y5 <- rbinom(n * K, size = 1, prob = mu_y5)
  DF$y5 <- rnorm(n_test * K, mean = eta_y5, sd = sigma3)
  
  ####Test data
  # we simulate random effects
  b5_test <- b_test[, c(5,6)]
  # linear predictor
  eta_y5_test <- as.vector(X5_test %*% betas5 + rowSums(Z5_test * b5_test[DF_test$id, ]))
  # mean of the binomial distribution
  #mu_y5_test <- plogis(eta_y5_test)
  # we simulate binomial longitudinal data
  #DF_test$y5 <- rbinom(n_test * K, size = 1, prob = mu_y5_test)
  
  DF_test$y5 <- rnorm(n_test * K, mean = eta_y5_test, sd = sigma3)
  
  #If using admin censoring:
  #AdminCens <- 10
  
  # we keep the longitudinal measurements before the event times: TRAIN
  #DF <- DF[DF$time <= AdminCens, ]
  
  # we keep the longitudinal measurements before the event times: TEST
  #DF_test <- DF_test[DF_test$time <= AdminCens, ]
  
  #Data set with just one obs per individual
  #DF.id <- DF[!duplicated(DF$id),]
  #DF_test.id <- DF_test[!duplicated(DF_test$id),]
  
  
  ### Here, we should program the DROPOUT!!!
  param1 <- -4
  param2 <- 0.4
  probs_vec <- numeric(n*K)
  kk <- 1
  for(i in 1:n){
    for(j in 1:K){
      if((kk-1)%%K==0){
        probs_vec[kk] <- 0
      } else{
        thet <- param1 + param2*DF[DF$id==i,]$time[j]
        probs_vec[kk] <- 1/(1+exp(-thet))
      }
      kk <- kk + 1
    }
  }
  
  DF$probs_drop <- probs_vec
  DF$dropout <- rbinom(n*K, size = 1, prob = probs_vec)
  
  DF_test$probs_dropout <- probs_vec
  DF_test$dropout <- rbinom(n*K, size = 1, prob = probs_vec)
  
  ###Putting NAs after dropout
  #DF_complete <- DF
  DF_miss <- DF %>%
    group_by(id) %>%
    group_modify(~ simulate_dropout(.x)) %>%
    ungroup()
  
  DF_test_miss <- DF_test %>%
    group_by(id) %>%
    group_modify(~ simulate_dropout(.x)) %>%
    ungroup()
  
  missing_ratio_train[count] <- counting_missing_data(DF_miss, n, K)
  missing_ratio_test[count] <- counting_missing_data(DF_test_miss, n, K)
  
  try(list_DF <- append(list_DF, list(DF)))
  try(list_DF_miss <- append(list_DF_miss, list(DF_miss)))
  try(list_DF_test <- append(list_DF_test, list(DF_test)))
  try(list_DF_test_miss <- append(list_DF_test_miss, list(DF_test_miss)))
  
  try(if(count==100){
    str100 <- "dataframes_MCAR_30.RData"
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

  ###OBS!!!!!
  ## DF_miss contiene el data frame con missing data!!
  
  ######################################
  #Fitting multivariate MM models
  # - We use snat_mvmer() from rstanarm
  # - We fit the true model plus two
  #.  (or 3) more models
  ######################################
  
  #checking time of computation using data frame with missings 
  # when n=300, k=10, we have 21mins, max rhat 1.67
  #      n=250, k=10: 17mins, max rhat 1.51
  #      n=200, k=10:  8.6 mins, max rhat 1.09
  try(true_model_miss <- stan_mvmer(
    formula = list(
      y1 ~ sex * time + (time | id),
      y3 ~ sqrt(time) + (sqrt(time) | id),
      y5 ~ treatment + time + (time | id)),
    data = list_DF_miss[[count]],
    family = list(gaussian, gaussian, gaussian),
    chains = 3, cores = 8, seed = 12345, iter = 1500))
  
  try(rhat_true_model_miss <- rhat(true_model_miss))
  
  #OJO!!! para que devuelva las predictions a pesar de estar pasandole NAs 
  #hay que pasarle el data frame sin los NAs
  
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
  
  #We fit with NO weights and pve=0.9
  try(mfpca1 <- MFPCA(m1, M = 3, 
                      uniExpansions = list(list(type="uFPCA", pve = 0.8),
                                           list(type ="uFPCA", pve=0.8),
                                           list(type="uFPCA", pve=0.8)),
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
  
  #plot(pred_mfpca1_test)
  
  #AQUI YA TENGO TODA LA INFO QUE QUIERO! las predicciones de esto último se 
  #hacen para cada uno de los puntos de la grid!! También para aquellos donde
  #hay NA! Por tanto para aquellos puntos estoy haciendo una predicción!
  #Y es esa predicción la que me tengo que quedar y entrar a comparar con 
  #la predicción que me puedan devolver los MM
  
  #Cuando tenemos el predict(MFPCAfit), el objeto que pred_mfpca1[[1]]@X es
  # una matriz n X max_length +1, ie cada fila es un individuo y cada columna
  # es uno de los puntos de la fixed time grid. Y esto para el outcome 
  # longitudinal 1 
  # Lo que tenemos aquí es una predicción para cada punto de la fixed grid
  # del outcome longitudinal 1.
  
  
  
  #####################################################
  ## Comparing predictions by means of standardized
  ## RMSE
  #####################################################
  
  
  
  
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
    strr <- "results_1_MCAR30_5feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr)
  })
  try(if(count==2){
    strr2 <- "results_2_MCAR30_5feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr2)
  })
  try(if(count==5){
    str1 <- "results_5_MCAR30_5feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str1)
  })
  try(if(count==10){
    str10 <- "results_10_MCAR30_5feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str10)
  })
  try(if(count==20){
    str20 <- "results_20_MCAR30_5feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str20)
  })
  try(if(count==30){
    str30 <- "results_30_MCAR30_5feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str30)
  })
  try(if(count==40){
    str40 <- "results_40_MCAR30_5feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss,file=str40)
  })
  try(if(count==50){
    str50 <- "results_50_MCAR30_5feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str50)
  })
  try(if(count==75){
    str75 <- "results_75_MCAR30_5feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str75)
  })
  try(if(count==100){
    str100 <- "results_100_MCAR30_5feb2025.RData"
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





###############################################################################
###############################################################################
##MCAR 50% of dropout
###############################################################################
###############################################################################
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
## I generate the data frames: MCAR 50% of dropout
############################################################
############################################################
for(count in 1:num_datasets){
  n <- 200 # number of subjects
  n_test <- 200
  K <- 10 # number of measurements per subject
  t_max <- 10 # maximum follow-up time
  min_sep <- 0.5
  
  # we construct a data frame with the design:
  # everyone has a baseline measurement, and then measurements at random 
  # follow-up times up to t_max 
  # !! for Machine Learning purposes: two consecutive measurements will have
  #more than 0.5 of
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
  
  #X2 <- model.matrix(~ treatment * time, data = DF )
  #Z2 <- model.matrix(~ time, data = DF)
  
  #X3 <- model.matrix(~ time, data = DF )
  #Z3 <- model.matrix(~ time, data = DF)
  
  X3 <- model.matrix(~ sqrt(time), data = DF )
  Z3 <- model.matrix(~ sqrt(time), data = DF)
  
  #X4 <- model.matrix(~ sex + time, data = DF)
  #Z4 <- model.matrix(~ time, data = DF)
  
  #X4 <- model.matrix(~ sex + time , data = DF)
  #Z4 <- model.matrix(~ time, data = DF)
  
  #X5 <- model.matrix(~ time, data = DF)
  #Z5 <- model.matrix(~ time, data = DF)
  
  X5 <- model.matrix(~ treatment + time, data = DF)
  Z5 <- model.matrix(~ time, data = DF)
  
  ####for test data set
  X1_test <- model.matrix(~ sex * time, data = DF_test)
  Z1_test <- model.matrix(~ time, data = DF_test)
  
  #X2_test <- model.matrix(~ treatment * time, data = DF_test )
  #Z2_test <- model.matrix(~ time, data = DF_test)
  
  #X3_test <- model.matrix(~ time, data = DF_test )
  #Z3_test <- model.matrix(~ time, data = DF_test)
  
  X3_test <- model.matrix(~ sqrt(time), data = DF_test )
  Z3_test <- model.matrix(~ sqrt(time), data = DF_test)
  
  #X4_test <- model.matrix(~ sex + time, data = DF_test)
  #Z4_test <- model.matrix(~ time, data = DF_test)
  
  #X4_test <- model.matrix(~ sex + time, data = DF_test)
  #Z4_test <- model.matrix(~ time, data = DF_test)
  
  #X5_test <- model.matrix(~ time, data = DF_test)
  #Z5_test <- model.matrix(~ time, data = DF_test)
  
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
  
  #Simulate second longitudinal outcome
  ###############################################
  #betas2 <- c(-1.8, -0.06, 0.5, 0.06) # fixed effects coefficients
  #betas2 <- c(10, -0.02, 2.222, 0.07) # fixed effects coefficients
  #sigma2 <- 0.25 # errors sd
  
  # we simulate random effects
  #b2 <- b[, c(3,4)]
  # linear predictor
  #eta_y2 <- as.vector(X2 %*% betas2 + rowSums(Z2 * b2[DF$id, ]))
  # we simulate normal longitudinal data
  #DF$y2 <- rnorm(n * K, mean = eta_y2, sd = sigma2)
  # we assume that values below 0 are not observed, and set equal to 0
  #DF$ind2 <- as.numeric(DF$y2 < -4)
  #DF$y2 <- pmax(DF$y2, -4)
  
  ####Test data
  # we simulate random effects
  #b2_test <- b_test[, c(3,4)]
  # linear predictor
  #eta_y2_test <- as.vector(X2_test %*% betas2 + rowSums(Z2_test * b2[DF_test$id, ]))
  # we simulate normal longitudinal data
  #DF_test$y2 <- rnorm(n_test * K, mean = eta_y2_test, sd = sigma2)
  # we assume that values below 0 are not observed, and set equal to 0
  #DF_test$ind2 <- as.numeric(DF_test$y2 < -4)
  #DF_test$y2 <- pmax(DF_test$y2, -4)
  
  #Simulate third longitudinal outcome
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
  
  #Simulate forth longitudinal outcome
  ###############################################
  #betas4 <- c(0.01, 0.5, -0.31416) # fixed effects coefficients
  #betas4 <- c(80, -4.25, -3.789) # fixed effects coefficients
  
  # we simulate random effects
  #b4 <- b[, c(7,8)]
  # linear predictor
  #eta_y4 <- as.vector(X4 %*% betas4 + rowSums(Z4 * b4[DF$id, ]))
  # mean of the binomial distribution
  #mu_y4 <- plogis(eta_y4)
  # we simulate binomial longitudinal data
  #DF$y4 <- rbinom(n * K, size = 1, prob = mu_y4)
  #DF$y4 <- rnorm(n_test * K, mean = eta_y4, sd = sigma3)
  
  #####Test data
  # we simulate random effects
  #b4_test <- b_test[, c(7,8)]
  # linear predictor
  #eta_y4_test <- as.vector(X4_test %*% betas4 + rowSums(Z4_test * b4_test[DF_test$id, ]))
  # mean of the binomial distribution
  #mu_y4_test <- plogis(eta_y4_test)
  # we simulate binomial longitudinal data
  #DF_test$y4 <- rbinom(n_test * K, size = 1, prob = mu_y4_test)
  #DF_test$y4 <- rnorm(n_test * K, mean = eta_y4_test, sd = sigma3)
  
  #Simulate fifth longitudinal outcome
  ###############################################
  betas5 <- c(1, 0.155, 0.12345) # fixed effects coefficients
  betas5 <- c(70, -5.44, -3.761289) # fixed effects coefficients
  
  # we simulate random effects
  b5 <- b[, c(5,6)]
  # linear predictor
  eta_y5 <- as.vector(X5 %*% betas5 + rowSums(Z5 * b5[DF$id, ]))
  # mean of the binomial distribution
  #mu_y5 <- plogis(eta_y5)
  # we simulate binomial longitudinal data
  #DF$y5 <- rbinom(n * K, size = 1, prob = mu_y5)
  DF$y5 <- rnorm(n_test * K, mean = eta_y5, sd = sigma3)
  
  ####Test data
  # we simulate random effects
  b5_test <- b_test[, c(5,6)]
  # linear predictor
  eta_y5_test <- as.vector(X5_test %*% betas5 + rowSums(Z5_test * b5_test[DF_test$id, ]))
  # mean of the binomial distribution
  #mu_y5_test <- plogis(eta_y5_test)
  # we simulate binomial longitudinal data
  #DF_test$y5 <- rbinom(n_test * K, size = 1, prob = mu_y5_test)
  
  DF_test$y5 <- rnorm(n_test * K, mean = eta_y5_test, sd = sigma3)
  
  #If using admin censoring:
  #AdminCens <- 10
  
  # we keep the longitudinal measurements before the event times: TRAIN
  #DF <- DF[DF$time <= AdminCens, ]
  
  # we keep the longitudinal measurements before the event times: TEST
  #DF_test <- DF_test[DF_test$time <= AdminCens, ]
  
  #Data set with just one obs per individual
  #DF.id <- DF[!duplicated(DF$id),]
  #DF_test.id <- DF_test[!duplicated(DF_test$id),]
  
  
  ### Here, we should program the DROPOUT!!!
  param1 <- -4
  param2 <- 0.675
  probs_vec <- numeric(n*K)
  kk <- 1
  for(i in 1:n){
    for(j in 1:K){
      if((kk-1)%%K==0){
        probs_vec[kk] <- 0
      } else{
        thet <- param1 + param2*DF[DF$id==i,]$time[j]
        probs_vec[kk] <- 1/(1+exp(-thet))
      }
      kk <- kk + 1
    }
  }
  
  DF$probs_drop <- probs_vec
  DF$dropout <- rbinom(n*K, size = 1, prob = probs_vec)
  
  DF_test$probs_dropout <- probs_vec
  DF_test$dropout <- rbinom(n*K, size = 1, prob = probs_vec)
  
  ###Putting NAs after dropout
  #DF_complete <- DF
  DF_miss <- DF %>%
    group_by(id) %>%
    group_modify(~ simulate_dropout(.x)) %>%
    ungroup()
  
  DF_test_miss <- DF_test %>%
    group_by(id) %>%
    group_modify(~ simulate_dropout(.x)) %>%
    ungroup()
  
  missing_ratio_train[count] <- counting_missing_data(DF_miss, n, K)
  missing_ratio_test[count] <- counting_missing_data(DF_test_miss, n, K)
  
  try(list_DF <- append(list_DF, list(DF)))
  try(list_DF_miss <- append(list_DF_miss, list(DF_miss)))
  try(list_DF_test <- append(list_DF_test, list(DF_test)))
  try(list_DF_test_miss <- append(list_DF_test_miss, list(DF_test_miss)))
  
  try(if(count==100){
    str100 <- "dataframes_MCAR_50.RData"
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
  
  ###OBS!!!!!
  ## DF_miss contiene el data frame con missing data!!
  
  ######################################
  #Fitting multivariate MM models
  # - We use snat_mvmer() from rstanarm
  # - We fit the true model plus two
  #.  (or 3) more models
  ######################################
  
  #checking time of computation using data frame with missings 
  # when n=300, k=10, we have 21mins, max rhat 1.67
  #      n=250, k=10: 17mins, max rhat 1.51
  #      n=200, k=10:  8.6 mins, max rhat 1.09
  try(true_model_miss <- stan_mvmer(
    formula = list(
      y1 ~ sex * time + (time | id),
      y3 ~ sqrt(time) + (sqrt(time) | id),
      y5 ~ treatment + time + (time | id)),
    data = list_DF_miss[[count]],
    family = list(gaussian, gaussian, gaussian),
    chains = 3, cores = 8, seed = 12345, iter = 1750))
  
  try(rhat_true_model_miss <- rhat(true_model_miss))
  
  #OJO!!! para que devuelva las predictions a pesar de estar pasandole NAs 
  #hay que pasarle el data frame sin los NAs
  
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
                           uniExpansions = list(list(type="uFPCA", pve = 0.65),
                                                list(type ="uFPCA", pve=0.65),
                                                list(type="uFPCA", pve=0.65)),
                           fit = TRUE))
  
  try(pred_mfpca1_test <- predict(mfpca1, scores = mfpca1_test$scores))
  
  try(pred_mfpca1_y1_test <- NULL)
  try(pred_mfpca1_y1_test <- pred_mfpca1_test[[1]]@X)
  try(pred_mfpca1_y2_test <- NULL)
  try(pred_mfpca1_y2_test <- pred_mfpca1_test[[2]]@X)
  try(pred_mfpca1_y3_test <- NULL)
  try(pred_mfpca1_y3_test <- pred_mfpca1_test[[3]]@X)
  
  #plot(pred_mfpca1_test)
  
  #AQUI YA TENGO TODA LA INFO QUE QUIERO! las predicciones de esto último se 
  #hacen para cada uno de los puntos de la grid!! También para aquellos donde
  #hay NA! Por tanto para aquellos puntos estoy haciendo una predicción!
  #Y es esa predicción la que me tengo que quedar y entrar a comparar con 
  #la predicción que me puedan devolver los MM
  
  #Cuando tenemos el predict(MFPCAfit), el objeto que pred_mfpca1[[1]]@X es
  # una matriz n X max_length +1, ie cada fila es un individuo y cada columna
  # es uno de los puntos de la fixed time grid. Y esto para el outcome 
  # longitudinal 1 
  # Lo que tenemos aquí es una predicción para cada punto de la fixed grid
  # del outcome longitudinal 1.
  
  
  
  #####################################################
  ## Comparing predictions by means of standardized
  ## RMSE
  #####################################################
  
  
  
  
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
    strr <- "results_1_MCAR50_6feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr)
  })
  try(if(count==2){
    strr2 <- "results_2_MCAR50_6feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr2)
  })
  try(if(count==5){
    str1 <- "results_5_MCAR50_6feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str1)
  })
  try(if(count==10){
    str10 <- "results_10_MCAR50_6feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str10)
  })
  try(if(count==20){
    str20 <- "results_20_MCAR50_6feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str20)
  })
  try(if(count==30){
    str30 <- "results_30_MCAR50_6feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str30)
  })
  try(if(count==40){
    str40 <- "results_40_MCAR50_6feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss,file=str40)
  })
  try(if(count==50){
    str50 <- "results_50_MCAR50_6feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str50)
  })
  try(if(count==75){
    str75 <- "results_75_MCAR50_6feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str75)
  })
  try(if(count==100){
    str100 <- "results_100_MCAR50_6feb2025.RData"
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


###############################################################################
###############################################################################
##MAR 30% of dropout
###############################################################################
###############################################################################
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
## I generate the data frames: MAR 30% of dropout
############################################################
############################################################
for(count in 1:num_datasets){
  n <- 200 # number of subjects
  n_test <- 200
  K <- 10 # number of measurements per subject
  t_max <- 10 # maximum follow-up time
  min_sep <- 0.5
  
  # we construct a data frame with the design:
  # everyone has a baseline measurement, and then measurements at random 
  # follow-up times up to t_max 
  # !! for Machine Learning purposes: two consecutive measurements will have
  #more than 0.5 of
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
  
  #X2 <- model.matrix(~ treatment * time, data = DF )
  #Z2 <- model.matrix(~ time, data = DF)
  
  #X3 <- model.matrix(~ time, data = DF )
  #Z3 <- model.matrix(~ time, data = DF)
  
  X3 <- model.matrix(~ sqrt(time), data = DF )
  Z3 <- model.matrix(~ sqrt(time), data = DF)
  
  #X4 <- model.matrix(~ sex + time, data = DF)
  #Z4 <- model.matrix(~ time, data = DF)
  
  #X4 <- model.matrix(~ sex + time , data = DF)
  #Z4 <- model.matrix(~ time, data = DF)
  
  #X5 <- model.matrix(~ time, data = DF)
  #Z5 <- model.matrix(~ time, data = DF)
  
  X5 <- model.matrix(~ treatment + time, data = DF)
  Z5 <- model.matrix(~ time, data = DF)
  
  ####for test data set
  X1_test <- model.matrix(~ sex * time, data = DF_test)
  Z1_test <- model.matrix(~ time, data = DF_test)
  
  #X2_test <- model.matrix(~ treatment * time, data = DF_test )
  #Z2_test <- model.matrix(~ time, data = DF_test)
  
  #X3_test <- model.matrix(~ time, data = DF_test )
  #Z3_test <- model.matrix(~ time, data = DF_test)
  
  X3_test <- model.matrix(~ sqrt(time), data = DF_test )
  Z3_test <- model.matrix(~ sqrt(time), data = DF_test)
  
  #X4_test <- model.matrix(~ sex + time, data = DF_test)
  #Z4_test <- model.matrix(~ time, data = DF_test)
  
  #X4_test <- model.matrix(~ sex + time, data = DF_test)
  #Z4_test <- model.matrix(~ time, data = DF_test)
  
  #X5_test <- model.matrix(~ time, data = DF_test)
  #Z5_test <- model.matrix(~ time, data = DF_test)
  
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
  
  #Simulate second longitudinal outcome
  ###############################################
  #betas2 <- c(-1.8, -0.06, 0.5, 0.06) # fixed effects coefficients
  #betas2 <- c(10, -0.02, 2.222, 0.07) # fixed effects coefficients
  #sigma2 <- 0.25 # errors sd
  
  # we simulate random effects
  #b2 <- b[, c(3,4)]
  # linear predictor
  #eta_y2 <- as.vector(X2 %*% betas2 + rowSums(Z2 * b2[DF$id, ]))
  # we simulate normal longitudinal data
  #DF$y2 <- rnorm(n * K, mean = eta_y2, sd = sigma2)
  # we assume that values below 0 are not observed, and set equal to 0
  #DF$ind2 <- as.numeric(DF$y2 < -4)
  #DF$y2 <- pmax(DF$y2, -4)
  
  ####Test data
  # we simulate random effects
  #b2_test <- b_test[, c(3,4)]
  # linear predictor
  #eta_y2_test <- as.vector(X2_test %*% betas2 + rowSums(Z2_test * b2[DF_test$id, ]))
  # we simulate normal longitudinal data
  #DF_test$y2 <- rnorm(n_test * K, mean = eta_y2_test, sd = sigma2)
  # we assume that values below 0 are not observed, and set equal to 0
  #DF_test$ind2 <- as.numeric(DF_test$y2 < -4)
  #DF_test$y2 <- pmax(DF_test$y2, -4)
  
  #Simulate third longitudinal outcome
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
  
  #Simulate forth longitudinal outcome
  ###############################################
  #betas4 <- c(0.01, 0.5, -0.31416) # fixed effects coefficients
  #betas4 <- c(80, -4.25, -3.789) # fixed effects coefficients
  
  # we simulate random effects
  #b4 <- b[, c(7,8)]
  # linear predictor
  #eta_y4 <- as.vector(X4 %*% betas4 + rowSums(Z4 * b4[DF$id, ]))
  # mean of the binomial distribution
  #mu_y4 <- plogis(eta_y4)
  # we simulate binomial longitudinal data
  #DF$y4 <- rbinom(n * K, size = 1, prob = mu_y4)
  #DF$y4 <- rnorm(n_test * K, mean = eta_y4, sd = sigma3)
  
  #####Test data
  # we simulate random effects
  #b4_test <- b_test[, c(7,8)]
  # linear predictor
  #eta_y4_test <- as.vector(X4_test %*% betas4 + rowSums(Z4_test * b4_test[DF_test$id, ]))
  # mean of the binomial distribution
  #mu_y4_test <- plogis(eta_y4_test)
  # we simulate binomial longitudinal data
  #DF_test$y4 <- rbinom(n_test * K, size = 1, prob = mu_y4_test)
  #DF_test$y4 <- rnorm(n_test * K, mean = eta_y4_test, sd = sigma3)
  
  #Simulate fifth longitudinal outcome
  ###############################################
  betas5 <- c(1, 0.155, 0.12345) # fixed effects coefficients
  betas5 <- c(70, -5.44, -3.761289) # fixed effects coefficients
  
  # we simulate random effects
  b5 <- b[, c(5,6)]
  # linear predictor
  eta_y5 <- as.vector(X5 %*% betas5 + rowSums(Z5 * b5[DF$id, ]))
  # mean of the binomial distribution
  #mu_y5 <- plogis(eta_y5)
  # we simulate binomial longitudinal data
  #DF$y5 <- rbinom(n * K, size = 1, prob = mu_y5)
  DF$y5 <- rnorm(n_test * K, mean = eta_y5, sd = sigma3)
  
  ####Test data
  # we simulate random effects
  b5_test <- b_test[, c(5,6)]
  # linear predictor
  eta_y5_test <- as.vector(X5_test %*% betas5 + rowSums(Z5_test * b5_test[DF_test$id, ]))
  # mean of the binomial distribution
  #mu_y5_test <- plogis(eta_y5_test)
  # we simulate binomial longitudinal data
  #DF_test$y5 <- rbinom(n_test * K, size = 1, prob = mu_y5_test)
  
  DF_test$y5 <- rnorm(n_test * K, mean = eta_y5_test, sd = sigma3)
  
  #If using admin censoring:
  #AdminCens <- 10
  
  # we keep the longitudinal measurements before the event times: TRAIN
  #DF <- DF[DF$time <= AdminCens, ]
  
  # we keep the longitudinal measurements before the event times: TEST
  #DF_test <- DF_test[DF_test$time <= AdminCens, ]
  
  #Data set with just one obs per individual
  #DF.id <- DF[!duplicated(DF$id),]
  #DF_test.id <- DF_test[!duplicated(DF_test$id),]
  
  
  ### MAR dropout with 30% intensity
  nu1 <- 15.5
  nu2 <- 7.5
  nu3 <- 42
  
  DF_miss <- DF %>%
    group_by(id) %>%                       # Group by each individual
    mutate(
      # Check dropout conditions on the previous time point
      dropout_point = which.max(lag(y1) > nu1 | lag(y3) > nu2 | lag(y5) < nu3),
      
      # Set NA to y1, y2, y3 after the dropout point (ignoring the first visit)
      y1 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y1),
      y3 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y3),
      y5 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y5)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point) 
 
  DF_test_miss <- DF_test %>%
    group_by(id) %>%                       # Group by each individual
    mutate(
      # Check dropout conditions on the previous time point
      dropout_point = which.max(lag(y1) > nu1 | lag(y3) > nu2 | lag(y5) < nu3),
      
      # Set NA to y1, y2, y3 after the dropout point (ignoring the first visit)
      y1 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y1),
      y3 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y3),
      y5 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y5)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point) 

  
  missing_ratio_train[count] <- counting_missing_data(DF_miss, n, K)
  missing_ratio_test[count] <- counting_missing_data(DF_test_miss, n, K)
  
  try(list_DF <- append(list_DF, list(DF)))
  try(list_DF_miss <- append(list_DF_miss, list(DF_miss)))
  try(list_DF_test <- append(list_DF_test, list(DF_test)))
  try(list_DF_test_miss <- append(list_DF_test_miss, list(DF_test_miss)))
  
  try(if(count==100){
    str100 <- "dataframes_MAR_30.RData"
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
  
  ###OBS!!!!!
  ## DF_miss contiene el data frame con missing data!!
  
  ######################################
  #Fitting multivariate MM models
  # - We use snat_mvmer() from rstanarm
  # - We fit the true model plus two
  #.  (or 3) more models
  ######################################
  
  #checking time of computation using data frame with missings 
  # when n=300, k=10, we have 21mins, max rhat 1.67
  #      n=250, k=10: 17mins, max rhat 1.51
  #      n=200, k=10:  8.6 mins, max rhat 1.09
  try(true_model_miss <- stan_mvmer(
    formula = list(
      y1 ~ sex * time + (time | id),
      y3 ~ sqrt(time) + (sqrt(time) | id),
      y5 ~ treatment + time + (time | id)),
    data = list_DF_miss[[count]],
    family = list(gaussian, gaussian, gaussian),
    chains = 3, cores = 8, seed = 12345, iter = 1750))
  
  try(rhat_true_model_miss <- rhat(true_model_miss))
  
  #OJO!!! para que devuelva las predictions a pesar de estar pasandole NAs 
  #hay que pasarle el data frame sin los NAs
  
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
  
  #We fit with NO weights and pve=0.8
  try(mfpca1 <- MFPCA(m1, M = 3, 
                      uniExpansions = list(list(type="uFPCA", pve = 0.8),
                                           list(type ="uFPCA", pve=0.8),
                                           list(type="uFPCA", pve=0.8)),
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
  
  #plot(pred_mfpca1_test)
  
  #AQUI YA TENGO TODA LA INFO QUE QUIERO! las predicciones de esto último se 
  #hacen para cada uno de los puntos de la grid!! También para aquellos donde
  #hay NA! Por tanto para aquellos puntos estoy haciendo una predicción!
  #Y es esa predicción la que me tengo que quedar y entrar a comparar con 
  #la predicción que me puedan devolver los MM
  
  #Cuando tenemos el predict(MFPCAfit), el objeto que pred_mfpca1[[1]]@X es
  # una matriz n X max_length +1, ie cada fila es un individuo y cada columna
  # es uno de los puntos de la fixed time grid. Y esto para el outcome 
  # longitudinal 1 
  # Lo que tenemos aquí es una predicción para cada punto de la fixed grid
  # del outcome longitudinal 1.
  
  
  
  #####################################################
  ## Comparing predictions by means of standardized
  ## RMSE
  #####################################################
  
  
  
  
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
    strr <- "results_1_MAR30_7feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr)
  })
  try(if(count==2){
    strr2 <- "results_2_MAR30_7feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr2)
  })
  try(if(count==5){
    str1 <- "results_5_MAR30_7feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str1)
  })
  try(if(count==10){
    str10 <- "results_10_MAR30_7feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str10)
  })
  try(if(count==20){
    str20 <- "results_20_MAR30_7feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str20)
  })
  try(if(count==30){
    str30 <- "results_30_MAR30_7feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str30)
  })
  try(if(count==40){
    str40 <- "results_40_MAR30_7feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss,file=str40)
  })
  try(if(count==50){
    str50 <- "results_50_MAR30_7feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str50)
  })
  try(if(count==75){
    str75 <- "results_75_MAR30_7feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_mo111del_miss, file=str75)
  })
  try(if(count==100){
    str100 <- "results_100_MAR30_7feb2025.RData"
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


###############################################################################
###############################################################################
##MAR 50% of dropout
###############################################################################
###############################################################################
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
## I generate the data frames: MAR 50% of dropout
############################################################
############################################################
for(count in 1:num_datasets){
  n <- 200 # number of subjects
  n_test <- 200
  K <- 10 # number of measurements per subject
  t_max <- 10 # maximum follow-up time
  min_sep <- 0.5
  
  # we construct a data frame with the design:
  # everyone has a baseline measurement, and then measurements at random 
  # follow-up times up to t_max 
  # !! for Machine Learning purposes: two consecutive measurements will have
  #more than 0.5 of
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
  
  #X2 <- model.matrix(~ treatment * time, data = DF )
  #Z2 <- model.matrix(~ time, data = DF)
  
  #X3 <- model.matrix(~ time, data = DF )
  #Z3 <- model.matrix(~ time, data = DF)
  
  X3 <- model.matrix(~ sqrt(time), data = DF )
  Z3 <- model.matrix(~ sqrt(time), data = DF)
  
  #X4 <- model.matrix(~ sex + time, data = DF)
  #Z4 <- model.matrix(~ time, data = DF)
  
  #X4 <- model.matrix(~ sex + time , data = DF)
  #Z4 <- model.matrix(~ time, data = DF)
  
  #X5 <- model.matrix(~ time, data = DF)
  #Z5 <- model.matrix(~ time, data = DF)
  
  X5 <- model.matrix(~ treatment + time, data = DF)
  Z5 <- model.matrix(~ time, data = DF)
  
  ####for test data set
  X1_test <- model.matrix(~ sex * time, data = DF_test)
  Z1_test <- model.matrix(~ time, data = DF_test)
  
  #X2_test <- model.matrix(~ treatment * time, data = DF_test )
  #Z2_test <- model.matrix(~ time, data = DF_test)
  
  #X3_test <- model.matrix(~ time, data = DF_test )
  #Z3_test <- model.matrix(~ time, data = DF_test)
  
  X3_test <- model.matrix(~ sqrt(time), data = DF_test )
  Z3_test <- model.matrix(~ sqrt(time), data = DF_test)
  
  #X4_test <- model.matrix(~ sex + time, data = DF_test)
  #Z4_test <- model.matrix(~ time, data = DF_test)
  
  #X4_test <- model.matrix(~ sex + time, data = DF_test)
  #Z4_test <- model.matrix(~ time, data = DF_test)
  
  #X5_test <- model.matrix(~ time, data = DF_test)
  #Z5_test <- model.matrix(~ time, data = DF_test)
  
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
  
  #Simulate second longitudinal outcome
  ###############################################
  #betas2 <- c(-1.8, -0.06, 0.5, 0.06) # fixed effects coefficients
  #betas2 <- c(10, -0.02, 2.222, 0.07) # fixed effects coefficients
  #sigma2 <- 0.25 # errors sd
  
  # we simulate random effects
  #b2 <- b[, c(3,4)]
  # linear predictor
  #eta_y2 <- as.vector(X2 %*% betas2 + rowSums(Z2 * b2[DF$id, ]))
  # we simulate normal longitudinal data
  #DF$y2 <- rnorm(n * K, mean = eta_y2, sd = sigma2)
  # we assume that values below 0 are not observed, and set equal to 0
  #DF$ind2 <- as.numeric(DF$y2 < -4)
  #DF$y2 <- pmax(DF$y2, -4)
  
  ####Test data
  # we simulate random effects
  #b2_test <- b_test[, c(3,4)]
  # linear predictor
  #eta_y2_test <- as.vector(X2_test %*% betas2 + rowSums(Z2_test * b2[DF_test$id, ]))
  # we simulate normal longitudinal data
  #DF_test$y2 <- rnorm(n_test * K, mean = eta_y2_test, sd = sigma2)
  # we assume that values below 0 are not observed, and set equal to 0
  #DF_test$ind2 <- as.numeric(DF_test$y2 < -4)
  #DF_test$y2 <- pmax(DF_test$y2, -4)
  
  #Simulate third longitudinal outcome
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
  
  #Simulate forth longitudinal outcome
  ###############################################
  #betas4 <- c(0.01, 0.5, -0.31416) # fixed effects coefficients
  #betas4 <- c(80, -4.25, -3.789) # fixed effects coefficients
  
  # we simulate random effects
  #b4 <- b[, c(7,8)]
  # linear predictor
  #eta_y4 <- as.vector(X4 %*% betas4 + rowSums(Z4 * b4[DF$id, ]))
  # mean of the binomial distribution
  #mu_y4 <- plogis(eta_y4)
  # we simulate binomial longitudinal data
  #DF$y4 <- rbinom(n * K, size = 1, prob = mu_y4)
  #DF$y4 <- rnorm(n_test * K, mean = eta_y4, sd = sigma3)
  
  #####Test data
  # we simulate random effects
  #b4_test <- b_test[, c(7,8)]
  # linear predictor
  #eta_y4_test <- as.vector(X4_test %*% betas4 + rowSums(Z4_test * b4_test[DF_test$id, ]))
  # mean of the binomial distribution
  #mu_y4_test <- plogis(eta_y4_test)
  # we simulate binomial longitudinal data
  #DF_test$y4 <- rbinom(n_test * K, size = 1, prob = mu_y4_test)
  #DF_test$y4 <- rnorm(n_test * K, mean = eta_y4_test, sd = sigma3)
  
  #Simulate fifth longitudinal outcome
  ###############################################
  betas5 <- c(1, 0.155, 0.12345) # fixed effects coefficients
  betas5 <- c(70, -5.44, -3.761289) # fixed effects coefficients
  
  # we simulate random effects
  b5 <- b[, c(5,6)]
  # linear predictor
  eta_y5 <- as.vector(X5 %*% betas5 + rowSums(Z5 * b5[DF$id, ]))
  # mean of the binomial distribution
  #mu_y5 <- plogis(eta_y5)
  # we simulate binomial longitudinal data
  #DF$y5 <- rbinom(n * K, size = 1, prob = mu_y5)
  DF$y5 <- rnorm(n_test * K, mean = eta_y5, sd = sigma3)
  
  ####Test data
  # we simulate random effects
  b5_test <- b_test[, c(5,6)]
  # linear predictor
  eta_y5_test <- as.vector(X5_test %*% betas5 + rowSums(Z5_test * b5_test[DF_test$id, ]))
  # mean of the binomial distribution
  #mu_y5_test <- plogis(eta_y5_test)
  # we simulate binomial longitudinal data
  #DF_test$y5 <- rbinom(n_test * K, size = 1, prob = mu_y5_test)
  
  DF_test$y5 <- rnorm(n_test * K, mean = eta_y5_test, sd = sigma3)
  
  #If using admin censoring:
  #AdminCens <- 10
  
  # we keep the longitudinal measurements before the event times: TRAIN
  #DF <- DF[DF$time <= AdminCens, ]
  
  # we keep the longitudinal measurements before the event times: TEST
  #DF_test <- DF_test[DF_test$time <= AdminCens, ]
  
  #Data set with just one obs per individual
  #DF.id <- DF[!duplicated(DF$id),]
  #DF_test.id <- DF_test[!duplicated(DF_test$id),]
  
  
  ### MAR dropout with 30% intensity
  nu1 <- 8
  nu2 <- 6.5
  nu3 <- 35
  
  DF_miss <- DF %>%
    group_by(id) %>%                       # Group by each individual
    mutate(
      # Check dropout conditions on the previous time point
      dropout_point = which.max(lag(y1) > nu1 | lag(y3) > nu2 | lag(y5) < nu3),
      
      # Set NA to y1, y2, y3 after the dropout point (ignoring the first visit)
      y1 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y1),
      y3 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y3),
      y5 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y5)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point) 
  
  DF_test_miss <- DF_test %>%
    group_by(id) %>%                       # Group by each individual
    mutate(
      # Check dropout conditions on the previous time point
      dropout_point = which.max(lag(y1) > nu1 | lag(y3) > nu2 | lag(y5) < nu3),
      
      # Set NA to y1, y2, y3 after the dropout point (ignoring the first visit)
      y1 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y1),
      y3 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y3),
      y5 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y5)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point) 
  
  
  missing_ratio_train[count] <- counting_missing_data(DF_miss, n, K)
  missing_ratio_test[count] <- counting_missing_data(DF_test_miss, n, K)
  
  try(list_DF <- append(list_DF, list(DF)))
  try(list_DF_miss <- append(list_DF_miss, list(DF_miss)))
  try(list_DF_test <- append(list_DF_test, list(DF_test)))
  try(list_DF_test_miss <- append(list_DF_test_miss, list(DF_test_miss)))
  
  try(if(count==100){
    str100 <- "dataframes_MAR_50.RData"
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
  
  ###OBS!!!!!
  ## DF_miss contiene el data frame con missing data!!
  
  ######################################
  #Fitting multivariate MM models
  # - We use snat_mvmer() from rstanarm
  # - We fit the true model plus two
  #.  (or 3) more models
  ######################################
  
  #checking time of computation using data frame with missings 
  # when n=300, k=10, we have 21mins, max rhat 1.67
  #      n=250, k=10: 17mins, max rhat 1.51
  #      n=200, k=10:  8.6 mins, max rhat 1.09
  try(true_model_miss <- stan_mvmer(
    formula = list(
      y1 ~ sex * time + (time | id),
      y3 ~ sqrt(time) + (sqrt(time) | id),
      y5 ~ treatment + time + (time | id)),
    data = list_DF_miss[[count]],
    family = list(gaussian, gaussian, gaussian),
    chains = 3, cores = 8, seed = 12345, iter = 1750))
  
  try(rhat_true_model_miss <- rhat(true_model_miss))
  
  #OJO!!! para que devuelva las predictions a pesar de estar pasandole NAs 
  #hay que pasarle el data frame sin los NAs
  
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
                           uniExpansions = list(list(type="uFPCA", pve = 0.65),
                                                list(type ="uFPCA", pve=0.65),
                                                list(type="uFPCA", pve=0.65)),
                           fit = TRUE))
  
  try(pred_mfpca1_test <- predict(mfpca1, scores = mfpca1_test$scores))
  
  try(pred_mfpca1_y1_test <- NULL)
  try(pred_mfpca1_y1_test <- pred_mfpca1_test[[1]]@X)
  try(pred_mfpca1_y2_test <- NULL)
  try(pred_mfpca1_y2_test <- pred_mfpca1_test[[2]]@X)
  try(pred_mfpca1_y3_test <- NULL)
  try(pred_mfpca1_y3_test <- pred_mfpca1_test[[3]]@X)
  
  #plot(pred_mfpca1_test)
  
  #AQUI YA TENGO TODA LA INFO QUE QUIERO! las predicciones de esto último se 
  #hacen para cada uno de los puntos de la grid!! También para aquellos donde
  #hay NA! Por tanto para aquellos puntos estoy haciendo una predicción!
  #Y es esa predicción la que me tengo que quedar y entrar a comparar con 
  #la predicción que me puedan devolver los MM
  
  #Cuando tenemos el predict(MFPCAfit), el objeto que pred_mfpca1[[1]]@X es
  # una matriz n X max_length +1, ie cada fila es un individuo y cada columna
  # es uno de los puntos de la fixed time grid. Y esto para el outcome 
  # longitudinal 1 
  # Lo que tenemos aquí es una predicción para cada punto de la fixed grid
  # del outcome longitudinal 1.
  
  
  
  #####################################################
  ## Comparing predictions by means of standardized
  ## RMSE
  #####################################################
  
  
  
  
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
    strr <- "results_1_MAR50_10feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr)
  })
  try(if(count==2){
    strr2 <- "results_2_MAR50_10feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr2)
  })
  try(if(count==5){
    str1 <- "results_5_MAR50_10feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str1)
  })
  try(if(count==10){
    str10 <- "results_10_MAR50_10feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str10)
  })
  try(if(count==20){
    str20 <- "results_20_MAR50_10feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str20)
  })
  try(if(count==30){
    str30 <- "results_30_MAR50_10feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str30)
  })
  try(if(count==40){
    str40 <- "results_40_MAR50_10feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss,file=str40)
  })
  try(if(count==50){
    str50 <- "results_50_MAR50_10feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str50)
  })
  try(if(count==75){
    str75 <- "results_75_MAR50_10feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_mo111del_miss, file=str75)
  })
  try(if(count==100){
    str100 <- "results_100_MAR50_10feb2025.RData"
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


###############################################################################
###############################################################################
## Threshold MAR 30% of dropout
###############################################################################
###############################################################################
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
  t_max <- 10 # maximum follow-up time
  min_sep <- 0.5
  
  # we construct a data frame with the design:
  # everyone has a baseline measurement, and then measurements at random
  # follow-up times up to t_max
  # !! for Machine Learning purposes: two consecutive measurements will have
  #more than 0.5 of
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
  
  #X2 <- model.matrix(~ treatment * time, data = DF )
  #Z2 <- model.matrix(~ time, data = DF)
  
  #X3 <- model.matrix(~ time, data = DF )
  #Z3 <- model.matrix(~ time, data = DF)
  
  X3 <- model.matrix(~ sqrt(time), data = DF )
  Z3 <- model.matrix(~ sqrt(time), data = DF)
  
  #X4 <- model.matrix(~ sex + time, data = DF)
  #Z4 <- model.matrix(~ time, data = DF)
  
  #X4 <- model.matrix(~ sex + time , data = DF)
  #Z4 <- model.matrix(~ time, data = DF)
  
  #X5 <- model.matrix(~ time, data = DF)
  #Z5 <- model.matrix(~ time, data = DF)
  
  X5 <- model.matrix(~ treatment + time, data = DF)
  Z5 <- model.matrix(~ time, data = DF)
  
  ####for test data set
  X1_test <- model.matrix(~ sex * time, data = DF_test)
  Z1_test <- model.matrix(~ time, data = DF_test)
  
  #X2_test <- model.matrix(~ treatment * time, data = DF_test )
  #Z2_test <- model.matrix(~ time, data = DF_test)
  
  #X3_test <- model.matrix(~ time, data = DF_test )
  #Z3_test <- model.matrix(~ time, data = DF_test)
  
  X3_test <- model.matrix(~ sqrt(time), data = DF_test )
  Z3_test <- model.matrix(~ sqrt(time), data = DF_test)
  
  #X4_test <- model.matrix(~ sex + time, data = DF_test)
  #Z4_test <- model.matrix(~ time, data = DF_test)
  
  #X4_test <- model.matrix(~ sex + time, data = DF_test)
  #Z4_test <- model.matrix(~ time, data = DF_test)
  
  #X5_test <- model.matrix(~ time, data = DF_test)
  #Z5_test <- model.matrix(~ time, data = DF_test)
  
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
  
  #Simulate second longitudinal outcome
  ###############################################
  #betas2 <- c(-1.8, -0.06, 0.5, 0.06) # fixed effects coefficients
  #betas2 <- c(10, -0.02, 2.222, 0.07) # fixed effects coefficients
  #sigma2 <- 0.25 # errors sd
  
  # we simulate random effects
  #b2 <- b[, c(3,4)]
  # linear predictor
  #eta_y2 <- as.vector(X2 %*% betas2 + rowSums(Z2 * b2[DF$id, ]))
  # we simulate normal longitudinal data
  #DF$y2 <- rnorm(n * K, mean = eta_y2, sd = sigma2)
  # we assume that values below 0 are not observed, and set equal to 0
  #DF$ind2 <- as.numeric(DF$y2 < -4)
  #DF$y2 <- pmax(DF$y2, -4)
  
  ####Test data
  # we simulate random effects
  #b2_test <- b_test[, c(3,4)]
  # linear predictor
  #eta_y2_test <- as.vector(X2_test %*% betas2 + rowSums(Z2_test * b2[DF_test$id, ]))
  # we simulate normal longitudinal data
  #DF_test$y2 <- rnorm(n_test * K, mean = eta_y2_test, sd = sigma2)
  # we assume that values below 0 are not observed, and set equal to 0
  #DF_test$ind2 <- as.numeric(DF_test$y2 < -4)
  #DF_test$y2 <- pmax(DF_test$y2, -4)
  
  #Simulate third longitudinal outcome
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
  
  #Simulate forth longitudinal outcome
  ###############################################
  #betas4 <- c(0.01, 0.5, -0.31416) # fixed effects coefficients
  #betas4 <- c(80, -4.25, -3.789) # fixed effects coefficients
  
  # we simulate random effects
  #b4 <- b[, c(7,8)]
  # linear predictor
  #eta_y4 <- as.vector(X4 %*% betas4 + rowSums(Z4 * b4[DF$id, ]))
  # mean of the binomial distribution
  #mu_y4 <- plogis(eta_y4)
  # we simulate binomial longitudinal data
  #DF$y4 <- rbinom(n * K, size = 1, prob = mu_y4)
  #DF$y4 <- rnorm(n_test * K, mean = eta_y4, sd = sigma3)
  
  #####Test data
  # we simulate random effects
  #b4_test <- b_test[, c(7,8)]
  # linear predictor
  #eta_y4_test <- as.vector(X4_test %*% betas4 + rowSums(Z4_test * b4_test[DF_test$id, ]))
  # mean of the binomial distribution
  #mu_y4_test <- plogis(eta_y4_test)
  # we simulate binomial longitudinal data
  #DF_test$y4 <- rbinom(n_test * K, size = 1, prob = mu_y4_test)
  #DF_test$y4 <- rnorm(n_test * K, mean = eta_y4_test, sd = sigma3)
  
  #Simulate fifth longitudinal outcome
  ###############################################
  betas5 <- c(1, 0.155, 0.12345) # fixed effects coefficients
  betas5 <- c(70, -5.44, -3.761289) # fixed effects coefficients
  
  # we simulate random effects
  b5 <- b[, c(5,6)]
  # linear predictor
  eta_y5 <- as.vector(X5 %*% betas5 + rowSums(Z5 * b5[DF$id, ]))
  # mean of the binomial distribution
  #mu_y5 <- plogis(eta_y5)
  # we simulate binomial longitudinal data
  #DF$y5 <- rbinom(n * K, size = 1, prob = mu_y5)
  DF$y5 <- rnorm(n_test * K, mean = eta_y5, sd = sigma3)
  
  ####Test data
  # we simulate random effects
  b5_test <- b_test[, c(5,6)]
  # linear predictor
  eta_y5_test <- as.vector(X5_test %*% betas5 + rowSums(Z5_test * b5_test[DF_test$id, ]))
  # mean of the binomial distribution
  #mu_y5_test <- plogis(eta_y5_test)
  # we simulate binomial longitudinal data
  #DF_test$y5 <- rbinom(n_test * K, size = 1, prob = mu_y5_test)
  
  DF_test$y5 <- rnorm(n_test * K, mean = eta_y5_test, sd = sigma3)
  
  #If using admin censoring:
  #AdminCens <- 10
  
  # we keep the longitudinal measurements before the event times: TRAIN
  #DF <- DF[DF$time <= AdminCens, ]
  
  # we keep the longitudinal measurements before the event times: TEST
  #DF_test <- DF_test[DF_test$time <= AdminCens, ]
  
  #Data set with just one obs per individual
  #DF.id <- DF[!duplicated(DF$id),]
  #DF_test.id <- DF_test[!duplicated(DF_test$id),]
  
  
  ### threshold MAR dropout with 30% intensity
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
  
  param1 <- -3.75
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
    str100 <- "dataframes_MAR_threshold_30.RData"
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
  
  ###OBS!!!!!
  ## DF_miss contiene el data frame con missing data!!
  
  ######################################
  #Fitting multivariate MM models
  # - We use snat_mvmer() from rstanarm
  # - We fit the true model plus two
  #.  (or 3) more models
  ######################################
  
  #checking time of computation using data frame with missings 
  # when n=300, k=10, we have 21mins, max rhat 1.67
  #      n=250, k=10: 17mins, max rhat 1.51
  #      n=200, k=10:  8.6 mins, max rhat 1.09
  try(true_model_miss <- stan_mvmer(
    formula = list(
      y1 ~ sex * time + (time | id),
      y3 ~ sqrt(time) + (sqrt(time) | id),
      y5 ~ treatment + time + (time | id)),
    data = list_DF_miss[[count]],
    family = list(gaussian, gaussian, gaussian),
    chains = 3, cores = 8, seed = 12345, iter = 1750))
  
  try(rhat_true_model_miss <- rhat(true_model_miss))
  
  #OJO!!! para que devuelva las predictions a pesar de estar pasandole NAs 
  #hay que pasarle el data frame sin los NAs
  
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
  
  #We fit with NO weights and pve=0.8
  try(mfpca1 <- MFPCA(m1, M = 3, 
                      uniExpansions = list(list(type="uFPCA", pve = 0.8),
                                           list(type ="uFPCA", pve=0.8),
                                           list(type="uFPCA", pve=0.8)),
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
  
  #plot(pred_mfpca1_test)
  
  #AQUI YA TENGO TODA LA INFO QUE QUIERO! las predicciones de esto último se 
  #hacen para cada uno de los puntos de la grid!! También para aquellos donde
  #hay NA! Por tanto para aquellos puntos estoy haciendo una predicción!
  #Y es esa predicción la que me tengo que quedar y entrar a comparar con 
  #la predicción que me puedan devolver los MM
  
  #Cuando tenemos el predict(MFPCAfit), el objeto que pred_mfpca1[[1]]@X es
  # una matriz n X max_length +1, ie cada fila es un individuo y cada columna
  # es uno de los puntos de la fixed time grid. Y esto para el outcome 
  # longitudinal 1 
  # Lo que tenemos aquí es una predicción para cada punto de la fixed grid
  # del outcome longitudinal 1.
  
  
  
  #####################################################
  ## Comparing predictions by means of standardized
  ## RMSE
  #####################################################
  
  
  
  
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
    strr <- "results_1_MAR30_thr_11feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr)
  })
  try(if(count==2){
    strr2 <- "results_2_MAR30_thr_11feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr2)
  })
  try(if(count==5){
    str1 <- "results_5_MAR30_thr_11feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str1)
  })
  try(if(count==10){
    str10 <- "results_10_MAR30_thr_11feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str10)
  })
  try(if(count==20){
    str20 <- "results_20_MAR30_thr_11feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str20)
  })
  try(if(count==30){
    str30 <- "results_30_MAR30_thr_11feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str30)
  })
  try(if(count==40){
    str40 <- "results_40_MAR30_thr_11feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss,file=str40)
  })
  try(if(count==50){
    str50 <- "results_50_MAR30_thr_11feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str50)
  })
  try(if(count==75){
    str75 <- "results_75_MAR30_thr_11feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_mo111del_miss, file=str75)
  })
  try(if(count==100){
    str100 <- "results_100_MAR30_thr_11feb2025.RData"
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


###############################################################################
###############################################################################
## Threshold MAR 50% of dropout
###############################################################################
###############################################################################
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
  t_max <- 10 # maximum follow-up time
  min_sep <- 0.5
  
  # we construct a data frame with the design:
  # everyone has a baseline measurement, and then measurements at random
  # follow-up times up to t_max
  # !! for Machine Learning purposes: two consecutive measurements will have
  #more than 0.5 of
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
  
  #X2 <- model.matrix(~ treatment * time, data = DF )
  #Z2 <- model.matrix(~ time, data = DF)
  
  #X3 <- model.matrix(~ time, data = DF )
  #Z3 <- model.matrix(~ time, data = DF)
  
  X3 <- model.matrix(~ sqrt(time), data = DF )
  Z3 <- model.matrix(~ sqrt(time), data = DF)
  
  #X4 <- model.matrix(~ sex + time, data = DF)
  #Z4 <- model.matrix(~ time, data = DF)
  
  #X4 <- model.matrix(~ sex + time , data = DF)
  #Z4 <- model.matrix(~ time, data = DF)
  
  #X5 <- model.matrix(~ time, data = DF)
  #Z5 <- model.matrix(~ time, data = DF)
  
  X5 <- model.matrix(~ treatment + time, data = DF)
  Z5 <- model.matrix(~ time, data = DF)
  
  ####for test data set
  X1_test <- model.matrix(~ sex * time, data = DF_test)
  Z1_test <- model.matrix(~ time, data = DF_test)
  
  #X2_test <- model.matrix(~ treatment * time, data = DF_test )
  #Z2_test <- model.matrix(~ time, data = DF_test)
  
  #X3_test <- model.matrix(~ time, data = DF_test )
  #Z3_test <- model.matrix(~ time, data = DF_test)
  
  X3_test <- model.matrix(~ sqrt(time), data = DF_test )
  Z3_test <- model.matrix(~ sqrt(time), data = DF_test)
  
  #X4_test <- model.matrix(~ sex + time, data = DF_test)
  #Z4_test <- model.matrix(~ time, data = DF_test)
  
  #X4_test <- model.matrix(~ sex + time, data = DF_test)
  #Z4_test <- model.matrix(~ time, data = DF_test)
  
  #X5_test <- model.matrix(~ time, data = DF_test)
  #Z5_test <- model.matrix(~ time, data = DF_test)
  
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
  
  #Simulate second longitudinal outcome
  ###############################################
  #betas2 <- c(-1.8, -0.06, 0.5, 0.06) # fixed effects coefficients
  #betas2 <- c(10, -0.02, 2.222, 0.07) # fixed effects coefficients
  #sigma2 <- 0.25 # errors sd
  
  # we simulate random effects
  #b2 <- b[, c(3,4)]
  # linear predictor
  #eta_y2 <- as.vector(X2 %*% betas2 + rowSums(Z2 * b2[DF$id, ]))
  # we simulate normal longitudinal data
  #DF$y2 <- rnorm(n * K, mean = eta_y2, sd = sigma2)
  # we assume that values below 0 are not observed, and set equal to 0
  #DF$ind2 <- as.numeric(DF$y2 < -4)
  #DF$y2 <- pmax(DF$y2, -4)
  
  ####Test data
  # we simulate random effects
  #b2_test <- b_test[, c(3,4)]
  # linear predictor
  #eta_y2_test <- as.vector(X2_test %*% betas2 + rowSums(Z2_test * b2[DF_test$id, ]))
  # we simulate normal longitudinal data
  #DF_test$y2 <- rnorm(n_test * K, mean = eta_y2_test, sd = sigma2)
  # we assume that values below 0 are not observed, and set equal to 0
  #DF_test$ind2 <- as.numeric(DF_test$y2 < -4)
  #DF_test$y2 <- pmax(DF_test$y2, -4)
  
  #Simulate third longitudinal outcome
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
  
  #Simulate forth longitudinal outcome
  ###############################################
  #betas4 <- c(0.01, 0.5, -0.31416) # fixed effects coefficients
  #betas4 <- c(80, -4.25, -3.789) # fixed effects coefficients
  
  # we simulate random effects
  #b4 <- b[, c(7,8)]
  # linear predictor
  #eta_y4 <- as.vector(X4 %*% betas4 + rowSums(Z4 * b4[DF$id, ]))
  # mean of the binomial distribution
  #mu_y4 <- plogis(eta_y4)
  # we simulate binomial longitudinal data
  #DF$y4 <- rbinom(n * K, size = 1, prob = mu_y4)
  #DF$y4 <- rnorm(n_test * K, mean = eta_y4, sd = sigma3)
  
  #####Test data
  # we simulate random effects
  #b4_test <- b_test[, c(7,8)]
  # linear predictor
  #eta_y4_test <- as.vector(X4_test %*% betas4 + rowSums(Z4_test * b4_test[DF_test$id, ]))
  # mean of the binomial distribution
  #mu_y4_test <- plogis(eta_y4_test)
  # we simulate binomial longitudinal data
  #DF_test$y4 <- rbinom(n_test * K, size = 1, prob = mu_y4_test)
  #DF_test$y4 <- rnorm(n_test * K, mean = eta_y4_test, sd = sigma3)
  
  #Simulate fifth longitudinal outcome
  ###############################################
  betas5 <- c(1, 0.155, 0.12345) # fixed effects coefficients
  betas5 <- c(70, -5.44, -3.761289) # fixed effects coefficients
  
  # we simulate random effects
  b5 <- b[, c(5,6)]
  # linear predictor
  eta_y5 <- as.vector(X5 %*% betas5 + rowSums(Z5 * b5[DF$id, ]))
  # mean of the binomial distribution
  #mu_y5 <- plogis(eta_y5)
  # we simulate binomial longitudinal data
  #DF$y5 <- rbinom(n * K, size = 1, prob = mu_y5)
  DF$y5 <- rnorm(n_test * K, mean = eta_y5, sd = sigma3)
  
  ####Test data
  # we simulate random effects
  b5_test <- b_test[, c(5,6)]
  # linear predictor
  eta_y5_test <- as.vector(X5_test %*% betas5 + rowSums(Z5_test * b5_test[DF_test$id, ]))
  # mean of the binomial distribution
  #mu_y5_test <- plogis(eta_y5_test)
  # we simulate binomial longitudinal data
  #DF_test$y5 <- rbinom(n_test * K, size = 1, prob = mu_y5_test)
  
  DF_test$y5 <- rnorm(n_test * K, mean = eta_y5_test, sd = sigma3)
  
  #If using admin censoring:
  #AdminCens <- 10
  
  # we keep the longitudinal measurements before the event times: TRAIN
  #DF <- DF[DF$time <= AdminCens, ]
  
  # we keep the longitudinal measurements before the event times: TEST
  #DF_test <- DF_test[DF_test$time <= AdminCens, ]
  
  #Data set with just one obs per individual
  #DF.id <- DF[!duplicated(DF$id),]
  #DF_test.id <- DF_test[!duplicated(DF_test$id),]
  
  
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
  
  ###OBS!!!!!
  ## DF_miss contiene el data frame con missing data!!
  
  ######################################
  #Fitting multivariate MM models
  # - We use snat_mvmer() from rstanarm
  # - We fit the true model plus two
  #.  (or 3) more models
  ######################################
  
  #checking time of computation using data frame with missings 
  # when n=300, k=10, we have 21mins, max rhat 1.67
  #      n=250, k=10: 17mins, max rhat 1.51
  #      n=200, k=10:  8.6 mins, max rhat 1.09
  try(true_model_miss <- stan_mvmer(
    formula = list(
      y1 ~ sex * time + (time | id),
      y3 ~ sqrt(time) + (sqrt(time) | id),
      y5 ~ treatment + time + (time | id)),
    data = list_DF_miss[[count]],
    family = list(gaussian, gaussian, gaussian),
    chains = 3, cores = 8, seed = 12345, iter = 1750))
  
  try(rhat_true_model_miss <- rhat(true_model_miss))
  
  #OJO!!! para que devuelva las predictions a pesar de estar pasandole NAs 
  #hay que pasarle el data frame sin los NAs
  
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
  
  #plot(pred_mfpca1_test)
  
  #AQUI YA TENGO TODA LA INFO QUE QUIERO! las predicciones de esto último se 
  #hacen para cada uno de los puntos de la grid!! También para aquellos donde
  #hay NA! Por tanto para aquellos puntos estoy haciendo una predicción!
  #Y es esa predicción la que me tengo que quedar y entrar a comparar con 
  #la predicción que me puedan devolver los MM
  
  #Cuando tenemos el predict(MFPCAfit), el objeto que pred_mfpca1[[1]]@X es
  # una matriz n X max_length +1, ie cada fila es un individuo y cada columna
  # es uno de los puntos de la fixed time grid. Y esto para el outcome 
  # longitudinal 1 
  # Lo que tenemos aquí es una predicción para cada punto de la fixed grid
  # del outcome longitudinal 1.
  
  
  
  #####################################################
  ## Comparing predictions by means of standardized
  ## RMSE
  #####################################################
  
  
  
  
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



###############################################################################
###############################################################################
## Increasing MAR 30% of dropout
###############################################################################
###############################################################################
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
## I generate the data frames: increasing MAR 30% of dropout
##.  - participant drops at t_ij with a prob determined by
##.    a logistic model with indicator Y_i(j-1) as a
##.    predictor
############################################################
############################################################
for(count in 1:num_datasets){
  n <- 200 # number of subjects
  n_test <- 200
  K <- 10 # number of measurements per subject
  t_max <- 10 # maximum follow-up time
  min_sep <- 0.5
  
  # we construct a data frame with the design:
  # everyone has a baseline measurement, and then measurements at random
  # follow-up times up to t_max
  # !! for Machine Learning purposes: two consecutive measurements will have
  #more than 0.5 of
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
  
  
  #Simulate second longitudinal outcome
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
  
  #Simulate third longitudinal outcome
  ###############################################
  betas5 <- c(1, 0.155, 0.12345) # fixed effects coefficients
  betas5 <- c(70, -5.44, -3.761289) # fixed effects coefficients
  
  # we simulate random effects
  b5 <- b[, c(5,6)]
  # linear predictor
  eta_y5 <- as.vector(X5 %*% betas5 + rowSums(Z5 * b5[DF$id, ]))
  # mean of the binomial distribution
  #mu_y5 <- plogis(eta_y5)
  # we simulate binomial longitudinal data
  #DF$y5 <- rbinom(n * K, size = 1, prob = mu_y5)
  DF$y5 <- rnorm(n_test * K, mean = eta_y5, sd = sigma3)
  
  ####Test data
  # we simulate random effects
  b5_test <- b_test[, c(5,6)]
  # linear predictor
  eta_y5_test <- as.vector(X5_test %*% betas5 + rowSums(Z5_test * b5_test[DF_test$id, ]))
  # mean of the binomial distribution
  #mu_y5_test <- plogis(eta_y5_test)
  # we simulate binomial longitudinal data
  #DF_test$y5 <- rbinom(n_test * K, size = 1, prob = mu_y5_test)
  
  DF_test$y5 <- rnorm(n_test * K, mean = eta_y5_test, sd = sigma3)
  
  #If using admin censoring:
  #AdminCens <- 10
  
  # we keep the longitudinal measurements before the event times: TRAIN
  #DF <- DF[DF$time <= AdminCens, ]
  
  # we keep the longitudinal measurements before the event times: TEST
  #DF_test <- DF_test[DF_test$time <= AdminCens, ]
  
  #Data set with just one obs per individual
  #DF.id <- DF[!duplicated(DF$id),]
  #DF_test.id <- DF_test[!duplicated(DF_test$id),]
  
  param1 <- -3.75
  param2 <- 0.00675
  probs_vec1 <- numeric(n*K)
  probs_vec2 <- numeric(n*K)
  probs_vec3 <- numeric(n*K)
  
  kk <- 1
  for(i in 1:n){
    for(j in 2:(K+1)){
      if((kk-1)%%K==0){
        probs_vec1[kk] <- 0
        probs_vec2[kk] <- 0
        probs_vec3[kk] <- 0
      } else{
        thet1 <- param1 + param2*DF[DF$id==i,]$y1[j-1]
        probs_vec1[kk] <- 1/(1+exp(-thet1))
        thet2 <- param1 + param2*DF[DF$id==i,]$y3[j-1]
        probs_vec2[kk] <- 1/(1+exp(-thet2))
        thet3 <- param1 + param2*DF[DF$id==i,]$y5[j-1]
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
      # checking dropout
      dropout_point = which.max(dropout_y1 == 1 | dropout_y3 == 1 | dropout_y5 == 1),
      
      # NA after dropout
      y1 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y1),
      y3 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y3),
      y5 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y5)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point)
  
  DF_test_miss <- DF_test %>%
    group_by(id) %>%                     # Group by each individual
    mutate(
      # checking dropout
      dropout_point = which.max(dropout_y1 == 1 | dropout_y3 == 1 | dropout_y5 == 1),
      
      # NA after dropout
      y1 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y1),
      y3 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y3),
      y5 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y5)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point)
  
  
  missing_ratio_train[count] <- counting_missing_data(DF_miss, n, K)
  missing_ratio_test[count] <- counting_missing_data(DF_test_miss, n, K)
  
  try(list_DF <- append(list_DF, list(DF)))
  try(list_DF_miss <- append(list_DF_miss, list(DF_miss)))
  try(list_DF_test <- append(list_DF_test, list(DF_test)))
  try(list_DF_test_miss <- append(list_DF_test_miss, list(DF_test_miss)))
  
  try(if(count==100){
    setwd("D:/La meva unitat/TFM/ResultsMMvsMFPCA")
    str100 <- "dataframes_increasing_MAR_30.RData"
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

load("D:/La meva unitat/TFM/ResultsMMvsMFPCA/dataframes_increasing_MAR_30.RData")

for(count in 1:num_datasets){
  
  ###OBS!!!!!
  ## DF_miss contiene el data frame con missing data!!
  
  ######################################
  #Fitting multivariate MM models
  # - We use snat_mvmer() from rstanarm
  # - We fit the true model plus two
  #.  (or 3) more models
  ######################################
  
  #checking time of computation using data frame with missings 
  # when n=300, k=10, we have 21mins, max rhat 1.67
  #      n=250, k=10: 17mins, max rhat 1.51
  #      n=200, k=10:  8.6 mins, max rhat 1.09
  try(true_model_miss <- stan_mvmer(
    formula = list(
      y1 ~ sex * time + (time | id),
      y3 ~ sqrt(time) + (sqrt(time) | id),
      y5 ~ treatment + time + (time | id)),
    data = list_DF_miss[[count]],
    family = list(gaussian, gaussian, gaussian),
    chains = 3, cores = 8, seed = 12345, iter = 1750))
  
  try(rhat_true_model_miss <- rhat(true_model_miss))
  
  #OJO!!! para que devuelva las predictions a pesar de estar pasandole NAs 
  #hay que pasarle el data frame sin los NAs
  
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
  
  #We fit with NO weights and pve=0.8
  try(mfpca1 <- MFPCA(m1, M = 3, 
                      uniExpansions = list(list(type="uFPCA", pve = 0.8),
                                           list(type ="uFPCA", pve=0.8),
                                           list(type="uFPCA", pve=0.8)),
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
  
  #plot(pred_mfpca1_test)
  
  #AQUI YA TENGO TODA LA INFO QUE QUIERO! las predicciones de esto último se 
  #hacen para cada uno de los puntos de la grid!! También para aquellos donde
  #hay NA! Por tanto para aquellos puntos estoy haciendo una predicción!
  #Y es esa predicción la que me tengo que quedar y entrar a comparar con 
  #la predicción que me puedan devolver los MM
  
  #Cuando tenemos el predict(MFPCAfit), el objeto que pred_mfpca1[[1]]@X es
  # una matriz n X max_length +1, ie cada fila es un individuo y cada columna
  # es uno de los puntos de la fixed time grid. Y esto para el outcome 
  # longitudinal 1 
  # Lo que tenemos aquí es una predicción para cada punto de la fixed grid
  # del outcome longitudinal 1.
  
  
  
  
  
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
    strr <- "results_1_mar30_incr_26feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr)
  })
  try(if(count==2){
    strr2 <- "results_2_mar30_incr_26feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr2)
  })
  try(if(count==5){
    str1 <- "results_5_mar30_incr_26feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str1)
  })
  try(if(count==10){
    str10 <- "results_10_mar30_incr_26feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str10)
  })
  try(if(count==20){
    str20 <- "results_20_mar30_incr_26feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str20)
  })
  try(if(count==30){
    str30 <- "results_30_mar30_incr_26feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str30)
  })
  try(if(count==40){
    str40 <- "results_40_mar30_incr_26feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss,file=str40)
  })
  try(if(count==50){
    str50 <- "results_50_mar30_incr_26feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str50)
  })
  try(if(count==75){
    str75 <- "results_75_mar30_incr_26feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_mo111del_miss, file=str75)
  })
  try(if(count==100){
    str100 <- "results_100_mar30_incr_26feb2025.RData"
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


###############################################################################
###############################################################################
## Increasing MAR 50% of dropout
###############################################################################
###############################################################################
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
## I generate the data frames: increasing MAR 50% of dropout
##.  - participant drops at t_ij with a prob determined by
##.    a logistic model with indicator Y_i(j-1) as a 
##.    predictor
############################################################
############################################################
for(count in 1:num_datasets){
  n <- 200 # number of subjects
  n_test <- 200
  K <- 10 # number of measurements per subject
  t_max <- 10 # maximum follow-up time
  min_sep <- 0.5
  
  # we construct a data frame with the design:
  # everyone has a baseline measurement, and then measurements at random 
  # follow-up times up to t_max 
  # !! for Machine Learning purposes: two consecutive measurements will have
  #more than 0.5 of
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
  
  
  #Simulate second longitudinal outcome
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
  
  #Simulate third longitudinal outcome
  ###############################################
  betas5 <- c(1, 0.155, 0.12345) # fixed effects coefficients
  betas5 <- c(70, -5.44, -3.761289) # fixed effects coefficients
  
  # we simulate random effects
  b5 <- b[, c(5,6)]
  # linear predictor
  eta_y5 <- as.vector(X5 %*% betas5 + rowSums(Z5 * b5[DF$id, ]))
  # mean of the binomial distribution
  #mu_y5 <- plogis(eta_y5)
  # we simulate binomial longitudinal data
  #DF$y5 <- rbinom(n * K, size = 1, prob = mu_y5)
  DF$y5 <- rnorm(n_test * K, mean = eta_y5, sd = sigma3)
  
  ####Test data
  # we simulate random effects
  b5_test <- b_test[, c(5,6)]
  # linear predictor
  eta_y5_test <- as.vector(X5_test %*% betas5 + rowSums(Z5_test * b5_test[DF_test$id, ]))
  # mean of the binomial distribution
  #mu_y5_test <- plogis(eta_y5_test)
  # we simulate binomial longitudinal data
  #DF_test$y5 <- rbinom(n_test * K, size = 1, prob = mu_y5_test)
  
  DF_test$y5 <- rnorm(n_test * K, mean = eta_y5_test, sd = sigma3)
  
  #If using admin censoring:
  #AdminCens <- 10
  
  # we keep the longitudinal measurements before the event times: TRAIN
  #DF <- DF[DF$time <= AdminCens, ]
  
  # we keep the longitudinal measurements before the event times: TEST
  #DF_test <- DF_test[DF_test$time <= AdminCens, ]
  
  #Data set with just one obs per individual
  #DF.id <- DF[!duplicated(DF$id),]
  #DF_test.id <- DF_test[!duplicated(DF_test$id),]
  
  param1 <- -3
  param2 <- 0.0099
  probs_vec1 <- numeric(n*K)
  probs_vec2 <- numeric(n*K)
  probs_vec3 <- numeric(n*K)
  
  kk <- 1
  for(i in 1:n){
    for(j in 2:(K+1)){
      if((kk-1)%%K==0){
        probs_vec1[kk] <- 0
        probs_vec2[kk] <- 0
        probs_vec3[kk] <- 0
      } else{
        thet1 <- param1 + param2*DF[DF$id==i,]$y1[j-1]
        probs_vec1[kk] <- 1/(1+exp(-thet1))
        thet2 <- param1 + param2*DF[DF$id==i,]$y3[j-1]
        probs_vec2[kk] <- 1/(1+exp(-thet2))
        thet3 <- param1 + param2*DF[DF$id==i,]$y5[j-1]
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
      # checking dropout
      dropout_point = which.max(dropout_y1 == 1 | dropout_y3 == 1 | dropout_y5 == 1),
      
      # NA after dropout
      y1 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y1),
      y3 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y3),
      y5 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y5)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point) 
  
  DF_test_miss <- DF_test %>%
    group_by(id) %>%                     # Group by each individual
    mutate(
      # checking dropout
      dropout_point = which.max(dropout_y1 == 1 | dropout_y3 == 1 | dropout_y5 == 1),
      
      # NA after dropout
      y1 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y1),
      y3 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y3),
      y5 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y5)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point) 
  
  
  missing_ratio_train[count] <- counting_missing_data(DF_miss, n, K)
  missing_ratio_test[count] <- counting_missing_data(DF_test_miss, n, K)
  
  try(list_DF <- append(list_DF, list(DF)))
  try(list_DF_miss <- append(list_DF_miss, list(DF_miss)))
  try(list_DF_test <- append(list_DF_test, list(DF_test)))
  try(list_DF_test_miss <- append(list_DF_test_miss, list(DF_test_miss)))
  
  try(if(count==100){
    setwd("D:/La meva unitat/TFM/ResultsMMvsMFPCA")
    str100 <- "dataframes_increasing_MAR_50.RData"
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
  
  ###OBS!!!!!
  ## DF_miss contiene el data frame con missing data!!
  
  ######################################
  #Fitting multivariate MM models
  # - We use snat_mvmer() from rstanarm
  # - We fit the true model plus two
  #.  (or 3) more models
  ######################################
  
  #checking time of computation using data frame with missings 
  # when n=300, k=10, we have 21mins, max rhat 1.67
  #      n=250, k=10: 17mins, max rhat 1.51
  #      n=200, k=10:  8.6 mins, max rhat 1.09
  try(true_model_miss <- stan_mvmer(
    formula = list(
      y1 ~ sex * time + (time | id),
      y3 ~ sqrt(time) + (sqrt(time) | id),
      y5 ~ treatment + time + (time | id)),
    data = list_DF_miss[[count]],
    family = list(gaussian, gaussian, gaussian),
    chains = 3, cores = 8, seed = 12345, iter = 1750))
  
  try(rhat_true_model_miss <- rhat(true_model_miss))
  
  #OJO!!! para que devuelva las predictions a pesar de estar pasandole NAs 
  #hay que pasarle el data frame sin los NAs
  
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
                           uniExpansions = list(list(type="uFPCA", pve = 0.65),
                                                list(type ="uFPCA", pve=0.65),
                                                list(type="uFPCA", pve=0.65)),
                           fit = TRUE))
  
  try(pred_mfpca1_test <- predict(mfpca1, scores = mfpca1_test$scores))
  
  try(pred_mfpca1_y1_test <- NULL)
  try(pred_mfpca1_y1_test <- pred_mfpca1_test[[1]]@X)
  try(pred_mfpca1_y2_test <- NULL)
  try(pred_mfpca1_y2_test <- pred_mfpca1_test[[2]]@X)
  try(pred_mfpca1_y3_test <- NULL)
  try(pred_mfpca1_y3_test <- pred_mfpca1_test[[3]]@X)
  
  #plot(pred_mfpca1_test)
  
  #AQUI YA TENGO TODA LA INFO QUE QUIERO! las predicciones de esto último se 
  #hacen para cada uno de los puntos de la grid!! También para aquellos donde
  #hay NA! Por tanto para aquellos puntos estoy haciendo una predicción!
  #Y es esa predicción la que me tengo que quedar y entrar a comparar con 
  #la predicción que me puedan devolver los MM
  
  #Cuando tenemos el predict(MFPCAfit), el objeto que pred_mfpca1[[1]]@X es
  # una matriz n X max_length +1, ie cada fila es un individuo y cada columna
  # es uno de los puntos de la fixed time grid. Y esto para el outcome 
  # longitudinal 1 
  # Lo que tenemos aquí es una predicción para cada punto de la fixed grid
  # del outcome longitudinal 1.
  
  
  
  #####################################################
  ## Comparing predictions by means of standardized
  ## RMSE
  #####################################################
  
  
  
  
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
    strr <- "results_1_MAR50_incr_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr)
  })
  try(if(count==2){
    strr2 <- "results_2_MAR50_incr_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr2)
  })
  try(if(count==5){
    str1 <- "results_5_MAR50_incr_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str1)
  })
  try(if(count==10){
    str10 <- "results_10_MAR50_incr_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str10)
  })
  try(if(count==20){
    str20 <- "results_20_MAR50_incr_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str20)
  })
  try(if(count==30){
    str30 <- "results_30_MAR50_incr_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str30)
  })
  try(if(count==40){
    str40 <- "results_40_MAR50_incr_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss,file=str40)
  })
  try(if(count==50){
    str50 <- "results_50_MAR50_incr_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str50)
  })
  try(if(count==75){
    str75 <- "results_75_MAR50_incr_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str75)
  })
  try(if(count==100){
    str100 <- "results_100_MAR50_incr_28feb2025.RData"
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


###############################################################################
###############################################################################
## Threshold MNAR 30% of dropout
###############################################################################
###############################################################################
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
## I generate the data frames: therhsold MNAR 30% of dropout
##.  - participant drops at t_ij with a prob determined by
##.    a logistic model with indicator Y_i(j)>nu as a
##.    predictor
############################################################
############################################################
for(count in 1:num_datasets){
  n <- 200 # number of subjects
  n_test <- 200
  K <- 10 # number of measurements per subject
  t_max <- 10 # maximum follow-up time
  min_sep <- 0.5
  
  # we construct a data frame with the design:
  # everyone has a baseline measurement, and then measurements at random
  # follow-up times up to t_max
  # !! for Machine Learning purposes: two consecutive measurements will have
  #more than 0.5 of
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
  
  
  #Simulate second longitudinal outcome
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
  
  #Simulate third longitudinal outcome
  ###############################################
  betas5 <- c(1, 0.155, 0.12345) # fixed effects coefficients
  betas5 <- c(70, -5.44, -3.761289) # fixed effects coefficients
  
  # we simulate random effects
  b5 <- b[, c(5,6)]
  # linear predictor
  eta_y5 <- as.vector(X5 %*% betas5 + rowSums(Z5 * b5[DF$id, ]))
  # mean of the binomial distribution
  #mu_y5 <- plogis(eta_y5)
  # we simulate binomial longitudinal data
  #DF$y5 <- rbinom(n * K, size = 1, prob = mu_y5)
  DF$y5 <- rnorm(n_test * K, mean = eta_y5, sd = sigma3)
  
  ####Test data
  # we simulate random effects
  b5_test <- b_test[, c(5,6)]
  # linear predictor
  eta_y5_test <- as.vector(X5_test %*% betas5 + rowSums(Z5_test * b5_test[DF_test$id, ]))
  # mean of the binomial distribution
  #mu_y5_test <- plogis(eta_y5_test)
  # we simulate binomial longitudinal data
  #DF_test$y5 <- rbinom(n_test * K, size = 1, prob = mu_y5_test)
  
  DF_test$y5 <- rnorm(n_test * K, mean = eta_y5_test, sd = sigma3)
  
  #If using admin censoring:
  #AdminCens <- 10
  
  # we keep the longitudinal measurements before the event times: TRAIN
  #DF <- DF[DF$time <= AdminCens, ]
  
  # we keep the longitudinal measurements before the event times: TEST
  #DF_test <- DF_test[DF_test$time <= AdminCens, ]
  
  #Data set with just one obs per individual
  #DF.id <- DF[!duplicated(DF$id),]
  #DF_test.id <- DF_test[!duplicated(DF_test$id),]
  
  
  ### threshold MAR dropout with 30% intensity
  nu1 <- 15.5
  nu2 <- 7.5
  nu3 <- 42
  
  DF <- DF %>%
    group_by(id) %>%                      
    mutate(
      # check conditions in the same obs (not the previous one)
      dropout_point = which.max(y1 > nu1 | y3 > nu2 | y5 < nu3),
      
      # crearting indicators
      y1_ind = ifelse(row_number() >= dropout_point & dropout_point > 1, 1, 0),
      y3_ind = ifelse(row_number() >= dropout_point & dropout_point > 1, 1, 0),
      y5_ind = ifelse(row_number() >= dropout_point & dropout_point > 1, 1, 0)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point)
  
  DF_test <- DF_test %>%
    group_by(id) %>%                      
    mutate(
      dropout_point = which.max(y1 > nu1 | y3 > nu2 | y5 < nu3),
      
      # Set NA to y1, y2, y3 after the dropout point (ignoring the first visit)
      y1_ind = ifelse(row_number() >= dropout_point & dropout_point > 1, 1, 0),
      y3_ind = ifelse(row_number() >= dropout_point & dropout_point > 1, 1, 0),
      y5_ind = ifelse(row_number() >= dropout_point & dropout_point > 1, 1, 0)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point)
  
  param1 <- -3.5
  param2 <- 0.7
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
      # checking dropout
      dropout_point = which.max(dropout_y1 == 1 | dropout_y3 == 1 | dropout_y5 == 1),
      
      # NA after dropout
      y1 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y1),
      y3 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y3),
      y5 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y5)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point)
  
  DF_test_miss <- DF_test %>%
    group_by(id) %>%                     # Group by each individual
    mutate(
      # checking dropout
      dropout_point = which.max(dropout_y1 == 1 | dropout_y3 == 1 | dropout_y5 == 1),
      
      # NA after dropout
      y1 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y1),
      y3 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y3),
      y5 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y5)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point)
  
  
  missing_ratio_train[count] <- counting_missing_data(DF_miss, n, K)
  missing_ratio_test[count] <- counting_missing_data(DF_test_miss, n, K)
  
  try(list_DF <- append(list_DF, list(DF)))
  try(list_DF_miss <- append(list_DF_miss, list(DF_miss)))
  try(list_DF_test <- append(list_DF_test, list(DF_test)))
  try(list_DF_test_miss <- append(list_DF_test_miss, list(DF_test_miss)))
  
  try(if(count==100){
    setwd("D:/La meva unitat/TFM/ResultsMMvsMFPCA")
    str100 <- "dataframes_MNAR_th_30.RData"
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
  
  ###OBS!!!!!
  ## DF_miss contiene el data frame con missing data!!
  
  ######################################
  #Fitting multivariate MM models
  # - We use snat_mvmer() from rstanarm
  # - We fit the true model plus two
  #.  (or 3) more models
  ######################################
  
  #checking time of computation using data frame with missings 
  # when n=300, k=10, we have 21mins, max rhat 1.67
  #      n=250, k=10: 17mins, max rhat 1.51
  #      n=200, k=10:  8.6 mins, max rhat 1.09
  try(true_model_miss <- stan_mvmer(
    formula = list(
      y1 ~ sex * time + (time | id),
      y3 ~ sqrt(time) + (sqrt(time) | id),
      y5 ~ treatment + time + (time | id)),
    data = list_DF_miss[[count]],
    family = list(gaussian, gaussian, gaussian),
    chains = 3, cores = 8, seed = 12345, iter = 1750))
  
  try(rhat_true_model_miss <- rhat(true_model_miss))
  
  #OJO!!! para que devuelva las predictions a pesar de estar pasandole NAs 
  #hay que pasarle el data frame sin los NAs
  
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
  
  #We fit with NO weights and pve=0.8
  try(mfpca1 <- MFPCA(m1, M = 3, 
                      uniExpansions = list(list(type="uFPCA", pve = 0.8),
                                           list(type ="uFPCA", pve=0.8),
                                           list(type="uFPCA", pve=0.8)),
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
  
  #plot(pred_mfpca1_test)
  
  #AQUI YA TENGO TODA LA INFO QUE QUIERO! las predicciones de esto último se 
  #hacen para cada uno de los puntos de la grid!! También para aquellos donde
  #hay NA! Por tanto para aquellos puntos estoy haciendo una predicción!
  #Y es esa predicción la que me tengo que quedar y entrar a comparar con 
  #la predicción que me puedan devolver los MM
  
  #Cuando tenemos el predict(MFPCAfit), el objeto que pred_mfpca1[[1]]@X es
  # una matriz n X max_length +1, ie cada fila es un individuo y cada columna
  # es uno de los puntos de la fixed time grid. Y esto para el outcome 
  # longitudinal 1 
  # Lo que tenemos aquí es una predicción para cada punto de la fixed grid
  # del outcome longitudinal 1.
  
  
  
  #####################################################
  ## Comparing predictions by means of standardized
  ## RMSE
  #####################################################
  
  
  
  
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
    strr <- "results_1_MNAR30_th_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr)
  })
  try(if(count==2){
    strr2 <- "results_2_MNAR30_th_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr2)
  })
  try(if(count==5){
    str1 <- "results_5_MNAR30_th_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str1)
  })
  try(if(count==10){
    str10 <- "results_10_MNAR30_th_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str10)
  })
  try(if(count==20){
    str20 <- "results_20_MNAR30_th_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str20)
  })
  try(if(count==30){
    str30 <- "results_30_MNAR30_th_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str30)
  })
  try(if(count==40){
    str40 <- "results_40_MNAR30_th_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss,file=str40)
  })
  try(if(count==50){
    str50 <- "results_50_MNAR30_th_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str50)
  })
  try(if(count==75){
    str75 <- "results_75_MNAR30_th_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str75)
  })
  try(if(count==100){
    str100 <- "results_100_MNAR30_th_28feb2025.RData"
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



##############################################################################
###############################################################################
## Threshold MNAR 50% of dropout
###############################################################################
###############################################################################
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
## I generate the data frames: therhsold MNAR 50% of dropout
##.  - participant drops at t_ij with a prob determined by
##.    a logistic model with indicator Y_i(j)>nu as a
##.    predictor
############################################################
############################################################
for(count in 1:num_datasets){
  n <- 200 # number of subjects
  n_test <- 200
  K <- 10 # number of measurements per subject
  t_max <- 10 # maximum follow-up time
  min_sep <- 0.5
  
  # we construct a data frame with the design:
  # everyone has a baseline measurement, and then measurements at random
  # follow-up times up to t_max
  # !! for Machine Learning purposes: two consecutive measurements will have
  #more than 0.5 of
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
  
  
  #Simulate second longitudinal outcome
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
  
  #Simulate third longitudinal outcome
  ###############################################
  betas5 <- c(1, 0.155, 0.12345) # fixed effects coefficients
  betas5 <- c(70, -5.44, -3.761289) # fixed effects coefficients
  
  # we simulate random effects
  b5 <- b[, c(5,6)]
  # linear predictor
  eta_y5 <- as.vector(X5 %*% betas5 + rowSums(Z5 * b5[DF$id, ]))
  # mean of the binomial distribution
  #mu_y5 <- plogis(eta_y5)
  # we simulate binomial longitudinal data
  #DF$y5 <- rbinom(n * K, size = 1, prob = mu_y5)
  DF$y5 <- rnorm(n_test * K, mean = eta_y5, sd = sigma3)
  
  ####Test data
  # we simulate random effects
  b5_test <- b_test[, c(5,6)]
  # linear predictor
  eta_y5_test <- as.vector(X5_test %*% betas5 + rowSums(Z5_test * b5_test[DF_test$id, ]))
  # mean of the binomial distribution
  #mu_y5_test <- plogis(eta_y5_test)
  # we simulate binomial longitudinal data
  #DF_test$y5 <- rbinom(n_test * K, size = 1, prob = mu_y5_test)
  
  DF_test$y5 <- rnorm(n_test * K, mean = eta_y5_test, sd = sigma3)
  
  #If using admin censoring:
  #AdminCens <- 10
  
  # we keep the longitudinal measurements before the event times: TRAIN
  #DF <- DF[DF$time <= AdminCens, ]
  
  # we keep the longitudinal measurements before the event times: TEST
  #DF_test <- DF_test[DF_test$time <= AdminCens, ]
  
  #Data set with just one obs per individual
  #DF.id <- DF[!duplicated(DF$id),]
  #DF_test.id <- DF_test[!duplicated(DF_test$id),]
  
  
  ### threshold MAR dropout with 30% intensity
  nu1 <- 15.5
  nu2 <- 7.5
  nu3 <- 42
  
  DF <- DF %>%
    group_by(id) %>%                      
    mutate(
      # check conditions in the same obs (not the previous one)
      dropout_point = which.max(y1 > nu1 | y3 > nu2 | y5 < nu3),
      
      # crearting indicators
      y1_ind = ifelse(row_number() >= dropout_point & dropout_point > 1, 1, 0),
      y3_ind = ifelse(row_number() >= dropout_point & dropout_point > 1, 1, 0),
      y5_ind = ifelse(row_number() >= dropout_point & dropout_point > 1, 1, 0)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point)
  
  DF_test <- DF_test %>%
    group_by(id) %>%                      
    mutate(
      dropout_point = which.max(y1 > nu1 | y3 > nu2 | y5 < nu3),
      
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
      # checking dropout
      dropout_point = which.max(dropout_y1 == 1 | dropout_y3 == 1 | dropout_y5 == 1),
      
      # NA after dropout
      y1 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y1),
      y3 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y3),
      y5 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y5)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point)
  
  DF_test_miss <- DF_test %>%
    group_by(id) %>%                     # Group by each individual
    mutate(
      # checking dropout
      dropout_point = which.max(dropout_y1 == 1 | dropout_y3 == 1 | dropout_y5 == 1),
      
      # NA after dropout
      y1 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y1),
      y3 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y3),
      y5 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y5)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point)
  
  
  missing_ratio_train[count] <- counting_missing_data(DF_miss, n, K)
  missing_ratio_test[count] <- counting_missing_data(DF_test_miss, n, K)
  
  try(list_DF <- append(list_DF, list(DF)))
  try(list_DF_miss <- append(list_DF_miss, list(DF_miss)))
  try(list_DF_test <- append(list_DF_test, list(DF_test)))
  try(list_DF_test_miss <- append(list_DF_test_miss, list(DF_test_miss)))
  
  try(if(count==100){
    setwd("D:/La meva unitat/TFM/ResultsMMvsMFPCA")
    str100 <- "dataframes_MNAR_th_50.RData"
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
  
  ###OBS!!!!!
  ## DF_miss contiene el data frame con missing data!!
  
  ######################################
  #Fitting multivariate MM models
  # - We use snat_mvmer() from rstanarm
  # - We fit the true model plus two
  #.  (or 3) more models
  ######################################
  
  #checking time of computation using data frame with missings 
  # when n=300, k=10, we have 21mins, max rhat 1.67
  #      n=250, k=10: 17mins, max rhat 1.51
  #      n=200, k=10:  8.6 mins, max rhat 1.09
  try(true_model_miss <- stan_mvmer(
    formula = list(
      y1 ~ sex * time + (time | id),
      y3 ~ sqrt(time) + (sqrt(time) | id),
      y5 ~ treatment + time + (time | id)),
    data = list_DF_miss[[count]],
    family = list(gaussian, gaussian, gaussian),
    chains = 3, cores = 8, seed = 12345, iter = 1750))
  
  try(rhat_true_model_miss <- rhat(true_model_miss))
  
  #OJO!!! para que devuelva las predictions a pesar de estar pasandole NAs 
  #hay que pasarle el data frame sin los NAs
  
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
                           uniExpansions = list(list(type="uFPCA", pve = 0.65),
                                                list(type ="uFPCA", pve=0.65),
                                                list(type="uFPCA", pve=0.65)),
                           fit = TRUE))
  
  try(pred_mfpca1_test <- predict(mfpca1, scores = mfpca1_test$scores))
  
  try(pred_mfpca1_y1_test <- NULL)
  try(pred_mfpca1_y1_test <- pred_mfpca1_test[[1]]@X)
  try(pred_mfpca1_y2_test <- NULL)
  try(pred_mfpca1_y2_test <- pred_mfpca1_test[[2]]@X)
  try(pred_mfpca1_y3_test <- NULL)
  try(pred_mfpca1_y3_test <- pred_mfpca1_test[[3]]@X)
  
  #plot(pred_mfpca1_test)
  
  #AQUI YA TENGO TODA LA INFO QUE QUIERO! las predicciones de esto último se 
  #hacen para cada uno de los puntos de la grid!! También para aquellos donde
  #hay NA! Por tanto para aquellos puntos estoy haciendo una predicción!
  #Y es esa predicción la que me tengo que quedar y entrar a comparar con 
  #la predicción que me puedan devolver los MM
  
  #Cuando tenemos el predict(MFPCAfit), el objeto que pred_mfpca1[[1]]@X es
  # una matriz n X max_length +1, ie cada fila es un individuo y cada columna
  # es uno de los puntos de la fixed time grid. Y esto para el outcome 
  # longitudinal 1 
  # Lo que tenemos aquí es una predicción para cada punto de la fixed grid
  # del outcome longitudinal 1.
  
  
  
  #####################################################
  ## Comparing predictions by means of standardized
  ## RMSE
  #####################################################
  
  
  
  
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
    strr <- "results_1_MNAR50_th_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr)
  })
  try(if(count==2){
    strr2 <- "results_2_MNAR50_th_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr2)
  })
  try(if(count==5){
    str1 <- "results_5_MNAR50_th_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str1)
  })
  try(if(count==10){
    str10 <- "results_10_MNAR50_th_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str10)
  })
  try(if(count==20){
    str20 <- "results_20_MNAR50_th_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str20)
  })
  try(if(count==30){
    str30 <- "results_30_MNAR50_th_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str30)
  })
  try(if(count==40){
    str40 <- "results_40_MNAR50_th_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss,file=str40)
  })
  try(if(count==50){
    str50 <- "results_50_MNAR50_th_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str50)
  })
  try(if(count==75){
    str75 <- "results_75_MNAR50_th_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str75)
  })
  try(if(count==100){
    str100 <- "results_100_MNAR50_th_28feb2025.RData"
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


###############################################################################
###############################################################################
## Increasing MNAR 30% of dropout
###############################################################################
###############################################################################
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
## I generate the data frames: increasing MNAR 30% of dropout
##.  - participant drops at t_ij with a prob determined by
##.    a logistic model with indicator Y_ij as a
##.    predictor
############################################################
############################################################
for(count in 1:num_datasets){
  n <- 200 # number of subjects
  n_test <- 200
  K <- 10 # number of measurements per subject
  t_max <- 10 # maximum follow-up time
  min_sep <- 0.5
  
  # we construct a data frame with the design:
  # everyone has a baseline measurement, and then measurements at random
  # follow-up times up to t_max
  # !! for Machine Learning purposes: two consecutive measurements will have
  #more than 0.5 of
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
  
  
  #Simulate second longitudinal outcome
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
  
  #Simulate third longitudinal outcome
  ###############################################
  betas5 <- c(1, 0.155, 0.12345) # fixed effects coefficients
  betas5 <- c(70, -5.44, -3.761289) # fixed effects coefficients
  
  # we simulate random effects
  b5 <- b[, c(5,6)]
  # linear predictor
  eta_y5 <- as.vector(X5 %*% betas5 + rowSums(Z5 * b5[DF$id, ]))
  # mean of the binomial distribution
  #mu_y5 <- plogis(eta_y5)
  # we simulate binomial longitudinal data
  #DF$y5 <- rbinom(n * K, size = 1, prob = mu_y5)
  DF$y5 <- rnorm(n_test * K, mean = eta_y5, sd = sigma3)
  
  ####Test data
  # we simulate random effects
  b5_test <- b_test[, c(5,6)]
  # linear predictor
  eta_y5_test <- as.vector(X5_test %*% betas5 + rowSums(Z5_test * b5_test[DF_test$id, ]))
  # mean of the binomial distribution
  #mu_y5_test <- plogis(eta_y5_test)
  # we simulate binomial longitudinal data
  #DF_test$y5 <- rbinom(n_test * K, size = 1, prob = mu_y5_test)
  
  DF_test$y5 <- rnorm(n_test * K, mean = eta_y5_test, sd = sigma3)
  
  #If using admin censoring:
  #AdminCens <- 10
  
  # we keep the longitudinal measurements before the event times: TRAIN
  #DF <- DF[DF$time <= AdminCens, ]
  
  # we keep the longitudinal measurements before the event times: TEST
  #DF_test <- DF_test[DF_test$time <= AdminCens, ]
  
  #Data set with just one obs per individual
  #DF.id <- DF[!duplicated(DF$id),]
  #DF_test.id <- DF_test[!duplicated(DF_test$id),]
  
  param1 <- -3.75
  param2 <- 0.0069
  probs_vec1 <- numeric(n*K)
  probs_vec2 <- numeric(n*K)
  probs_vec3 <- numeric(n*K)
  
  kk <- 1
  for(i in 1:n){
    for(j in 1:K){
      if((kk-1)%%K==0){
        probs_vec1[kk] <- 0
        probs_vec2[kk] <- 0
        probs_vec3[kk] <- 0
      } else{
        thet1 <- param1 + param2*DF[DF$id==i,]$y1[j]
        probs_vec1[kk] <- 1/(1+exp(-thet1))
        thet2 <- param1 + param2*DF[DF$id==i,]$y3[j]
        probs_vec2[kk] <- 1/(1+exp(-thet2))
        thet3 <- param1 + param2*DF[DF$id==i,]$y5[j]
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
      # checking dropout
      dropout_point = which.max(dropout_y1 == 1 | dropout_y3 == 1 | dropout_y5 == 1),
      
      # NA after dropout
      y1 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y1),
      y3 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y3),
      y5 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y5)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point)
  
  DF_test_miss <- DF_test %>%
    group_by(id) %>%                     # Group by each individual
    mutate(
      # checking dropout
      dropout_point = which.max(dropout_y1 == 1 | dropout_y3 == 1 | dropout_y5 == 1),
      
      # NA after dropout
      y1 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y1),
      y3 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y3),
      y5 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y5)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point)
  
  
  missing_ratio_train[count] <- counting_missing_data(DF_miss, n, K)
  missing_ratio_test[count] <- counting_missing_data(DF_test_miss, n, K)
  
  try(list_DF <- append(list_DF, list(DF)))
  try(list_DF_miss <- append(list_DF_miss, list(DF_miss)))
  try(list_DF_test <- append(list_DF_test, list(DF_test)))
  try(list_DF_test_miss <- append(list_DF_test_miss, list(DF_test_miss)))
  
  try(if(count==100){
    setwd("D:/La meva unitat/TFM/ResultsMMvsMFPCA")
    str100 <- "dataframes_increasing_MNAR_30.RData"
    save(list_DF, list_DF_miss, list_DF_test, list_DF_test_miss,
         missing_ratio_train, missing_ratio_test, file=str100)
  })
  
  print(count)
}


###################################################
####################################################
## Fitting MM and MFPCA models and doing predictions
####################################################
####################################################

for(count in 1:num_datasets){
  
  ###OBS!!!!!
  ## DF_miss contiene el data frame con missing data!!
  
  ######################################
  #Fitting multivariate MM models
  # - We use snat_mvmer() from rstanarm
  # - We fit the true model plus two
  #.  (or 3) more models
  ######################################
  
  #checking time of computation using data frame with missings 
  # when n=300, k=10, we have 21mins, max rhat 1.67
  #      n=250, k=10: 17mins, max rhat 1.51
  #      n=200, k=10:  8.6 mins, max rhat 1.09
  try(true_model_miss <- stan_mvmer(
    formula = list(
      y1 ~ sex * time + (time | id),
      y3 ~ sqrt(time) + (sqrt(time) | id),
      y5 ~ treatment + time + (time | id)),
    data = list_DF_miss[[count]],
    family = list(gaussian, gaussian, gaussian),
    chains = 3, cores = 8, seed = 12345, iter = 1750))
  
  try(rhat_true_model_miss <- rhat(true_model_miss))
  
  #OJO!!! para que devuelva las predictions a pesar de estar pasandole NAs 
  #hay que pasarle el data frame sin los NAs
  
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
  
  #We fit with NO weights and pve=0.8
  try(mfpca1 <- MFPCA(m1, M = 3, 
                      uniExpansions = list(list(type="uFPCA", pve = 0.8),
                                           list(type ="uFPCA", pve=0.8),
                                           list(type="uFPCA", pve=0.8)),
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
  
  #plot(pred_mfpca1_test)
  
  #AQUI YA TENGO TODA LA INFO QUE QUIERO! las predicciones de esto último se 
  #hacen para cada uno de los puntos de la grid!! También para aquellos donde
  #hay NA! Por tanto para aquellos puntos estoy haciendo una predicción!
  #Y es esa predicción la que me tengo que quedar y entrar a comparar con 
  #la predicción que me puedan devolver los MM
  
  #Cuando tenemos el predict(MFPCAfit), el objeto que pred_mfpca1[[1]]@X es
  # una matriz n X max_length +1, ie cada fila es un individuo y cada columna
  # es uno de los puntos de la fixed time grid. Y esto para el outcome 
  # longitudinal 1 
  # Lo que tenemos aquí es una predicción para cada punto de la fixed grid
  # del outcome longitudinal 1.
  
  
  
  #####################################################
  ## Comparing predictions by means of standardized
  ## RMSE
  #####################################################
  
  
  
  
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
    strr <- "results_1_MNAR30_incr_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr)
  })
  try(if(count==2){
    strr2 <- "results_2_MNAR30_incr_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr2)
  })
  try(if(count==5){
    str1 <- "results_5_MNAR30_incr_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str1)
  })
  try(if(count==10){
    str10 <- "results_10_MNAR30_incr_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str10)
  })
  try(if(count==20){
    str20 <- "results_20_MNAR30_incr_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str20)
  })
  try(if(count==30){
    str30 <- "results_30_MNAR30_incr_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str30)
  })
  try(if(count==40){
    str40 <- "results_40_MNAR30_incr_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss,file=str40)
  })
  try(if(count==50){
    str50 <- "results_50_MNAR30_incr_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str50)
  })
  try(if(count==75){
    str75 <- "results_75_MNAR30_incr_28feb2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str75)
  })
  try(if(count==100){
    str100 <- "results_100_MNAR30_incr_28feb2025.RData"
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


################################
###############################
###simuls remaining: 
## 1) INCREASING MNAR 50%
############################
##############################

###############################################################################
###############################################################################
## Increasing MNAR 50% of dropout
###############################################################################
###############################################################################
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
## I generate the data frames: increasing MNAR 50% of dropout
##.  - participant drops at t_ij with a prob determined by
##.    a logistic model with indicator Y_ij as a
##.    predictor
############################################################
############################################################
for(count in 1:num_datasets){
  n <- 200 # number of subjects
  n_test <- 200
  K <- 10 # number of measurements per subject
  t_max <- 10 # maximum follow-up time
  min_sep <- 0.5
  
  # we construct a data frame with the design:
  # everyone has a baseline measurement, and then measurements at random
  # follow-up times up to t_max
  # !! for Machine Learning purposes: two consecutive measurements will have
  #more than 0.5 of
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
  
  
  #Simulate second longitudinal outcome
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
  
  #Simulate third longitudinal outcome
  ###############################################
  betas5 <- c(1, 0.155, 0.12345) # fixed effects coefficients
  betas5 <- c(70, -5.44, -3.761289) # fixed effects coefficients
  
  # we simulate random effects
  b5 <- b[, c(5,6)]
  # linear predictor
  eta_y5 <- as.vector(X5 %*% betas5 + rowSums(Z5 * b5[DF$id, ]))
  # mean of the binomial distribution
  #mu_y5 <- plogis(eta_y5)
  # we simulate binomial longitudinal data
  #DF$y5 <- rbinom(n * K, size = 1, prob = mu_y5)
  DF$y5 <- rnorm(n_test * K, mean = eta_y5, sd = sigma3)
  
  ####Test data
  # we simulate random effects
  b5_test <- b_test[, c(5,6)]
  # linear predictor
  eta_y5_test <- as.vector(X5_test %*% betas5 + rowSums(Z5_test * b5_test[DF_test$id, ]))
  # mean of the binomial distribution
  #mu_y5_test <- plogis(eta_y5_test)
  # we simulate binomial longitudinal data
  #DF_test$y5 <- rbinom(n_test * K, size = 1, prob = mu_y5_test)
  
  DF_test$y5 <- rnorm(n_test * K, mean = eta_y5_test, sd = sigma3)
  
  #If using admin censoring:
  #AdminCens <- 10
  
  # we keep the longitudinal measurements before the event times: TRAIN
  #DF <- DF[DF$time <= AdminCens, ]
  
  # we keep the longitudinal measurements before the event times: TEST
  #DF_test <- DF_test[DF_test$time <= AdminCens, ]
  
  #Data set with just one obs per individual
  #DF.id <- DF[!duplicated(DF$id),]
  #DF_test.id <- DF_test[!duplicated(DF_test$id),]
  
  param1 <- -3
  param2 <- 0.01
  probs_vec1 <- numeric(n*K)
  probs_vec2 <- numeric(n*K)
  probs_vec3 <- numeric(n*K)
  
  kk <- 1
  for(i in 1:n){
    for(j in 1:K){
      if((kk-1)%%K==0){
        probs_vec1[kk] <- 0
        probs_vec2[kk] <- 0
        probs_vec3[kk] <- 0
      } else{
        thet1 <- param1 + param2*DF[DF$id==i,]$y1[j]
        probs_vec1[kk] <- 1/(1+exp(-thet1))
        thet2 <- param1 + param2*DF[DF$id==i,]$y3[j]
        probs_vec2[kk] <- 1/(1+exp(-thet2))
        thet3 <- param1 + param2*DF[DF$id==i,]$y5[j]
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
      # checking dropout
      dropout_point = which.max(dropout_y1 == 1 | dropout_y3 == 1 | dropout_y5 == 1),
      
      # NA after dropout
      y1 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y1),
      y3 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y3),
      y5 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y5)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point)
  
  DF_test_miss <- DF_test %>%
    group_by(id) %>%                     # Group by each individual
    mutate(
      # checking dropout
      dropout_point = which.max(dropout_y1 == 1 | dropout_y3 == 1 | dropout_y5 == 1),
      
      # NA after dropout
      y1 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y1),
      y3 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y3),
      y5 = ifelse(row_number() >= dropout_point & dropout_point > 1, NA, y5)
    ) %>%
    ungroup() %>%
    dplyr::select(-dropout_point)
  
  
  missing_ratio_train[count] <- counting_missing_data(DF_miss, n, K)
  missing_ratio_test[count] <- counting_missing_data(DF_test_miss, n, K)
  
  try(list_DF <- append(list_DF, list(DF)))
  try(list_DF_miss <- append(list_DF_miss, list(DF_miss)))
  try(list_DF_test <- append(list_DF_test, list(DF_test)))
  try(list_DF_test_miss <- append(list_DF_test_miss, list(DF_test_miss)))
  
  try(if(count==100){
    setwd("D:/La meva unitat/TFM/ResultsMMvsMFPCA")
    str100 <- "dataframes_increasing_MNAR_50.RData"
    save(list_DF, list_DF_miss, list_DF_test, list_DF_test_miss,
         missing_ratio_train, missing_ratio_test, file=str100)
  })
  
  print(count)
}


###################################################
####################################################
## Fitting MM and MFPCA models and doing predictions
####################################################
####################################################

for(count in 1:num_datasets){
  
  ###OBS!!!!!
  ## DF_miss contiene el data frame con missing data!!
  
  ######################################
  #Fitting multivariate MM models
  # - We use snat_mvmer() from rstanarm
  # - We fit the true model plus two
  #.  (or 3) more models
  ######################################
  
  #checking time of computation using data frame with missings 
  # when n=300, k=10, we have 21mins, max rhat 1.67
  #      n=250, k=10: 17mins, max rhat 1.51
  #      n=200, k=10:  8.6 mins, max rhat 1.09
  try(true_model_miss <- stan_mvmer(
    formula = list(
      y1 ~ sex * time + (time | id),
      y3 ~ sqrt(time) + (sqrt(time) | id),
      y5 ~ treatment + time + (time | id)),
    data = list_DF_miss[[count]],
    family = list(gaussian, gaussian, gaussian),
    chains = 3, cores = 8, seed = 12345, iter = 1750))
  
  try(rhat_true_model_miss <- rhat(true_model_miss))
  
  #OJO!!! para que devuelva las predictions a pesar de estar pasandole NAs 
  #hay que pasarle el data frame sin los NAs
  
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
                           uniExpansions = list(list(type="uFPCA", pve = 0.65),
                                                list(type ="uFPCA", pve=0.65),
                                                list(type="uFPCA", pve=0.65)),
                           fit = TRUE))
  
  try(pred_mfpca1_test <- predict(mfpca1, scores = mfpca1_test$scores))
  
  try(pred_mfpca1_y1_test <- NULL)
  try(pred_mfpca1_y1_test <- pred_mfpca1_test[[1]]@X)
  try(pred_mfpca1_y2_test <- NULL)
  try(pred_mfpca1_y2_test <- pred_mfpca1_test[[2]]@X)
  try(pred_mfpca1_y3_test <- NULL)
  try(pred_mfpca1_y3_test <- pred_mfpca1_test[[3]]@X)
  
  #plot(pred_mfpca1_test)
  
  #AQUI YA TENGO TODA LA INFO QUE QUIERO! las predicciones de esto último se 
  #hacen para cada uno de los puntos de la grid!! También para aquellos donde
  #hay NA! Por tanto para aquellos puntos estoy haciendo una predicción!
  #Y es esa predicción la que me tengo que quedar y entrar a comparar con 
  #la predicción que me puedan devolver los MM
  
  #Cuando tenemos el predict(MFPCAfit), el objeto que pred_mfpca1[[1]]@X es
  # una matriz n X max_length +1, ie cada fila es un individuo y cada columna
  # es uno de los puntos de la fixed time grid. Y esto para el outcome 
  # longitudinal 1 
  # Lo que tenemos aquí es una predicción para cada punto de la fixed grid
  # del outcome longitudinal 1.
  
  
  
  #####################################################
  ## Comparing predictions by means of standardized
  ## RMSE
  #####################################################
  
  
  
  
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
    strr <- "results_1_MNAR50_incr_03mar2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr)
  })
  try(if(count==2){
    strr2 <- "results_2_MNAR50_incr_03mar2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=strr2)
  })
  try(if(count==5){
    str1 <- "results_5_MNAR50_incr_03mar2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str1)
  })
  try(if(count==10){
    str10 <- "results_10_MNAR50_incr_03mar2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str10)
  })
  try(if(count==20){
    str20 <- "results_20_MNAR50_incr_03mar2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str20)
  })
  try(if(count==30){
    str30 <- "results_30_MNAR50_incr_03mar2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str30)
  })
  try(if(count==40){
    str40 <- "results_40_MNAR50_incr_03mar2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss,file=str40)
  })
  try(if(count==50){
    str50 <- "results_50_MNAR50_incr_03mar2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str50)
  })
  try(if(count==75){
    str75 <- "results_75_MNAR50_incr_03mar2025.RData"
    save(list_pred_mfpca1_y1, list_pred_mfpca1_y2, list_pred_mfpca1_y3,
         list_pred_mfpca1_y1_test, list_pred_mfpca1_y2_test, list_pred_mfpca1_y3_test,
         list_preds_tmod_y1_miss, list_preds_tmod_y1_miss2, list_preds_tmod_y1_miss_test,
         list_preds_tmod_y2_miss, list_preds_tmod_y2_miss2, list_preds_tmod_y2_miss_test,
         list_preds_tmod_y3_miss, list_preds_tmod_y3_miss2, list_preds_tmod_y3_miss_test,
         list_rhat_true_model_miss, file=str75)
  })
  try(if(count==100){
    str100 <- "results_100_MNAR50_incr_03mar2025.RData"
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

