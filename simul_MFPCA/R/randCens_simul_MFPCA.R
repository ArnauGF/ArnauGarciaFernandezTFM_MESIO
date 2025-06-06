################################################################################
### Random censoring 30%
################################################################################

library(JMbayes2)
library(MASS)
library(lattice)
library(dplyr)
library(tidyr)

pow2 <- function(x){
  x^2
}

## Loading random cens 30% datasets
#load("D:/La meva unitat/TFM/MFPCA_simuls/R/datasets/datasets_randCens_MFPCA_29apr2025.RData")
load("G:/TFM/SLinJointModels/MFPCA_simuls/R/datasets_randCens_MFPCA_29apr2025.RData")

repl <- 100
list_rhats_MJM <- list()
perc_cens_train <- numeric(repl)
perc_cens_test <- numeric(repl)
list_full_rhats_MJM <- list()

ibs_train <- numeric(repl)
ibs_test <- numeric(repl)
n_risk_train <- numeric(repl)
n_risk_test <- numeric(repl)
n_cens_train <- numeric(repl)
n_cens_test <- numeric(repl)
n_event_train <- numeric(repl)
n_event_test <- numeric(repl)
list_brier_model_train <- list()
list_w_model_train <- list()
list_brier_model_test <- list()
list_w_model_test <- list()
dSL_cv_IBS <- eSL_cv_IBS <- dSL_test_IBS <- eSL_test_IBS <- numeric(repl)
disSL_ibs <- numeric(repl)
IBS_univ <- IBS_w  <- matrix(nrow = repl, ncol = 4)
IBS_univ_test  <- IBS_w_test <- matrix(nrow = repl, ncol = 4)


##metrics for the second interval:
ibs_train_2 <- numeric(repl)
ibs_test_2 <- numeric(repl)
n_risk_train_2 <- numeric(repl)
n_risk_test_2 <- numeric(repl)
n_cens_train_2 <- numeric(repl)
n_cens_test_2 <- numeric(repl)
n_event_train_2 <- numeric(repl)
n_event_test_2 <- numeric(repl)
list_brier_model_train_2 <- list()
list_w_model_train_2 <- list()
list_brier_model_test_2 <- list()
list_w_model_test_2 <- list()
dSL_cv_IBS_2 <- eSL_cv_IBS_2 <- dSL_test_IBS_2 <- eSL_test_IBS_2 <- numeric(repl)
disSL_ibs_2 <- numeric(repl)
IBS_univ_2 <- IBS_w_2  <- matrix(nrow = repl, ncol = 4)
IBS_univ_test_2  <- IBS_w_test_2 <- matrix(nrow = repl, ncol = 4)


for(count in 1:repl){
  n <- 175 # number of subjects
  n_test <- 175
  K <- 10 # number of measurements per subject
  
  # we set DF ans DF_test to the corresponding element of the list
  DF <- list_DF_randCens[[count]]
  DF_test <- list_DF_test_randCens[[count]]
  
  DF.id <- DF[!duplicated(DF$id),]
  DF_test.id <- DF_test[!duplicated(DF_test$id),] 
  
  ## Checking % of censoring
  perc_cens_train[count] <- sum(DF.id$event==0)/n
  perc_cens_test[count] <- sum(DF_test.id$event==0)/n_test
  
  ## Survival model
  try(CoxFit <- coxph(Surv(Time, event) ~ sex, data = DF.id))
  
  ## Longitudinal models
  try(LM1 <- lme(y1 ~ sex + time, data = DF, random = ~ time | id))
  try(LM2 <- lme(y2 ~ time, data = DF, random = ~ time | id))
  try(LM3 <- lme(y3 ~ treat + pow2(time), data = DF, random = ~ time | id))
  try(LM4 <- lme(y4 ~ sex + treat + time, data = DF, random = ~ time | id))
  
  
  #Fitting the multivariate JM
  try(multiJM <- jm(CoxFit, list(LM1, LM2, LM3, LM4), time_var = "time",
                    n_iter = 13000L, n_burnin = 300L, n_thin = 5L))
  
  
  try(rhats <- numeric())
  try(for(i in 1:(length(multiJM$statistics$Rhat)-2)){
    rhats <- c(rhats, multiJM$statistics$Rhat[[i]][,1]) 
  })
  try(rhats <- as.numeric(rhats))
  
  #we save convergence metrics of the mJM
  try(list_rhats_MJM <- append(list_rhats_MJM, list(rhats)))
  try(list_full_rhats_MJM <- append(list_full_rhats_MJM, 
                                    list(multiJM$statistics$Rhat)))
  
  
  # It seems that those bs_gammas (the ones related with B-splines used to
  # fit the besaline hazard) are too big (more than 2). Could be that we are
  # assuming ctt baseline hazard function?
  
  ## First interval:
  t0 <- 6
  dt <- 1.5
  #Second interval:
  t0_2 <- 4
  dt_2 <- 1.5
  
  ## IBS im (6,7.5]:
  try(brier_score_multi_train <- tvBrier(multiJM, newdata = DF, Tstart = t0, Dt = dt, 
                                         integrated = TRUE, type_weights = "IPCW"))
  try(ibs_train[count] <- brier_score_multi_train$Brier)
  try(n_risk_train[count] <- brier_score_multi_train$nr)
  try(n_event_train[count] <- brier_score_multi_train$nint)
  try(n_cens_train[count] <- brier_score_multi_train$ncens)
  
  ### IBS in (4, 5.5]:
  try(brier_score_multi_train_2 <- tvBrier(multiJM, newdata = DF, Tstart = t0_2, Dt = dt_2, 
                                           integrated = TRUE, type_weights = "IPCW"))
  try(ibs_train_2[count] <- brier_score_multi_train_2$Brier)
  try(n_risk_train_2[count] <- brier_score_multi_train_2$nr)
  try(n_event_train_2[count] <- brier_score_multi_train_2$nint)
  try(n_cens_train_2[count] <- brier_score_multi_train_2$ncens)
  
  ## IBS in test in (6,7.5]:
  try(brier_score_multi_test <- tvBrier(multiJM, newdata = DF_test, Tstart = t0, Dt = dt, 
                                        integrated = TRUE, type_weights = "IPCW"))
  try(ibs_test[count] <- brier_score_multi_test$Brier)
  try(n_risk_test[count] <- brier_score_multi_test$nr)
  try(n_event_test[count] <- brier_score_multi_test$nint)
  try(n_cens_test[count] <- brier_score_multi_test$ncens)
  
  ## IBS in test in (4,5.5]:
  try(brier_score_multi_test_2 <- tvBrier(multiJM, newdata = DF_test, Tstart = t0_2, Dt = dt_2, 
                                          integrated = TRUE, type_weights = "IPCW"))
  try(ibs_test_2[count] <- brier_score_multi_test_2$Brier)
  try(n_risk_test_2[count] <- brier_score_multi_test_2$nr)
  try(n_event_test_2[count] <- brier_score_multi_test_2$nint)
  try(n_cens_test_2[count] <- brier_score_multi_test_2$ncens)
  
  #SuperLearning with the library of models built with the univariate JM
  
  try(CVdats <- create_folds(DF, V = 3, id_var = "id"))
  
  fit_models <- function (data) {
    library("JMbayes2")
    pow2 <- function(x){
      x^2
    }
    data_id <- data[!duplicated(data$id), ]
    CoxFit <- coxph(Surv(Time, event) ~ sex, data = data_id)
    LM1 <- lme(y1 ~ sex + time, data = data, random = ~ time | id)
    LM2 <- lme(y2 ~ time, data = data, random = ~ time | id)
    LM3 <- lme(y3 ~ treat + pow2(time), data = data, random = ~ time | id)
    LM4 <- lme(y4 ~ sex + treat + time, data = data, random = ~ time | id)
    JM1 <- jm(CoxFit, LM1, time_var = "time")
    JM2 <- jm(CoxFit, LM2, time_var = "time")
    JM3 <- jm(CoxFit, LM3, time_var = "time")
    JM4 <- jm(CoxFit, LM4, time_var = "time")
    out <- list(M1 = JM1, M2 = JM2, M3 = JM3, M4 = JM4)
    class(out) <- "jmList"
    out
  }

  try(cl <- parallel::makeCluster(4L))
  try(Models_folds <- parallel::parLapply(cl, CVdats$training, fit_models))
  try(parallel::stopCluster(cl))
  
  
  ## IBS for CV data in (6,7.5]:
  try(Brier_weights <- tvBrier(Models_folds, newdata = CVdats$testing, 
                               integrated = TRUE, Tstart = t0, Dt = dt,
                               type_weights = "IPCW"))
  
  ## IBS for CV data in (4,5.5]:
  try(Brier_weights_2 <- tvBrier(Models_folds, newdata = CVdats$testing, 
                                 integrated = TRUE, Tstart = t0_2, Dt = dt_2,
                                 type_weights = "IPCW"))
  
  #Now with testing data
  #We fit the models in the whole training data set and test in testing data
  try(Models <- fit_models(DF))
  
  ## IBS in (6,7.5]
  try(bw <- Brier_weights$weights)
  try(Brier_weights_test <- tvBrier(Models, newdata = DF_test, model_weights = bw, 
                                    Tstart = t0, Dt = dt, integrated = TRUE,
                                    type_weights = "IPCW"))
  
  ## IBS in (4,5.5]
  try(bw_2 <- Brier_weights_2$weights)
  try(Brier_weights_test_2 <- tvBrier(Models, newdata = DF_test, model_weights = bw_2, 
                                      Tstart = t0_2, Dt = dt_2, integrated = TRUE,
                                      type_weights = "IPCW"))
  
  ## IBS in (6,7.5]:
  disSL_ibs[count] <- 0
  try(disSL_ibs[count] <- which.min(Brier_weights$Brier_per_model))
  if(disSL_ibs[count] == 1){
    try(Brier_dSL_test <- tvBrier(Models$M1, newdata = DF_test, 
                                  Tstart = t0, Dt = dt, integrated = TRUE,
                                  type_weights = "IPCW"))
  } else if(disSL_ibs[count] == 2){
    try(Brier_dSL_test <- tvBrier(Models$M2, newdata = DF_test, 
                                  Tstart = t0, Dt = dt, integrated = TRUE,
                                  type_weights = "IPCW"))
  } else if(disSL_ibs[count] == 3){
    try(Brier_dSL_test <- tvBrier(Models$M3, newdata = DF_test, 
                                  Tstart = t0, Dt = dt, integrated = TRUE,
                                  type_weights = "IPCW"))
  } else if(disSL_ibs[count] == 4){
    try(Brier_dSL_test <- tvBrier(Models$M4, newdata = DF_test, 
                                  Tstart = t0, Dt = dt, integrated = TRUE,
                                  type_weights = "IPCW"))
  } 
  
  ## IBS in (4, 5.5]:
  disSL_ibs_2[count] <- 0
  try(disSL_ibs_2[count] <- which.min(Brier_weights_2$Brier_per_model))
  if(disSL_ibs_2[count] == 1){
    try(Brier_dSL_test_2 <- tvBrier(Models$M1, newdata = DF_test, 
                                    Tstart = t0_2, Dt = dt_2, integrated = TRUE,
                                    type_weights = "IPCW"))
  } else if(disSL_ibs_2[count] == 2){
    try(Brier_dSL_test_2 <- tvBrier(Models$M2, newdata = DF_test, 
                                    Tstart = t0_2, Dt = dt_2, integrated = TRUE,
                                    type_weights = "IPCW"))
  } else if(disSL_ibs_2[count] == 3){
    try(Brier_dSL_test_2 <- tvBrier(Models$M3, newdata = DF_test, 
                                    Tstart = t0_2, Dt = dt_2, integrated = TRUE,
                                    type_weights = "IPCW"))
  } else if(disSL_ibs_2[count] == 4){
    try(Brier_dSL_test_2 <- tvBrier(Models$M4, newdata = DF_test, 
                                    Tstart = t0_2, Dt = dt_2, integrated = TRUE,
                                    type_weights = "IPCW"))
  } 
  
  ########################
  #Save the desired metrics
  ###########################
  
  ## saving metrics for the interval (6,7.5]:
  try(IBS_univ[count, ] <- Brier_weights$Brier_per_model)

  try(IBS_w[count, ] <- Brier_weights$weights)

  try(dSL_cv_IBS[count] <- min(Brier_weights$Brier_per_model))

  try(eSL_cv_IBS[count] <- Brier_weights$Brier)

  try(eSL_test_IBS[count] <- Brier_weights_test$Brier)

  try(dSL_test_IBS[count] <- Brier_dSL_test$Brier)

  
  ## saving metrics for the interval (4,5.5]:
  try(IBS_univ_2[count, ] <- Brier_weights_2$Brier_per_model)

  try(IBS_w_2[count, ] <- Brier_weights_2$weights)

  try(dSL_cv_IBS_2[count] <- min(Brier_weights_2$Brier_per_model))

  try(eSL_cv_IBS_2[count] <- Brier_weights_2$Brier)

  try(eSL_test_IBS_2[count] <- Brier_weights_test_2$Brier)

  try(dSL_test_IBS_2[count] <- Brier_dSL_test_2$Brier)

  
  
  try(if(count==10){
    str10 <- "G:/TFM/SLinJointModels/MFPCA_simuls/results/repl10_randCens_30_2apr2025.RData"
    save(perc_cens_test, perc_cens_train,
         list_rhats_MJM, list_full_rhats_MJM, 
         n_risk_train, n_risk_test, n_event_train, n_event_test,
         n_cens_train, n_cens_test,
         ibs_train, ibs_test, 
         IBS_univ, IBS_w, dSL_cv_IBS, eSL_test_IBS, eSL_cv_IBS,
         dSL_test_IBS, disSL_ibs, 
         n_risk_train_2, n_risk_test_2, n_event_train_2, n_event_test_2,
         n_cens_train_2, n_cens_test_2,
         ibs_train_2, ibs_test_2, 
         IBS_univ_2, IBS_w_2, dSL_cv_IBS_2, eSL_test_IBS_2, eSL_cv_IBS_2,
         dSL_test_IBS_2, disSL_ibs_2, 
         file=str10)
  })
  #try(if(count==25){
  #  str25 <- "D:/La meva unitat/TFM/MFPCA_simuls/results/repl25_randCens_30_2apr2025.RData"
  #  save(checkTimes_test, checkTimes, perc_cens_test, perc_cens_train,
  #       list_rhats_MJM, list_full_rhats_MJM, 
  #       n_risk_train, n_risk_test, n_event_train, n_event_test,
  #       n_cens_train, n_cens_test,
  #       ibs_train, ibs_test, 
  #       IBS_univ, IBS_w, dSL_cv_IBS, eSL_test_IBS, eSL_cv_IBS,
  #       dSL_test_IBS, disSL_ibs, 
  #       epce_train, epce_test,
  #       EPCE_univ, EPCE_w, dSL_cv_EPCE, eSL_cv_EPCE, dSL_test_EPCE, eSL_test_EPCE,
  #       disSL_epce,
  #       n_risk_train_2, n_risk_test_2, n_event_train_2, n_event_test_2,
  #       n_cens_train_2, n_cens_test_2,
  #       ibs_train_2, ibs_test_2, 
  #       IBS_univ_2, IBS_w_2, dSL_cv_IBS_2, eSL_test_IBS_2, eSL_cv_IBS_2,
  #       dSL_test_IBS_2, disSL_ibs_2, 
  #       epce_train_2, epce_test_2,
  #       EPCE_univ_2, EPCE_w_2, dSL_cv_EPCE_2, eSL_cv_EPCE_2, dSL_test_EPCE_2, eSL_test_EPCE_2,
  #       disSL_epce_2,
  #       file=str25)
  #})
  try(if(count==50){
    str50 <- "G:/TFM/SLinJointModels/MFPCA_simuls/results/repl50_randCens_30_2apr2025.RData"
    save(perc_cens_test, perc_cens_train,
         list_rhats_MJM, list_full_rhats_MJM, 
         n_risk_train, n_risk_test, n_event_train, n_event_test,
         n_cens_train, n_cens_test,
         ibs_train, ibs_test, 
         IBS_univ, IBS_w, dSL_cv_IBS, eSL_test_IBS, eSL_cv_IBS,
         dSL_test_IBS, disSL_ibs, 
         n_risk_train_2, n_risk_test_2, n_event_train_2, n_event_test_2,
         n_cens_train_2, n_cens_test_2,
         ibs_train_2, ibs_test_2, 
         IBS_univ_2, IBS_w_2, dSL_cv_IBS_2, eSL_test_IBS_2, eSL_cv_IBS_2,
         dSL_test_IBS_2, disSL_ibs_2, 
         file=str50)
  })
  try(if(count==75){
    str75 <- "G:/TFM/SLinJointModels/MFPCA_simuls/results/repl75_randCens_30_2apr2025.RData"
    save(perc_cens_test, perc_cens_train,
         list_rhats_MJM, list_full_rhats_MJM, 
         n_risk_train, n_risk_test, n_event_train, n_event_test,
         n_cens_train, n_cens_test,
         ibs_train, ibs_test, 
         IBS_univ, IBS_w, dSL_cv_IBS, eSL_test_IBS, eSL_cv_IBS,
         dSL_test_IBS, disSL_ibs, 
         n_risk_train_2, n_risk_test_2, n_event_train_2, n_event_test_2,
         n_cens_train_2, n_cens_test_2,
         ibs_train_2, ibs_test_2, 
         IBS_univ_2, IBS_w_2, dSL_cv_IBS_2, eSL_test_IBS_2, eSL_cv_IBS_2,
         dSL_test_IBS_2, disSL_ibs_2, 
         file=str75)
  })
  
  if(count == repl){
    strr <- "G:/TFM/SLinJointModels/MFPCA_simuls/results/repl100_randCens_30_2apr2025.RData"
    save(perc_cens_test, perc_cens_train,
         list_rhats_MJM, list_full_rhats_MJM, 
         n_risk_train, n_risk_test, n_event_train, n_event_test,
         n_cens_train, n_cens_test,
         ibs_train, ibs_test, 
         IBS_univ, IBS_w, dSL_cv_IBS, eSL_test_IBS, eSL_cv_IBS,
         dSL_test_IBS, disSL_ibs, 
         n_risk_train_2, n_risk_test_2, n_event_train_2, n_event_test_2,
         n_cens_train_2, n_cens_test_2,
         ibs_train_2, ibs_test_2, 
         IBS_univ_2, IBS_w_2, dSL_cv_IBS_2, eSL_test_IBS_2, eSL_cv_IBS_2,
         dSL_test_IBS_2, disSL_ibs_2, 
         file=strr)
  }
  
  print(count)
}

