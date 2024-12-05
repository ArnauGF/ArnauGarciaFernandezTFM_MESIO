##############################################################################
##############################################################################
## Analysis of the dataset used with ML models using SL
## See how the data sets are generated in simul_sce2.R file
#############################################################################
##############################################################################
library(JMbayes2)

num_dat <- 100
IBS_multi  <- IBS_multi_test <- numeric(num_dat)
dSL_cv_IBS <- eSL_cv_IBS <- dSL_test_IBS <- eSL_test_IBS <- numeric(num_dat)
IBS_univ  <- IBS_w <- matrix(nrow = num_dat, ncol = 5)
IBS_univ_test <- IBS_w_test <- matrix(nrow = num_dat, ncol = 5)
censoring_train <- censoring_test <- numeric(num_dat)
disSL_ibs <- numeric(num_dat)
count <- 0

###Read the data frames:
load("SceII_L5_LLMs_100df.RData")

### we save it in a drive folder:
setwd("D:/La meva unitat/TFM/code")

for(count in 1:num_dat){
  #we print the current iteration in order to see the progress
  print(count)
  #We take the current data sets
  df_name <- paste0("DF_", count)
  #df.id_name <- paste0("DF.id_", count)
  df_test_name <- paste0("DF_test_", count)
  #df_test.id_name <- paste0("DF_test.id_", count)
  
  DF <- get(df_name)
  DF.id <- DF[!duplicated(DF$id),]
  DF_test <- get(df_test_name)
  DF_test.id <- DF_test[!duplicated(DF_test$id),]
  
  
  try(CoxFit <- coxph(Surv(Time, event) ~ sex, data = DF.id))
  
  try(LM1 <- lme(y1 ~ sex*time, data = DF, random = ~ time | id))
  try(LM2 <- lme(y2 ~ treatment*time, data = DF, random = ~ time | id))
  try(LM3 <- lme(y3 ~ time, data = DF, random = ~ time | id))
  #try(LM4 <- mixed_model(y4 ~ sex + time, data = DF, random = ~ time || id, family = binomial()))
  try(LM4 <- lme(y4 ~ sex + time, data = DF, random = ~ time | id))
  #try(LM5 <- mixed_model(y5 ~ time, data = DF, random = ~ time || id, family = binomial()))
  try(LM5 <- lme(y5 ~ treatment + time, data = DF, random = ~ time | id))
  
  try(multiJM <- jm(CoxFit, list(LM1, LM2, LM3, LM4, LM5), time_var = "time",
                    n_iter = 12000L, n_burnin = 2000L, n_thin = 5L))
  
  t0 <- 2.5
  dt <- 1
  try(brier_score_multi_IPCW <- tvBrier(multiJM, newdata = DF, Tstart = t0, Dt = dt, 
                                        integrated = TRUE, type_weights = "IPCW"))
  try(brier_score_multi_test_IPCW <- tvBrier(multiJM, newdata = DF_test, Tstart = t0, Dt = dt, 
                                             integrated = TRUE, type_weights = "IPCW"))
  
  try(CVdats <- create_folds(DF, V = 3, id_var = "id"))
  
  fit_models <- function (data) {
    library("JMbayes2")
    data_id <- data[!duplicated(data$id), ]
    CoxFit <- coxph(Surv(Time, event) ~ sex, data = data_id)
    LM1 <- lme(y1 ~ sex*time, data = data, random = ~ time | id)
    LM2 <- lme(y2 ~ treatment*time, data = data, random = ~ time | id)
    LM3 <- lme(y3 ~ time, data = data, random = ~ time | id)
    LM4 <- lme(y4 ~ sex + time, data = data, random = ~ time | id)
    LM5 <- lme(y5 ~ treatment + time, data = data, random = ~ time | id)
    JM1 <- jm(CoxFit, LM1, time_var = "time")
    JM2 <- jm(CoxFit, LM2, time_var = "time")
    JM3 <- jm(CoxFit, LM3, time_var = "time")
    JM4 <- jm(CoxFit, LM4, time_var = "time")
    JM5 <- jm(CoxFit, LM5, time_var = "time")
    out <- list(M1 = JM1, M2 = JM2, M3 = JM3, M4 = JM4, M5 = JM5)
    class(out) <- "jmList"
    out
  }
  
  try(cl <- parallel::makeCluster(5L))
  try(Models_folds <- parallel::parLapply(cl, CVdats$training, fit_models))
  try(parallel::stopCluster(cl))
  
  try(Brier_sl_IPCW <- tvBrier(Models_folds, newdata = CVdats$testing, 
                               integrated = TRUE, Tstart = t0, Dt = dt,
                               type_weights = "IPCW"))
  
  try(Models <- fit_models(DF))
  cat("Full models fitted\n")
  
  try(bw <- Brier_sl_IPCW$weights)
  try(Brier_SL_IPCW_test <- tvBrier(Models, newdata = DF_test, model_weights = bw, 
                                    Tstart = t0, Dt = dt, integrated = TRUE,
                                    type_weights = "IPCW"))
  
  disSL_ibs[count] <- 0
  try(disSL_ibs[count] <- which.min(Brier_sl_IPCW$Brier_per_model))
  
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
  } else if(disSL_ibs[count] == 5){
    try(Brier_dSL_test <- tvBrier(Models$M5, newdata = DF_test, 
                                  Tstart = t0, Dt = dt, integrated = TRUE,
                                  type_weights = "IPCW"))
  }
  
  ################################################################
  ##### Saving results:
  ################################################################
  try(IBS_multi[count] <- brier_score_multi_IPCW$Brier)
  
  try(IBS_univ[count, ] <- Brier_sl_IPCW$Brier_per_model)
  
  try(IBS_w[count, ] <- Brier_sl_IPCW$weights)
  
  try(dSL_cv_IBS[count] <- min(Brier_sl_IPCW$Brier_per_model))
  
  try(eSL_cv_IBS[count] <- Brier_sl_IPCW$Brier)
  
  try(eSL_test_IBS[count] <- Brier_SL_IPCW_test$Brier)
  
  try(IBS_multi_test[count] <- brier_score_multi_test_IPCW$Brier)
  
  try(dSL_test_IBS[count] <- Brier_dSL_test$Brier)
  
  
  
  try(censoring_train[count] <- sum(DF.id$event==0))
  try(censoring_test[count] <- sum(DF_test.id$event==0))
  
  
  if(count==2){
    strr <- "results2_sce2_2dec_5llms.RData"
    save(IBS_multi, IBS_univ, IBS_w, dSL_cv_IBS,
         eSL_cv_IBS, eSL_test_IBS, 
         IBS_multi_test, censoring_train,
         censoring_test, dSL_test_IBS, disSL_ibs, file=strr)
  }
  if(count==5){
    str1 <- "results5_sce2_2dec_5llms.RData"
    save(IBS_multi, IBS_univ, IBS_w, dSL_cv_IBS,
         eSL_cv_IBS, eSL_test_IBS, 
         IBS_multi_test, censoring_train,
         censoring_test, dSL_test_IBS, disSL_ibs, file=str1)
  }
  try(if(count==10){
    str10 <- "results10_sce2_2dec_5llms.RData"
    save(IBS_multi, IBS_univ, IBS_w, dSL_cv_IBS,
         eSL_cv_IBS, eSL_test_IBS, 
         IBS_multi_test, censoring_train,
         censoring_test, dSL_test_IBS, disSL_ibs, file=str10)
  })
  try(if(count==20){
    str20 <- "results20_sce2_2dec_5llms.RData"
    save(IBS_multi, IBS_univ, IBS_w, dSL_cv_IBS,
         eSL_cv_IBS, eSL_test_IBS, 
         IBS_multi_test, censoring_train,
         censoring_test, dSL_test_IBS, disSL_ibs, file=str20)
  })
  try(if(count==30){
    str30 <- "results30_sce2_2dec_5llms.RData"
    save(IBS_multi, IBS_univ, IBS_w, dSL_cv_IBS,
         eSL_cv_IBS, eSL_test_IBS, 
         IBS_multi_test, censoring_train,
         censoring_test, dSL_test_IBS, disSL_ibs, file=str30)
  })
  try(if(count==40){
    str40 <- "results40_sce2_2dec_5llms.RData"
    save(IBS_multi, IBS_univ, IBS_w, dSL_cv_IBS,
         eSL_cv_IBS, eSL_test_IBS, 
         IBS_multi_test, censoring_train,
         censoring_test, dSL_test_IBS, disSL_ibs, file=str30)
  })
  try(if(count==50){
    str50 <- "results50_sce2_2dec_5llms.RData"
    save(IBS_multi, IBS_univ, IBS_w, dSL_cv_IBS,
         eSL_cv_IBS, eSL_test_IBS, 
         IBS_multi_test, censoring_train,
         censoring_test, dSL_test_IBS, disSL_ibs, file=str50)
  })
  try(if(count==75){
    str75 <- "results75_sce2_2dec_5llms.RData"
    save(IBS_multi, IBS_univ, IBS_w, dSL_cv_IBS,
         eSL_cv_IBS, eSL_test_IBS, 
         IBS_multi_test, censoring_train,
         censoring_test, dSL_test_IBS, disSL_ibs, file=str75)
  })
  try(if(count==100){
    str100 <- "results100_sce2_2dec_5llms.RData"
    save(IBS_multi, IBS_univ, IBS_w, dSL_cv_IBS,
         eSL_cv_IBS, eSL_test_IBS, 
         IBS_multi_test, censoring_train,
         censoring_test, dSL_test_IBS, disSL_ibs, file=str100)
  })
}
