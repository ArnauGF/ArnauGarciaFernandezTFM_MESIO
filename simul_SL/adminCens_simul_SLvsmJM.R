####################################################################
####################################################################
## We do several replications
## With L=4, the last outcome being a GLMM (binary longitudinal outcome)
## type 1 cens
####################################################################
####################################################################
repl <- 100
list_rhats_MJM <- list()
checkTimes <- numeric(repl)
checkTimes_test <- numeric(repl)
perc_cens_train <- numeric(repl)
perc_cens_test <- numeric(repl)
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
list_full_rhats_MJM <- list()
IBS_multi <- IBS_multi_test  <- numeric(repl)
dSL_cv_IBS <- eSL_cv_IBS <- dSL_test_IBS <- eSL_test_IBS <- numeric(repl)
disSL_ibs <- numeric(repl)
IBS_univ <- IBS_w  <- matrix(nrow = repl, ncol = 4)
IBS_univ_test  <- IBS_w_test <- matrix(nrow = repl, ncol = 4)


epce_train <- epce_test <- numeric(repl)
dSL_cv_EPCE <- eSL_cv_EPCE <- dSL_test_EPCE <- eSL_test_EPCE <- numeric(repl)
EPCE_univ <- EPCE_w <- matrix(nrow = repl, ncol = 4)
EPCE_univ_test <- EPCE_w_test <- matrix(nrow = repl, ncol = 4)
disSL_epce <- numeric(repl)


for(count in 1:repl){
  n <- 250 # number of subjects
  n_test <- 250
  K <- 10 # number of measurements per subject
  
  # we set DF ans DF_test to the corresponding element of the list
  DF <- list_DF[[count]]
  DF_test <- list_DF_test[[count]]
  
  ##saving true times
  trueTimes <- DF[!duplicated(DF$id),]$Time
  checkTimes[count] <- sum(trueTimes==150)/n
  
  ##saving true times
  trueTimes_test <- DF_test[!duplicated(DF_test$id),]$Time
  checkTimes_test[count] <- sum(trueTimes_test==150)/n_test
  
  #simulating censoring; Type I
  Ctimes <- 7.5
  Time <- pmin(trueTimes, Ctimes)
  event <- as.numeric(trueTimes <= Ctimes) # event indicator
  # we keep the longitudinal measurements before the event times
  DF$Time <- Time[DF$id]
  DF$event <- event[DF$id]
  DF <- DF[DF$time <= DF$Time, ]
  
  ## Checking % of censoring
  perc_cens_train[count] <- sum(event==0)/n
  
  ## Testing data:
  Time <- pmin(trueTimes_test, Ctimes)
  event <- as.numeric(trueTimes_test <= Ctimes) # event indicator
  # we keep the longitudinal measurements before the event times
  DF_test$Time <- Time[DF_test$id]
  DF_test$event <- event[DF_test$id]
  DF_test <- DF_test[DF_test$time <= DF_test$Time, ]
  
  ## Checking % of censoring
  perc_cens_test[count] <- sum(event==0)/n_test
  
  
  DF.id <- DF[!duplicated(DF$id),]
  DF_test.id <- DF_test[!duplicated(DF_test$id),] 
  ## Survival model
  try(CoxFit <- coxph(Surv(Time, event) ~ sex, data = DF.id))
  
  ## Longitudinal models
  try(LM1 <- lme(y1 ~ sex + time, data = DF, random = ~ time | id))
  try(LM2 <- lme(y2 ~ time, data = DF, random = ~ time | id))
  try(LM3 <- lme(y3 ~ treat + pow2(time), data = DF, random = ~ time | id))
  #try(GLM4 <- mixed_model(y4 ~ sex + treat + time, data = DF,
  #                        random = ~ time || id, family = binomial()))
  
  #we fit the GLMM with just random intercepts
  try(GLM4 <- mixed_model(y4 ~ sex + treat + time, data = DF,
                          random = ~ 1 | id, family = binomial()))
  
  
  #Fitting the multivariate JM
  try(multiJM <- jm(CoxFit, list(LM1, LM2, LM3, GLM4), time_var = "time",
                    n_iter = 13000L, n_burnin = 300L, n_thin = 5L))
  
  
  # If we need to specify the knots we should put:
  # control = list(knots = list(c(4,5,6,7,8,9,10)))))
  # control = list(knots = list(c(4,5,5.5,6,6.5,7,8,9)))
  
  ## OBS: by adding which_independent = "all" we say the model matrix of var-cov
  ## is block diagonal, and it simplifies the problem
  
  try(rhats <- numeric())
  try(for(i in 1:(length(multiJM$statistics$Rhat)-2)){
    rhats <- c(rhats, multiJM$statistics$Rhat[[i]][,1]) 
  })
  try(rhats <- as.numeric(rhats))
  
  
  # It seems that those bs_gammas (the ones related with B-splines used to
  # fit the besaline hazard) are too big (more than 2). Could be that we are
  # assuming ctt baseline hazard function?
  
  ## We compute the metrics, first IBS:
  t0 <- 6
  dt <- 1.5
  try(brier_score_multi_train <- tvBrier(multiJM, newdata = DF, Tstart = t0, Dt = dt, 
                                         integrated = TRUE, type_weights = "IPCW"))
  try(ibs_train[count] <- brier_score_multi_train$Brier)
  try(n_risk_train[count] <- brier_score_multi_train$nr)
  try(n_event_train[count] <- brier_score_multi_train$nint)
  try(n_cens_train[count] <- brier_score_multi_train$ncens)
  
  try(EPCE_score_multi_train <- tvEPCE(multiJM, newdata = DF, Tstart = t0, Dt = dt))
  try(epce_train[count] <- EPCE_score_multi_train$EPCE)
  
  
  
  try(brier_score_multi_test <- tvBrier(multiJM, newdata = DF_test, Tstart = t0, Dt = dt, 
                                        integrated = TRUE, type_weights = "IPCW"))
  try(ibs_test[count] <- brier_score_multi_test$Brier)
  try(n_risk_test[count] <- brier_score_multi_test$nr)
  try(n_event_test[count] <- brier_score_multi_test$nint)
  try(n_cens_test[count] <- brier_score_multi_test$ncens)
  
  try(EPCE_score_multi_test <- tvEPCE(multiJM, newdata = DF_test, Tstart = t0, Dt = dt))
  try(epce_test[count] <- EPCE_score_multi_test$EPCE)
  
  try(list_rhats_MJM <- append(list_rhats_MJM, list(rhats)))
  try(list_full_rhats_MJM <- append(list_full_rhats_MJM, 
                                    list(multiJM$statistics$Rhat)))
  
  
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
    LM4 <- mixed_model(y4 ~ sex + treat + time, data = data,
                       random = ~ time | id, family = binomial())
    JM1 <- jm(CoxFit, LM1, time_var = "time")
    JM2 <- jm(CoxFit, LM2, time_var = "time")
    JM3 <- jm(CoxFit, LM3, time_var = "time")
    JM4 <- jm(CoxFit, LM4, time_var = "time")
    out <- list(M1 = JM1, M2 = JM2, M3 = JM3, M4 = JM4)
    class(out) <- "jmList"
    out
  }
  ## obs: although we have used just random intercepts for the binary outcome,
  ## when fitting the univariate joint model, we cannot use just one long 
  ## outcome and just radom intercepts. Then, we have used for the univariate
  ## model also random slope.
  
  try(cl <- parallel::makeCluster(4L))
  try(Models_folds <- parallel::parLapply(cl, CVdats$training, fit_models))
  try(parallel::stopCluster(cl))
  
  
  
  try(Brier_weights <- tvBrier(Models_folds, newdata = CVdats$testing, 
                               integrated = TRUE, Tstart = t0, Dt = dt,
                               type_weights = "IPCW"))
  
  try(EPCE_weights <- tvEPCE(Models_folds, newdata = CVdats$testing, 
                             Tstart = t0, Dt = dt))
  
  #Now with testing data
  #We fit the models in the whole training data set and test in testing data
  try(Models <- fit_models(DF))
  
  try(bw <- Brier_weights$weights)
  try(Brier_weights_test <- tvBrier(Models, newdata = DF_test, model_weights = bw, 
                                    Tstart = t0, Dt = dt, integrated = TRUE,
                                    type_weights = "IPCW"))
  
  try(ew <- EPCE_weights$weights)
  try(EPCE_weights_test <- tvEPCE(Models, newdata = DF_test, model_weights = ew,
                                  Tstart = t0, Dt = dt))
  
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
  
  disSL_epce[count] <- 0
  try(disSL_epce[count] <- which.min(EPCE_weights$EPCE_per_model))
  if(disSL_epce[count] == 1){
    try(EPCE_dSL_test <- tvEPCE(Models$M1, newdata = DF_test,
                                Tstart = t0, Dt = dt))
  } else if(disSL_epce[count] == 2){
    try(EPCE_dSL_test <- tvEPCE(Models$M2, newdata = DF_test,
                                Tstart = t0, Dt = dt))
  } else if(disSL_epce[count] == 3){
    try(EPCE_dSL_test <- tvEPCE(Models$M3, newdata = DF_test,
                                Tstart = t0, Dt = dt))
  } else if(disSL_epce[count] == 4){
    try(EPCE_dSL_test <- tvEPCE(Models$M4, newdata = DF_test,
                                Tstart = t0, Dt = dt))
  } 
  
  ########################
  #Save the desired metrics
  ###########################
  
  try(IBS_univ[count, ] <- Brier_weights$Brier_per_model)
  try(EPCE_univ[count, ] <- EPCE_weights$EPCE_per_model)
  try(IBS_w[count, ] <- Brier_weights$weights)
  try(EPCE_w[count, ] <- EPCE_weights$weights)
  try(dSL_cv_IBS[count] <- min(Brier_weights$Brier_per_model))
  try(dSL_cv_EPCE[count] <- min(EPCE_weights$EPCE_per_model))
  try(eSL_cv_IBS[count] <- Brier_weights$Brier)
  try(eSL_cv_EPCE[count] <- EPCE_weights$EPCE)
  try(eSL_test_IBS[count] <- Brier_weights_test$Brier)
  try(eSL_test_EPCE[count] <- EPCE_weights_test$EPCE)
  try(dSL_test_IBS[count] <- Brier_dSL_test$Brier)
  try(dSL_test_EPCE[count] <- EPCE_dSL_test$EPCE)
  
  
  try(if(count==25){
    setwd("D:/La meva unitat/TFM/results_SLvsmJM")
    str25 <- "repl25_adminCens_30_16mar2025.RData"
    save(checkTimes_test, checkTimes, perc_cens_test, perc_cens_train,
         list_rhats_MJM, list_full_rhats_MJM, ibs_train, ibs_test, 
         n_risk_train, n_risk_test, n_event_train, n_event_test,
         n_cens_train, n_cens_test, list_w_model_train, list_w_model_test,
         list_brier_model_train, list_brier_model_test,
         IBS_multi, IBS_univ, IBS_w, dSL_cv_IBS, eSL_test_IBS, eSL_cv_IBS,
         IBS_multi_test, dSL_test_IBS, disSL_ibs, 
         EPCE_univ, EPCE_w, dSL_cv_EPCE, eSL_cv_EPCE, dSL_test_EPCE, eSL_test_EPCE,
         epce_train, epce_test, disSL_epce,
         file=str25)
  })
  try(if(count==50){
    setwd("D:/La meva unitat/TFM/results_SLvsmJM")
    str50 <- "repl50_adminCens_30_16mar2025.RData"
    save(checkTimes_test, checkTimes, perc_cens_test, perc_cens_train,
         list_rhats_MJM, list_full_rhats_MJM, ibs_train, ibs_test, 
         n_risk_train, n_risk_test, n_event_train, n_event_test,
         n_cens_train, n_cens_test, list_w_model_train, list_w_model_test,
         list_brier_model_train, list_brier_model_test,
         IBS_multi, IBS_univ, IBS_w, dSL_cv_IBS, eSL_test_IBS, eSL_cv_IBS,
         IBS_multi_test, dSL_test_IBS, disSL_ibs, 
         EPCE_univ, EPCE_w, dSL_cv_EPCE, eSL_cv_EPCE, dSL_test_EPCE, eSL_test_EPCE,
         epce_train, epce_test, disSL_epce,
         file=str50)
  })
  try(if(count==75){
    setwd("D:/La meva unitat/TFM/results_SLvsmJM")
    str75 <- "repl75_adminCens_30_16mar2025.RData"
    save(checkTimes_test, checkTimes, perc_cens_test, perc_cens_train,
         list_rhats_MJM, list_full_rhats_MJM, ibs_train, ibs_test, 
         n_risk_train, n_risk_test, n_event_train, n_event_test,
         n_cens_train, n_cens_test, list_w_model_train, list_w_model_test,
         list_brier_model_train, list_brier_model_test,
         IBS_multi, IBS_univ, IBS_w, dSL_cv_IBS, eSL_test_IBS, eSL_cv_IBS,
         IBS_multi_test, dSL_test_IBS, disSL_ibs, 
         EPCE_univ, EPCE_w, dSL_cv_EPCE, eSL_cv_EPCE, dSL_test_EPCE, eSL_test_EPCE,
         epce_train, epce_test, disSL_epce,
         file=str75)
  })
  
  if(count == repl){
    setwd("D:/La meva unitat/TFM/results_SLvsmJM")
    strr <- "repl100_adminCens_30_16mar2025.RData"
    save(checkTimes_test, checkTimes, perc_cens_test, perc_cens_train,
         list_rhats_MJM, list_full_rhats_MJM, ibs_train, ibs_test, 
         n_risk_train, n_risk_test, n_event_train, n_event_test,
         n_cens_train, n_cens_test, list_w_model_train, list_w_model_test,
         list_brier_model_train, list_brier_model_test,
         IBS_multi, IBS_univ, IBS_w, dSL_cv_IBS, eSL_test_IBS, eSL_cv_IBS,
         IBS_multi_test, dSL_test_IBS, disSL_ibs, 
         EPCE_univ, EPCE_w, dSL_cv_EPCE, eSL_cv_EPCE, dSL_test_EPCE, eSL_test_EPCE,
         epce_train, epce_test, disSL_epce,
         file=strr)
  }
  paste0(count, " Admin Cens")
}