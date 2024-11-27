## Let's see the IBS IPCW when using the data set simulated for the deepsurv:
library(JMbayes2)

try(CoxFit <- coxph(Surv(Time, event) ~ sex, data = DF.id))

try(LM1 <- lme(y1 ~ sex*time, data = DF, random = ~ time | id))
try(LM2 <- lme(y2 ~ sex + time, data = DF, random = ~ time | id))
try(LM3 <- lme(y3 ~ time, data = DF, random = ~ time | id))
try(LM4 <- mixed_model(y4 ~ sex + time, data = DF,
                       random = ~ time || id, family = binomial()))
try(LM5 <- mixed_model(y5 ~ time, data = DF,
                       random = ~ time || id, family = binomial()))

try(multiJM <- jm(CoxFit, list(LM1, LM2, LM3, LM4, LM5), time_var = "time",
                  n_iter = 12000L, n_burnin = 2000L, n_thin = 5L))

t0 <- 4
dt <- 1.5
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
  LM2 <- lme(y2 ~ sex + time, data = data, random = ~ time | id)
  LM3 <- lme(y3 ~ time, data = data, random = ~ time | id)
  LM4 <- mixed_model(y4 ~ sex + time, data = data,
                     random = ~ time || id, family = binomial())
  LM5 <- mixed_model(y5 ~ time, data = data,
                     random = ~ time || id, family = binomial())
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