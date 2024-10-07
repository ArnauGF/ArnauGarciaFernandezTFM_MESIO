library(JMbayes2)
pbc2 <- pbc2
pbc2.id <- pbc2.id

## Survival model
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)

## Longitudinal models
LM1 <- lme(log(serBilir) ~ year, data = pbc2, random = ~ year | id)
LM2 <- lme(prothrombin ~ year * sex, data = pbc2, random = ~ year | id)
LM3 <- mixed_model(ascites ~ year + sex, data = pbc2,
                   random = ~ year || id, family = binomial())
LM4 <- lme(albumin ~ year, data = pbc2, random = ~ year | id)
LM5 <- mixed_model(hepatomegaly ~ year, data = pbc2,
                   random = ~ year || id, family = binomial())

##########################################################################
#Fitting the multivariate JM
#see https://drizopoulos.github.io/JMbayes2/articles/JMbayes2.html#multivariate
multiJM <- jm(CoxFit, list(LM1, LM2, LM3, LM4, LM5), time_var = "year",
                which_independent = "all",
                n_iter = 12000L, n_burnin = 2000L, n_thin = 5L)
summary(multiJM)
#We compute a prediction for Patients 25 and 93. 
t0 <- 5
dt <- 3
ND <- pbc2[pbc2$id %in% c(25, 93), ]
ND <- ND[ND$year < t0, ]
ND$event <- 0
ND$years <- t0

#In the following example, we calculate predictions from time t0 to time 12
predSurv <- predict(multiJM, newdata = ND, process = "event",
                    times = seq(t0, 12, length.out = 51),
                    return_newdata = TRUE)
plot(predSurv)

#calculate the integrated Brier score as an overall measure of predictive 
#performance in (t0, t0+ Dt]=(5, 5+3]
brier_score_multi <- tvBrier(multiJM, newdata = pbc2, Tstart = t0, Dt = dt, 
                             integrated = TRUE)
EPCE_score_multi <- tvEPCE(multiJM, newdata = pbc2, Tstart = t0, Dt = dt)

#calculate the integrated Brier score as an overall measure of predictive 
#performance in (t0, t0+ Dt]=(7, 7+3]
brier_score_multi_2 <- tvBrier(multiJM, newdata = pbc2, Tstart = 7, Dt = 3, 
                             integrated = TRUE)
EPCE_score_multi_2 <- tvEPCE(multiJM, newdata = pbc2, Tstart = 7, Dt = 3)

##########################################################################
#SuperLearning with the library of models built with the univariate JM

CVdats <- create_folds(pbc2, V = 3, id_var = "id")

fit_models <- function (data) {
  library("JMbayes2")
  data$status2 <- as.numeric(data$status != "alive")
  data_id <- data[!duplicated(data$id), ]
  CoxFit <- coxph(Surv(years, status2) ~ sex, data = data_id)
  LM1 <- lme(log(serBilir) ~ year, data = data, random = ~ year | id)
  LM2 <- lme(prothrombin ~ year * sex, data = data, random = ~ year | id)
  LM3 <- mixed_model(ascites ~ year + sex, data = data,
                     random = ~ year || id, family = binomial())
  LM4 <- lme(albumin ~ year, data = data, random = ~ year | id)
  LM5 <- mixed_model(hepatomegaly ~ year, data = data,
                     random = ~ year || id, family = binomial())
  JM1 <- jm(CoxFit, LM1, time_var = "year")
  JM2 <- jm(CoxFit, LM2, time_var = "year")
  JM3 <- jm(CoxFit, LM3, time_var = "year")
  JM4 <- jm(CoxFit, LM4, time_var = "year")
  JM5 <- jm(CoxFit, LM5, time_var = "year")
  out <- list(M1 = JM1, M2 = JM2, M3 = JM3, M4 = JM4, M5 = JM5)
  class(out) <- "jmList"
  out
}

cl <- parallel::makeCluster(5L)
Models_folds <- parallel::parLapply(cl, CVdats$training, fit_models)
parallel::stopCluster(cl)

#computing Brier weights
Brier_weights <- tvBrier(Models_folds, newdata = CVdats$testing, 
                         integrated = TRUE, Tstart = t0, Dt = dt)
#computing Brier weights in another time window
Brier_weights_2 <- tvBrier(Models_folds, newdata = CVdats$testing, 
                         integrated = TRUE, Tstart = 7, Dt = 3)
#computing EPCE weights
EPCE_weights <- tvEPCE(Models_folds, newdata = CVdats$testing, 
                       Tstart = t0, Dt = dt)
#computing EPCE weights in another time window
EPCE_weights_2 <- tvEPCE(Models_folds, newdata = CVdats$testing, 
                       Tstart = 7, Dt = 3)

#Now we use the weights with the whole data set:
Models <- fit_models(pbc2)

#we use Brier weights to compute IBS with full data:
bw <- Brier_weights$weights
brier_score_SL_full <- tvBrier(Models, newdata = pbc2, model_weights = bw, 
                               Tstart = t0, Dt = dt, integrated = TRUE)
#we use Brier weights to compute IBS with full data (TIME WINDOW 2):
bw2 <- Brier_weights_2$weights
brier_score_SL_full_2 <- tvBrier(Models, newdata = pbc2, model_weights = bw2, 
                               Tstart = 7, Dt = 3, integrated = TRUE)

#we use EPCE weights to compute EPCE with full data:
ew <- EPCE_weights$weights
EPCE_score_SL_full <- tvEPCE(Models, newdata = pbc2, model_weights = ew,
                             Tstart = t0, Dt = dt)
#we use EPCE weights to compute EPCE with full data (NEW TIME WINDOW):
ew2 <- EPCE_weights_2$weights
EPCE_score_SL_full_2 <- tvEPCE(Models, newdata = pbc2, model_weights = ew2,
                             Tstart = 7, Dt = 3)


load("CaseStudy.RData")

