###############################################################################
###############################################################################
#SMIULATION STUDY
#we follow (https://drizopoulos.github.io/JMbayes2/articles/Non_Gaussian_Mixed_Models.html)
###############################################################################
###############################################################################
library(survival)
library(JMbayes2)
#We have L=5 longitudinal outcomes, 3 LMM, 2 GLMM. We simulate then JM data


set.seed(12345)
n <- 6*30 # number of subjects
K <- 15 # number of measurements per subject
t_max <- 20 # maximum follow-up time

# we construct a data frame with the design:
# everyone has a baseline measurement, and then measurements at random 
# follow-up times up to t_max
DF <- data.frame(id = rep(seq_len(n), each = K),
                 time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                 sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

# design matrices for the fixed and random effects for longitudinal submodels
X1 <- model.matrix(~ sex * time, data = DF)
Z1 <- model.matrix(~ time, data = DF)

X2 <- model.matrix(~ sex + time, data = DF )
Z2 <- model.matrix(~ time, data = DF)

X3 <- model.matrix(~ time, data = DF )
Z3 <- model.matrix(~ time, data = DF)

X4 <- model.matrix(~ sex + time, data = DF)
Z4 <- model.matrix(~ time, data = DF)

X5 <- model.matrix(~ time, data = DF)
Z5 <- model.matrix(~ time, data = DF)

#Simulate first longitudinal outcome
###############################################
betas1 <- c(-2.2, -0.25, 1.24, -0.05) # fixed effects coefficients
sigma1 <- 0.125 # errors sd
D11_1 <- 1.0 # variance of random intercepts
D22_1 <- 0.5 # variance of random slopes

# we simulate random effects
b1 <- cbind(rnorm(n, sd = sqrt(D11_1)), rnorm(n, sd = sqrt(D22_1)))
# linear predictor
eta_y1 <- as.vector(X1 %*% betas1 + rowSums(Z1 * b1[DF$id, ]))
# we simulate normal longitudinal data
DF$y1 <- rnorm(n * K, mean = eta_y1, sd = sigma1)
# we assume that values below 4 are not observed, and set equal to 0
DF$ind1 <- as.numeric(DF$y1 < -4)
DF$y1 <- pmax(DF$y1, -4)

#Simulate second longitudinal outcome
###############################################
betas2 <- c(-1.8, -0.06, 0.5) # fixed effects coefficients
sigma2 <- 0.25 # errors sd
D11_2 <- 1.2 # variance of random intercepts
D22_2 <- 0.25 # variance of random slopes

# we simulate random effects
b2 <- cbind(rnorm(n, sd = sqrt(D11_2)), rnorm(n, sd = sqrt(D22_2)))
# linear predictor
eta_y2 <- as.vector(X2 %*% betas2 + rowSums(Z2 * b2[DF$id, ]))
# we simulate normal longitudinal data
DF$y2 <- rnorm(n * K, mean = eta_y2, sd = sigma2)
# we assume that values below 0 are not observed, and set equal to 0
DF$ind2 <- as.numeric(DF$y2 < -4)
DF$y2 <- pmax(DF$y2, -4)

#Simulate third longitudinal outcome
###############################################
betas3 <- c(-2.5, 0.0333) # fixed effects coefficients
sigma3 <- 0.25 # errors sd
D11_3 <- 0.15 # variance of random intercepts
D22_3 <- 0.05 # variance of random slopes


# we simulate random effects
b3 <- cbind(rnorm(n, sd = sqrt(D11_3)), rnorm(n, sd = sqrt(D22_3)))
# linear predictor
eta_y3 <- as.vector(X3 %*% betas3 + rowSums(Z3 * b3[DF$id, ]))
# we simulate normal longitudinal data
DF$y3 <- rnorm(n * K, mean = eta_y3, sd = sigma3)
# we assume that values below 0 are not observed, and set equal to 0
#DF$ind3 <- as.numeric(DF$y3 < 0)
#DF$y3 <- pmax(DF$y3, 0)


#Simulate forth longitudinal outcome
###############################################
betas4 <- c(0.01, 0.5, -0.31416) # fixed effects coefficients
D11_4 <- 0.212 # variance of random intercepts
D22_4 <- 0.0125 # variance of random slopes


# we simulate random effects
b4 <- cbind(rnorm(n, sd = sqrt(D11_4)), rnorm(n, sd = sqrt(D22_4)))
# linear predictor
eta_y4 <- as.vector(X4 %*% betas4 + rowSums(Z4 * b4[DF$id, ]))
# mean of the binomial distribution
mu_y4 <- plogis(eta_y4)
# we simulate binomial longitudinal data
DF$y4 <- rbinom(n * K, size = 1, prob = mu_y4)


#Simulate fifth longitudinal outcome
###############################################
betas5 <- c(1, 0.155) # fixed effects coefficients
D11_5 <- 0.0121 # variance of random intercepts
D22_5 <- 0.0754 # variance of random slopes


# we simulate random effects
b5 <- cbind(rnorm(n, sd = sqrt(D11_5)), rnorm(n, sd = sqrt(D22_5)))
# linear predictor
eta_y5 <- as.vector(X5 %*% betas5 + rowSums(Z5 * b5[DF$id, ]))
# mean of the binomial distribution
mu_y5 <- plogis(eta_y5)
# we simulate binomial longitudinal data
DF$y5 <- rbinom(n * K, size = 1, prob = mu_y5)

#OBS: other glmm distr (beta, gamma, negative binom, etc) can be used


#Simulate event times
########################################
upp_Cens <- 25 # fixed Type I censoring time
shape_wb <- 4.5 # shape Weibull
alpha <- c(0.8, 2.777, 3, 5, 0.74) # association coefficients
gammas <- c("(Intercept)" = -15, "sex" = 0.5)
W <- model.matrix(~ sex, data = DF[!duplicated(DF$id), ])
# linear predictor for the survival model
eta_t <- as.vector(W %*% gammas)
# to simulate event times we use inverse transform sampling
# (https://en.wikipedia.org/wiki/Inverse_transform_sampling). Namely, we want 
# to find t, such that S(t) = u, where S(.) is the survival function, and u a 
# number from the Unif(0, 1) distribution. The function below calculates 
# log(u) - log(S(t)), and for a given u, we want to find t for which it equals
# zero. We do that below using the uniroot() function
invS <- function (t, i) {
  # i denotes the subject
  sex_i <- W[i, 2L]
  # h() is the hazard function and we assume a Weibull baseline hazard
  h <- function (s) {
    X1_at_s <- cbind(1, sex_i, s, sex_i * s)
    Z1_at_s <- cbind(1, s)
    X2_at_s <- cbind(1, sex_i, s)
    Z2_at_s <- cbind(1, s)
    X3_at_s <- cbind(1, s)
    Z3_at_s <- cbind(1, s)
    X4_at_s <- cbind(1, sex_i, s)
    Z4_at_s <- cbind(1, s)
    X5_at_s <- cbind(1, s)
    Z5_at_s <- cbind(1, s)
    # the linear predictor from the mixed model evaluated at time s
    f1 <- as.vector(X1_at_s %*% betas1 +
                     rowSums(Z1_at_s * b1[rep(i, nrow(Z1_at_s)), ]))
    f2 <- as.vector(X2_at_s %*% betas2 +
                      rowSums(Z2_at_s * b2[rep(i, nrow(Z2_at_s)), ]))
    f3 <- as.vector(X3_at_s %*% betas3 +
                      rowSums(Z3_at_s * b3[rep(i, nrow(Z3_at_s)), ]))
    f4 <- as.vector(X4_at_s %*% betas4 +
                      rowSums(Z4_at_s * b4[rep(i, nrow(Z4_at_s)), ]))
    f5 <- as.vector(X5_at_s %*% betas5 +
                      rowSums(Z5_at_s * b5[rep(i, nrow(Z5_at_s)), ]))
    exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] 
          + f1*alpha[1] + f2*alpha[2] + f3*alpha[3] + f4*alpha[4] + f5*alpha[5])
  }
  # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
  integrate(h, lower = 0, upper = t)$value + log(u[i])
}

invS1 <- function (t, i) {
  # i denotes the subject
  sex_i <- W[i, 2L]
  # h() is the hazard function and we assume a Weibull baseline hazard
  h1 <- function (s) {
    X1_at_s <- cbind(1, sex_i, s, sex_i * s)
    Z1_at_s <- cbind(1, s)
    X2_at_s <- cbind(1, sex_i, s)
    Z2_at_s <- cbind(1, s)
    X3_at_s <- cbind(1, s)
    Z3_at_s <- cbind(1, s)
    X4_at_s <- cbind(1, sex_i, s)
    Z4_at_s <- cbind(1, s)
    X5_at_s <- cbind(1, s)
    Z5_at_s <- cbind(1, s)
    # the linear predictor from the mixed model evaluated at time s
    f1 <- as.vector(X1_at_s %*% betas1 +
                      rowSums(Z1_at_s * b1[rep(i, nrow(Z1_at_s)), ]))
    exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] 
        + f1*alpha[1])
  }
  # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
  integrate(h1, lower = 0, upper = t)$value + log(u[i])
}

invS2 <- function (t, i) {
  # i denotes the subject
  sex_i <- W[i, 2L]
  # h() is the hazard function and we assume a Weibull baseline hazard
  h2 <- function (s) {
    X1_at_s <- cbind(1, sex_i, s, sex_i * s)
    Z1_at_s <- cbind(1, s)
    X2_at_s <- cbind(1, sex_i, s)
    Z2_at_s <- cbind(1, s)
    X3_at_s <- cbind(1, s)
    Z3_at_s <- cbind(1, s)
    X4_at_s <- cbind(1, sex_i, s)
    Z4_at_s <- cbind(1, s)
    X5_at_s <- cbind(1, s)
    Z5_at_s <- cbind(1, s)
    # the linear predictor from the mixed model evaluated at time s
    f2 <- as.vector(X2_at_s %*% betas2 +
                      rowSums(Z2_at_s * b2[rep(i, nrow(Z2_at_s)), ]))
    exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] 
        + f2*alpha[2])
  }
  # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
  integrate(h2, lower = 0, upper = t)$value + log(u[i])
}

invS3 <- function (t, i) {
  # i denotes the subject
  sex_i <- W[i, 2L]
  # h() is the hazard function and we assume a Weibull baseline hazard
  h3 <- function (s) {
    X1_at_s <- cbind(1, sex_i, s, sex_i * s)
    Z1_at_s <- cbind(1, s)
    X2_at_s <- cbind(1, sex_i, s)
    Z2_at_s <- cbind(1, s)
    X3_at_s <- cbind(1, s)
    Z3_at_s <- cbind(1, s)
    X4_at_s <- cbind(1, sex_i, s)
    Z4_at_s <- cbind(1, s)
    X5_at_s <- cbind(1, s)
    Z5_at_s <- cbind(1, s)
    # the linear predictor from the mixed model evaluated at time s
    f3 <- as.vector(X3_at_s %*% betas3 +
                      rowSums(Z3_at_s * b3[rep(i, nrow(Z3_at_s)), ]))
    exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] 
        + f3*alpha[3])
  }
  # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
  integrate(h3, lower = 0, upper = t)$value + log(u[i])
}

invS4 <- function (t, i) {
  # i denotes the subject
  sex_i <- W[i, 2L]
  # h() is the hazard function and we assume a Weibull baseline hazard
  h4 <- function (s) {
    X1_at_s <- cbind(1, sex_i, s, sex_i * s)
    Z1_at_s <- cbind(1, s)
    X2_at_s <- cbind(1, sex_i, s)
    Z2_at_s <- cbind(1, s)
    X3_at_s <- cbind(1, s)
    Z3_at_s <- cbind(1, s)
    X4_at_s <- cbind(1, sex_i, s)
    Z4_at_s <- cbind(1, s)
    X5_at_s <- cbind(1, s)
    Z5_at_s <- cbind(1, s)
    # the linear predictor from the mixed model evaluated at time s
    f4 <- as.vector(X4_at_s %*% betas4 +
                      rowSums(Z4_at_s * b4[rep(i, nrow(Z4_at_s)), ]))
    exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] 
        + f4*alpha[4])
  }
  # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
  integrate(h4, lower = 0, upper = t)$value + log(u[i])
}

invS5 <- function (t, i) {
  # i denotes the subject
  sex_i <- W[i, 2L]
  # h() is the hazard function and we assume a Weibull baseline hazard
  h5 <- function (s) {
    X1_at_s <- cbind(1, sex_i, s, sex_i * s)
    Z1_at_s <- cbind(1, s)
    X2_at_s <- cbind(1, sex_i, s)
    Z2_at_s <- cbind(1, s)
    X3_at_s <- cbind(1, s)
    Z3_at_s <- cbind(1, s)
    X4_at_s <- cbind(1, sex_i, s)
    Z4_at_s <- cbind(1, s)
    X5_at_s <- cbind(1, s)
    Z5_at_s <- cbind(1, s)
    # the linear predictor from the mixed model evaluated at time s
    f5 <- as.vector(X5_at_s %*% betas5 +
                      rowSums(Z5_at_s * b5[rep(i, nrow(Z5_at_s)), ]))
    exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] 
        + f5*alpha[5])
  }
  # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
  integrate(h5, lower = 0, upper = t)$value + log(u[i])
}
# we simulate the event times
u <- runif(n)
trueTimes <- numeric(n)
for (i in seq_len(n/6)) {
  Up <- 100
  Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
  trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
}
for (i in (n/6):(2*(n/6))) {
  Up <- 100
  Root <- try(uniroot(invS1, interval = c(1e-05, Up), i = i)$root, TRUE)
  trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
}
for (i in (2*(n/6)):(3*(n/6))) {
  Up <- 100
  Root <- try(uniroot(invS2, interval = c(1e-05, Up), i = i)$root, TRUE)
  trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
}
for (i in (3*(n/6)):(4*(n/6))) {
  Up <- 100
  Root <- try(uniroot(invS3, interval = c(1e-05, Up), i = i)$root, TRUE)
  trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
}
for (i in (4*(n/6)):(5*(n/6))) {
  Up <- 100
  Root <- try(uniroot(invS4, interval = c(1e-05, Up), i = i)$root, TRUE)
  trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
}
for (i in (5*(n/6)):(6*(n/6))) {
  Up <- 100
  Root <- try(uniroot(invS5, interval = c(1e-05, Up), i = i)$root, TRUE)
  trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
}
# we use fixed Type I right censoring denoting the end of the trial.
Ctimes <- upp_Cens
Time <- pmin(trueTimes, Ctimes)
event <- as.numeric(trueTimes <= Ctimes) # event indicator

# we keep the longitudinal measurements before the event times
DF$Time <- Time[DF$id]
DF$event <- event[DF$id]
DF <- DF[DF$time <= DF$Time, ]

#Now we have the data frame with the simulated data





################################################################################
###############################################################################
#Applying the analysis to the simulated data frame
################################################################################
###############################################################################
DF.id <- DF[!duplicated(DF$id),]
sum(DF.id$event==0)

## Survival model
CoxFit <- coxph(Surv(Time, event) ~ sex, data = DF.id)

## Longitudinal models
LM1 <- lme(y1 ~ sex*time, data = DF, random = ~ time | id)
LM2 <- lme(y2 ~ sex + time, data = DF, random = ~ time | id)
LM3 <- lme(y3 ~ time, data = DF, random = ~ time | id)
LM4 <- mixed_model(y4 ~ sex + time, data = DF,
                   random = ~ time || id, family = binomial())
LM5 <- mixed_model(y5 ~ time, data = DF,
                   random = ~ time || id, family = binomial())

##########################################################################
#Fitting the multivariate JM
#see https://drizopoulos.github.io/JMbayes2/articles/JMbayes2.html#multivariate
multiJM <- jm(CoxFit, list(LM1, LM2, LM3, LM4, LM5), time_var = "time",
              which_independent = "all",
              n_iter = 12000L, n_burnin = 2000L, n_thin = 5L)
summary(multiJM)


#calculate the integrated Brier score as an overall measure of predictive 
#performance in (t0, t0+ Dt]=(12, 12+6]
t0 <- 15
dt <- 8
brier_score_multi <- tvBrier(multiJM, newdata = DF, Tstart = t0, Dt = dt, 
                             integrated = TRUE)
EPCE_score_multi <- tvEPCE(multiJM, newdata = DF, Tstart = t0, Dt = dt)

#calculate the integrated Brier score as an overall measure of predictive 
#performance in (t0, t0+ Dt]=(7, 7+3]
brier_score_multi_2 <- tvBrier(multiJM, newdata = DF, Tstart = 7, Dt = 3, 
                               integrated = TRUE)
EPCE_score_multi_2 <- tvEPCE(multiJM, newdata = DF, Tstart = 7, Dt = 3)

##########################################################################
#SuperLearning with the library of models built with the univariate JM

CVdats <- create_folds(DF, V = 3, id_var = "id")

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
brier_score_SL_full <- tvBrier(Models, newdata = DF, model_weights = bw, 
                               Tstart = t0, Dt = dt, integrated = TRUE)
#we use Brier weights to compute IBS with full data (TIME WINDOW 2):
bw2 <- Brier_weights_2$weights
brier_score_SL_full_2 <- tvBrier(Models, newdata = DF, model_weights = bw2, 
                                 Tstart = 7, Dt = 3, integrated = TRUE)

#we use EPCE weights to compute EPCE with full data:
ew <- EPCE_weights$weights
EPCE_score_SL_full <- tvEPCE(Models, newdata = DF, model_weights = ew,
                             Tstart = t0, Dt = dt)
#we use EPCE weights to compute EPCE with full data (NEW TIME WINDOW):
ew2 <- EPCE_weights_2$weights
EPCE_score_SL_full_2 <- tvEPCE(Models, newdata = DF, model_weights = ew2,
                               Tstart = 7, Dt = 3)



##############################################################################
##############################################################################
##############################################################################
#Let us do the simulation of num_dat datasets in a loop
##############################################################################
##############################################################################
setwd("G:/TFM/SLinJointModels")

library(JMbayes2)
num_dat <- 150
IBS_multi <- EPCE_multi <- IBS_multi_test <- EPCE_multi_test <- numeric(num_dat)
dSL_cv_IBS <- eSL_cv_IBS <- dSL_test_IBS <- eSL_test_IBS <- numeric(num_dat)
dSL_cv_EPCE <- eSL_cv_EPCE <- dSL_test_EPCE <- eSL_test_EPCE <- numeric(num_dat)
IBS_univ <- EPCE_univ <- IBS_w <- EPCE_w <- matrix(nrow = num_dat, ncol = 5)
IBS_univ_test <- EPCE_univ_test <- IBS_w_test <- 
  EPCE_w_test <- matrix(nrow = num_dat, ncol = 5)
censoring_train <- censoring_test <- numeric(num_dat)
count <- 0

for(count in 1:num_dat){
  ############################
  #We generate the dataset
  ###########################
  n <- 6*50 # number of subjects
  n_test <- 6*50
  K <- 15 # number of measurements per subject
  t_max <- 20 # maximum follow-up time
  
  # we construct a data frame with the design:
  # everyone has a baseline measurement, and then measurements at random 
  # follow-up times up to t_max
  DF <- data.frame(id = rep(seq_len(n), each = K),
                   time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                   sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))
  DF_test <- data.frame(id = rep(seq_len(n_test), each = K),
                   time = c(replicate(n_test, c(0, sort(runif(K - 1, 0, t_max))))),
                   sex = rep(gl(2, n_test/2, labels = c("male", "female")), each = K))
  
  
  # design matrices for the fixed and random effects for longitudinal submodels
  X1 <- model.matrix(~ sex * time, data = DF)
  Z1 <- model.matrix(~ time, data = DF)
  
  X2 <- model.matrix(~ sex + time, data = DF )
  Z2 <- model.matrix(~ time, data = DF)
  
  X3 <- model.matrix(~ time, data = DF )
  Z3 <- model.matrix(~ time, data = DF)
  
  X4 <- model.matrix(~ sex + time, data = DF)
  Z4 <- model.matrix(~ time, data = DF)
  
  X5 <- model.matrix(~ time, data = DF)
  Z5 <- model.matrix(~ time, data = DF)
  
  ####for test data set
  X1_test <- model.matrix(~ sex * time, data = DF_test)
  Z1_test <- model.matrix(~ time, data = DF_test)
  
  X2_test <- model.matrix(~ sex + time, data = DF_test )
  Z2_test <- model.matrix(~ time, data = DF_test)
  
  X3_test <- model.matrix(~ time, data = DF_test )
  Z3_test <- model.matrix(~ time, data = DF_test)
  
  X4_test <- model.matrix(~ sex + time, data = DF_test)
  Z4_test <- model.matrix(~ time, data = DF_test)
  
  X5_test <- model.matrix(~ time, data = DF_test)
  Z5_test <- model.matrix(~ time, data = DF_test)
  
  #Simulate first longitudinal outcome
  ###############################################
  betas1 <- c(-2.2, -0.25, 1.24, -0.05) # fixed effects coefficients
  sigma1 <- 0.125 # errors sd
  D11_1 <- 1.0 # variance of random intercepts
  D22_1 <- 0.5 # variance of random slopes
  
  # we simulate random effects
  b1 <- cbind(rnorm(n, sd = sqrt(D11_1)), rnorm(n, sd = sqrt(D22_1)))
  # linear predictor
  eta_y1 <- as.vector(X1 %*% betas1 + rowSums(Z1 * b1[DF$id, ]))
  # we simulate normal longitudinal data
  DF$y1 <- rnorm(n * K, mean = eta_y1, sd = sigma1)
  # we assume that values below 4 are not observed, and set equal to 0
  DF$ind1 <- as.numeric(DF$y1 < -4)
  DF$y1 <- pmax(DF$y1, -4)
  
  #Test data:
  # we simulate random effects
  b1_test <- cbind(rnorm(n_test, sd = sqrt(D11_1)), rnorm(n_test, sd = sqrt(D22_1)))
  # linear predictor
  eta_y1_test <- as.vector(X1_test %*% betas1 + rowSums(Z1_test * b1_test[DF_test$id, ]))
  # we simulate normal longitudinal data
  DF_test$y1 <- rnorm(n_test * K, mean = eta_y1_test, sd = sigma1)
  # we assume that values below 4 are not observed, and set equal to 0
  DF_test$ind1 <- as.numeric(DF_test$y1 < -4)
  DF_test$y1 <- pmax(DF_test$y1, -4)
  
  #Simulate second longitudinal outcome
  ###############################################
  betas2 <- c(-1.8, -0.06, 0.5) # fixed effects coefficients
  sigma2 <- 0.25 # errors sd
  D11_2 <- 1.2 # variance of random intercepts
  D22_2 <- 0.25 # variance of random slopes
  
  # we simulate random effects
  b2 <- cbind(rnorm(n, sd = sqrt(D11_2)), rnorm(n, sd = sqrt(D22_2)))
  # linear predictor
  eta_y2 <- as.vector(X2 %*% betas2 + rowSums(Z2 * b2[DF$id, ]))
  # we simulate normal longitudinal data
  DF$y2 <- rnorm(n * K, mean = eta_y2, sd = sigma2)
  # we assume that values below 0 are not observed, and set equal to 0
  DF$ind2 <- as.numeric(DF$y2 < -4)
  DF$y2 <- pmax(DF$y2, -4)
  
  ####Test data
  # we simulate random effects
  b2_test <- cbind(rnorm(n_test, sd = sqrt(D11_2)), rnorm(n_test, sd = sqrt(D22_2)))
  # linear predictor
  eta_y2_test <- as.vector(X2_test %*% betas2 + rowSums(Z2_test * b2[DF_test$id, ]))
  # we simulate normal longitudinal data
  DF_test$y2 <- rnorm(n_test * K, mean = eta_y2_test, sd = sigma2)
  # we assume that values below 0 are not observed, and set equal to 0
  DF_test$ind2 <- as.numeric(DF_test$y2 < -4)
  DF_test$y2 <- pmax(DF_test$y2, -4)
  
  #Simulate third longitudinal outcome
  ###############################################
  betas3 <- c(-2.5, 0.0333) # fixed effects coefficients
  sigma3 <- 0.25 # errors sd
  D11_3 <- 0.15 # variance of random intercepts
  D22_3 <- 0.05 # variance of random slopes
  
  
  # we simulate random effects
  b3 <- cbind(rnorm(n, sd = sqrt(D11_3)), rnorm(n, sd = sqrt(D22_3)))
  # linear predictor
  eta_y3 <- as.vector(X3 %*% betas3 + rowSums(Z3 * b3[DF$id, ]))
  # we simulate normal longitudinal data
  DF$y3 <- rnorm(n * K, mean = eta_y3, sd = sigma3)
  # we assume that values below 0 are not observed, and set equal to 0
  #DF$ind3 <- as.numeric(DF$y3 < 0)
  #DF$y3 <- pmax(DF$y3, 0)
  
  ########Test data
  # we simulate random effects
  b3_test <- cbind(rnorm(n_test, sd = sqrt(D11_3)), rnorm(n_test, sd = sqrt(D22_3)))
  # linear predictor
  eta_y3_test <- as.vector(X3_test %*% betas3 + rowSums(Z3_test * b3_test[DF_test$id, ]))
  # we simulate normal longitudinal data
  DF_test$y3 <- rnorm(n_test * K, mean = eta_y3_test, sd = sigma3)
  
  #Simulate forth longitudinal outcome
  ###############################################
  betas4 <- c(0.01, 0.5, -0.31416) # fixed effects coefficients
  D11_4 <- 0.212 # variance of random intercepts
  D22_4 <- 0.0125 # variance of random slopes
  
  
  # we simulate random effects
  b4 <- cbind(rnorm(n, sd = sqrt(D11_4)), rnorm(n, sd = sqrt(D22_4)))
  # linear predictor
  eta_y4 <- as.vector(X4 %*% betas4 + rowSums(Z4 * b4[DF$id, ]))
  # mean of the binomial distribution
  mu_y4 <- plogis(eta_y4)
  # we simulate binomial longitudinal data
  DF$y4 <- rbinom(n * K, size = 1, prob = mu_y4)
  
  #####Test data
  # we simulate random effects
  b4_test <- cbind(rnorm(n_test, sd = sqrt(D11_4)), rnorm(n_test, sd = sqrt(D22_4)))
  # linear predictor
  eta_y4_test <- as.vector(X4_test %*% betas4 + rowSums(Z4_test * b4_test[DF_test$id, ]))
  # mean of the binomial distribution
  mu_y4_test <- plogis(eta_y4_test)
  # we simulate binomial longitudinal data
  DF_test$y4 <- rbinom(n_test * K, size = 1, prob = mu_y4_test)
  
  
  #Simulate fifth longitudinal outcome
  ###############################################
  betas5 <- c(1, 0.155) # fixed effects coefficients
  D11_5 <- 0.0121 # variance of random intercepts
  D22_5 <- 0.0754 # variance of random slopes
  
  
  # we simulate random effects
  b5 <- cbind(rnorm(n, sd = sqrt(D11_5)), rnorm(n, sd = sqrt(D22_5)))
  # linear predictor
  eta_y5 <- as.vector(X5 %*% betas5 + rowSums(Z5 * b5[DF$id, ]))
  # mean of the binomial distribution
  mu_y5 <- plogis(eta_y5)
  # we simulate binomial longitudinal data
  DF$y5 <- rbinom(n * K, size = 1, prob = mu_y5)
  
  ####Test data
  # we simulate random effects
  b5_test <- cbind(rnorm(n_test, sd = sqrt(D11_5)), rnorm(n_test, sd = sqrt(D22_5)))
  # linear predictor
  eta_y5_test <- as.vector(X5_test %*% betas5 + rowSums(Z5_test * b5_test[DF_test$id, ]))
  # mean of the binomial distribution
  mu_y5_test <- plogis(eta_y5_test)
  # we simulate binomial longitudinal data
  DF_test$y5 <- rbinom(n_test * K, size = 1, prob = mu_y5_test)
  
  #OBS: other glmm distr (beta, gamma, negative binom, etc) can be used
  
  
  #Simulate event times
  ########################################
  upp_Cens <- 25 # fixed Type I censoring time
  shape_wb <- 4.5 # shape Weibull
  alpha <- c(0.8, 2.777, 3, 5, 0.74) # association coefficients
  gammas <- c("(Intercept)" = -15, "sex" = 0.5)
  W <- model.matrix(~ sex, data = DF[!duplicated(DF$id), ])
  # linear predictor for the survival model
  eta_t <- as.vector(W %*% gammas)
  # to simulate event times we use inverse transform sampling
  # (https://en.wikipedia.org/wiki/Inverse_transform_sampling). Namely, we want 
  # to find t, such that S(t) = u, where S(.) is the survival function, and u a 
  # number from the Unif(0, 1) distribution. The function below calculates 
  # log(u) - log(S(t)), and for a given u, we want to find t for which it equals
  # zero. We do that below using the uniroot() function
  invS <- function (t, i) {
    # i denotes the subject
    sex_i <- W[i, 2L]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h <- function (s) {
      X1_at_s <- cbind(1, sex_i, s, sex_i * s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, sex_i, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, s)
      Z3_at_s <- cbind(1, s)
      X4_at_s <- cbind(1, sex_i, s)
      Z4_at_s <- cbind(1, s)
      X5_at_s <- cbind(1, s)
      Z5_at_s <- cbind(1, s)
      # the linear predictor from the mixed model evaluated at time s
      f1 <- as.vector(X1_at_s %*% betas1 +
                        rowSums(Z1_at_s * b1[rep(i, nrow(Z1_at_s)), ]))
      f2 <- as.vector(X2_at_s %*% betas2 +
                        rowSums(Z2_at_s * b2[rep(i, nrow(Z2_at_s)), ]))
      f3 <- as.vector(X3_at_s %*% betas3 +
                        rowSums(Z3_at_s * b3[rep(i, nrow(Z3_at_s)), ]))
      f4 <- as.vector(X4_at_s %*% betas4 +
                        rowSums(Z4_at_s * b4[rep(i, nrow(Z4_at_s)), ]))
      f5 <- as.vector(X5_at_s %*% betas5 +
                        rowSums(Z5_at_s * b5[rep(i, nrow(Z5_at_s)), ]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] 
          + f1*alpha[1] + f2*alpha[2] + f3*alpha[3] + f4*alpha[4] + f5*alpha[5])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h, lower = 0, upper = t)$value + log(u[i])
  }
  invS1 <- function (t, i) {
    # i denotes the subject
    sex_i <- W[i, 2L]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h1 <- function (s) {
      X1_at_s <- cbind(1, sex_i, s, sex_i * s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, sex_i, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, s)
      Z3_at_s <- cbind(1, s)
      X4_at_s <- cbind(1, sex_i, s)
      Z4_at_s <- cbind(1, s)
      X5_at_s <- cbind(1, s)
      Z5_at_s <- cbind(1, s)
      # the linear predictor from the mixed model evaluated at time s
      f1 <- as.vector(X1_at_s %*% betas1 +
                        rowSums(Z1_at_s * b1[rep(i, nrow(Z1_at_s)), ]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] 
          + f1*alpha[1])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h1, lower = 0, upper = t)$value + log(u[i])
  }
  
  invS2 <- function (t, i) {
    # i denotes the subject
    sex_i <- W[i, 2L]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h2 <- function (s) {
      X1_at_s <- cbind(1, sex_i, s, sex_i * s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, sex_i, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, s)
      Z3_at_s <- cbind(1, s)
      X4_at_s <- cbind(1, sex_i, s)
      Z4_at_s <- cbind(1, s)
      X5_at_s <- cbind(1, s)
      Z5_at_s <- cbind(1, s)
      # the linear predictor from the mixed model evaluated at time s
      f2 <- as.vector(X2_at_s %*% betas2 +
                        rowSums(Z2_at_s * b2[rep(i, nrow(Z2_at_s)), ]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] 
          + f2*alpha[2])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h2, lower = 0, upper = t)$value + log(u[i])
  }
  
  invS3 <- function (t, i) {
    # i denotes the subject
    sex_i <- W[i, 2L]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h3 <- function (s) {
      X1_at_s <- cbind(1, sex_i, s, sex_i * s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, sex_i, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, s)
      Z3_at_s <- cbind(1, s)
      X4_at_s <- cbind(1, sex_i, s)
      Z4_at_s <- cbind(1, s)
      X5_at_s <- cbind(1, s)
      Z5_at_s <- cbind(1, s)
      # the linear predictor from the mixed model evaluated at time s
      f3 <- as.vector(X3_at_s %*% betas3 +
                        rowSums(Z3_at_s * b3[rep(i, nrow(Z3_at_s)), ]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] 
          + f3*alpha[3])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h3, lower = 0, upper = t)$value + log(u[i])
  }
  
  invS4 <- function (t, i) {
    # i denotes the subject
    sex_i <- W[i, 2L]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h4 <- function (s) {
      X1_at_s <- cbind(1, sex_i, s, sex_i * s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, sex_i, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, s)
      Z3_at_s <- cbind(1, s)
      X4_at_s <- cbind(1, sex_i, s)
      Z4_at_s <- cbind(1, s)
      X5_at_s <- cbind(1, s)
      Z5_at_s <- cbind(1, s)
      # the linear predictor from the mixed model evaluated at time s
      f4 <- as.vector(X4_at_s %*% betas4 +
                        rowSums(Z4_at_s * b4[rep(i, nrow(Z4_at_s)), ]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] 
          + f4*alpha[4])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h4, lower = 0, upper = t)$value + log(u[i])
  }
  
  invS5 <- function (t, i) {
    # i denotes the subject
    sex_i <- W[i, 2L]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h5 <- function (s) {
      X1_at_s <- cbind(1, sex_i, s, sex_i * s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, sex_i, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, s)
      Z3_at_s <- cbind(1, s)
      X4_at_s <- cbind(1, sex_i, s)
      Z4_at_s <- cbind(1, s)
      X5_at_s <- cbind(1, s)
      Z5_at_s <- cbind(1, s)
      # the linear predictor from the mixed model evaluated at time s
      f5 <- as.vector(X5_at_s %*% betas5 +
                        rowSums(Z5_at_s * b5[rep(i, nrow(Z5_at_s)), ]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] 
          + f5*alpha[5])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h5, lower = 0, upper = t)$value + log(u[i])
  }
  # we simulate the event times
  u <- runif(n)
  trueTimes <- numeric(n)
  for (i in seq_len(n/6)) {
    Up <- 100
    Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  for (i in (n/6):(2*(n/6))) {
    Up <- 100
    Root <- try(uniroot(invS1, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  for (i in (2*(n/6)):(3*(n/6))) {
    Up <- 100
    Root <- try(uniroot(invS2, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  for (i in (3*(n/6)):(4*(n/6))) {
    Up <- 100
    Root <- try(uniroot(invS3, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  for (i in (4*(n/6)):(5*(n/6))) {
    Up <- 100
    Root <- try(uniroot(invS4, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  for (i in (5*(n/6)):(6*(n/6))) {
    Up <- 100
    Root <- try(uniroot(invS5, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  # we use fixed Type I right censoring denoting the end of the trial.
  Ctimes <- upp_Cens
  Time <- pmin(trueTimes, Ctimes)
  event <- as.numeric(trueTimes <= Ctimes) # event indicator
  
  # we keep the longitudinal measurements before the event times
  DF$Time <- Time[DF$id]
  DF$event <- event[DF$id]
  DF <- DF[DF$time <= DF$Time, ]
  
  
  #######Test data##################
  W_test <- model.matrix(~ sex, data = DF_test[!duplicated(DF_test$id), ])
  # linear predictor for the survival model
  eta_t_test <- as.vector(W_test %*% gammas)
  # to simulate event times we use inverse transform sampling
  # (https://en.wikipedia.org/wiki/Inverse_transform_sampling). Namely, we want 
  # to find t, such that S(t) = u, where S(.) is the survival function, and u a 
  # number from the Unif(0, 1) distribution. The function below calculates 
  # log(u) - log(S(t)), and for a given u, we want to find t for which it equals
  # zero. We do that below using the uniroot() function
  invS <- function (t, i) {
    # i denotes the subject
    sex_i <- W_test[i, 2L]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h <- function (s) {
      X1_at_s <- cbind(1, sex_i, s, sex_i * s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, sex_i, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, s)
      Z3_at_s <- cbind(1, s)
      X4_at_s <- cbind(1, sex_i, s)
      Z4_at_s <- cbind(1, s)
      X5_at_s <- cbind(1, s)
      Z5_at_s <- cbind(1, s)
      # the linear predictor from the mixed model evaluated at time s
      f1 <- as.vector(X1_at_s %*% betas1 +
                        rowSums(Z1_at_s * b1_test[rep(i, nrow(Z1_at_s)), ]))
      f2 <- as.vector(X2_at_s %*% betas2 +
                        rowSums(Z2_at_s * b2_test[rep(i, nrow(Z2_at_s)), ]))
      f3 <- as.vector(X3_at_s %*% betas3 +
                        rowSums(Z3_at_s * b3_test[rep(i, nrow(Z3_at_s)), ]))
      f4 <- as.vector(X4_at_s %*% betas4 +
                        rowSums(Z4_at_s * b4_test[rep(i, nrow(Z4_at_s)), ]))
      f5 <- as.vector(X5_at_s %*% betas5 +
                        rowSums(Z5_at_s * b5_test[rep(i, nrow(Z5_at_s)), ]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t_test[i] 
          + f1*alpha[1] + f2*alpha[2] + f3*alpha[3] + f4*alpha[4] + f5*alpha[5])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h, lower = 0, upper = t)$value + log(u[i])
  }
  
  invS1 <- function (t, i) {
    # i denotes the subject
    sex_i <- W_test[i, 2L]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h1 <- function (s) {
      X1_at_s <- cbind(1, sex_i, s, sex_i * s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, sex_i, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, s)
      Z3_at_s <- cbind(1, s)
      X4_at_s <- cbind(1, sex_i, s)
      Z4_at_s <- cbind(1, s)
      X5_at_s <- cbind(1, s)
      Z5_at_s <- cbind(1, s)
      # the linear predictor from the mixed model evaluated at time s
      f1 <- as.vector(X1_at_s %*% betas1 +
                        rowSums(Z1_at_s * b1_test[rep(i, nrow(Z1_at_s)), ]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t_test[i] 
          + f1*alpha[1])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h1, lower = 0, upper = t)$value + log(u[i])
  }
  
  invS2 <- function (t, i) {
    # i denotes the subject
    sex_i <- W_test[i, 2L]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h2 <- function (s) {
      X1_at_s <- cbind(1, sex_i, s, sex_i * s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, sex_i, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, s)
      Z3_at_s <- cbind(1, s)
      X4_at_s <- cbind(1, sex_i, s)
      Z4_at_s <- cbind(1, s)
      X5_at_s <- cbind(1, s)
      Z5_at_s <- cbind(1, s)
      # the linear predictor from the mixed model evaluated at time s
      f2 <- as.vector(X2_at_s %*% betas2 +
                        rowSums(Z2_at_s * b2_test[rep(i, nrow(Z2_at_s)), ]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t_test[i] 
          + f2*alpha[2])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h2, lower = 0, upper = t)$value + log(u[i])
  }
  
  invS3 <- function (t, i) {
    # i denotes the subject
    sex_i <- W_test[i, 2L]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h3 <- function (s) {
      X1_at_s <- cbind(1, sex_i, s, sex_i * s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, sex_i, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, s)
      Z3_at_s <- cbind(1, s)
      X4_at_s <- cbind(1, sex_i, s)
      Z4_at_s <- cbind(1, s)
      X5_at_s <- cbind(1, s)
      Z5_at_s <- cbind(1, s)
      # the linear predictor from the mixed model evaluated at time s
      f3 <- as.vector(X3_at_s %*% betas3 +
                        rowSums(Z3_at_s * b3_test[rep(i, nrow(Z3_at_s)), ]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t_test[i] 
          + f3*alpha[3])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h3, lower = 0, upper = t)$value + log(u[i])
  }
  
  invS4 <- function (t, i) {
    # i denotes the subject
    sex_i <- W_test[i, 2L]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h4 <- function (s) {
      X1_at_s <- cbind(1, sex_i, s, sex_i * s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, sex_i, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, s)
      Z3_at_s <- cbind(1, s)
      X4_at_s <- cbind(1, sex_i, s)
      Z4_at_s <- cbind(1, s)
      X5_at_s <- cbind(1, s)
      Z5_at_s <- cbind(1, s)
      # the linear predictor from the mixed model evaluated at time s
      f4 <- as.vector(X4_at_s %*% betas4 +
                        rowSums(Z4_at_s * b4_test[rep(i, nrow(Z4_at_s)), ]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t_test[i] 
          + f4*alpha[4])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h4, lower = 0, upper = t)$value + log(u[i])
  }
  
  invS5 <- function (t, i) {
    # i denotes the subject
    sex_i <- W_test[i, 2L]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h5 <- function (s) {
      X1_at_s <- cbind(1, sex_i, s, sex_i * s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, sex_i, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, s)
      Z3_at_s <- cbind(1, s)
      X4_at_s <- cbind(1, sex_i, s)
      Z4_at_s <- cbind(1, s)
      X5_at_s <- cbind(1, s)
      Z5_at_s <- cbind(1, s)
      # the linear predictor from the mixed model evaluated at time s
      f5 <- as.vector(X5_at_s %*% betas5 +
                        rowSums(Z5_at_s * b5_test[rep(i, nrow(Z5_at_s)), ]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t_test[i] 
          + f5*alpha[5])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h5, lower = 0, upper = t)$value + log(u[i])
  }
  
  # we simulate the event times
  u <- runif(n)
  trueTimes <- numeric(n)
  for (i in seq_len(n/6)) {
    Up <- 100
    Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  for (i in (n/6):(2*(n/6))) {
    Up <- 100
    Root <- try(uniroot(invS1, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  for (i in (2*(n/6)):(3*(n/6))) {
    Up <- 100
    Root <- try(uniroot(invS2, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  for (i in (3*(n/6)):(4*(n/6))) {
    Up <- 100
    Root <- try(uniroot(invS3, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  for (i in (4*(n/6)):(5*(n/6))) {
    Up <- 100
    Root <- try(uniroot(invS4, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  for (i in (5*(n/6)):(6*(n/6))) {
    Up <- 100
    Root <- try(uniroot(invS5, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  
  # we use fixed Type I right censoring denoting the end of the trial.
  Ctimes <- upp_Cens
  Time <- pmin(trueTimes, Ctimes)
  event <- as.numeric(trueTimes <= Ctimes) # event indicator
  
  # we keep the longitudinal measurements before the event times
  DF_test$Time <- Time[DF_test$id]
  DF_test$event <- event[DF_test$id]
  DF_test <- DF_test[DF_test$time <= DF_test$Time, ]
  
  ###################################################
  #We fit the models and save the desired metrics
  ###################################################
  DF.id <- DF[!duplicated(DF$id),]
  DF_test.id <- DF_test[!duplicated(DF_test$id),] 
  ## Survival model
  try(CoxFit <- coxph(Surv(Time, event) ~ sex, data = DF.id))
  
  ## Longitudinal models
  try(LM1 <- lme(y1 ~ sex*time, data = DF, random = ~ time | id))
  try(LM2 <- lme(y2 ~ sex + time, data = DF, random = ~ time | id))
  try(LM3 <- lme(y3 ~ time, data = DF, random = ~ time | id))
  try(LM4 <- mixed_model(y4 ~ sex + time, data = DF,
                     random = ~ time || id, family = binomial()))
  try(LM5 <- mixed_model(y5 ~ time, data = DF,
                     random = ~ time || id, family = binomial()))
  

  #Fitting the multivariate JM
  try(multiJM <- jm(CoxFit, list(LM1, LM2, LM3, LM4, LM5), time_var = "time",
                which_independent = "all",
                n_iter = 12000L, n_burnin = 2000L, n_thin = 5L))
  cat("Multi JM fitted\n")
  
  #calculate the integrated Brier score as an overall measure of predictive 
  #performance in (t0, t0+ Dt]=(5, 5+3]
  t0 <- 12
  dt <- 6
  try(brier_score_multi <- tvBrier(multiJM, newdata = DF, Tstart = t0, Dt = dt, 
                               integrated = TRUE))
  try(EPCE_score_multi <- tvEPCE(multiJM, newdata = DF, Tstart = t0, Dt = dt))
  cat("Brier and EPCE multi computed\n")
  #for the test data
  try(brier_score_multi_test <- tvBrier(multiJM, newdata = DF_test, Tstart = t0, Dt = dt, 
                               integrated = TRUE))
  try(EPCE_score_multi_test <- tvEPCE(multiJM, newdata = DF_test, Tstart = t0, Dt = dt))
  cat("Brier and EPCE multi test computed\n")
  
  #SuperLearning with the library of models built with the univariate JM
  
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
  
  cat("Cross-validated models fitted\n")
  
  #computing Brier weights
  try(Brier_weights <- tvBrier(Models_folds, newdata = CVdats$testing, 
                           integrated = TRUE, Tstart = t0, Dt = dt))
  #computing EPCE weights
  try(EPCE_weights <- tvEPCE(Models_folds, newdata = CVdats$testing, 
                         Tstart = t0, Dt = dt))
  cat("IBS and EPCE for SL computed\n")
  
  #Now with testing data
  #We fir the models in the whole training data set and test in testing data
  try(Models <- fit_models(DF))
  cat("Full models fitted\n")

  try(bw <- Brier_weights$weights)
  try(Brier_weights_test <- tvBrier(Models, newdata = DF_test, model_weights = bw, 
                                 Tstart = t0, Dt = dt, integrated = TRUE))
  try(ew <- EPCE_weights$weights)
  try(EPCE_weights_test <- tvEPCE(Models, newdata = DF_test, model_weights = ew,
                               Tstart = t0, Dt = dt))
  cat("Brier and EPCE SL in test data computed\n")
  
  ########################
  #Save the desired metrics
  ###########################
  try(IBS_multi[count] <- brier_score_multi$Brier)
  try(EPCE_multi[count] <- EPCE_score_multi$EPCE)
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
  try(IBS_multi_test[count] <- brier_score_multi_test$Brier)
  try(EPCE_multi_test[count] <- EPCE_score_multi_test$EPCE)
  
  try(censoring_train[count] <- sum(DF.id$event==0))
  try(censoring_test[count] <- sum(DF_test.id$event==0))
  
  #hay que guardar solo los objetos que necesitamos yh no todo!!!!!!
  #si no pesa demasiado
  if(count==1){
    strr <- "results1iter_Scenario2.RData"
    save(IBS_multi, EPCE_multi, IBS_univ, EPCE_univ, IBS_w, EPCE_w, dSL_cv_IBS,
         dSL_cv_EPCE, eSL_cv_IBS, eSL_cv_EPCE, eSL_test_IBS,  eSL_test_EPCE, 
         IBS_multi_test, EPCE_multi_test, censoring_train,
         censoring_test, file=strr)
  }
  if(count==2){
    str1 <- "results2_Scenario2.RData"
    save(IBS_multi, EPCE_multi, IBS_univ, EPCE_univ, IBS_w, EPCE_w, dSL_cv_IBS,
         dSL_cv_EPCE, eSL_cv_IBS, eSL_cv_EPCE, eSL_test_IBS,  eSL_test_EPCE, 
         IBS_multi_test, EPCE_multi_test, censoring_train,
         censoring_test, file=str1)
  }
  if(count==10){
    str10 <- "results10_Scenario2.RData"
    save(IBS_multi, EPCE_multi, IBS_univ, EPCE_univ, IBS_w, EPCE_w, dSL_cv_IBS,
         dSL_cv_EPCE, eSL_cv_IBS, eSL_cv_EPCE, eSL_test_IBS,  eSL_test_EPCE, 
         IBS_multi_test, EPCE_multi_test, censoring_train,
         censoring_test, file=str10)
  }
  if(count==20){
    str20 <- "results20_sce2.RData"
    save(IBS_multi, EPCE_multi, IBS_univ, EPCE_univ, IBS_w, EPCE_w, dSL_cv_IBS,
         dSL_cv_EPCE, eSL_cv_IBS, eSL_cv_EPCE, eSL_test_IBS,  eSL_test_EPCE, 
         IBS_multi_test, EPCE_multi_test, censoring_train,
         censoring_test, file=str20)
  }
  if(count==30){
    str30 <- "results30_sce2.RData"
    save(IBS_multi, EPCE_multi, IBS_univ, EPCE_univ, IBS_w, EPCE_w, dSL_cv_IBS,
         dSL_cv_EPCE, eSL_cv_IBS, eSL_cv_EPCE, eSL_test_IBS,  eSL_test_EPCE, 
         IBS_multi_test, EPCE_multi_test, censoring_train,
         censoring_test, file=str30)
  }
  if(count==40){
    str40 <- "results40_sce2.RData"
    save(IBS_multi, EPCE_multi, IBS_univ, EPCE_univ, IBS_w, EPCE_w, dSL_cv_IBS,
         dSL_cv_EPCE, eSL_cv_IBS, eSL_cv_EPCE, eSL_test_IBS,  eSL_test_EPCE, 
         IBS_multi_test, EPCE_multi_test, censoring_train,
         censoring_test, file=str40)
  }
  if(count==50){
    str50 <- "results50_sce2.RData"
    save(IBS_multi, EPCE_multi, IBS_univ, EPCE_univ, IBS_w, EPCE_w, dSL_cv_IBS,
         dSL_cv_EPCE, eSL_cv_IBS, eSL_cv_EPCE, eSL_test_IBS,  eSL_test_EPCE, 
         IBS_multi_test, EPCE_multi_test, censoring_train,
         censoring_test, file=str50)
  }
  if(count==60){
    str60 <- "results60_sce2.RData"
    save(IBS_multi, EPCE_multi, IBS_univ, EPCE_univ, IBS_w, EPCE_w, dSL_cv_IBS,
         dSL_cv_EPCE, eSL_cv_IBS, eSL_cv_EPCE, eSL_test_IBS,  eSL_test_EPCE, 
         IBS_multi_test, EPCE_multi_test, censoring_train,
         censoring_test, file=str60)
  }
  if(count==70){
    str70 <- "results70_sce2.RData"
    save(IBS_multi, EPCE_multi, IBS_univ, EPCE_univ, IBS_w, EPCE_w, dSL_cv_IBS,
         dSL_cv_EPCE, eSL_cv_IBS, eSL_cv_EPCE, eSL_test_IBS,  eSL_test_EPCE, 
         IBS_multi_test, EPCE_multi_test, censoring_train,
         censoring_test, file=str70)
  }
  if(count==100){
    str100 <- "results100_sce2.RData"
    save(IBS_multi, EPCE_multi, IBS_univ, EPCE_univ, IBS_w, EPCE_w, dSL_cv_IBS,
         dSL_cv_EPCE, eSL_cv_IBS, eSL_cv_EPCE, eSL_test_IBS,  eSL_test_EPCE, 
         IBS_multi_test, EPCE_multi_test, censoring_train,
         censoring_test, file=str100)
  }
  if(count==120){
    str120 <- "results120_sce.RData"
    save(IBS_multi, EPCE_multi, IBS_univ, EPCE_univ, IBS_w, EPCE_w, dSL_cv_IBS,
         dSL_cv_EPCE, eSL_cv_IBS, eSL_cv_EPCE, eSL_test_IBS,  eSL_test_EPCE, 
         IBS_multi_test, EPCE_multi_test, censoring_train,
         censoring_test, file=str120)
  }
  if(count==150){
    str150 <- "results150_sce2.RData"
    save(IBS_multi, EPCE_multi, IBS_univ, EPCE_univ, IBS_w, EPCE_w, dSL_cv_IBS,
         dSL_cv_EPCE, eSL_cv_IBS, eSL_cv_EPCE, eSL_test_IBS,  eSL_test_EPCE, 
         IBS_multi_test, EPCE_multi_test, censoring_train,
         censoring_test, file=str150)
  }
  if(count==175){
    str175 <- "results175_sce2.RData"
    save(IBS_multi, EPCE_multi, IBS_univ, EPCE_univ, IBS_w, EPCE_w, dSL_cv_IBS,
         dSL_cv_EPCE, eSL_cv_IBS, eSL_cv_EPCE, eSL_test_IBS,  eSL_test_EPCE, 
         IBS_multi_test, EPCE_multi_test, censoring_train,
         censoring_test, file=str175)
  }
  if(count==199){
    str200 <- "results200.RData"
    save(IBS_multi, EPCE_multi, IBS_univ, EPCE_univ, IBS_w, EPCE_w, dSL_cv_IBS,
         dSL_cv_EPCE, eSL_cv_IBS, eSL_cv_EPCE, eSL_test_IBS,  eSL_test_EPCE, 
         IBS_multi_test, EPCE_multi_test, censoring_train,
         censoring_test, file=str200)
  }
}



