#############################################################################
#############################################################################
## Building main simulations (MJM vs SL-bases predictions) form the scratch.
## Just L=2, block Diagonal matrix (no correlated outcomes), baseline
## hazard function ctt (so exponential distr)
## Also longitudinal time observations will be a fixed grid 1 to 10 
## (instead of U(0,20))
#############################################################################
#############################################################################
##Libraries:
library(JMbayes2)
library(MASS)
library(lattice)
#############################################################################


##################################
### l=1
##################################

### We simulate the first longitudinal outcome, LMM (we want positive slope)

n <- 250 # number of subjects
n_test <- 250
K <- 10 # number of measurements per subject
t_max <- 20 # maximum follow-up time

set.seed(1234)
sigmaa <- matrix(c(runif(1,0.25,0.75), runif(7, -0.0005, 0.0005), 
                   runif(1,0,0.15), runif(7, -0.0005, 0.0005),
                   runif(1,0.25,0.75), runif(7, -0.0005, 0.0005),
                   runif(1,0,0.15), runif(7, -0.0005, 0.0005),
                   runif(1,0.25,0.75), runif(7, -0.0005, 0.0005),
                   runif(1,0,0.025), runif(7, -0.0005, 0.0005),
                   runif(1,0,0.05))
                 ,nrow = 7)
treat <- rep(as.factor(sample(c("A", "B"), size = n, replace = TRUE))
             ,each = K)

DF <- data.frame(id = rep(seq_len(n), each = K),
                 sex = rep(gl(2, n/2, labels = c("male", "female")), each = K),
                 treat = treat,
                 time = c(replicate(n, c(0, 1:(K-1)))))
DF_test <- data.frame(id = rep(seq_len(n_test), each = K),
                      sex = rep(gl(2, n_test/2, labels = c("male", "female")), each = K),
                      treat = treat,
                      time = c(replicate(n, c(0, 1:(K-1)))))

#correlated random effects
b <-  mvrnorm(n = n, mu = rep(0,7), Sigma = sigmaa%*%t(sigmaa))
b_test <-  mvrnorm(n = n, mu = rep(0,7), Sigma = sigmaa%*%t(sigmaa))

## Now we define the model

##Train:
X1 <- model.matrix(~ sex + time, data = DF)
Z1 <- model.matrix(~ time, data = DF)

##Test:
X1_test <- model.matrix(~ sex + time, data = DF_test)
Z1_test <- model.matrix(~ time, data = DF_test)

##Paramaters specification:
#b1 <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(0.067, 0.0025, 0.0000178, 0.0008988), 
#                                                 nrow = 2, ncol = 2))
b1 <- b[, c(1,2)]

betas1 <- c(0.123, 0.098765, 1.5)
sigma1 <- 0.125 #error term sd

## Longitudinal profile simulation:
eta_y1 <- as.vector(X1 %*% betas1 + rowSums(Z1 * b1[DF$id, ]))
DF$y1 <- rnorm(n * K, mean = eta_y1, sd = sigma1)


#b1_test <- mvrnorm(n = n_test, mu = c(0,0), Sigma = matrix(c(0.067, 0.0025, 0.0000178, 0.0008988), 
#                                                 nrow = 2, ncol = 2))
b1_test <- b_test[, c(1,2)]

eta_y1_test <- as.vector(X1_test %*% betas1 + rowSums(Z1_test * b1_test[DF_test$id, ]))
DF_test$y1 <- rnorm(n * K, mean = eta_y1_test, sd = sigma1)

# Let us see the longitudinal profiles
xyplot(y1 ~ time, groups = id,
       data = DF,
       type = "l" ,xlab="Time",ylab="Y1")

xyplot(y1 ~ time, groups = id,
       data = DF_test,
       type = "l" ,xlab="Time",ylab="Y1")


##obs: Having smaller random slopes we get the longitudinal profiles in an 
## smaller range of values, also big parameters for binary baseline covariates
## implies a jump in the longitudinal profiles


##################################
### l=2
##################################

X2 <- model.matrix(~ time, data = DF)
Z2 <- model.matrix(~ time, data = DF)
X2_test <- model.matrix(~ time, data = DF_test)
Z2_test <- model.matrix(~ time, data = DF_test)

##Paramaters specification:
#b2 <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(0.067, 0.00025, 0.00178, 0.0008988), 
#                                                 nrow = 2, ncol = 2))
b2 <- b[, c(3,4)]

betas2 <- c(15.641, -1.25918)
sigma2 <- 0.125 #error term sd

##Longitudinal profile simulation
eta_y2 <- as.vector(X2 %*% betas2 + rowSums(Z2 * b2[DF$id, ]))
DF$y2 <- rnorm(n * K, mean = eta_y2, sd = sigma2)

#b2_test <- mvrnorm(n = n_test, mu = c(0,0), Sigma = matrix(c(0.067, 0.00025, 0.00178, 0.0008988), 
#                                                           nrow = 2, ncol = 2))
b2_test <- b_test[, c(3,4)]

eta_y2_test <- as.vector(X2_test %*% betas2 + rowSums(Z2_test * b2_test[DF_test$id, ]))
DF_test$y2 <- rnorm(n * K, mean = eta_y2_test, sd = sigma2)

# Let us see the longitudinal profiles
xyplot(y2 ~ time, groups = id,
       data = DF,
       type = "l" ,xlab="Time",ylab="Y2")

xyplot(y2 ~ time, groups = id,
       data = DF_test,
       type = "l" ,xlab="Time",ylab="Y2")


##################################
### l=3
##################################
#We define a function to compute the square of a given number
pow2 <- function(x){
  x^2
}

X3 <- model.matrix(~ treat + pow2(time), data = DF)
Z3 <- model.matrix(~ pow2(time), data = DF)
X3_test <- model.matrix(~ treat + pow2(time), data = DF_test)
Z3_test <- model.matrix(~ pow2(time), data = DF_test)

##Paramaters specification:
#b3 <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(0.009123, 0.000025, 0.0005716, 0.00008988), 
#                                                 nrow = 2, ncol = 2))
b3 <- b[, c(5,6)]
betas3 <- c(1.314159, 0.125, 0.18181)
sigma3 <- 0.125 #error term sd

##Longitudinal profile simulation
eta_y3 <- as.vector(X3 %*% betas3 + rowSums(Z3 * b3[DF$id, ]))
DF$y3 <- rnorm(n * K, mean = eta_y3, sd = sigma2)


#b3_test <- mvrnorm(n = n_test, mu = c(0,0), Sigma = matrix(c(0.009123, 0.000025, 0.0005716, 0.00008988), 
#                                                           nrow = 2, ncol = 2))
b3_test <- b_test[, c(5,6)]
eta_y3_test <- as.vector(X3_test %*% betas3 + rowSums(Z3_test * b3_test[DF_test$id, ]))
DF_test$y3 <- rnorm(n * K, mean = eta_y3_test, sd = sigma3)

# Let us see the longitudinal profiles
xyplot(y3 ~ time, groups = id,
       data = DF,
       type = "l" ,xlab="Time",ylab="Y3")

xyplot(y3 ~ time, groups = id,
       data = DF_test,
       type = "l" ,xlab="Time",ylab="Y3")


##################################
### l=4
## now we want to use a glmm to
## model a binary longitudinal
## outcome
##################################

X4 <- model.matrix(~ sex + treat + time, data = DF)
Z4 <- model.matrix(~ 1, data = DF)
X4_test <- model.matrix(~ sex + treat + time, data = DF_test)
Z4_test <- model.matrix(~ 1, data = DF_test)

##Paramaters specification:
#b4 <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(0.009123, 0.000025, 0.0005716, 0.00008988), 
#                                                 nrow = 2, ncol = 2))
b4 <- b[, 7]
betas4 <- c(-2.33, -0.07152, 0.42313, 1)
sigma4 <- 0.125 #error term sd

##Longitudinal profile simulation
eta_y4 <- as.vector(X4 %*% betas4 + rowSums(Z4 * b4[DF$id]))
DF$eta_y4 <- eta_y4
DF$y4 <- rbinom(n*K, size = 1, prob = plogis(eta_y4))


#b4_test <- mvrnorm(n = n_test, mu = c(0,0), Sigma = matrix(c(0.009123, 0.000025, 0.0005716, 0.00008988), 
#                                                           nrow = 2, ncol = 2))
b4_test <- b_test[, 7]
eta_y4_test <- as.vector(X4_test %*% betas4 + rowSums(Z4_test * b4_test[DF_test$id]))
DF_test$eta_y4 <- eta_y4_test
DF_test$y4 <- rbinom(n_test*K, size = 1, prob = plogis(eta_y4_test))

# Let us see the longitudinal profiles
xyplot(y4 ~ time, groups = id,
       data = DF,
       type = "l" ,xlab="Time",ylab="Y4")

xyplot(y4 ~ time, groups = id,
       data = DF_test,
       type = "l" ,xlab="Time",ylab="Y4")


xyplot(eta_y4 ~ time, groups = id,
       data = DF,
       type = "l" ,xlab="Time",ylab="eta Y4")

xyplot(eta_y4 ~ time, groups = id,
       data = DF_test,
       type = "l" ,xlab="Time",ylab="eta Y4")

xyplot(plogis(eta_y4) ~ time, groups = id,
       data = DF,
       type = "l" ,xlab="Time",ylab="eta Y4")

xyplot(plogis(eta_y4) ~ time, groups = id,
       data = DF_test,
       type = "l" ,xlab="Time",ylab="eta Y4")

##################################
### Time-to-event simulation
##################################

# We assume baseline hazard to be ctt (exponential distr). We use the inverse 
# transform sampling method to simulate survival times.

#parameters specification:
shape_wb <- 1.2 # we assume exp distr and ctt baseline hazard
alpha <- c(0.1, -0.2, 0.312) # association coefficients
alpha <- c(0.0879, -0.43, 0.2, -0.01)
#alpha <- 2*alpha
gammas <- c("(Intercept)" = 0.25, "sex" = -0.5)
#gammas <- c("(Intercept)" = 0.15, "sex" = -0.25)
#model matrix and linear predictor:
W <- model.matrix(~ sex, data = DF[!duplicated(DF$id), ])
treat_ind <- DF[!duplicated(DF$id), ]$treat
eta_t <- as.vector(W %*% gammas)

invS <- function (t, i) {
  # i denotes the subject
  sex_i <- W[i, 2L]
  treat_i <- treat_ind[i]
  # h() is the hazard function and we assume a Weibull baseline hazard
  h <- function (s) {
    X1_at_s <- cbind(1, sex_i, s)
    Z1_at_s <- cbind(1, s)
    X2_at_s <- cbind(1, s)
    Z2_at_s <- cbind(1, s)
    X3_at_s <- cbind(1, treat_i, pow2(s))
    Z3_at_s <- cbind(1,pow2(s))
    X4_at_s <- cbind(1, sex_i, treat_i, s)
    Z4_at_s <- cbind(1)
    # the linear predictor from the mixed model evaluated at time s
    f1 <- as.vector(X1_at_s %*% betas1 +
                      rowSums(Z1_at_s * b1[rep(i, nrow(Z1_at_s)), ]))
    f2 <- as.vector(X2_at_s %*% betas2 +
                      rowSums(Z2_at_s * b2[rep(i, nrow(Z2_at_s)), ]))
    f3 <- as.vector(X3_at_s %*% betas3 +
                      rowSums(Z3_at_s * b3[rep(i, nrow(Z3_at_s)), ]))
    f4 <- as.vector(X4_at_s %*% betas4 +
                      rowSums(Z4_at_s * b4[rep(i, nrow(Z4_at_s))]))
    exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] + f1*alpha[1] + 
          f2*alpha[2] + f3*alpha[3] + f4*alpha[4])
  }
  # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
  integrate(h, lower = 0, upper = t)$value + log(u[i])
}

# we simulate the event times
u <- runif(n)
trueTimes <- numeric(n)
for (i in seq_len(n)) {
  Up <- 100
  Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
  trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
}


## Testing data:
#model matrix and linear predictor:
W_test <- model.matrix(~ sex, data = DF_test[!duplicated(DF_test$id), ])
treat_ind_test <- DF_test[!duplicated(DF_test$id), ]$treat
eta_t_test <- as.vector(W_test %*% gammas)

invS <- function (t, i) {
  # i denotes the subject
  sex_i <- W_test[i, 2L]
  treat_i <- treat_ind_test[i]
  # h() is the hazard function and we assume a Weibull baseline hazard
  h <- function (s) {
    X1_at_s <- cbind(1, sex_i, s)
    Z1_at_s <- cbind(1, s)
    X2_at_s <- cbind(1, s)
    Z2_at_s <- cbind(1, s)
    X3_at_s <- cbind(1, treat_i, pow2(s))
    Z3_at_s <- cbind(1,pow2(s))
    X4_at_s <- cbind(1, sex_i, treat_i, s)
    Z4_at_s <- cbind(1)
    # the linear predictor from the mixed model evaluated at time s
    f1 <- as.vector(X1_at_s %*% betas1 +
                      rowSums(Z1_at_s * b1_test[rep(i, nrow(Z1_at_s)), ]))
    f2 <- as.vector(X2_at_s %*% betas2 +
                      rowSums(Z2_at_s * b2_test[rep(i, nrow(Z2_at_s)), ]))
    f3 <- as.vector(X3_at_s %*% betas3 +
                      rowSums(Z3_at_s * b3_test[rep(i, nrow(Z3_at_s)), ]))
    f4 <- as.vector(X4_at_s %*% betas4 +
                      rowSums(Z4_at_s * b4_test[rep(i, nrow(Z4_at_s))]))
    exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t_test[i] + f1*alpha[1] + 
          f2*alpha[2] + f3*alpha[3] + f4*alpha[4])
  }
  # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
  integrate(h, lower = 0, upper = t)$value + log(u[i])
}

# we simulate the event times
u <- runif(n)
trueTimes_test <- numeric(n)
for (i in seq_len(n)) {
  Up <- 100
  Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
  trueTimes_test[i] <- if (!inherits(Root, "try-error")) Root else 150
}



##saving true times
#checkTimes[count] <- trueTimes

#save Times in data frame
Time <- trueTimes
DF$Time <- Time[DF$id]

##saving true times
#checkTimes_test[count] <- trueTimes_test

#save Times in data frame
Time <- trueTimes_test
DF_test$Time <- Time[DF_test$id]


sum(trueTimes==150)/250
sum(trueTimes_test==150)/250


mean(trueTimes[trueTimes<150])
mean(trueTimes_test[trueTimes_test<150])

hist(trueTimes[trueTimes<150])
hist(trueTimes_test[trueTimes_test<150])


#once we have the time-to-event simulation well-calibrated we can see how the 
# model fitting goes

##############################
## Fitting models
##############################

#simulating censoring; Type I
Ctimes <- 6.7
Time <- pmin(trueTimes, Ctimes)
event <- as.numeric(trueTimes <= Ctimes) # event indicator
# we keep the longitudinal measurements before the event times
DF$Time <- Time[DF$id]
DF$event <- event[DF$id]
DF <- DF[DF$time <= DF$Time, ]

## Checking % of censoring
sum(event==0)/n

## Testing data:
Time <- pmin(trueTimes_test, Ctimes)
event <- as.numeric(trueTimes_test <= Ctimes) # event indicator
# we keep the longitudinal measurements before the event times
DF_test$Time <- Time[DF_test$id]
DF_test$event <- event[DF_test$id]
DF_test <- DF_test[DF_test$time <= DF_test$Time, ]

## Checking % of censoring
sum(event==0)/n_test

###################################################
#We fit the models and save the desired metrics
###################################################
DF.id <- DF[!duplicated(DF$id),]
DF_test.id <- DF_test[!duplicated(DF_test$id),] 
## Survival model
try(CoxFit <- coxph(Surv(Time, event) ~ sex, data = DF.id))

## Longitudinal models
try(LM1 <- lme(y1 ~ sex + time, data = DF, random = ~ time | id))
try(LM2 <- lme(y2 ~ time, data = DF, random = ~ time | id))
try(LM3 <- lme(y3 ~ treat + pow2(time), data = DF, random = ~ time | id))
#try(GLM4 <- mixed_model(y4 ~ sex + treat + time, data = DF,
#                       random = ~ time | id, family = binomial()))
#we fit the GLMM with just random intercepts
try(GLM4 <- mixed_model(y4 ~ sex + treat + time, data = DF,
                        random = ~ 1 | id, family = binomial()))

#Fitting the multivariate JM
try(multiJM <- jm(CoxFit, list(LM1, LM2, LM3, GLM4), time_var = "time",
                  n_iter = 13000L, n_burnin = 3000L, n_thin = 5L))


uniJM <- jm(CoxFit, GLM4, time_var = "time", 
            n_iter = 13000L, n_burnin = 3000L, n_thin = 5L)

## OBS: by adding which_independent = "all" we say the model matrix of var-cov
## is block diagonal, and it simplifies the problem

try(rhats <- numeric())
try(for(i in 1:(length(multiJM$statistics$Rhat)-2)){
  rhats <- c(rhats, multiJM$statistics$Rhat[[i]][,1]) 
})
try(rhats <- as.numeric(rhats))
rhats

# It seems that those bs_gammas (the ones related with B-splines used to
# fit the besaline hazard) are too big (more than 2). Could be that we are
# assuming ctt baseline hazard function?

## We compute the metrics, first IBS:
t0 <- 4.5
dt <- 1
try(brier_score_multi_train <- tvBrier(multiJM, newdata = DF, Tstart = t0, Dt = dt, 
                                      integrated = TRUE))
try(brier_score_multi_test <- tvBrier(multiJM, newdata = DF_test, Tstart = t0, Dt = dt, 
                                      integrated = TRUE))

##with IPCW weigths
try(brier_score_multi_train_IPCW <- tvBrier(multiJM, newdata = DF, Tstart = t0, Dt = dt, 
                                       integrated = TRUE, type_weights = "IPCW"))
try(brier_score_multi_test_IPCW <- tvBrier(multiJM, newdata = DF_test, Tstart = t0, Dt = dt, 
                                      integrated = TRUE, type_weights = "IPCW"))






####################################################################
####################################################################
## We do several replications
####################################################################
####################################################################
repl <- 5
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

set.seed(1234)
sigmaa <- matrix(c(runif(1,0.25,0.75), runif(6, -0.0005, 0.0005), 
                   runif(1,0,0.25), runif(6, -0.0005, 0.0005),
                   runif(1,0.25,0.75), runif(6, -0.0005, 0.0005),
                   runif(1,0,0.25), runif(6, -0.0005, 0.0005),
                   runif(1,0.25,0.75), runif(6, -0.0005, 0.0005),
                   runif(1,0,0.015))
                 ,nrow = 6)
treat <- rep(as.factor(sample(c("A", "B"), size = n, replace = TRUE))
             ,each = K)

for(count in 1:repl){
  n <- 250 # number of subjects
  n_test <- 250
  K <- 10 # number of measurements per subject
  t_max <- 20 # maximum follow-up time
  
  DF <- data.frame(id = rep(seq_len(n), each = K),
                   sex = rep(gl(2, n/2, labels = c("male", "female")), each = K),
                   treat = treat,
                   time = c(replicate(n, c(0, 1:(K-1)))))
  DF_test <- data.frame(id = rep(seq_len(n_test), each = K),
                        sex = rep(gl(2, n_test/2, labels = c("male", "female")), each = K),
                        treat = treat,
                        time = c(replicate(n, c(0, 1:(K-1)))))
  
  #correlated random effects
  b <-  mvrnorm(n = n, mu = rep(0,6), Sigma = sigmaa%*%t(sigmaa))
  b_test <-  mvrnorm(n = n, mu = rep(0,6), Sigma = sigmaa%*%t(sigmaa))
  
  ## Now we define the model
  
  ##Train:
  X1 <- model.matrix(~ sex + time, data = DF)
  Z1 <- model.matrix(~ time, data = DF)
  
  ##Test:
  X1_test <- model.matrix(~ sex + time, data = DF_test)
  Z1_test <- model.matrix(~ time, data = DF_test)
  
  ##Paramaters specification:
  #b1 <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(0.067, 0.0025, 0.0000178, 0.0008988), 
  #                                                 nrow = 2, ncol = 2))
  b1 <- b[,c(1,2)]
  betas1 <- c(0.123, 0.098765, 1.5)
  sigma1 <- 0.125 #error term sd
  
  ## Longitudinal profile simulation:
  eta_y1 <- as.vector(X1 %*% betas1 + rowSums(Z1 * b1[DF$id, ]))
  DF$y1 <- rnorm(n * K, mean = eta_y1, sd = sigma1)
  
  
  #b1_test <- mvrnorm(n = n_test, mu = c(0,0), Sigma = matrix(c(0.067, 0.0025, 0.0000178, 0.0008988), 
  #                                                           nrow = 2, ncol = 2))
  b1_test <- b_test[,c(1,2)]
  eta_y1_test <- as.vector(X1_test %*% betas1 + rowSums(Z1_test * b1_test[DF_test$id, ]))
  DF_test$y1 <- rnorm(n * K, mean = eta_y1_test, sd = sigma1)
  
  X2 <- model.matrix(~ time, data = DF)
  Z2 <- model.matrix(~ time, data = DF)
  X2_test <- model.matrix(~ time, data = DF_test)
  Z2_test <- model.matrix(~ time, data = DF_test)
  
  ##Paramaters specification:
  #b2 <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(0.067, 0.00025, 0.00178, 0.0008988), 
  #                                                 nrow = 2, ncol = 2))
  b2 <- b[, c(3,4)]
  betas2 <- c(15.641, -1.25918)
  sigma2 <- 0.125 #error term sd
  
  ##Longitudinal profile simulation
  eta_y2 <- as.vector(X2 %*% betas2 + rowSums(Z2 * b2[DF$id, ]))
  DF$y2 <- rnorm(n * K, mean = eta_y2, sd = sigma2)
  
  #b2_test <- mvrnorm(n = n_test, mu = c(0,0), Sigma = matrix(c(0.067, 0.00025, 0.00178, 0.0008988), 
  #                                                           nrow = 2, ncol = 2))
  b2_test <- b_test[, c(3,4)]
  eta_y2_test <- as.vector(X2_test %*% betas2 + rowSums(Z2_test * b2_test[DF_test$id, ]))
  DF_test$y2 <- rnorm(n * K, mean = eta_y2_test, sd = sigma2)
  
  #We define a function to compute the square of a given number
  pow2 <- function(x){
    x^2
  }
  
  X3 <- model.matrix(~ treat + pow2(time), data = DF)
  Z3 <- model.matrix(~ pow2(time), data = DF)
  X3_test <- model.matrix(~ treat + pow2(time), data = DF_test)
  Z3_test <- model.matrix(~ pow2(time), data = DF_test)
  
  ##Paramaters specification:
  #b3 <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(0.009123, 0.000025, 0.0005716, 0.00008988), 
  #                                                 nrow = 2, ncol = 2))
  b3 <- b[, c(5,6)]
  betas3 <- c(1.314159, 0.125,0.08181)
  sigma3 <- 0.125 #error term sd
  
  ##Longitudinal profile simulation
  eta_y3 <- as.vector(X3 %*% betas3 + rowSums(Z3 * b3[DF$id, ]))
  DF$y3 <- rnorm(n * K, mean = eta_y3, sd = sigma2)
  
  
  #b3_test <- mvrnorm(n = n_test, mu = c(0,0), Sigma = matrix(c(0.009123, 0.000025, 0.0005716, 0.00008988), 
  #                                                           nrow = 2, ncol = 2))
  b3_test <- b_test[, c(5,6)]
  eta_y3_test <- as.vector(X3_test %*% betas3 + rowSums(Z3_test * b3_test[DF_test$id, ]))
  DF_test$y3 <- rnorm(n * K, mean = eta_y3_test, sd = sigma3)
  
  ##################################
  ### Time-to-event simulation
  ##################################
  
  # We assume baseline hazard to be ctt (exponential distr). We use the inverse 
  # transform sampling method to simulate survival times.
  
  #parameters specification:
  shape_wb <- 1.2 # we assume exp distr and ctt baseline hazard
  #alpha <- c(0.1, -0.2, 0.312) # association coefficients
  alpha <- c(0.0879, -0.43, 0.2)
  #alpha <- 2*alpha
  gammas <- c("(Intercept)" = 0.25, "sex" = -0.5)
  #model matrix and linear predictor:
  W <- model.matrix(~ sex, data = DF[!duplicated(DF$id), ])
  treat_ind <- DF[!duplicated(DF$id), ]$treat
  eta_t <- as.vector(W %*% gammas)
  
  invS <- function (t, i) {
    # i denotes the subject
    sex_i <- W[i, 2L]
    treat_i <- treat_ind[i]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h <- function (s) {
      X1_at_s <- cbind(1, sex_i, s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, treat_i, pow2(s))
      Z3_at_s <- cbind(1,pow2(s))
      # the linear predictor from the mixed model evaluated at time s
      f1 <- as.vector(X1_at_s %*% betas1 +
                        rowSums(Z1_at_s * b1[rep(i, nrow(Z1_at_s)), ]))
      f2 <- as.vector(X2_at_s %*% betas2 +
                        rowSums(Z2_at_s * b2[rep(i, nrow(Z2_at_s)), ]))
      f3 <- as.vector(X3_at_s %*% betas3 +
                        rowSums(Z3_at_s * b3[rep(i, nrow(Z3_at_s)), ]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] + f1*alpha[1] + 
            f2*alpha[2] + f3*alpha[3])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h, lower = 0, upper = t)$value + log(u[i])
  }
  
  # we simulate the event times
  u <- runif(n)
  trueTimes <- numeric(n)
  for (i in seq_len(n)) {
    Up <- 100
    Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  
  
  ## Testing data:
  #model matrix and linear predictor:
  W_test <- model.matrix(~ sex, data = DF_test[!duplicated(DF_test$id), ])
  treat_ind_test <- DF_test[!duplicated(DF_test$id), ]$treat
  eta_t_test <- as.vector(W_test %*% gammas)
  
  invS <- function (t, i) {
    # i denotes the subject
    sex_i <- W_test[i, 2L]
    treat_i <- treat_ind_test[i]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h <- function (s) {
      X1_at_s <- cbind(1, sex_i, s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, treat_i, pow2(s))
      Z3_at_s <- cbind(1,pow2(s))
      # the linear predictor from the mixed model evaluated at time s
      f1 <- as.vector(X1_at_s %*% betas1 +
                        rowSums(Z1_at_s * b1_test[rep(i, nrow(Z1_at_s)), ]))
      f2 <- as.vector(X2_at_s %*% betas2 +
                        rowSums(Z2_at_s * b2_test[rep(i, nrow(Z2_at_s)), ]))
      f3 <- as.vector(X3_at_s %*% betas3 +
                        rowSums(Z3_at_s * b3_test[rep(i, nrow(Z3_at_s)), ]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t_test[i] + f1*alpha[1] + 
            f2*alpha[2] + f3*alpha[3])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h, lower = 0, upper = t)$value + log(u[i])
  }
  
  # we simulate the event times
  u <- runif(n)
  trueTimes_test <- numeric(n)
  for (i in seq_len(n)) {
    Up <- 100
    Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes_test[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  
  
  ##saving true times
  checkTimes[count] <- sum(trueTimes==150)/n
  
  #save Times in data frame
  Time <- trueTimes
  DF$Time <- Time[DF$id]
  
  ##saving true times
  checkTimes_test[count] <- sum(trueTimes_test==150)/n_test
  
  #save Times in data frame
  Time <- trueTimes_test
  DF_test$Time <- Time[DF_test$id]
  
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
  
  
  #Fitting the multivariate JM
  try(multiJM <- jm(CoxFit, list(LM1, LM2, LM3), time_var = "time",
                    n_iter = 13000L, n_burnin = 3000L, n_thin = 4L))
      
      
  # If we need to specify the knots we should put:
  # control = list(knots = list(c(4,5,6,7,8,9,10)))))
  
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
  t0 <- 4.5
  dt <- 1
  try(brier_score_multi_train <- tvBrier(multiJM, newdata = DF, Tstart = t0, Dt = dt, 
                                         integrated = TRUE))
  try(ibs_train[count] <- brier_score_multi_train$Brier)
  try(n_risk_train[count] <- brier_score_multi_train$nr)
  try(n_event_train[count] <- brier_score_multi_train$nint)
  try(n_cens_train[count] <- brier_score_multi_train$ncens)
  try(list_brier_model_train <- append(list_brier_model_train, 
                                 list(brier_score_multi_train$Brier_per_model)))
  try(list_w_model_train <- append(list_w_model_train, 
                             list(brier_score_multi_train$weights)))
  
  
  try(brier_score_multi_test <- tvBrier(multiJM, newdata = DF_test, Tstart = t0, Dt = dt, 
                                        integrated = TRUE))
  try(ibs_test[count] <- brier_score_multi_test$Brier)
  try(n_risk_test[count] <- brier_score_multi_test$nr)
  try(n_event_test[count] <- brier_score_multi_test$nint)
  try(n_cens_test[count] <- brier_score_multi_test$ncens)
  try(list_brier_model_test <- append(list_brier_model_test, 
                                 list(brier_score_multi_test$Brier_per_model)))
  try(list_w_model_test <- append(list_w_model_test, 
                             list(brier_score_multi_test$weights)))
  
  
  try(list_rhats_MJM <- append(list_rhats_MJM, list(rhats)))
  try(list_full_rhats_MJM <- append(list_full_rhats_MJM, 
                                    list(multiJM$statistics$Rhat)))
  
  if(count == repl){
    strr <- "rebuilding_simul_3outcomes_correl_typ1_21feb2025.RData"
    save(checkTimes_test, checkTimes, perc_cens_test, perc_cens_train,
         list_rhats_MJM, ibs_train, ibs_test, 
         n_risk_train, n_risk_test, n_event_train, n_event_test,
         n_cens_train, n_cens_test, list_w_model_train, list_w_model_test,
         list_brier_model_train, list_brier_model_test, 
         file=strr)
  }
  print(count)
}

#### results from above loop
load("rebuilding_simul_attempt3_typ1_20feb2025.RData")

df_comparing <- data.frame(IBS_train = ibs_train,
                           IBS_test = ibs_test,
                           perc_cens_train = perc_cens_train,
                           perc_cens_test = perc_cens_test,
                           n_event_train = n_event_train,
                           n_event_test = n_event_test,
                           n_risk_train = n_risk_train,
                           n_risk_test = n_risk_test)

mean(df_comparing$IBS_test - df_comparing$IBS_train)

list_rhats_MJM

list_full_rhats_MJM[[1]]


### Trace plots:

multiJM$statistics$Rhat

multiJM$statistics$Mean

ggtraceplot(multiJM, parm = "bs_gammas")

ggtraceplot(multiJM, parm = "alphas")







####################################################################
####################################################################
## We do several replications
## With L=4, the last outcome being a GLMM (binary longitudinal outcome)
####################################################################
####################################################################
repl <- 30
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



set.seed(1234)
sigmaa <- matrix(c(runif(1,0.25,0.75), runif(7, -0.0005, 0.0005), 
                   runif(1,0,0.15), runif(7, -0.0005, 0.0005),
                   runif(1,0.25,0.75), runif(7, -0.0005, 0.0005),
                   runif(1,0,0.15), runif(7, -0.0005, 0.0005),
                   runif(1,0.25,0.75), runif(7, -0.0005, 0.0005),
                   runif(1,0,0.025), runif(7, -0.0005, 0.0005),
                   runif(1,0,0.05))
                 ,nrow = 7)
treat <- rep(as.factor(sample(c("A", "B"), size = n, replace = TRUE))
             ,each = K)

for(count in 1:repl){
  n <- 250 # number of subjects
  n_test <- 250
  K <- 10 # number of measurements per subject
  t_max <- 20 # maximum follow-up time
  
  DF <- data.frame(id = rep(seq_len(n), each = K),
                   sex = rep(gl(2, n/2, labels = c("male", "female")), each = K),
                   treat = treat,
                   time = c(replicate(n, c(0, 1:(K-1)))))
  DF_test <- data.frame(id = rep(seq_len(n_test), each = K),
                        sex = rep(gl(2, n_test/2, labels = c("male", "female")), each = K),
                        treat = treat,
                        time = c(replicate(n, c(0, 1:(K-1)))))
  
  #correlated random effects
  b <-  mvrnorm(n = n, mu = rep(0,7), Sigma = sigmaa%*%t(sigmaa))
  b_test <-  mvrnorm(n = n, mu = rep(0,7), Sigma = sigmaa%*%t(sigmaa))
  
  ## Now we define the model
  
  ##Train:
  X1 <- model.matrix(~ sex + time, data = DF)
  Z1 <- model.matrix(~ time, data = DF)
  
  ##Test:
  X1_test <- model.matrix(~ sex + time, data = DF_test)
  Z1_test <- model.matrix(~ time, data = DF_test)
  
  ##Paramaters specification:
  #b1 <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(0.067, 0.0025, 0.0000178, 0.0008988), 
  #                                                 nrow = 2, ncol = 2))
  b1 <- b[,c(1,2)]
  betas1 <- c(0.123, 0.098765, 1.5)
  sigma1 <- 0.125 #error term sd
  
  ## Longitudinal profile simulation:
  eta_y1 <- as.vector(X1 %*% betas1 + rowSums(Z1 * b1[DF$id, ]))
  DF$y1 <- rnorm(n * K, mean = eta_y1, sd = sigma1)
  
  
  #b1_test <- mvrnorm(n = n_test, mu = c(0,0), Sigma = matrix(c(0.067, 0.0025, 0.0000178, 0.0008988), 
  #                                                           nrow = 2, ncol = 2))
  b1_test <- b_test[,c(1,2)]
  eta_y1_test <- as.vector(X1_test %*% betas1 + rowSums(Z1_test * b1_test[DF_test$id, ]))
  DF_test$y1 <- rnorm(n * K, mean = eta_y1_test, sd = sigma1)
  
  X2 <- model.matrix(~ time, data = DF)
  Z2 <- model.matrix(~ time, data = DF)
  X2_test <- model.matrix(~ time, data = DF_test)
  Z2_test <- model.matrix(~ time, data = DF_test)
  
  ##Paramaters specification:
  #b2 <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(0.067, 0.00025, 0.00178, 0.0008988), 
  #                                                 nrow = 2, ncol = 2))
  b2 <- b[, c(3,4)]
  betas2 <- c(15.641, -1.25918)
  sigma2 <- 0.125 #error term sd
  
  ##Longitudinal profile simulation
  eta_y2 <- as.vector(X2 %*% betas2 + rowSums(Z2 * b2[DF$id, ]))
  DF$y2 <- rnorm(n * K, mean = eta_y2, sd = sigma2)
  
  #b2_test <- mvrnorm(n = n_test, mu = c(0,0), Sigma = matrix(c(0.067, 0.00025, 0.00178, 0.0008988), 
  #                                                           nrow = 2, ncol = 2))
  b2_test <- b_test[, c(3,4)]
  eta_y2_test <- as.vector(X2_test %*% betas2 + rowSums(Z2_test * b2_test[DF_test$id, ]))
  DF_test$y2 <- rnorm(n * K, mean = eta_y2_test, sd = sigma2)
  
  #We define a function to compute the square of a given number
  pow2 <- function(x){
    x^2
  }
  
  X3 <- model.matrix(~ treat + pow2(time), data = DF)
  Z3 <- model.matrix(~ pow2(time), data = DF)
  X3_test <- model.matrix(~ treat + pow2(time), data = DF_test)
  Z3_test <- model.matrix(~ pow2(time), data = DF_test)
  
  ##Paramaters specification:
  #b3 <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(0.009123, 0.000025, 0.0005716, 0.00008988), 
  #                                                 nrow = 2, ncol = 2))
  b3 <- b[, c(5,6)]
  betas3 <- c(1.314159, 0.125,0.08181)
  sigma3 <- 0.125 #error term sd
  
  ##Longitudinal profile simulation
  eta_y3 <- as.vector(X3 %*% betas3 + rowSums(Z3 * b3[DF$id, ]))
  DF$y3 <- rnorm(n * K, mean = eta_y3, sd = sigma2)
  
  
  #b3_test <- mvrnorm(n = n_test, mu = c(0,0), Sigma = matrix(c(0.009123, 0.000025, 0.0005716, 0.00008988), 
  #                                                           nrow = 2, ncol = 2))
  b3_test <- b_test[, c(5,6)]
  eta_y3_test <- as.vector(X3_test %*% betas3 + rowSums(Z3_test * b3_test[DF_test$id, ]))
  DF_test$y3 <- rnorm(n * K, mean = eta_y3_test, sd = sigma3)
  
  ## y4, glmm to model binary long outcome
  X4 <- model.matrix(~ sex + treat + time, data = DF)
  Z4 <- model.matrix(~ 1, data = DF)
  X4_test <- model.matrix(~ sex + treat + time, data = DF_test)
  Z4_test <- model.matrix(~ 1, data = DF_test)
  
  ##Paramaters specification:
  #b4 <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(0.009123, 0.000025, 0.0005716, 0.00008988), 
  #                                                 nrow = 2, ncol = 2))
  b4 <- b[, 7]
  betas4 <- c(-2.33, -0.07152, 0.42313, 1)
  sigma4 <- 0.125 #error term sd
  
  ##Longitudinal profile simulation
  eta_y4 <- as.vector(X4 %*% betas4 + rowSums(Z4 * b4[DF$id]))
  DF$y4 <- rbinom(n*K, size = 1, prob = plogis(eta_y4))
  
  
  #b4_test <- mvrnorm(n = n_test, mu = c(0,0), Sigma = matrix(c(0.009123, 0.000025, 0.0005716, 0.00008988), 
  #                                                           nrow = 2, ncol = 2))
  b4_test <- b_test[, 7]
  eta_y4_test <- as.vector(X4_test %*% betas4 + rowSums(Z4_test * b4_test[DF_test$id]))
  DF_test$y4 <- rbinom(n_test*K, size = 1, prob = plogis(eta_y4_test))
  
  ##################################
  ### Time-to-event simulation
  ##################################
  
  # We assume baseline hazard to be ctt (exponential distr). We use the inverse 
  # transform sampling method to simulate survival times.
  
  #parameters specification:
  shape_wb <- 1.2 # we assume exp distr and ctt baseline hazard
  #alpha <- c(0.1, -0.2, 0.312) # association coefficients
  #alpha <- c(0.0879, -0.43, 0.2, -0.001)
  alpha <- c(0.0879, -0.43, 0.2, -0.0001)
  ##obs: when alpha=0 we know we have convergence
  #alpha <- 2*alpha
  gammas <- c("(Intercept)" = 0.25, "sex" = -0.5)
  #gammas <- c("(Intercept)" = 0.15, "sex" = -0.25)
  #model matrix and linear predictor:
  W <- model.matrix(~ sex, data = DF[!duplicated(DF$id), ])
  treat_ind <- DF[!duplicated(DF$id), ]$treat
  eta_t <- as.vector(W %*% gammas)
  
  invS <- function (t, i) {
    # i denotes the subject
    sex_i <- W[i, 2L]
    treat_i <- treat_ind[i]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h <- function (s) {
      X1_at_s <- cbind(1, sex_i, s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, treat_i, pow2(s))
      Z3_at_s <- cbind(1,pow2(s))
      X4_at_s <- cbind(1, sex_i, treat_i, s)
      Z4_at_s <- cbind(1)
      # the linear predictor from the mixed model evaluated at time s
      f1 <- as.vector(X1_at_s %*% betas1 +
                        rowSums(Z1_at_s * b1[rep(i, nrow(Z1_at_s)), ]))
      f2 <- as.vector(X2_at_s %*% betas2 +
                        rowSums(Z2_at_s * b2[rep(i, nrow(Z2_at_s)), ]))
      f3 <- as.vector(X3_at_s %*% betas3 +
                        rowSums(Z3_at_s * b3[rep(i, nrow(Z3_at_s)), ]))
      f4 <- as.vector(X4_at_s %*% betas4 +
                        rowSums(Z4_at_s * b4[rep(i, nrow(Z4_at_s)) ]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] + f1*alpha[1] + 
            f2*alpha[2] + f3*alpha[3] + f4*alpha[4])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h, lower = 0, upper = t)$value + log(u[i])
  }
  
  # we simulate the event times
  u <- runif(n)
  trueTimes <- numeric(n)
  for (i in seq_len(n)) {
    Up <- 100
    Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  
  
  ## Testing data:
  #model matrix and linear predictor:
  W_test <- model.matrix(~ sex, data = DF_test[!duplicated(DF_test$id), ])
  treat_ind_test <- DF_test[!duplicated(DF_test$id), ]$treat
  eta_t_test <- as.vector(W_test %*% gammas)
  
  invS <- function (t, i) {
    # i denotes the subject
    sex_i <- W_test[i, 2L]
    treat_i <- treat_ind_test[i]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h <- function (s) {
      X1_at_s <- cbind(1, sex_i, s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, treat_i, pow2(s))
      Z3_at_s <- cbind(1,pow2(s))
      X4_at_s <- cbind(1, sex_i, treat_i, s)
      Z4_at_s <- cbind(1)
      # the linear predictor from the mixed model evaluated at time s
      f1 <- as.vector(X1_at_s %*% betas1 +
                        rowSums(Z1_at_s * b1_test[rep(i, nrow(Z1_at_s)), ]))
      f2 <- as.vector(X2_at_s %*% betas2 +
                        rowSums(Z2_at_s * b2_test[rep(i, nrow(Z2_at_s)), ]))
      f3 <- as.vector(X3_at_s %*% betas3 +
                        rowSums(Z3_at_s * b3_test[rep(i, nrow(Z3_at_s)), ]))
      f4 <- as.vector(X4_at_s %*% betas4 +
                        rowSums(Z4_at_s * b4_test[rep(i, nrow(Z4_at_s))]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t_test[i] + f1*alpha[1] + 
            f2*alpha[2] + f3*alpha[3] + f4*alpha[4])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h, lower = 0, upper = t)$value + log(u[i])
  }
  
  # we simulate the event times
  u <- runif(n)
  trueTimes_test <- numeric(n)
  for (i in seq_len(n)) {
    Up <- 100
    Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes_test[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  
  ##saving true times
  checkTimes[count] <- sum(trueTimes==150)/n
  
  #save Times in data frame
  Time <- trueTimes
  DF$Time <- Time[DF$id]
  
  ##saving true times
  checkTimes_test[count] <- sum(trueTimes_test==150)/n_test
  
  #save Times in data frame
  Time <- trueTimes_test
  DF_test$Time <- Time[DF_test$id]
  
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
  t0 <- 4.5
  dt <- 1
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


  
  

  
  if(count == repl){
    strr <- "ibsipcw_epce_simul_4outcomes_correl_typ1_SL_06mar2025.RData"
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
  print(count)
}

## We have convergene problems for the l=4 longitudinal outcome
## we should analyze the univariate joint model for this outcome
## and try to reach good convergence there before going to the 
## multivariate joint model

## Things for fixing it:
## 1) go just for random intercepts
## 2) make alpha even smaller (for alpha=0 we have convergence because
##.   we only take into account the rest of longitudinal outcomes)

df_comparing_ibs <- data.frame(IBS_mJM_train = ibs_train,
                           IBS_mJM_test = ibs_test,
                           IBS_eSL_train = eSL_cv_IBS,
                           IBS_eSL_test = eSL_test_IBS,
                           IBS_dSL_train = dSL_cv_IBS,
                           IBS_dSL_test = dSL_test_IBS,
                           perc_cens_train = perc_cens_train,
                           perc_cens_test = perc_cens_test,
                           n_event_train = n_event_train,
                           n_event_test = n_event_test,
                           n_risk_train = n_risk_train,
                           n_risk_test = n_risk_test)

df_comparing_ibs <- df_comparing_ibs[1:25,]

df_comparing_epce <- data.frame(EPCE_mJM_train = epce_train,
                           EPCE_mJM_test = epce_test,
                           EPCE_eSL_train = eSL_cv_EPCE,
                           EPCE_eSL_test = eSL_test_EPCE,
                           EPCE_dSL_train = dSL_cv_EPCE,
                           EPCE_dSL_test = dSL_test_EPCE,
                           perc_cens_train = perc_cens_train,
                           perc_cens_test = perc_cens_test,
                           n_event_train = n_event_train,
                           n_event_test = n_event_test,
                           n_risk_train = n_risk_train,
                           n_risk_test = n_risk_test)

df_comparing_epce <- df_comparing_epce[1:25,]

df_comparing3_pooled <- df_comparing3[-c(15,21,30),]

full <- rbind(df_comparing3_pooled, df_comparing2)


mean(df_comparing2$IBS_mJM_test - df_comparing2$IBS_mJM_train)
mean(df_comparing2$IBS_eSL_test - df_comparing2$IBS_eSL_train)
mean(df_comparing2$IBS_dSL_test - df_comparing2$IBS_dSL_train)

list_rhats_MJM

list_full_rhats_MJM[[1]]

#boxplot IBS
df <- with(df_comparing2, data.frame(IBS_mJM_train, IBS_mJM_test, IBS_eSL_train, IBS_eSL_test,  
                 IBS_dSL_train, IBS_dSL_test))
library(dplyr)
library(tidyr)
library(ggplot2)


#boxplot IBS
df <- with(df_comparing_ibs, data.frame("dSL CV Data" = IBS_dSL_train,
                 "eSL CV Data" = IBS_eSL_train,
                 "Multivariate model CV Data" = IBS_mJM_train,
                 "dSL Test Data" = IBS_dSL_test,
                 "eSL Test Data" = IBS_eSL_test,
                 "Multivariate model Test Data" = IBS_mJM_test))
df_long <- df %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

df_long$variable <- factor(df_long$variable,
                           levels = c("dSL.CV.Data",
                                      "eSL.CV.Data",
                                      "Multivariate.model.CV.Data",
                                      "dSL.Test.Data",
                                      "eSL.Test.Data",
                                      "Multivariate.model.Test.Data"))

ploty_2 <- ggplot(df_long, aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue") +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 20)) +
  xlab("") + ylab("IBS") +
  scale_x_discrete(labels = c("dSL.CV.Data" = "dSL CV Data",
                              "eSL.CV.Data" = "eSL CV Data",
                              "Multivariate.model.CV.Data" = "Multivariate JM\n CV Data",
                              "dSL.Test.Data" = "dSL Test Data",
                              "eSL.Test.Data" = "eSL Test Data",
                              "Multivariate.model.Test.Data" = "Multivariate JM\n Test Data"))+
  geom_vline(xintercept = (3+4)/2, linetype = "dashed", color = "grey", linewidth = 0.5)+
  ggtitle("Integrated Brier Score: (4.5,5.5]")


## EPCE plot

df2 <- with(df_comparing_epce, data.frame("dSL CV Data" = EPCE_dSL_train,
                                         "eSL CV Data" = EPCE_eSL_train,
                                         "Multivariate model CV Data" = EPCE_mJM_train,
                                         "dSL Test Data" = EPCE_dSL_test,
                                         "eSL Test Data" = EPCE_eSL_test,
                                         "Multivariate model Test Data" = EPCE_mJM_test))
df_long2 <- df2 %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

df_long2$variable <- factor(df_long2$variable,
                            levels = c("dSL.CV.Data",
                                       "eSL.CV.Data",
                                       "Multivariate.model.CV.Data",
                                       "dSL.Test.Data",
                                       "eSL.Test.Data",
                                       "Multivariate.model.Test.Data"))

plottt_2 <- ggplot(df_long2, aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue") +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 20)) +
  xlab("") + ylab("EPCE") +
  scale_x_discrete(labels = c("dSL.CV.Data" = "dSL CV Data",
                              "eSL.CV.Data" = "eSL CV Data",
                              "Multivariate.model.CV.Data" = "Multivariate JM\n CV Data",
                              "dSL.Test.Data" = "dSL Test Data",
                              "eSL.Test.Data" = "eSL Test Data",
                              "Multivariate.model.Test.Data" = "Multivariate JM\n Test Data"))+
  geom_vline(xintercept = (3+4)/2, linetype = "dashed", color = "grey", linewidth = 1)+
  ggtitle("EPCE: (4.5,5.5]")



################################################################################
### Now with random censoring
################################################################################

repl <- 20
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



set.seed(1234)
sigmaa <- matrix(c(runif(1,0.25,0.75), runif(7, -0.0005, 0.0005), 
                   runif(1,0,0.15), runif(7, -0.0005, 0.0005),
                   runif(1,0.25,0.75), runif(7, -0.0005, 0.0005),
                   runif(1,0,0.15), runif(7, -0.0005, 0.0005),
                   runif(1,0.25,0.75), runif(7, -0.0005, 0.0005),
                   runif(1,0,0.025), runif(7, -0.0005, 0.0005),
                   runif(1,0,0.05))
                 ,nrow = 7)
treat <- rep(as.factor(sample(c("A", "B"), size = n, replace = TRUE))
             ,each = K)


for(count in 1:repl){
  n <- 250 # number of subjects
  n_test <- 250
  K <- 10 # number of measurements per subject
  t_max <- 20 # maximum follow-up time
  
  DF <- data.frame(id = rep(seq_len(n), each = K),
                   sex = rep(gl(2, n/2, labels = c("male", "female")), each = K),
                   treat = treat,
                   time = c(replicate(n, c(0, 1:(K-1)))))
  DF_test <- data.frame(id = rep(seq_len(n_test), each = K),
                        sex = rep(gl(2, n_test/2, labels = c("male", "female")), each = K),
                        treat = treat,
                        time = c(replicate(n, c(0, 1:(K-1)))))
  
  #correlated random effects
  b <-  mvrnorm(n = n, mu = rep(0,7), Sigma = sigmaa%*%t(sigmaa))
  b_test <-  mvrnorm(n = n, mu = rep(0,7), Sigma = sigmaa%*%t(sigmaa))
  
  ## Now we define the model
  
  ##Train:
  X1 <- model.matrix(~ sex + time, data = DF)
  Z1 <- model.matrix(~ time, data = DF)
  
  ##Test:
  X1_test <- model.matrix(~ sex + time, data = DF_test)
  Z1_test <- model.matrix(~ time, data = DF_test)
  
  ##Paramaters specification:
  #b1 <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(0.067, 0.0025, 0.0000178, 0.0008988), 
  #                                                 nrow = 2, ncol = 2))
  b1 <- b[,c(1,2)]
  betas1 <- c(0.123, 0.098765, 1.5)
  sigma1 <- 0.125 #error term sd
  
  ## Longitudinal profile simulation:
  eta_y1 <- as.vector(X1 %*% betas1 + rowSums(Z1 * b1[DF$id, ]))
  DF$y1 <- rnorm(n * K, mean = eta_y1, sd = sigma1)
  
  
  #b1_test <- mvrnorm(n = n_test, mu = c(0,0), Sigma = matrix(c(0.067, 0.0025, 0.0000178, 0.0008988), 
  #                                                           nrow = 2, ncol = 2))
  b1_test <- b_test[,c(1,2)]
  eta_y1_test <- as.vector(X1_test %*% betas1 + rowSums(Z1_test * b1_test[DF_test$id, ]))
  DF_test$y1 <- rnorm(n * K, mean = eta_y1_test, sd = sigma1)
  
  X2 <- model.matrix(~ time, data = DF)
  Z2 <- model.matrix(~ time, data = DF)
  X2_test <- model.matrix(~ time, data = DF_test)
  Z2_test <- model.matrix(~ time, data = DF_test)
  
  ##Paramaters specification:
  #b2 <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(0.067, 0.00025, 0.00178, 0.0008988), 
  #                                                 nrow = 2, ncol = 2))
  b2 <- b[, c(3,4)]
  betas2 <- c(15.641, -1.25918)
  sigma2 <- 0.125 #error term sd
  
  ##Longitudinal profile simulation
  eta_y2 <- as.vector(X2 %*% betas2 + rowSums(Z2 * b2[DF$id, ]))
  DF$y2 <- rnorm(n * K, mean = eta_y2, sd = sigma2)
  
  #b2_test <- mvrnorm(n = n_test, mu = c(0,0), Sigma = matrix(c(0.067, 0.00025, 0.00178, 0.0008988), 
  #                                                           nrow = 2, ncol = 2))
  b2_test <- b_test[, c(3,4)]
  eta_y2_test <- as.vector(X2_test %*% betas2 + rowSums(Z2_test * b2_test[DF_test$id, ]))
  DF_test$y2 <- rnorm(n * K, mean = eta_y2_test, sd = sigma2)
  
  #We define a function to compute the square of a given number
  pow2 <- function(x){
    x^2
  }
  
  X3 <- model.matrix(~ treat + pow2(time), data = DF)
  Z3 <- model.matrix(~ pow2(time), data = DF)
  X3_test <- model.matrix(~ treat + pow2(time), data = DF_test)
  Z3_test <- model.matrix(~ pow2(time), data = DF_test)
  
  ##Paramaters specification:
  #b3 <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(0.009123, 0.000025, 0.0005716, 0.00008988), 
  #                                                 nrow = 2, ncol = 2))
  b3 <- b[, c(5,6)]
  betas3 <- c(1.314159, 0.125,0.08181)
  sigma3 <- 0.125 #error term sd
  
  ##Longitudinal profile simulation
  eta_y3 <- as.vector(X3 %*% betas3 + rowSums(Z3 * b3[DF$id, ]))
  DF$y3 <- rnorm(n * K, mean = eta_y3, sd = sigma2)
  
  
  #b3_test <- mvrnorm(n = n_test, mu = c(0,0), Sigma = matrix(c(0.009123, 0.000025, 0.0005716, 0.00008988), 
  #                                                           nrow = 2, ncol = 2))
  b3_test <- b_test[, c(5,6)]
  eta_y3_test <- as.vector(X3_test %*% betas3 + rowSums(Z3_test * b3_test[DF_test$id, ]))
  DF_test$y3 <- rnorm(n * K, mean = eta_y3_test, sd = sigma3)
  
  ## y4, glmm to model binary long outcome
  X4 <- model.matrix(~ sex + treat + time, data = DF)
  Z4 <- model.matrix(~ 1, data = DF)
  X4_test <- model.matrix(~ sex + treat + time, data = DF_test)
  Z4_test <- model.matrix(~ 1, data = DF_test)
  
  ##Paramaters specification:
  #b4 <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(0.009123, 0.000025, 0.0005716, 0.00008988), 
  #                                                 nrow = 2, ncol = 2))
  b4 <- b[, 7]
  betas4 <- c(-2.33, -0.07152, 0.42313, 1)
  sigma4 <- 0.125 #error term sd
  
  ##Longitudinal profile simulation
  eta_y4 <- as.vector(X4 %*% betas4 + rowSums(Z4 * b4[DF$id]))
  DF$y4 <- rbinom(n*K, size = 1, prob = plogis(eta_y4))
  
  
  #b4_test <- mvrnorm(n = n_test, mu = c(0,0), Sigma = matrix(c(0.009123, 0.000025, 0.0005716, 0.00008988), 
  #                                                           nrow = 2, ncol = 2))
  b4_test <- b_test[, 7]
  eta_y4_test <- as.vector(X4_test %*% betas4 + rowSums(Z4_test * b4_test[DF_test$id]))
  DF_test$y4 <- rbinom(n_test*K, size = 1, prob = plogis(eta_y4_test))
  
  ##################################
  ### Time-to-event simulation
  ##################################
  
  # We assume baseline hazard to be ctt (exponential distr). We use the inverse 
  # transform sampling method to simulate survival times.
  
  #parameters specification:
  shape_wb <- 1.2 # we assume exp distr and ctt baseline hazard
  #alpha <- c(0.1, -0.2, 0.312) # association coefficients
  #alpha <- c(0.0879, -0.43, 0.2, -0.001)
  alpha <- c(0.0879, -0.43, 0.2, -0.0001)
  ##obs: when alpha=0 we know we have convergence
  #alpha <- 2*alpha
  gammas <- c("(Intercept)" = 0.25, "sex" = -0.5)
  #gammas <- c("(Intercept)" = 0.15, "sex" = -0.25)
  #model matrix and linear predictor:
  W <- model.matrix(~ sex, data = DF[!duplicated(DF$id), ])
  treat_ind <- DF[!duplicated(DF$id), ]$treat
  eta_t <- as.vector(W %*% gammas)
  
  invS <- function (t, i) {
    # i denotes the subject
    sex_i <- W[i, 2L]
    treat_i <- treat_ind[i]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h <- function (s) {
      X1_at_s <- cbind(1, sex_i, s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, treat_i, pow2(s))
      Z3_at_s <- cbind(1,pow2(s))
      X4_at_s <- cbind(1, sex_i, treat_i, s)
      Z4_at_s <- cbind(1)
      # the linear predictor from the mixed model evaluated at time s
      f1 <- as.vector(X1_at_s %*% betas1 +
                        rowSums(Z1_at_s * b1[rep(i, nrow(Z1_at_s)), ]))
      f2 <- as.vector(X2_at_s %*% betas2 +
                        rowSums(Z2_at_s * b2[rep(i, nrow(Z2_at_s)), ]))
      f3 <- as.vector(X3_at_s %*% betas3 +
                        rowSums(Z3_at_s * b3[rep(i, nrow(Z3_at_s)), ]))
      f4 <- as.vector(X4_at_s %*% betas4 +
                        rowSums(Z4_at_s * b4[rep(i, nrow(Z4_at_s)) ]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] + f1*alpha[1] + 
            f2*alpha[2] + f3*alpha[3] + f4*alpha[4])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h, lower = 0, upper = t)$value + log(u[i])
  }
  
  # we simulate the event times
  u <- runif(n)
  trueTimes <- numeric(n)
  for (i in seq_len(n)) {
    Up <- 100
    Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  
  
  ## Testing data:
  #model matrix and linear predictor:
  W_test <- model.matrix(~ sex, data = DF_test[!duplicated(DF_test$id), ])
  treat_ind_test <- DF_test[!duplicated(DF_test$id), ]$treat
  eta_t_test <- as.vector(W_test %*% gammas)
  
  invS <- function (t, i) {
    # i denotes the subject
    sex_i <- W_test[i, 2L]
    treat_i <- treat_ind_test[i]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h <- function (s) {
      X1_at_s <- cbind(1, sex_i, s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, treat_i, pow2(s))
      Z3_at_s <- cbind(1,pow2(s))
      X4_at_s <- cbind(1, sex_i, treat_i, s)
      Z4_at_s <- cbind(1)
      # the linear predictor from the mixed model evaluated at time s
      f1 <- as.vector(X1_at_s %*% betas1 +
                        rowSums(Z1_at_s * b1_test[rep(i, nrow(Z1_at_s)), ]))
      f2 <- as.vector(X2_at_s %*% betas2 +
                        rowSums(Z2_at_s * b2_test[rep(i, nrow(Z2_at_s)), ]))
      f3 <- as.vector(X3_at_s %*% betas3 +
                        rowSums(Z3_at_s * b3_test[rep(i, nrow(Z3_at_s)), ]))
      f4 <- as.vector(X4_at_s %*% betas4 +
                        rowSums(Z4_at_s * b4_test[rep(i, nrow(Z4_at_s))]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t_test[i] + f1*alpha[1] + 
            f2*alpha[2] + f3*alpha[3] + f4*alpha[4])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h, lower = 0, upper = t)$value + log(u[i])
  }
  
  # we simulate the event times
  u <- runif(n)
  trueTimes_test <- numeric(n)
  for (i in seq_len(n)) {
    Up <- 100
    Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes_test[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  
  ##saving true times
  checkTimes[count] <- sum(trueTimes==150)/n
  
  #save Times in data frame
  Time <- trueTimes
  DF$Time <- Time[DF$id]
  
  ##saving true times
  checkTimes_test[count] <- sum(trueTimes_test==150)/n_test
  
  #save Times in data frame
  Time <- trueTimes_test
  DF_test$Time <- Time[DF_test$id]
  
  #simulating censoring; Random censoring
  Ctimes <- rexp(n, 1/18.5)
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
  
  
  
  
  
  
  if(count == repl){
    strr <- "ibsipcw_epce_simul_4outcomes_correl_randCens_SL_06mar2025.RData"
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
  print(count)
}


df_comparing_ibs_randcens2 <- data.frame(IBS_mJM_train = ibs_train,
                               IBS_mJM_test = ibs_test,
                               IBS_eSL_train = eSL_cv_IBS,
                               IBS_eSL_test = eSL_test_IBS,
                               IBS_dSL_train = dSL_cv_IBS,
                               IBS_dSL_test = dSL_test_IBS,
                               perc_cens_train = perc_cens_train,
                               perc_cens_test = perc_cens_test,
                               n_event_train = n_event_train,
                               n_event_test = n_event_test,
                               n_risk_train = n_risk_train,
                               n_risk_test = n_risk_test)



df_comparing_epce_randcens2 <- data.frame(EPCE_mJM_train = epce_train,
                                EPCE_mJM_test = epce_test,
                                EPCE_eSL_train = eSL_cv_EPCE,
                                EPCE_eSL_test = eSL_test_EPCE,
                                EPCE_dSL_train = dSL_cv_EPCE,
                                EPCE_dSL_test = dSL_test_EPCE,
                                perc_cens_train = perc_cens_train,
                                perc_cens_test = perc_cens_test,
                                n_event_train = n_event_train,
                                n_event_test = n_event_test,
                                n_risk_train = n_risk_train,
                                n_risk_test = n_risk_test)


##obs: it seems that the convergence of bs_gamma parameters is even worse!
## also alpha parameters are not converging in this case

df_comp_full_ibs <- rbind(df_comparing_ibs_randcens2, df_comparing_ibs_randcens)
df_comp_full_ibs <- df_comp_full_ibs[-25,]


library(dplyr)
library(tidyr)
library(ggplot2)


#boxplot IBS
df <- with(df_comp_full_ibs, data.frame("dSL CV Data" = IBS_dSL_train,
                                        "eSL CV Data" = IBS_eSL_train,
                                        "Multivariate model CV Data" = IBS_mJM_train,
                                        "dSL Test Data" = IBS_dSL_test,
                                        "eSL Test Data" = IBS_eSL_test,
                                        "Multivariate model Test Data" = IBS_mJM_test))
df_long <- df %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

df_long$variable <- factor(df_long$variable,
                           levels = c("dSL.CV.Data",
                                      "eSL.CV.Data",
                                      "Multivariate.model.CV.Data",
                                      "dSL.Test.Data",
                                      "eSL.Test.Data",
                                      "Multivariate.model.Test.Data"))

ploty_2 <- ggplot(df_long, aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue") +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 20)) +
  xlab("") + ylab("IBS") +
  scale_x_discrete(labels = c("dSL.CV.Data" = "dSL CV Data",
                              "eSL.CV.Data" = "eSL CV Data",
                              "Multivariate.model.CV.Data" = "Multivariate JM\n CV Data",
                              "dSL.Test.Data" = "dSL Test Data",
                              "eSL.Test.Data" = "eSL Test Data",
                              "Multivariate.model.Test.Data" = "Multivariate JM\n Test Data"))+
  geom_vline(xintercept = (3+4)/2, linetype = "dashed", color = "grey", linewidth = 0.5)+
  ggtitle("Integrated Brier Score: (6,7.5]")


## EPCE plot

df_comp_full_epce <- rbind(df_comparing_epce_randcens2, df_comparing_epce_randcens)
df_comp_full_epce <- df_comp_full_epce[-25,]


df2 <- with(df_comp_full_epce, data.frame("dSL CV Data" = EPCE_dSL_train,
                                          "eSL CV Data" = EPCE_eSL_train,
                                          "Multivariate model CV Data" = EPCE_mJM_train,
                                          "dSL Test Data" = EPCE_dSL_test,
                                          "eSL Test Data" = EPCE_eSL_test,
                                          "Multivariate model Test Data" = EPCE_mJM_test))
df_long2 <- df2 %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

df_long2$variable <- factor(df_long2$variable,
                            levels = c("dSL.CV.Data",
                                       "eSL.CV.Data",
                                       "Multivariate.model.CV.Data",
                                       "dSL.Test.Data",
                                       "eSL.Test.Data",
                                       "Multivariate.model.Test.Data"))

plottt_2 <- ggplot(df_long2, aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue") +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 20)) +
  xlab("") + ylab("EPCE") +
  scale_x_discrete(labels = c("dSL.CV.Data" = "dSL CV Data",
                              "eSL.CV.Data" = "eSL CV Data",
                              "Multivariate.model.CV.Data" = "Multivariate JM\n CV Data",
                              "dSL.Test.Data" = "dSL Test Data",
                              "eSL.Test.Data" = "eSL Test Data",
                              "Multivariate.model.Test.Data" = "Multivariate JM\n Test Data"))+
  geom_vline(xintercept = (3+4)/2, linetype = "dashed", color = "grey", linewidth = 1)+
  ggtitle("EPCE: (6,7.5]")


