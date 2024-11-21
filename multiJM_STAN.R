library(rstan)
library(survival)
library(lme4)
##############################################################################
##############################################################################
## We try to fit multi JM using Stan, based on the code in rstanarm R paxkage
##############################################################################
##############################################################################

##############################################################################
##############################################################################
##############################################################################
##############################################################################
## R/stan_jm.R
##############################################################################
##############################################################################
##############################################################################
##############################################################################

##Load aux functions:
source("aux_functions_JM_STAN.R")

## Stan JM function
stan_jm <- function(formulaLong, dataLong, formulaEvent, dataEvent, time_var, 
                    id_var, family = gaussian, assoc = "etavalue", 
                    lag_assoc = 0, grp_assoc, scale_assoc = NULL, epsilon = 1E-5,
                    basehaz = c("bs", "weibull", "piecewise"), basehaz_ops,
                    qnodes = 15, init = "prefit", weights,	
                    priorLong = normal(autoscale=TRUE), priorLong_intercept = normal(autoscale=TRUE), 
                    priorLong_aux = cauchy(0, 5, autoscale=TRUE), priorEvent = normal(autoscale=TRUE), 
                    priorEvent_intercept = normal(autoscale=TRUE), priorEvent_aux = cauchy(autoscale=TRUE),
                    priorEvent_assoc = normal(autoscale=TRUE), prior_covariance = lkj(autoscale=TRUE), 
                    prior_PD = FALSE, algorithm = c("sampling", "meanfield", "fullrank"), 
                    adapt_delta = NULL, max_treedepth = 10L, QR = FALSE, 
                    sparse = FALSE, ...) {
  
  #-----------------------------
  # Pre-processing of arguments
  #-----------------------------  
  
  # Set seed if specified
  dots <- list(...)
  if ("seed" %in% names(dots))
    set.seed(dots$seed)
  
  algorithm <- match.arg(algorithm)
  basehaz   <- match.arg(basehaz)
  
  if (missing(basehaz_ops)) basehaz_ops <- NULL
  if (missing(weights))     weights     <- NULL
  if (missing(id_var))      id_var      <- NULL
  if (missing(time_var))    time_var    <- NULL
  if (missing(grp_assoc))   grp_assoc   <- NULL
  
  if (!is.null(weights)) 
    stop("'weights' are not yet implemented.")
  if (QR)               
    stop("'QR' decomposition is not yet implemented.")
  if (sparse)
    stop("'sparse' option is not yet implemented.")
  
  if (is.null(time_var))
    stop("'time_var' must be specified.")
  
  # Formula
  formulaLong <- validate_arg(formulaLong, "formula"); M <- length(formulaLong)
  if (M > 3L)
    stop("'stan_jm' is currently limited to a maximum of 3 longitudinal outcomes.")
  
  # Data
  dataLong <- validate_arg(dataLong, "data.frame", validate_length = M)  
  dataEvent <- as.data.frame(dataEvent)
  
  # Family
  ok_family_classes <- c("function", "family", "character")
  ok_families <- c("binomial", "gaussian", "Gamma", 
                   "inverse.gaussian", "poisson", "neg_binomial_2")
  family <- validate_arg(family, ok_family_classes, validate_length = M)
  family <- lapply(family, validate_famlink, ok_families)
  
  # Assoc
  ok_assoc_classes <- c("NULL", "character")
  assoc <- validate_arg(assoc, ok_assoc_classes, validate_length = M)
  
  # Is priorLong* already a list?
  priorLong <- broadcast_prior(priorLong, M)
  priorLong_intercept <- broadcast_prior(priorLong_intercept, M)
  priorLong_aux <- broadcast_prior(priorLong_aux, M)
  
  #-----------
  # Fit model
  #-----------
  
  stanfit <- stan_jm.fit(formulaLong = formulaLong, dataLong = dataLong, 
                         formulaEvent = formulaEvent, dataEvent = dataEvent, 
                         time_var = time_var, id_var = id_var, family = family,
                         assoc = assoc, lag_assoc = lag_assoc, grp_assoc = grp_assoc, 
                         epsilon = epsilon, basehaz = basehaz, basehaz_ops = basehaz_ops, 
                         qnodes = qnodes, init = init, weights = weights, scale_assoc = scale_assoc,
                         priorLong = priorLong, 
                         priorLong_intercept = priorLong_intercept, 
                         priorLong_aux = priorLong_aux, 
                         priorEvent = priorEvent, 
                         priorEvent_intercept = priorEvent_intercept, 
                         priorEvent_aux = priorEvent_aux, 
                         priorEvent_assoc = priorEvent_assoc, 
                         prior_covariance = prior_covariance, prior_PD = prior_PD, 
                         algorithm = algorithm, adapt_delta = adapt_delta, 
                         max_treedepth = max_treedepth, QR = QR, sparse = sparse, ...)
  if (algorithm != "optimizing" && !is(stanfit, "stanfit")) return(stanfit)
  y_mod <- attr(stanfit, "y_mod")
  e_mod <- attr(stanfit, "e_mod")
  a_mod <- attr(stanfit, "a_mod")
  cnms  <- attr(stanfit, "cnms")
  flevels <- attr(stanfit, "flevels")
  assoc <- attr(stanfit, "assoc")
  scale_assoc <- attr(stanfit, "scale_assoc")
  id_var <- attr(stanfit, "id_var")
  basehaz    <- attr(stanfit, "basehaz")
  grp_stuff  <- attr(stanfit, "grp_stuff")
  prior_info <- attr(stanfit, "prior_info")
  stanfit <- drop_attributes(stanfit, "y_mod", "e_mod", "a_mod", "cnms", 
                             "flevels", "assoc", "id_var", "basehaz", 
                             "grp_stuff", "prior_info","scale_assoc")
  
  terms <- c(fetch(y_mod, "terms"), list(terms(e_mod$mod)))
  n_yobs <- fetch_(y_mod, "x", "N")
  n_grps <- sapply(flevels, n_distinct)
  n_subjects <- e_mod$Npat
  
  fit <- nlist(stanfit, formula = c(formulaLong, formulaEvent), family,
               id_var, time_var, weights, scale_assoc, qnodes, basehaz, assoc,
               M, cnms, flevels, n_grps, n_subjects, n_yobs, epsilon,
               algorithm, terms, glmod = y_mod, survmod = e_mod, 
               assocmod = a_mod, grp_stuff, dataLong, dataEvent,
               prior.info = prior_info, stan_function = "stan_jm", 
               call = match.call(expand.dots = TRUE))
  
  out <- stanmvreg(fit)
  return(out)
}


###############################################################################
###############################################################################
### Using https://github.com/sambrilleman/2018-StanCon-Notebook/blob/master/notebook.pdf
##############################################################################
###############################################################################
library(rstan)
standata <- readRDS("Stan/data/standata.rds")
staninit <- readRDS("Stan/data/staninit.rds")
mod1 <- stan(
    file = "Stan/jm_SamBrilleman.stan", 
    data = standata,
    init = function() staninit,
    chains = 2, seed = 12345)

## let's see some results related with the model
print(mod1, pars = "y2_alpha")

library(tidyr)
posterior <- as.data.frame(mod1)
names(posterior)

library(bayesplot)
mcmc_areas(posterior,
           pars = c("y1_gamma"),
           prob = 0.8)

## the above code works, now let's try it adding the indicators to the 
# longitudinal outcomes

#we change the init list
staninit2 <- append(staninit, list(pind = 0.5))

mod2 <- stan(
  file = "Stan/jm_multi_simple_indicators.stan", 
  data = standata,
  init = function() staninit2,
  chains = 2, seed = 12345)

print(mod2, pars = "post_prob_z1")
print(mod2, pars = "post_prob_z2")

posterior2 <- as.data.frame(mod2)
names(posterior2)


mcmc_areas(posterior2,
           pars = c("post_prob_z1[1]"),
           prob = 0.8)

mcmc_areas(posterior2,
           pars = c("post_prob_z2[1]"),
           prob = 0.8)

## with these we have obtained weights 0.3589744 for the first longi outcome
# and 0.6410256 for the second one
#we now fit the same model with JMbayes2 package and apply SL to see differences
#between weigths
library(JMbayes2)
pbcLong <- readRDS("Data/pbcLong.rds")
pbcSurv <- readRDS("Data/pbcSurv.rds")

#we gather all the info into the longitudinal data set to be able to work with
#the functions we have:
library(dplyr)
pbc_full <- left_join(pbcLong, pbcSurv, by = "id")
pbc_full$age.y <- NULL
pbc_full$sex.y <- NULL
pbc_full$trt.y <- NULL
pbc_full <- pbc_full %>%
  rename(age = age.x, trt = trt.x, sex = sex.x)

CVdats <- create_folds(pbc_full, V = 3, id_var = "id")

fit_models <- function (dataLong) {
  library("JMbayes2")
  data_id <- dataLong[!duplicated(dataLong$id), ]
  CoxFit <- coxph(Surv(futimeYears, death) ~ sex + trt, data = data_id)
  LM1 <- lme(logBili ~ year, data = dataLong, random = ~ year | id)
  LM2 <- lme(albumin ~ year, data = dataLong, random = ~ year | id)
  JM1 <- jm(CoxFit, LM1, time_var = "year")
  JM2 <- jm(CoxFit, LM2, time_var = "year")
  out <- list(M1 = JM1, M2 = JM2)
  class(out) <- "jmList"
  out
}

cl <- parallel::makeCluster(5L)
Models_folds <- parallel::parLapply(cl, CVdats$training, fit_models)
parallel::stopCluster(cl)

# time window to compute the scores
t0 <- 5
dt <- 5
#computing Brier weights
Brier_weights <- tvBrier(Models_folds, newdata = CVdats$testing, 
                         integrated = TRUE, Tstart = t0, Dt = dt)
# EPCE weights
EPCE_weights <- tvEPCE(Models_folds, newdata = CVdats$testing, 
                       Tstart = t0, Dt = dt)

# computing with full data and taking into account the BMA indicators weights:

#fitting models with whole dataset
Models <- fit_models(pbc_full)

#we use Brier weights to compute IBS with full data:
bw <- Brier_weights$weights
brier_score_SL_full <- tvBrier(Models, newdata = pbc_full, model_weights = bw, 
                               Tstart = t0, Dt = dt, integrated = TRUE)

#we use EPCE weights to compute EPCE with full data:
ew <- EPCE_weights$weights
EPCE_score_SL_full <- tvEPCE(Models, newdata = pbc_full, model_weights = ew,
                             Tstart = t0, Dt = dt)

#we use the Indicator BMA weights:
BMA_IND_w <- c(0.3589744, 0.6410256)
brier_score_full_BMA <- tvBrier(Models, newdata = pbc_full, model_weights = BMA_IND_w, 
                               Tstart = t0, Dt = dt, integrated = TRUE)
EPCE_score_full_BMA <- tvEPCE(Models, newdata = pbc_full, model_weights = BMA_IND_w,
                             Tstart = t0, Dt = dt)

## we can see that in this case SL is giving slightly better results than BMA indicators
