#setwd("/Users/arnaugarcia/Desktop/TFM_super")
library(ggplot2)
library(tidyr)



load("results20_sce2BIS.RData")
cens_tst <- numeric(80)
cens_tr <- numeric(80)
cens_tst[1:20] <- censoring_test[1:20]
cens_tr[1:20] <- censoring_train[1:20]

#with ibs
IBS_m <- numeric(80)
dSL_IBS <- numeric(80)
eSL_tr_IBS <- numeric(80)
eSL_tst_IBS <- numeric(80)
IBS_m_tst <- numeric(80)

IBS_m[1:20] <- IBS_multi[1:20]
dSL_IBS[1:20] <- dSL_cv_IBS[1:20]
eSL_tr_IBS[1:20] <- eSL_cv_IBS[1:20]
IBS_m_tst[1:20] <- IBS_multi_test[1:20]
eSL_tst_IBS[1:20] <- eSL_test_IBS[1:20]

#with epce
EPCE_m <- numeric(80)
dSL_EPCE <- numeric(80)
eSL_tr_EPCE <- numeric(80)
eSL_tst_EPCE <- numeric(80)
EPCE_m_tst <- numeric(80)

EPCE_m[1:20] <- EPCE_multi[1:20]
dSL_EPCE[1:20] <- dSL_cv_EPCE[1:20]
eSL_tr_EPCE[1:20] <- eSL_cv_EPCE[1:20]
EPCE_m_tst[1:20] <- EPCE_multi_test[1:20]
eSL_tst_EPCE[1:20] <- eSL_test_EPCE[1:20]

load("results60_sce2.RData")

cens_tst[21:80] <- censoring_test[1:60]
cens_tr[21:80] <- censoring_train[1:60]

IBS_m[21:80] <- IBS_multi[1:60]
dSL_IBS[21:80] <- dSL_cv_IBS[1:60]
eSL_tr_IBS[21:80] <- eSL_cv_IBS[1:60]
IBS_m_tst[21:80] <- IBS_multi_test[1:60]
eSL_tst_IBS[21:80] <- eSL_test_IBS[1:60]

IBS_multi <- IBS_m
IBS_multi_t <- IBS_m_tst
eSL_cv_IBS <- eSL_tr_IBS
eSL_t_IBS <- eSL_tst_IBS
dSL_cv_IBS <- dSL_IBS

EPCE_m[21:80] <- EPCE_multi[1:60]
dSL_EPCE[21:80] <- dSL_cv_EPCE[1:60]
eSL_tr_EPCE[21:80] <- eSL_cv_EPCE[1:60]
EPCE_m_tst[21:80] <- EPCE_multi_test[1:60]
eSL_tst_EPCE[21:80] <- eSL_test_EPCE[1:60]

EPCE_multi <- EPCE_m
EPCE_multi_t <- EPCE_m_tst
eSL_cv_EPCE <- eSL_tr_EPCE
eSL_t_EPCE <- eSL_tst_EPCE
dSL_cv_EPCE <- dSL_EPCE

#boxplot IBS
df <- data.frame(IBS_multi, IBS_multi_t, eSL_cv_IBS, eSL_t_IBS,  
                 dSL_cv_IBS)
df_long <- df %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

ploty <- ggplot(df_long, aes(x = variable, y = value)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("") + ylab("IBS") +
  ggtitle("Integrated Brier Score: (12,18]")

pdf("IBS_10iters_allmodels_generator.pdf",width=12)
ploty
dev.off()


#boxplot EPCE
df2 <- data.frame(EPCE_multi, EPCE_multi_t, eSL_cv_EPCE, eSL_t_EPCE,
                  dSL_cv_EPCE)
df_long2 <- df2 %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

plottt <- ggplot(df_long2, aes(x = variable, y = value)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("") + ylab("EPCE") +
  ggtitle("EPCE: (12,18]")

pdf("EPCE_20iters_multi_generator.pdf",width=12)
plottt
dev.off()


(sum(censoring_test[1:20])/20)/300
(sum(censoring_train[1:20])/20)/300

# on avg 40% of censoring

