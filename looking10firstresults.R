#setwd("/Users/arnaugarcia/Desktop/TFM_super")
#load("results20.RData")
setwd("G:/TFM")
library(ggplot2)
library(tidyr)

#with ibs
IBS_multi <- IBS_multi[1:10]
dSL_cv_IBS <- dSL_cv_IBS[1:10]
eSL_cv_IBS <- eSL_cv_IBS[1:10]
IBS_multi_t <- IBS_multi_test[1:10]
eSL_t_IBS <- eSL_test_IBS[1:10]


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

#with epce
EPCE_multi <- EPCE_multi[1:10]
EPCE_multi_t <- EPCE_multi_test[1:10]
dSL_cv_EPCE <- dSL_cv_EPCE[1:10]
eSL_cv_EPCE <- eSL_cv_EPCE[1:10]
eSL_t_EPCE <- eSL_test_EPCE[1:10]

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


(sum(censoring_test[1:10])/10)/300
(sum(censoring_train[1:10])/10)/300
#on avg 34% of censoring

# 40% of censoring

