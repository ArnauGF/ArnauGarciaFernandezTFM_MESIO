##############################################################################
##############################################################################
## Plotting results for SL simulation study
##############################################################################
##############################################################################
library(dplyr)
library(tidyr)
library(ggplot2)

## - We define a function that generates boxplots for the results obtained
##.  for both metrics, IBS and EPCE
plot_boxs <- function(string, sce){
  #Load the RData file with results
  load(string)

  
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
  
  
  #boxplot IBS
  df <- with(df_comparing_ibs, data.frame("dSL CV Data" = IBS_dSL_train,
                                          "eSL CV Data" = IBS_eSL_train,
                                          "Multivariate model CV Data" = IBS_mJM_train,
                                          "dSL Test Data" = IBS_dSL_test,
                                          "eSL Test Data" = IBS_eSL_test,
                                          "Multivariate model Test Data" = IBS_mJM_test))
  df_long <- df %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
  
  titl_ibs <- paste0("IBS (IPCW); ", sce)
  
  df_long$variable <- factor(df_long$variable,
                             levels = c("dSL.CV.Data",
                                        "eSL.CV.Data",
                                        "Multivariate.model.CV.Data",
                                        "dSL.Test.Data",
                                        "eSL.Test.Data",
                                        "Multivariate.model.Test.Data"))
  
  plot_ibs <- ggplot(df_long, aes(x = variable, y = value)) +
    geom_boxplot(fill = "lightblue") +
    theme_minimal() +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 5),
          plot.title = element_text(size = 14)) +
    xlab("") + ylab("IBS") +
    scale_x_discrete(labels = c("dSL.CV.Data" = "dSL CV Data",
                                "eSL.CV.Data" = "eSL CV Data",
                                "Multivariate.model.CV.Data" = "Multivariate JM\n CV Data",
                                "dSL.Test.Data" = "dSL Test Data",
                                "eSL.Test.Data" = "eSL Test Data",
                                "Multivariate.model.Test.Data" = "Multivariate JM\n Test Data"))+
    geom_vline(xintercept = (3+4)/2, linetype = "dashed", color = "grey", linewidth = 0.5)+
    ggtitle(titl_ibs)
  
  
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
  
  titl_epce <- paste0("EPCE; ", sce)
  
  plot_epce <- ggplot(df_long2, aes(x = variable, y = value)) +
    geom_boxplot(fill = "lightblue") +
    theme_minimal() +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 5),
          plot.title = element_text(size = 14)) +
    xlab("") + ylab("EPCE") +
    scale_x_discrete(labels = c("dSL.CV.Data" = "dSL CV Data",
                                "eSL.CV.Data" = "eSL CV Data",
                                "Multivariate.model.CV.Data" = "Multivariate JM\n CV Data",
                                "dSL.Test.Data" = "dSL Test Data",
                                "eSL.Test.Data" = "eSL Test Data",
                                "Multivariate.model.Test.Data" = "Multivariate JM\n Test Data"))+
    geom_vline(xintercept = (3+4)/2, linetype = "dashed", color = "grey", linewidth = 0.5)+
    ggtitle(titl_epce)
  
  return(list(plot_ibs, plot_epce))
}


## an expl on how to use the function above with preliminar results:
box_randCens <- plot_boxs("/Users/arnaugarcia/Desktop/TFM_super/results_SL/ibsipcw_epce_simul_4outcomes_correl_randCens_SL_08mar2025.RData",
                          "(6, 7.5]; Random Censoring")
box_infoCens <- plot_boxs("/Users/arnaugarcia/Desktop/TFM_super/results_SL/ibsipcw_epce_simul_4outcomes_correl_infoCens_SL_08mar2025.RData",
                          "(6, 7.5]; Informative Censoring")

box_adminCens <- plot_boxs("/Users/arnaugarcia/Desktop/TFM_super/results_SL/ibsipcw_epce_simul_4outcomes_correl_typ1_SL_06mar2025.RData",
                           "(4,5.5]; Admin Censoring")

## we can use ggarrenge() function to put boxplots together:
library(ggpubr)

# Create a blank plot with a red dashed horizontal line (simulating a separator between rows)
separator <- ggplot() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1) + 
  theme_void() 

boxs_ibs <- ggarrange(box_adminCens[[1]], box_randCens[[1]], box_infoCens[[1]],
                      nrow = 1, ncol = 3)

boxs_epce <- ggarrange(box_adminCens[[2]], box_randCens[[2]], box_infoCens[[2]],
                       nrow = 1, ncol = 3)


boxplots_prev <- ggarrange(boxs_ibs, 
                           separator,
                           boxs_epce, 
                           ncol = 1,
                           heights = c(2, 0.1, 2))

setwd("/Users/arnaugarcia/Desktop/TFM_super")
pdf("boxplots_SL.pdf",width=12)
boxplots_prev
dev.off()