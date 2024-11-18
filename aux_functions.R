# we code a couple of functions for passing the mcmc from the JM object
# where the mcmc samples are separated in the three chains, and are giving by a list
# to a normal data.frame

#the first c
JMtodataframe_params <- function(x){
  require(JMbayes2)
  
  ####### BS gammas
  my_list <- x$mcmc$bs_gammas
  gamma_bs_post <- do.call(rbind, my_list)
  gamma_bs_post <- cbind(ID = rownames(gamma_bs_post), gamma_bs_post)
  rownames(gamma_bs_post) <- NULL
  gamma_bs_post <- as.data.frame(gamma_bs_post)
  
  ##### gammas
  my_list <- x$mcmc$gammas
  gamma_post <- do.call(rbind, my_list)
  gamma_post <- cbind(ID = rownames(gamma_post), gamma_post)
  rownames(gamma_post) <- NULL
  gamma_post <- as.data.frame(gamma_post)
  
  #### tau_bs_gammas
  my_list <- x$mcmc$tau_bs_gammas
  gamma_bs_tau_post <- do.call(rbind, my_list)
  gamma_bs_tau_post <- cbind(ID = rownames(gamma_bs_tau_post), gamma_bs_tau_post)
  rownames(gamma_bs_tau_post) <- NULL
  gamma_bs_tau_post <- as.data.frame(gamma_bs_tau_post)
  
  ### alphas
  my_list <- x$mcmc$alphas
  alphas_post <- do.call(rbind, my_list)
  alphas_post <- cbind(ID = rownames(alphas_post), alphas_post)
  rownames(alphas_post) <- NULL
  alphas_post <- as.data.frame(alphas_post)
  
  ### W_std_gammas
  my_list <- x$mcmc$W_std_gammas
  W_std_gammas_post <- do.call(rbind, my_list)
  W_std_gammas_post <- cbind(ID = rownames(W_std_gammas_post), W_std_gammas_post)
  rownames(W_std_gammas_post) <- NULL
  W_std_gammas_post <- as.data.frame(W_std_gammas_post)
  
  ### Wlong_std_alphas
  my_list <- x$mcmc$Wlong_std_alphas
  Wlong_std_alphas_post <- do.call(rbind, my_list)
  Wlong_std_alphas_post <- cbind(ID = rownames(Wlong_std_alphas_post), Wlong_std_alphas_post)
  rownames(Wlong_std_alphas_post) <- NULL
  Wlong_std_alphas_post <- as.data.frame(Wlong_std_alphas_post)
  
  ### D
  my_list <- x$mcmc$D
  D_post <- do.call(rbind, my_list)
  D_post <- cbind(ID = rownames(D_post), D_post)
  rownames(D_post) <- NULL
  D_post <- as.data.frame(D_post)
  
  ### betas1
  my_list <- x$mcmc$betas1
  betas1_post <- do.call(rbind, my_list)
  betas1_post <- cbind(ID = rownames(betas1_post), betas1_post)
  rownames(betas1_post) <- NULL
  betas1_post <- as.data.frame(betas1_post)
  
  ### sigmas
  my_list <- x$mcmc$sigmas
  sigmas_post <- do.call(rbind, my_list)
  sigmas_post <- cbind(ID = rownames(sigmas_post), sigmas_post)
  rownames(sigmas_post) <- NULL
  sigmas_post <- as.data.frame(sigmas_post)
  
  ### alphaF
  my_list <- x$mcmc$alphaF
  alphaF_post <- do.call(rbind, my_list)
  alphaF_post <- cbind(ID = rownames(alphaF_post), alphaF_post)
  rownames(alphaF_post) <- NULL
  alphaF_post <- as.data.frame(alphaF_post)
  
  ### sigmaF
  my_list <- x$mcmc$sigmaF
  sigmaF_post <- do.call(rbind, my_list)
  sigmaF_post <- cbind(ID = rownames(sigmaF_post), sigmaF_post)
  rownames(sigmaF_post) <- NULL
  sigmaF_post <- as.data.frame(sigmaF_post)
  
  
  ### joint all data frames
  posterior <- cbind(gamma_bs_post, gamma_post, gamma_bs_tau_post, alphas_post,
                     W_std_gammas_post, Wlong_std_alphas_post, D_post, 
                     betas1_post, sigmas_post, alphaF_post, 
                     sigmaF_post)
  return(posterior)
}


#function for the mcmc of the random effects, IT IS NOT WELL coded!!!
# we have to take into account that we have n individuals, and for each
# r random effects!
JMtodataframe_randeffects <- function(x){
  require(JMbayes2)
  
  ### b
  my_list <- x$mcmc$b
  b_post <- do.call(rbind, my_list)
  b_post <- cbind(ID = rownames(b_post), b_post)
  rownames(b_post) <- NULL
  b_post <- as.data.frame(b_post)
  
  ### frailty
  my_list <- Models$M1$mcmc$frailty
  frailty_post <- do.call(rbind, my_list)
  frailty_post <- cbind(ID = rownames(frailty_post), frailty_post)
  rownames(frailty_post) <- NULL
  frailty_post <- as.data.frame(frailty_post)
  
  posterior <- cbind(frailty_post, b_post)
  
  return(posterior)
}


## we try to use density estimation to compute the mode of the posterior distr
### problems: 
## - we need it for a lot of variables, and we need the argmax{} for all these
# variables at the same time. Thus, this is involving hard multivariate calculus
## and some packages in R are not working for this amount of dimensions and observations.

posterior_df_M1 <- JMtodataframe_params(Models$M1)
post_b_M1 <- JMtodataframe_randeffects(Models$M1)

my_list <- Models$M1$mcmc$bs_gammas
gamma_bs_post <- do.call(rbind, my_list)
gamma_bs_post <- cbind(ID = rownames(gamma_bs_post), gamma_bs_post)
rownames(gamma_bs_post) <- NULL
gamma_bs_post <- as.data.frame(gamma_bs_post)
head(gamma_bs_post)
library(ks)
kde_res <- kde(gamma_bs_post[,1:2])

grid_points <- kde_res$eval.points
density_values <- kde_res$estimate

max_idx <- which(density_values == max(density_values), arr.ind = TRUE)

# Get the corresponding mode values (grid points at max density)
mode_x <- grid_points[[1]][max_idx[1]]
mode_y <- grid_points[[2]][max_idx[2]]

cat("Mode of the KDE is at (x, y):", mode_x, mode_y, "\n")


plot(kde_res)
abline(v=mode_estimate, col="tomato", lty = 2)

my_list <- Models$M1$mcmc$b
b_post <- do.call(rbind, my_list)
b_post <- cbind(ID = rownames(b_post), b_post)
rownames(b_post) <- NULL
b_post <- as.data.frame(b_post)
