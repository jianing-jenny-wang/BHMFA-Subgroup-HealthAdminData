#####################################
######### Post Processing ###########
#####################################
library(ggplot2)
library(tidyr)

color1 <- rgb(255,192,203, max = 255, alpha = 90, names = "lt.pink")
color2 <- rgb(158,162,163,max = 255, alpha = 90, names = "lt.gray")

##################################
### Summary of True Prevalence ###
##################################
summary_truePrev <- function(true_samples){
 nsim <- dim(true_samples)[1] 
 K <- dim(TruePrev_all_simu)[2]
 summary_true_prev_across_simu <- apply(true_samples, 2, FUN = function(x){
   mean_x <- mean(x)
   sd_x <- sd(x)
   median_x <- median(x)
   lb95_x <- quantile(x, probs = 0.025)
   ub95_x <- quantile(x, probs = 0.975)
   min_x <- min(x)
   max_x <- max(x)
   res_vec <- c(mean_x, sd_x, median_x, lb95_x, ub95_x, min_x, max_x)
   return(res_vec)
 })
 summary_true_prev_across_simu <- t(summary_true_prev_across_simu)
 colnames(summary_true_prev_across_simu) <- c("Mean_TruePrev_k", "SD_TruePrev_k","Median_TruePrev_k",
                                              "LB95_TruePrev_k", "UB95_TruePrev_k", "Min_TruePrev_k", "Max_TruePrev_k")
 summary_true_prev_across_simu <- data.frame(K = seq(1:K), summary_true_prev_across_simu)
 summary_true_prev <- colMeans(summary_true_prev_across_simu[,c(2:ncol(summary_true_prev_across_simu))])
 names(summary_true_prev) <- c("Ave_Mean_TruePrev_k", "Ave_SD_TruePrev_k","Ave_Median_TruePrev_k",
                                  "Ave_LB95_TruePrev_k", "Ave_UB95_TruePrev_k", "Ave_Min_TruePrev_k", "Ave_Max_TruePrev_k")
 summary_true_prev <- round(summary_true_prev, digits = 4)
 return(summary_true_prev)
 }



#####################
### Target size N ###
#####################

## Draw trace plot of total N target ##
trace_postNtarget <- function(samples, title = "(C1) Number of population"){
  plot(jitter(samples), xlab = "iteration", ylab = "N", main = title, type = "l")
}

## Draw posterior distribution of total N target ##
plot_postNtarget <- function(samples, xlimit, TotalN_obs, TotalN_target, names_of_chain = "(C1)"){
  hist(samples, nclass = 50, col = "gray95", main = paste(names_of_chain, "Posterior of Total N target", sep = " "), xlab = "Posterior of total target N",
       las = 1, xlim = xlimit)
  abline(v = quantile(samples, prob = 0.25), col="gray40", lwd = 3) # Lower bound of samples
  abline(v = quantile(samples, prob = 0.75), col="gray40", lwd = 3) # Upper bound of samples
  abline(v = quantile(samples, prob = 0.025), col="gray75", lwd = 3) # Lower bound of samples
  abline(v = quantile(samples, prob = 0.975), col="gray75", lwd = 3) # Upper bound of samples
  abline(v = mean(samples), col="darkorange", lwd = 3) # Mean of samples
  abline(v = median(samples), col="red", lwd = 3) # Median of samples
  # abline(v = TotalN_obs, col="black", lwd = 4) # Observed (Total)
  abline(v = TotalN_target, col = "blue", lwd = 4) # Truth N target (Total)
  legend("topright", c("Observed","95% CrI", "50% CrI", "Median N", "Mean N", "True Total N"), fill=c("black", "gray75", "gray40", "red", "darkorange", "blue"), 
         xpd=TRUE, cex=0.8, bty='n')
}

## Compute summary statistics of post samples N ##
stat_fun <- function(x) {
  mean <- mean(x)
  med <- median(x)
  Qt2.5 <- quantile(x, prob = 0.025)
  Qt97.5 <- quantile(x, prob = 0.975)
  res <- c(mean, med, Qt2.5, Qt97.5)
  names(res) <- c("Mean", "Median", "Qt 2.5%", "Qt 97.5%")
  return(res)
}

##################
### Prevalence ###
##################

## Save samples from lambda when lambda was generated individually ##
save_prev_g <- function(samples, names_in_chain) {
  index_prev_g <- grep("^.*prev_g*.*",names_in_chain)
  chain_c_prev_g <- samples[,index_prev_g] # G columns
  return(chain_c_prev_g)
}

plot_prev <- function(G, samples, truevalue, names_of_chain, names_of_g){
  par(mfrow=c(2,2)) 
  for(g in 1:G){
    hist(samples[,g], col = "gray", main = paste(names_of_chain, " Post Prev. for ", names_of_g[g], sep = ""), xlab = "Prevalence", las = 1, cex=1)
    abline(v = quantile(samples[,g], prob = 0.25), col="black", lwd = 3) # Lower bound of samples
    abline(v = quantile(samples[,g], prob = 0.75), col="black", lwd = 3) # Upper bound of samples
    abline(v = quantile(samples[,g], prob = 0.025), col="gray50", lwd = 3) # Lower bound of samples
    abline(v = quantile(samples[,g], prob = 0.975), col="gray50", lwd = 3) # Upper bound of samples
    abline(v = mean(samples[,g]), col="darkorange", lwd = 3) # Mean of samples
    abline(v = median(samples[,g]), col="red", lwd = 3) # Median of samples
    abline(v = truevalue[g], col = "blue", lwd = 4) # Truth prevalence = N/P
  }
  # legend("topright", c("(95%)50% CrI", "Med-Prev","Mean-Prev", "True-Prev"), fill=c("gray75", "red", "darkorange", "blue"), 
  #        xpd=TRUE, cex=0.7, bty='n')
}


###############################
### logit list effect mu_jg ###
###############################
## Obtain post samples of logit_list_effect ##
save_mu_jg <- function(samples, names_in_chain) {
  index_mu_jg <- grep("^.*mu_jg*.*",names_in_chain)
  chain_c_mu_jg <- samples[,index_mu_jg] # J*G columns
  return(chain_c_mu_jg)
}
## Compute summary statistics of post samples logit_list_effect ##
stat_fun_mu_jg <- function(J = J, chain = 1, samples, true_mu_jg = dtListEffects_one_simu){
  mean <- colMeans(samples)
  median <- apply(samples,2,median)
  Qt2.5 <- apply(samples, 2, FUN = function(x) quantile(x, prob = 0.025))
  Qt97.5 <- apply(samples, 2, FUN = function(x) quantile(x, prob = 0.975))
  res <- data.frame(Chain = rep(chain, length(mean)), List = seq(1:J), Mean = mean, Median = median, Qt2.5 = Qt2.5, Qt97.5 = Qt97.5, true_logitlisteff = true_logitlisteff)
  names(res) <- c("Chain", "List", "Mean", "Median", "Qt 2.5%", "Qt 97.5%","True logit lists effect")
  row.names(res) <- NULL
  return(res)     
}

## Plot the comparison between posterior and prior distributions
plot_listeffect <- function(J, post_listeffect, prior_listeffect, truevalue, names_of_chain){
  if(J %% 2 == 0){
  par(mfrow=c(J/2,2))  
  }
  if(J %% 2 != 0){
  par(mfrow=c((J+1)/2,2))   
  }

  for(j in 1:J){
    post <- hist(post_listeffect[,j], plot = FALSE)
    prior <- hist(prior_listeffect[,j], plot = FALSE)
    xlim <- range(c(post$breaks, prior$breaks))
    ylim <- max(c(post$count, prior$count)) 
    plot(post, col = color2, xlim = xlim, ylim = c(0,ylim), main = paste(names_of_chain, " logit list effect j =", j, sep = ""), xlab = "logit list effect", las = 1)
    plot(prior, add = TRUE, col = color1)
    abline(v = truevalue[j], col = "blue", lwd = 3) # true value of logit list effect
    legend("topright", c("Post","Prior"), fill=c(color2,color1), 
           xpd=TRUE, cex=1.2, bty='n')
  }
  
}


###############################
### logit list effect mu_j ####
##### In homogeneous model ####
###############################
## Obtain post samples of logit_list_effect ##
save_mu_j <- function(samples, names_in_chain) {
  index_mu_j <- grep("^.*mu_j*.*",names_in_chain)
  chain_c_mu_j <- samples[,index_mu_j] # G columns
  return(chain_c_mu_j)
}



###########################
## Factor Loading Matrix ##
###########################
save_loading_mat <- function(samples, names_in_chain) {
  index_loading <- grep("^.*loading_mat_upper_tri*.*",names_in_chain)
  chain_c_loading <- samples[,index_loading] # G columns
  return(chain_c_loading)
}

## Plot the comparison between posterior and prior distributions
plot_loading_mat <- function(J, D, post_loading_mat, truevalue, names_of_chain){
  if(J*D %% 2 == 0){
    par(mfrow=c(J*D/2,2))  
  }
  if(J*D %% 2 != 0){
    par(mfrow=c((J*D+1)/2,2))   
  }
  # post_loading_mat has dim = nsamples * D*J
  # at column-wise, it will go through 1:J for D=1, then go through 1:J for D=2, etc
  # Therefore, when D = 1, post_loading_mat has J columns
  # when D > 1, post_loading_mat has J*D columns
  if(D == 1){
    for(j in 1:J){
      post <- hist(post_loading_mat[,j], plot = FALSE)
      xlim <- range(c(post$breaks))
      ylim <- max(c(post$count)) 
      plot(post, col = color2, xlim = xlim, ylim = c(0,ylim), main = paste(names_of_chain, "loading matrix row j =", j, sep = ""), xlab = "loading matrix", las = 1)
      abline(v = truevalue[j,], col = "blue", lwd = 3) # true value of loading matrix
      legend("topright", c("Post"), fill=c(color2), 
             xpd=TRUE, cex=1.2, bty='n')
    }
  }
  if(D > 1){
    truevalue_vec <- as.vector(truevalue)
    element_index_j <- seq(1:J)
    element_index_d <- seq(1:D)
    permutation_jd <- expand.grid(element_index_j,element_index_d)
    element_index <- paste("[J = ", permutation_jd[,1], ", D = ", permutation_jd[,2], "]", sep = "")
    
    for(jd in 1:(J*D)){
      post <- hist(post_loading_mat[,jd], plot = FALSE)
      xlim <- range(c(post$breaks))
      ylim <- max(c(post$count)) 
      plot(post, col = color2, xlim = xlim, ylim = c(0,ylim), 
           main = paste(names_of_chain, " loading matrix element ", element_index[jd], sep = ""), 
           xlab = "loading matrix element value", las = 1)
      abline(v = truevalue_vec[jd], col = "blue", lwd = 3) # true value of loading matrix
      legend("topright", c("Post"), fill=c(color2), 
             xpd=TRUE, cex=1.2, bty='n')
      print(jd)
    }
  }
  par(mfrow=c(1,1)) 
}


##############
## Sigma_mu ##
##############
# ## Draw posterior distribution of sigma_mu ##
plot_postSigma_mu <- function(samples, TrueVals, xlimit, names_of_chain = "(C1)", model){
  if(model == "corr"){
  maintitle <- paste(names_of_chain, "Posterior of Sigma_mu", sep = " ")
  xlab <- "Posterior of Sigma_mu"
  }else{  maintitle <- paste(names_of_chain, "Posterior of computed SD of mu_jg", sep = " ")
          xlab <- "Posterior samples"}
  sample_hist <- hist(samples, xlim = xlimit, main = maintitle, xlab = xlab,
                   las = 1, nclass = 50)
  abline(v = quantile(samples, prob = 0.25), col="gray40", lwd = 3) # Lower bound of samples
  abline(v = quantile(samples, prob = 0.75), col="gray40", lwd = 3) # Upper bound of samples
  abline(v = quantile(samples, prob = 0.025), col="gray75", lwd = 3) # Lower bound of samples
  abline(v = quantile(samples, prob = 0.975), col="gray75", lwd = 3) # Upper bound of samples
  abline(v = mean(samples), col="darkorange", lwd = 3) # Mean of samples
  abline(v = median(samples), col="red", lwd = 3) # Median of samples
  abline(v = TrueVals, col = "blue", lwd = 4)
  legend("topright", c("95% CrI", "50% CrI", "Median", "Mean","TrueVals"),
         fill=c("gray75", "gray40", "red", "darkorange","blue"),
  xpd=TRUE, cex=0.8, bty='n')
}


# # If compare between two distributions 
# c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
# c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
# 
# ## Draw posterior distribution of sigma_mu ##
# plot_postSigma_mu <- function(samples, TrueVals, xlimit, names_of_chain = "(C1)"){
#   sample_hist <- hist(samples, xlim = xlimit,
#                    las = 1, nclass = 50, plot = FALSE)
#   true_hist <- hist(TrueVals, xlim = xlimit,
#                     las = 1, nclass = 50, plot = FALSE) # Truth values
#   plot(sample_hist, xlim = xlimit, main = paste(names_of_chain, "Posterior of Sigma_mu", sep = " "), xlab = "Posterior of Sigma_mu", col = c1)
#   plot(true_hist, xlim = xlimit, main = paste(names_of_chain, "Posterior of Sigma_mu", sep = " "), xlab = "Posterior of Sigma_mu", col = c2, add = TRUE)
#   abline(v = quantile(samples, prob = 0.25), col="gray40", lwd = 3) # Lower bound of samples
#   abline(v = quantile(samples, prob = 0.75), col="gray40", lwd = 3) # Upper bound of samples
#   abline(v = quantile(samples, prob = 0.025), col="gray75", lwd = 3) # Lower bound of samples
#   abline(v = quantile(samples, prob = 0.975), col="gray75", lwd = 3) # Upper bound of samples
#   abline(v = mean(samples), col="darkorange", lwd = 3) # Mean of samples
#   abline(v = median(samples), col="red", lwd = 3) # Median of samples
#   legend("topright", c("95% CrI", "50% CrI", "Median", "Mean","TrueVals","Post-Samples"),
#          fill=c("gray75", "gray40", "red", "darkorange","lightpink","lightblue"),
#   xpd=TRUE, cex=0.8, bty='n')
# }



##############
## Sigma_theta ##
##############
# ## Draw posterior distribution of sigma_mu ##
plot_postSigma_theta <- function(samples, TrueVals, xlimit, names_of_chain = "(C1)"){
  maintitle <- paste(names_of_chain, "Posterior of Sigma_theta", sep = " ")
  xlab <- "Posterior of Sigma_theta"
  sample_hist <- hist(samples, xlim = xlimit, main = maintitle, xlab = xlab,
                      las = 1, nclass = 50)
  abline(v = quantile(samples, prob = 0.25), col="gray40", lwd = 3) # Lower bound of samples
  abline(v = quantile(samples, prob = 0.75), col="gray40", lwd = 3) # Upper bound of samples
  abline(v = quantile(samples, prob = 0.025), col="gray75", lwd = 3) # Lower bound of samples
  abline(v = quantile(samples, prob = 0.975), col="gray75", lwd = 3) # Upper bound of samples
  abline(v = mean(samples), col="darkorange", lwd = 3) # Mean of samples
  abline(v = median(samples), col="red", lwd = 3) # Median of samples
  abline(v = TrueVals, col = "blue", lwd = 4)
  legend("topright", c("95% CrI", "50% CrI", "Median", "Mean","TrueVals"),
         fill=c("gray75", "gray40", "red", "darkorange","blue"),
         xpd=TRUE, cex=0.8, bty='n')
}


########################################
### Post process to compute measures ###
########################################
### Compute the Mean Diff/Relative Bias, RMSE, and MAE across nsim ###
Comp.Bias.RMSE.MAE.vec <- function(est,truevalue){
  if(length(truevalue) == 1){
    truevalue <- rep(truevalue, times = length(est))
  }
  if(length(est) == length(truevalue)){
    MSE <- mean(est - truevalue)^2  
    RMSE <- sqrt(MSE)
    MAE <- median(abs(est - truevalue))
    mean.diff.bias <- mean(est-truevalue)
    var.rel.bias <- sd(est/truevalue)
    mean.rel.bias <- mean(est/truevalue)
    Qt.2.5.rel.bias <- quantile(est/truevalue, probs = c(0.025, 0.975))[1]
    Qt.97.5.rel.bias <- quantile(est/truevalue, probs = c(0.025, 0.975))[2]
    T.Test.Diff.bias <- t.test(est-truevalue)
    T.Test.P.Value.Diff.bias <- T.Test.Diff.bias$p.value
    Bayesian.P.Value <- sum(est>=truevalue)/length(est) # the probability that the replicated data could be more extreme than the observed data (BDA P146)
    comp.bias <- c(MSE, RMSE, MAE, mean.diff.bias,
                    mean.rel.bias, var.rel.bias, 
                    Qt.2.5.rel.bias, Qt.97.5.rel.bias,
                    T.Test.P.Value.Diff.bias,
                    Bayesian.P.Value)
    names(comp.bias) <- c("MSE","RMSE","MAE", "Mean(EstVar-TrueVar)",
                          "Mean.Rel.Bias", "SD.Rel.Bias",
                          "Qt2.5.Rel.Bias", "Qt97.5.Rel.Bias", "T.Test.P.Value.Diff.Bias",
                          "Bayesian.P.Value")
  }
  if(length(est) != length(truevalue)){
    message("estimated values must be the same length of number of true values")
  }
  return(comp.bias)
}

