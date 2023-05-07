# Jianing Wang
# Boston University Dissertation
# Chapter 3


#########################################
####### Run other candidate model #######
#### Automated Log-linear CRC model #####
#########################################

if (!require("MASS")) install.packages("MASS") else (require("MASS", quietly = TRUE)) 
if (!require("dplyr")) install.packages("dplyr") else (require("dplyr", quietly = TRUE)) 
if (!require("tidyr")) install.packages("tidyr") else (require("tidyr", quietly = TRUE)) 
if (!require("data.table")) install.packages("data.table") else (require("data.table", quietly = TRUE)) 
if (!require("reshape2")) install.packages("reshape2") else (require("reshape2", quietly = TRUE)) 
if (!require("ggplot2")) install.packages("ggplot2") else (require("ggplot2", quietly = TRUE)) 
if (!require("stringr")) install.packages("stringr") else (require("stringr", quietly = TRUE)) 
if (!require("LaplacesDemon")) install.packages("LaplacesDemon") else (require("LaplacesDemon", quietly = TRUE)) 


## Path to functions
path_to_funcs <- ".../Simulation/"
source(paste(path_to_funcs, "data_simu_functions.r", sep = ""))
source(paste(path_to_funcs, "postprocess_functions.r", sep = ""))
source(paste0(path_to_funcs, "modelfit_functions.r"))

## Save data
path_to_data <- ".../Simulation/MapData/"


#### Choose the Scenario ####
scenario <- 1.2 # 1.1 or 1.2
J <- 3
G <- 4
## If including or excluding the extremely small data source (large variance will include small list)
Large_Variance_PDetect <- "Large" # Large # Small

if(scenario == 1.1){
  Prevalence_scenario <- "SamePrev" 
  case <- paste("Simu", scenario, "_", Prevalence_scenario, "_", Large_Variance_PDetect, "Var", "_LogLinearCRC", sep = "")
}else{
  Prevalence_scenario <- "VaryPrev" 
  case <- paste("Simu", scenario, "_", Prevalence_scenario, "_", Large_Variance_PDetect, "Var", "_LogLinearCRC", sep = "")
}


#### Which Case ####

generate_plots <- FALSE # if generate plots through batch jobs: mix-trace plot

J <- 3
G <- 4

#### Set Working Directory ####
if(case == "Simu1.1_SamePrev_SmallVar_LogLinearCRC"){
  path_to_input <- ".../Simulation/GenDt/Simu1.1_SamePrev_J3G4/mu_jg_SmallVar/CorrBHMFA/"
  path_to_corrBHM_output <- ".../Simulation/MCMC_Out/Simu1.1_SamePrev_J3G4/mu_jg_SmallVar/CorrBHMFA/"
  path_to_output <- ".../Simulation/MCMC_Out/Simu1.1_SamePrev_J3G4/mu_jg_SmallVar/LoglinearCRC/"
}
if(case == "Simu1.1_SamePrev_LargeVar_LogLinearCRC"){
  path_to_input <- ".../Simulation/GenDt/Simu1.1_SamePrev_J3G4/mu_jg_LargeVar/CorrBHMFA/"
  path_to_corrBHM_output <- ".../Simulation/MCMC_Out/Simu1.1_SamePrev_J3G4/mu_jg_LargeVar/CorrBHMFA/"
  path_to_output <- ".../Simulation/MCMC_Out/Simu1.1_SamePrev_J3G4/mu_jg_LargeVar/LoglinearCRC/"
}
if(case == "Simu1.2_VaryPrev_SmallVar_LogLinearCRC"){
  path_to_input <- ".../Simulation/GenDt/Simu1.2_VaryPrev_J3G4/mu_jg_SmallVar/CorrBHMFA/"
  path_to_corrBHM_output <- ".../Simulation/MCMC_Out/Simu1.2_VaryPrev_J3G4/mu_jg_SmallVar/CorrBHMFA/"
  path_to_output <- ".../Simulation/MCMC_Out/Simu1.2_VaryPrev_J3G4/mu_jg_SmallVar/LoglinearCRC/"
}
if(case == "Simu1.2_VaryPrev_LargeVar_LogLinearCRC"){
  path_to_input <- ".../Simulation/GenDt/Simu1.2_VaryPrev_J3G4/mu_jg_LargeVar/CorrBHMFA/"
  path_to_corrBHM_output <- ".../Simulation/MCMC_Out/Simu1.2_VaryPrev_J3G4/mu_jg_LargeVar/CorrBHMFA/"
  path_to_output <- ".../Simulation/MCMC_Out/Simu1.2_VaryPrev_J3G4/mu_jg_LargeVar/LoglinearCRC/"
}

path_to_funcs <- ".../Simulation/"

#### Load functions ####
source(paste0(path_to_funcs, "data_simu_functions.r"))
source(paste0(path_to_funcs, "modelfit_functions.r"))
source(paste0(path_to_funcs, "postprocess_functions.r"))

### The log linear CRC model will only run on those had generated corrBHMFA (true) model results ###
### Import which simulation ran through in corrBHMFA ###
TrueVals <- list.files(path = path_to_corrBHM_output, pattern = 'TrueVals_(simu\\d+).Rdata$')
# Find number of simulations
nsim <- length(TrueVals)
# Find which simulation runs through
simu_ran_through <- c()
for(task_id in 1:nsim){
  which_ran_through <- str_extract_all(TrueVals[task_id],"\\(?[0-9]+\\)?")
  simu_ran_through <- c(simu_ran_through,which_ran_through)
}
simu_ran_through <- sort(as.numeric(do.call(c,simu_ran_through)))

############################
### Fit Log-Linear Model ###
############################

### For each generated data, fit log-linear model with Negative Binomial Distribution assumption ###
rel.bias <- c()
Error <- c()
AbsoluteError <- c()
SquaredError <- c()
rel.bias_prev <- matrix(0, nrow = G, ncol = nsim)
Error_prev <- matrix(0, nrow = G, ncol = nsim)
AbsoluteError_prev <- matrix(0, nrow = G, ncol = nsim)
SquaredError_prev <- matrix(0, nrow = G, ncol = nsim)

for(i in 1:nsim){
  task_id <- simu_ran_through[i]
  ### Read the data  
  dtYfullTable_one_simu <- readRDS(paste0(path_to_input, "dtYfullTable_simu", task_id, ".RData"))
  dtNtarget_one_simu <- readRDS(paste0(path_to_input, "dtNtarget_simu", task_id, ".RData"))
  dtBasicSetUp_one_simu <- readRDS(paste0(path_to_input, "BasicSetUp.RData"))
  P <- dtBasicSetUp_one_simu$P
  J <- dtBasicSetUp_one_simu$J
  G <- dtBasicSetUp_one_simu$G
  
  ### Create the contingency table 
  if(J == 3){
    empty_grid_list <- expand.grid(List1 = c(0,1), List2 = c(0,1), List3 = c(0,1), G = seq(1:G)) 
    empty_grid_list <- empty_grid_list[-which(empty_grid_list$List1 == 0 & empty_grid_list$List2 == 0 & empty_grid_list$List3 == 0),]
    yObs_by_list_grp <- aggregate(dtYfullTable_one_simu$obs_label, by=list(List1 = dtYfullTable_one_simu$X1,
                                                                           List2 = dtYfullTable_one_simu$X2,
                                                                           List3 = dtYfullTable_one_simu$X3,
                                                                           G = dtYfullTable_one_simu$grp_label), FUN=sum)
    # Get complete yObs table by list and location
    yObs_contingency_final <- merge(x = empty_grid_list, y = yObs_by_list_grp, by.x = c("List1","List2","List3", "G"),
                                    by.y = c("List1","List2","List3", "G"), 
                                    all.x = TRUE)
    yObs_contingency_final$x <- ifelse(is.na(yObs_contingency_final$x),0, yObs_contingency_final$x)
    colnames(yObs_contingency_final) <- c("List1","List2","List3","G","n_obs")
    yObs_contingency_final <- with(yObs_contingency_final, yObs_contingency_final[order(G, List1, List2, List3),])
  }
  
  ### Separate into Group-specific contingency table ###
  yObs_contingency_by_grp_ls <- split(yObs_contingency_final, f = yObs_contingency_final$G)    
  
  ## Poisson vs NB Models
  tryCatch(
    expr = {if(J == 3){
      est_total_by_grp <- sapply(1:G, FUN = function(g){
        dt <- as.data.frame( yObs_contingency_by_grp_ls[g])
        colnames(dt) <- c("List1","List2","List3","G","n_obs")
        dt$List1 <- as.character(dt$List1)
        dt$List2 <- as.character(dt$List2)
        dt$List3 <- as.character(dt$List3)
        if(sum(dt$n_obs)>0){
          ## Poisson Models
          mod_poi0 <- glm(n_obs ~ List1 + List2 + List3, data = dt, family = poisson(link = "log"))
          mod_poi1 <- glm(n_obs ~ List1 + List2 + List3 + List1:List2, data = dt, family = poisson(link = "log"))
          mod_poi2 <- glm(n_obs ~ List1 + List2 + List3 + List1:List3, data = dt, family = poisson(link = "log"))
          mod_poi3 <- glm(n_obs ~ List1 + List2 + List3 + List2:List3, data = dt, family = poisson(link = "log"))
          mod_poi4 <- glm(n_obs ~ List1 + List2 + List3 + List2:List3 + List1:List2, data = dt, family = poisson(link = "log"))
          mod_poi5 <- glm(n_obs ~ List1 + List2 + List3 + List2:List3 + List1:List3, data = dt, family = poisson(link = "log"))
          mod_poi6 <- glm(n_obs ~ List1 + List2 + List3 + List1:List3 + List1:List2, data = dt, family = poisson(link = "log"))
          mod_poi7 <- glm(n_obs ~ List1 + List2 + List3 + List1:List3 + List1:List2 + List2:List3, data = dt, family = poisson(link = "log"))
          aic_vec_poi <- c(mod_poi0$aic, mod_poi1$aic, mod_poi2$aic, mod_poi3$aic, mod_poi4$aic, mod_poi5$aic, mod_poi6$aic, mod_poi7$aic)
          coef_vec_poi <- c(mod_poi0$coefficients[1], mod_poi1$coefficients[1], mod_poi2$coefficients[1], 
                            mod_poi3$coefficients[1], mod_poi4$coefficients[1], mod_poi5$coefficients[1], 
                            mod_poi6$coefficients[1], mod_poi7$coefficients[1])
          which_best_aic_poi <- which(aic_vec_poi == min(aic_vec_poi))
          ## NB Models
          mod_nb0 <- glm.nb(n_obs ~ List1 + List2 + List3, data = dt)
          mod_nb1 <- glm.nb(n_obs ~ List1 + List2 + List3 + List1:List2, data = dt)
          mod_nb2 <- glm.nb(n_obs ~ List1 + List2 + List3 + List1:List3, data = dt)
          mod_nb3 <- glm.nb(n_obs ~ List1 + List2 + List3 + List2:List3, data = dt)
          mod_nb4 <- glm.nb(n_obs ~ List1 + List2 + List3 + List2:List3 + List1:List2, data = dt)
          mod_nb5 <- glm.nb(n_obs ~ List1 + List2 + List3 + List2:List3 + List1:List3, data = dt)
          mod_nb6 <- glm.nb(n_obs ~ List1 + List2 + List3 + List1:List3 + List1:List2, data = dt)
          # mod_nb7 <- glm.nb(n_obs ~ List1 + List2 + List3 + List1:List3 + List1:List2 + List2:List3, data = dt) # not fitting, saturated
          aic_vec_nb <- c(mod_nb0$aic, mod_nb1$aic, mod_nb2$aic, mod_nb3$aic, mod_nb4$aic, mod_nb5$aic, mod_nb6$aic)
          coef_vec_nb <- c(mod_nb0$coefficients[1], mod_nb1$coefficients[1], mod_nb2$coefficients[1], 
                           mod_nb3$coefficients[1], mod_nb4$coefficients[1], mod_nb5$coefficients[1], 
                           mod_nb6$coefficients[1])
          which_best_aic_nb <- which(aic_vec_nb == min(aic_vec_nb))
          
          if(min(aic_vec_poi) <= min(aic_vec_nb)){
            intercept_g <- coef_vec_poi[which_best_aic_poi]
            est_unknown_g <- exp(intercept_g)
          }
          else{intercept_g <- coef_vec_nb[which_best_aic_nb]
          est_unknown_g <- exp(intercept_g)}
        }
        if(sum(dt$n_obs) == 0){
          est_unknown_g <- 0
        }
        known_g <- sum(dt$n_obs)
        est_total_g <- round(est_unknown_g + known_g, digits = 0)
        return(est_total_g)
      })
    }
    },
    error = function(e){
      print(paste("Error in model fit", i, sep = " "))
    }
  )
  
  ### Compute relative bias on total Ntarget
  rel.bias_task_id <- sum(est_total_by_grp)/sum(dtNtarget_one_simu)
  rel.bias <- c(rel.bias, rel.bias_task_id)
  
  ### Compute RMSE/MAE on total Ntarget
  Error_task_id <- sum(est_total_by_grp) - sum(dtNtarget_one_simu)
  Error <- c(Error, Error_task_id)
  SquaredError_task_id <- (Error_task_id)^2
  AbsoluteError <- c(AbsoluteError, abs(Error_task_id))
  SquaredError <- c(SquaredError, SquaredError_task_id)
  
  ### Compute relative bias on area-specific prevalence
  est_prev_by_grp <- est_total_by_grp/P
  true_prev_by_grp <- dtNtarget_one_simu/P
  rel.bias_prev_task_id <- est_prev_by_grp/true_prev_by_grp
  rel.bias_prev[,i] <- rel.bias_prev_task_id
  
  ### Compute RMSE/MAE on area-specific prevalence
  Error_prev_task_id <- est_prev_by_grp - true_prev_by_grp
  Error_prev[,i] <- Error_prev_task_id
  AbsoluteError_prev[,i] <- abs(Error_prev_task_id)
  SquaredError_prev_task_id <- (Error_prev_task_id)^2
  SquaredError_prev[,i] <- SquaredError_prev_task_id
  
  print(i)
} # end of the loop for task_id

### Effectively how many generate finite estimates
stable_nsim <- length(rel.bias[which(rel.bias < 100000)])
stable_nsim
### Which failed to run through ###
fail_runs <- which(rel.bias >= 100000)


### Compute Mean of the relative bias for total Ntarget
rel.bias_final <- rel.bias[which(rel.bias < 100000)]
Mean.rel.bias <- mean(rel.bias_final)
ttest_pvalue <- t.test(Error, mu = 0)$p.value
### Compute RMSE for total Ntarget
SquaredError_final <- SquaredError[which(SquaredError < 100000)]
RMSE <- sqrt(mean(SquaredError_final))
### Compute MAE for total Ntarget
AbsoluteError_final <- AbsoluteError[which(AbsoluteError < 100000)]
MAE <- median(AbsoluteError_final)


### Summary for area-specific prevalence
rel.bias_prev_final <- apply(rel.bias_prev, c(1,2), FUN = function(x) {
  if(x > 5) {x <- NA}
  else{x <- x}
  return(x)})
### Number of times across nsim for each area to produce unstable estimates
freq_infinite_est_by_grp <- apply(rel.bias_prev_final, 1, FUN = function(x) sum(is.na(x)))
freq_infinite_est_by_grp_df <- data.frame(G = 1:G, freq_infinite_est_by_grp = freq_infinite_est_by_grp)
### Number of areas that never have unstable estimates
length(freq_infinite_est_by_grp_df$G[freq_infinite_est_by_grp_df$freq_infinite_est_by_grp== 0])
### Compute Mean of the relative bias for area-specific prevalence
Mean.rel.bias_prev <- apply(rel.bias_prev_final, 1, FUN = function(x){mean(x, na.rm = TRUE)})
ttest_pvalue_prev <- apply(Error_prev, 1, FUN = function(x) t.test(x, mu = 0)$p.value)

AbsoluteError_prev_final <- apply(AbsoluteError_prev, c(1,2), FUN = function(x) {
  if(x > 1) {x <- NA}
  else{x <- x}
  return(x)})
### Compute Mean of the relative bias for area-specific prevalence
MAE_prev <- apply(AbsoluteError_prev_final, 1, FUN = function(x){median(x, na.rm = TRUE)})

SquaredError_prev_final <- apply(AbsoluteError_prev, c(1,2), FUN = function(x) {
  if(x > 1) {x <- NA}
  else{x <- x}
  return(x)})
### Compute Mean of the relative bias for area-specific prevalence
RMSE_prev <- apply(SquaredError_prev_final, 1, FUN = function(x){sqrt(mean(x, na.rm = TRUE))})


summary_stats_Ntarget <- data.frame(case = case,
                                    stable_nsim = stable_nsim,
                                    Mean.rel.bias_Ntarget = Mean.rel.bias,
                                    ttest_pvalue_Ntarget = ttest_pvalue,
                                    RMSE_Ntarget = RMSE,
                                    MAE_Ntarget = MAE
)
summary_stats_Ntarget


summary_stats_prev <- data.frame(Mean.rel.bias_prev = Mean.rel.bias_prev,
                                 ttest_pvalue_prev = ttest_pvalue_prev,
                                 MAE_prev = MAE_prev,
                                 RMSE_prev = RMSE_prev
)
summary_stats_prev
