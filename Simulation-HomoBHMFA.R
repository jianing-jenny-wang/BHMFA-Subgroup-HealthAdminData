###########################################################
######################## Fit Model ########################
######## Model BHMFA with homogeneous list effect #########
###########################################################

if (!require("dplyr")) install.packages("dplyr") else (require("dplyr", quietly = TRUE)) 
if (!require("reshape2")) install.packages("reshape2") else (require("reshape2", quietly = TRUE)) 
if (!require("stringr")) install.packages("stringr") else (require("stringr", quietly = TRUE)) 
if (!require("nimble")) install.packages("nimble") else (require("nimble", quietly = TRUE)) 
if (!require("coda")) install.packages("coda") else (require("coda", quietly = TRUE)) 
# if (!require("bayesplot")) install.packages("bayesplot") else (require("bayesplot", quietly = TRUE)) 
if (!require("LaplacesDemon")) install.packages("LaplacesDemon") else (require("LaplacesDemon", quietly = TRUE)) 
if (!require("ggmcmc")) install.packages("ggmcmc") else (require("ggmcmc", quietly = TRUE)) 
if (!require("jagsUI")) install.packages("jagsUI") else (require("jagsUI", quietly = TRUE)) # parallel computing
if (!require("parallel")) install.packages("parallel") else (require("parallel", quietly = TRUE)) # parallel computing




#### Which Case ####
#### Choose the Scenario ####
scenario <- 1.2 # 1.1 or 1.2

## If including or excluding the extremely small data source (large variance will include small list)
Large_Variance_PDetect <- "Large" # Large # Small

if(scenario == 1.1){
  Prevalence_scenario <- "SamePrev" 
  case <- paste("Simu", scenario, "_", Prevalence_scenario, "_", Large_Variance_PDetect, "Var", "_CorrBHMFA", sep = "")
}else{
  Prevalence_scenario <- "VaryPrev" 
  case <- paste("Simu", scenario, "_", Prevalence_scenario, "_", Large_Variance_PDetect, "Var", "_CorrBHMFA", sep = "")
}

#### Run with R-intel ####
run_with_intel <- "yes" 

#### Parallel the Chains ####
parallel_run <- "yes"

#### Set Working Directory ####
if(case == "Simu1.1_SamePrev_SmallVar_HomoBHMFA"){
  path_to_input <- ".../Simulation/GenDt/GenDatasets/Simu1.1_SamePrev_J3G4/mu_jg_SmallVar/CorrBHMFA/"
  path_to_output <- ".../Simulation/MCMC_Out/Simu1.1_SamePrev_J3G4/mu_jg_SmallVar/HomoBHMFA/"
}
if(case == "Simu1.1_SamePrev_LargeVar_HomoBHMFA"){
  path_to_input <- ".../Simulation/GenDt/GenDatasets/Simu1.1_SamePrev_J3G4/mu_jg_LargeVar/CorrBHMFA/"
  path_to_output <- ".../Simulation/MCMC_Out/Simu1.1_SamePrev_J3G4/mu_jg_LargeVar/HomoBHMFA/"
}
if(case == "Simu1.2_VaryPrev_SmallVar_HomoBHMFA"){
  path_to_input <- ".../Simulation/GenDt/GenDatasets/Simu1.2_VaryPrev_J3G4/mu_jg_SmallVar/CorrBHMFA/"
  path_to_output <- ".../Simulation/MCMC_Out/Simu1.2_VaryPrev_J3G4/mu_jg_SmallVar/HomoBHMFA/"
}
if(case == "Simu1.2_VaryPrev_LargeVar_HomoBHMFA"){
  path_to_input <- ".../Simulation/GenDt/GenDatasets/Simu1.2_VaryPrev_J3G4/mu_jg_LargeVar/CorrBHMFA/"
  path_to_output <- ".../Simulation/MCMC_Out/Simu1.2_VaryPrev_J3G4/mu_jg_LargeVar/HomoBHMFA/"
}

path_to_funcs <- ".../Simulation/"

#### Load functions ####
source(paste0(path_to_funcs, "data_simu_functions.r"))
source(paste0(path_to_funcs, "modelfit_functions.r"))
source(paste0(path_to_funcs, "postprocess_functions.r"))


#### Borrow the job ID to be seed ####
task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))
# task_id <- 40

#### Import the Datasets ####
### Import a list of simulated dataset ###
dtYfullTable_one_simu <- readRDS(paste0(path_to_input, "dtYfullTable_simu", task_id, ".RData"))
dtNtarget_one_simu <- readRDS(paste0(path_to_input, "dtNtarget_simu", task_id, ".RData"))
dtyObs_by_grp_one_simu <- readRDS(paste0(path_to_input, "dtyObs_by_grp_simu", task_id, ".RData"))
dtyObs_by_list_grp_simu <- readRDS(paste0(path_to_input, "dtyObs_by_list_grp_simu", task_id, ".RData"))
dtBasicSetUp_one_simu <- readRDS(paste0(path_to_input, "BasicSetUp.RData"))
dtPDetect_by_list_grp_val <- readRDS(paste0(path_to_input, "dtPDetect_by_list_grp_simu", task_id, ".RData"))
names_of_g <- dtBasicSetUp_one_simu$names_of_g
########## FINISH DATA IMPORT ##########

#### Data Manipulation ####
### Rename the data and Basic information ###
single_dt <- dtYfullTable_one_simu # one simulated dataset
G <- dtBasicSetUp_one_simu$G # number of groups
J <- dtBasicSetUp_one_simu$J # number of lists
D <- dtBasicSetUp_one_simu$D # number of latent factors
P <- dtBasicSetUp_one_simu$P # number of normal people by group
N_target <- dtNtarget_one_simu # number of target people by group
N_obs <- dtyObs_by_grp_one_simu$n_obs # number of observed people by group
TotalN_target <- sum(N_target) # total number of target population from all groups
TotalN_obs <- sum(N_obs) # total number of observed target population from all groups
Prev <- N_target/P # true prevalence by group
## Compute theoretical/empirical list effects in real generated data for posterior check
empirical_mu_jg_mat <- logit(dtPDetect_by_list_grp_val) # G*J in generated data 
empirical_pdetect_mat <- invlogit(empirical_mu_jg_mat)
theoretical_mu_jg_mat <- logit(dtBasicSetUp_one_simu$invlogit_mu_jg_mat) # J*G in theoretical simu design
theoretical_pdetect_mat <- invlogit(theoretical_mu_jg_mat)

# Create the observable dataset for one single simulated dataset
dtYObs <- subset(dtYfullTable_one_simu, obs_label == 1)

#### Data Augmentation ####
prev_ub <- 0.20
M_aug <- round(P*prev_ub,digits = 0)
aug_size <- sum(M_aug) # The augmented size (DA size)
pseudo_prev <- N_target/M_aug # Effectively what is the true value of prevalence after DA
logit_pseudo_prev <- logit(pseudo_prev) # Effectively what is the true value of lambda after DA

all_zeros_size <- M_aug - N_obs # N of additional zero rows by location
dt_all_zeros <- list()
for(g in 1:G){
  dt_all_zeros_grp_g <- array(0, dim = c(all_zeros_size[g],ncol(dtYObs)-1))
  dt_all_zeros_grp_g <- cbind(dt_all_zeros_grp_g, rep(g, nrow(dt_all_zeros_grp_g)))
  dt_all_zeros_grp_g <- as.data.frame(dt_all_zeros_grp_g)
  colnames(dt_all_zeros_grp_g) <- colnames(dtYObs)
  dt_all_zeros[[g]] <- dt_all_zeros_grp_g
}
dt_all_zeros <- do.call(rbind,dt_all_zeros)
# Augmentation to make all locations have untarget + target people
dtYObs_aug <- rbind(dtYObs, dt_all_zeros)
# Order by observed first, and then unobserved, nrow = sum(P)
dtYObs_aug <- arrange(dtYObs_aug, desc(obs_label))
dtYObs_yaug <- dtYObs_aug[,seq(J)] # keep the capture history columns 1:J
dtYObs_yaug$grp_label <- dtYObs_aug$grp_label # add the group labeling column to indicate the membership of the groups G
# dtYObs_yaug is the dataset to be passed into the model
########## FINISH DATA MANIPULATION ##########



###################################
### Set up parallel computation ###
###################################

# First you create a cluster using the total amount of cores you have but one to make sure your computer 
# can go on working:
nbcores <- as.numeric(Sys.getenv("NSLOTS")) 
my_cluster <- makeCluster(nbcores)

# Then you wrap your workflow in a function to be run in parallel:
workflow <- function(seed, config_param){
  
  if (!require("nimble")) install.packages("nimble") else (require("nimble", quietly = TRUE)) 
  if (!require("matrixsampling")) install.packages("matrixsampling") else (require("matrixsampling", quietly = TRUE)) 
  
  ################################################
  #### MODEL - Assuming Exchangeable RE Model ####
  ################################################
  ### Create a NIMBLE model ###
  ## Define the model code, all nodes, variables and relationship, and all constants
  Mod_Code <- nimbleCode( {
    ##################
    ## Hyper-priors ##
    ##################
    
    # variance parameter of the theta #
    ## Option 1: Variance ~ Inverse Gamma, precision ~ Gamma
    # prec_theta ~ dgamma(shape = 0.5, scale = 1/0.5) # inverse variance of theta prior
    # sigma_theta <- 1/prec_theta
    ## Option 2: SD ~ Half-Cauchy (t distribution with df = 1)
    sigma_theta ~ T(dt(mu = 0, tau = sd_HC_sigma_theta, df = 1), 0, )  # half-t truncated at 0 and set no upper limit
    
    ############
    ## Priors ##
    ############
    # Prevalence #
    for(g in 1:G){
      prev_g[g] ~ dunif(0,1)  # for each group g, assign uniform prior for the prevalence
    }
    
    # mu_jg = mu_j, FE across lists, the same for groups within a list #
    # ignoring the group heterogeneity within a list #
    for(j in 1:J){ 
      mu_j[j] ~ dnorm(mean = mean_mu, sd = sigma_mu)  # for each list j, all G grps have the same value mu_j
    }
    
    # Factor Loading Matrix - alpha FE #
    for(r in 1:J){
      for(c in 1:D){
        loading_mat[r,c] ~ dnorm(0, sd = sigma_loading)
        loading_mat_upper_tri[r,c] <- loading_mat[r,c]*mask_loading_mat[r,c]
      }
    }
    
    ################
    ## Likelihood ##
    ################
    # Prevalence Level #
    
    for (i in 1:M){ # M = augmented large population
      
      z[i] ~ dbern(prev_g[grp_index[i]]) # inclusion indicator to be a target subject, same location's subjects share the same pseudo prev
      
      if(D == 1){
        theta[i] ~ dnorm(0, sd = sigma_theta) # vector of individual heterogeneity effect is a random (latent) effect   
        loading_latent_factor[i,1:J] <- theta[i] * t(loading_mat_upper_tri[1:J,1]) # dim = 1 * trans(J*D)
      }
      
      for(j in 1:J){
        # detection prob model #
        logit(p.detect[i,j]) <- mu_j[j] + loading_latent_factor[i,j] # detection probability for list j and additional individual i's heterogeneity
        
        # compute effective p #
        p.eff[i,j] <- z[i] * p.detect[i,j]  # being a potential target * being observed
        y[i,j] ~ dbern(p.eff[i,j]) # individual capture histories  
      } # j
    } # i
    
    #########################
    ## Derive the quantity ##
    #########################
    # <1> Number of target people
    N <- sum(z[1:M]) 
  })
  ########## FINISH BAYESIAN MODEL WRITTING ##########
  
  run_stepwise_mcmc <- "yes" # "yes", if run stepwise, then the initial values cannot be specified with nested list.
  
  ## Use only when parallel computation ##
  J <- config_param$Constants$J
  D <- config_param$Constants$D
  G <- config_param$Constants$G
  aug_size <- config_param$Constants$M
  
  # modelfit_functions.r, write here because this has been wrapped into a function
  genLoading_mat_init <- function(seed, J, D){
    set.seed(seed = seed)
    loading_mat_init <- matrix(rnorm(n = J*D, mean = 0, sd = 1), nrow = J, ncol = D)
    loading_mat_init <- round(loading_mat_init, digits = 3)
    loading_mat_init[upper.tri(loading_mat_init)] <- 0  
    if(D == 1){
      loading_mat_init[1,1] <- 1
    }
    return(loading_mat_init)
  }
  
  ## Create initial values for individual latent factors ##
  genTheta_init <- function(seed, aug_size = aug_size, D){
    if(D == 1){
      set.seed(seed)
      theta <- rnorm(n = aug_size, mean = 0, sd = 0.1)
    }
    if(D > 1){
      set.seed(seed)
      theta <-  matrix(rnorm(n = aug_size*D, mean = 0, sd = 0.1), nrow = aug_size, ncol = D)
    }
    return(theta)
  }
  
  ### Configuration of MCMC ###
  ## Initials
  set.seed(123) # ensure reproducibility while drawing randomly initial values.
  
  ## Initials
  ## lambda and prev_g both refer to pseudo prevalence
  if(run_stepwise_mcmc == "yes"){
    Mod_Inits <- list(prev_g = runif(G, 0.05, 0.2), 
                      mu_j = rep(0,J), 
                      sigma_theta = runif(1, 1, 5),
                      loading_mat = genLoading_mat_init(seed = 234, J = J, D = D),
                      theta = genTheta_init(seed = 234, aug_size = aug_size, D = D),
                      z = rep(1, aug_size)
    )
  }
  
  #############################
  ### Run the MCMC sampling ###
  #############################
  ## Run MCMC step by step ##
  ## Build model ##
  Mod_HomoBHMFA <- nimbleModel(code = Mod_Code, name = "Mod_HomoBHMFA",
                               constants = Mod_Consts,
                               data = Mod_Data, inits = Mod_Inits,
                               calculate = FALSE # disable the calculation of all deterministic nodes and log-likelihood
  )
  ## Compile the nimble model ##
  cMod_HomoBHMFA <- compileNimble(Mod_HomoBHMFA)
  # ## (1) Configuration
  Mod_HomoBHMFAConf <- configureMCMC(Mod_HomoBHMFA, 
                                     useConjugacy = FALSE, # because the model is anyway not conjugacy, no need to check, speed up
                                     monitors = c("N","prev_g", "mu_j", "loading_mat_upper_tri"),
                                     enableWAIC = nimbleOptions('enableWAIC' = TRUE),
                                     print = TRUE)
  
  ## Note about enableWAIC
  # A logical argument, specifying whether to enable WAIC calculations for the resulting MCMC algorithm. 
  # Defaults to the value of nimbleOptions('MCMCenableWAIC'), which in turn defaults to FALSE. 
  # Setting nimbleOptions('enableWAIC' = TRUE) will ensure that WAIC is enabled for all calls to configureMCMC and buildMCMC.
  
  # ## (2) Build MCMC algorithm for the model
  Mod_HomoBHMFAMCMC <- buildMCMC(Mod_HomoBHMFAConf,
                                 useConjugacy = FALSE # disable the search for conjugate samplers
  )
  # ## Test for debugging
  # Mod_HomoBHMFAMCMC$run(1)
  # ## (3) Compile C++ after configuration
  C_mcmc_Mod_HomoBHMFA <- compileNimble(Mod_HomoBHMFAMCMC, project = Mod_HomoBHMFA) 
  # ## (4) Run MCMC sampling
  samples_Mod_HomoBHMFA <- runMCMC(C_mcmc_Mod_HomoBHMFA, 
                                   niter = config_param$niter, 
                                   nburnin = config_param$nburnin, 
                                   nchains = config_param$nchains,  # comment this line because we set three chains outside
                                   thin = config_param$nthin, 
                                   setSeed = seed,  # control the state of R the random number generator for reproducibility
                                   # progressBar = TRUE,
                                   summary = TRUE,
                                   WAIC = TRUE,
                                   samplesAsCodaMCMC = TRUE)
  ########## END OF RUNNING MODELS ##########
  
  
  return(samples_Mod_HomoBHMFA)
  
} ######## END OF PARALLEL COMPUTATION WRAP ##########


### Now we run the code using parLapply(), which uses cluster nodes to execute our workflow:

## Configuration of MCMC sampling
nchains <- 3
niter <- 500000 
nburnin <- 400000 
nthin <- 100

## Create a mask matrix to obtain loading factor matrix with upper tri = 0 ## 
mask_loading_factor_mat <- matrix(1, nrow = J, ncol = D)
mask_loading_factor_mat[upper.tri(mask_loading_factor_mat)] <- 0

### Configuration of MCMC ###
### Constants, data, initials, parameter monitors ###
## Constants ##
Mod_Consts <- list(M = aug_size,
                   J = J,
                   G = G,
                   D = D,
                   mean_mu = 0, # mean for homogeneous(across G groups) list effect mu_j 
                   sigma_mu = 10, # sd for homogeneous(across G groups) list effect mu_j 
                   sd_HC_sigma_mu = 10, # standard deviation of half-cauchy prior distribution for variance parameter sigma_mu for all j and g
                   grp_index =  dtYObs_yaug$grp_label, # location label for each individual
                   sd_HC_sigma_theta = 10, # standard deviation of half-cauchy prior distribution for variance parameter sigma_theta for i
                   sigma_loading = 10, # standard deviation of prior for alpha_jd
                   mask_loading_mat = mask_loading_factor_mat  # TRUE indicates upper triangulation of the factor loading matrix
)

## Data, pass in y ##
Mod_Data <- list(y = dtYObs_yaug[,seq(J)]) 

## Put together into a list ##
Mod_HomoBHMFA_config <- list(nchains = nchains,
                             niter = niter,
                             nburnin = nburnin,
                             nthin = nthin,
                             Constants = Mod_Consts,
                             Data = Mod_Data)
########## FINISH MODEL PIECES ##########
## Run parallel versions of lapply
output <- parLapply(cl = my_cluster, 
                    X = 9, 
                    fun = workflow, 
                    config_param = Mod_HomoBHMFA_config ## Configuration, Data and Constants ##
)

# Close the cluster with stopCluster() 
# so that processes do not continue to run in the background and slow down other processes:
stopCluster(my_cluster)

########## END WITH PARALLEL COMPUTATION ###########



########## Post-Process MCMC Samples ############
# length(output) = 1
samples_Mod_HomoBHMFA <- output[[1]] # length = 1


#############################################
### Save the outputs and plot the results ###
#############################################
chain1 <- samples_Mod_HomoBHMFA$samples$chain1
chain2 <- samples_Mod_HomoBHMFA$samples$chain2
chain3 <- samples_Mod_HomoBHMFA$samples$chain3
# each has dim = nsamples*(1 for N + dim for loading_mat_upper_tri[1, 1], [2, 1], [3, 1], 
#                          + dim for mu_j[1], [2], [3], ...[J],
#                          + prev_g[1]...[G])

names_chain1 <- colnames(chain1)
names_chain2 <- colnames(chain2)
names_chain3 <- colnames(chain3)

nsamples <- (niter - nburnin)/nthin

### Note: Below functions are designed for run_stepwise_mcmc == "yes" ###

###################
## Target size N ##
###################
samp1_N <- chain1[,'N']
samp2_N <- chain2[,'N']
samp3_N <- chain3[,'N']
ls_samp_N <- list(chain1 = samp1_N,
                  chain2 = samp2_N,
                  chain3 = samp3_N)

## Combine chains ##
samp_combo_N <- c(samp1_N, samp2_N, samp3_N)

################
## Prevalence ##
################
# Keep only lambda parameters' samp
samp1_prev_g <- save_prev_g(samples = chain1, names_in_chain = names_chain1)
samp2_prev_g <- save_prev_g(samples = chain2, names_in_chain = names_chain2)
samp3_prev_g <- save_prev_g(samples = chain3, names_in_chain = names_chain3)
ls_samp_prev_g <- list(chain1 = samp1_prev_g,
                       chain2 = samp2_prev_g,
                       chain3 = samp3_prev_g)

## Effective prevalence ##
samp1_effprev <- samp1_prev_g * prev_ub
samp2_effprev <- samp2_prev_g * prev_ub
samp3_effprev <- samp3_prev_g * prev_ub

## Combine chains for G groups
samp_combo_effprev <- as.data.frame(rbind(samp1_effprev, samp2_effprev, samp3_effprev))


####################
## Lists' Effects ##
####################
samp1_mu_j <- save_mu_j(samples = chain1, names_in_chain = names_chain1) # dim = nsamples*(J*G)
samp2_mu_j <- save_mu_j(samples = chain2, names_in_chain = names_chain2) # dim = nsamples*(J*G)
samp3_mu_j <- save_mu_j(samples = chain3, names_in_chain = names_chain3) # dim = nsamples*(J*G)
ls_mu_j <- list(chain1 = samp1_mu_j,
                chain2 = samp2_mu_j,
                chain3 = samp3_mu_j)

## Combine chains for G groups
samp_combo_mu_j <- as.data.frame(rbind(samp1_mu_j, samp2_mu_j, samp3_mu_j))
colnames(samp_combo_mu_j) <- paste("mu_j", 1:J, sep = "")


####################################
## Loading factors upper triangle ##
####################################
samp1_loading_upper_tri <- save_loading_mat(samples = chain1, names_in_chain = names_chain1)
samp2_loading_upper_tri <- save_loading_mat(samples = chain2, names_in_chain = names_chain2)
samp3_loading_upper_tri <- save_loading_mat(samples = chain3, names_in_chain = names_chain3)

## Combine chains for G groups
samp_combo_loading <- as.data.frame(rbind(samp1_loading_upper_tri, samp2_loading_upper_tri, samp3_loading_upper_tri))
colnames(samp_combo_loading) <- paste("loading_factor_j", 1:J, sep = "")

################################## 
### Save the samples from MCMC ###
##################################
## Save the output samples N ##
chain_index <- rep(1:nchains,each = nsamples)
combodt_samp_N <- data.frame(Chains = chain_index,
                             Combo_samp_N = samp_combo_N)
colnames(combodt_samp_N) <- c("Chains" ,"Ntarget")

### Save the computed prevalence ###
chain_index <- rep(1:nchains,each = nsamples)
combodt_samp_prev <- data.frame(Chains = chain_index,
                                Combo_samp_prev = samp_combo_effprev)
colnames(combodt_samp_prev) <- c("Chains", paste("effprev_g",seq(1,G),sep=""))

### Save the computed mu_j ###
chain_index <- rep(1:nchains,each = nsamples)
combodt_samp_mu_j <- data.frame(Chains = chain_index,
                                Combo_samp_mu_j = samp_combo_mu_j)
colnames(combodt_samp_mu_j) <- c("Chains", paste("mu_j", 1:J, sep = ""))

### Save the computed loading factors ###
chain_index <- rep(1:nchains,each = nsamples)
combodt_samp_loading_factor <- data.frame(Chains = chain_index,
                                          Combo_samp_loading_factor = samp_combo_loading)
colnames(combodt_samp_loading_factor) <- c("Chains", paste("loading_factor_j", 1:J, sep = ""))


## Extract the WAIC
mcmc.out.WAIC <- samples_Mod_HomoBHMFA$WAIC

# ## Extract all chains summary
mcmc.out.summary <- samples_Mod_HomoBHMFA$summary$all.chains

## Save the true values
TrueVals <- list(TotalN_target = TotalN_target,
                 TotalN_obs = TotalN_obs,
                 N_target = N_target,
                 Prev = Prev,
                 PDetect_by_list_grp_val = dtPDetect_by_list_grp_val,
                 Empirical_mu_jg_mat = empirical_mu_jg_mat,
                 Empirical_pdetect_jg_mat = empirical_pdetect_mat,
                 Theoretical_mu_jg_mat = theoretical_mu_jg_mat,
                 Theoretical_pdetect_jg_mat = theoretical_pdetect_mat,
                 alpha = dtBasicSetUp_one_simu$alpha
)

## Save MCMC Setup 
NIMBLE_Setup <- list(case = case,
                     Fix_FirstElement_Alpha = fix_first_element_in_alpha,
                     Parallel_Run = parallel_run,
                     nchains = nchains,
                     niter = niter,
                     nburnin = nburnin,
                     nthin = nthin,
                     nsamples = nsamples,
                     mean_mu = Mod_Consts$mean_mu, # prior for homogeneous mu_j
                     sigma_mu = Mod_Consts$sigma_mu, # prior for homogeneous mu_j
                     sd_mean_mu = Mod_Consts$sd_mean_mu, # prior value for list-specific fixed effect variance
                     sd_HC_sigma_mu = Mod_Consts$sd_HC_sigma_mu, # prior value for group-specific RE half-cauchy variance
                     sd_HC_sigma_theta = Mod_Consts$sd_HC_sigma_theta, # standard deviation of half-cauchy prior distribution for variance parameter sigma_theta for i
                     sigma_loading = Mod_Consts$sigma_loading # sigma of alpha
)

## Save data
saveRDS(samples_Mod_HomoBHMFA$samples, file = paste0(path_to_output, "mcmc.out_HomoBHMFA_samples", "_simu", task_id, ".Rdata"))
saveRDS(combodt_samp_N, file = paste0(path_to_output, "mcmc.out_HomoBHMFA_N", "_simu", task_id, ".Rdata"))
saveRDS(combodt_samp_prev, file = paste0(path_to_output, "mcmc.out_HomoBHMFA_prev", "_simu", task_id, ".Rdata"))
saveRDS(combodt_samp_mu_j, file = paste0(path_to_output, "mcmc.out_HomoBHMFA_mu_j", "_simu", task_id, ".Rdata"))
saveRDS(combodt_samp_loading_factor, file = paste0(path_to_output, "mcmc.out_HomoBHMFA_loadings", "_simu", task_id, ".Rdata"))
saveRDS(mcmc.out.WAIC, file = paste0(path_to_output, "mcmc.out_HomoBHMFA_WAIC", "_simu", task_id, ".Rdata"))
saveRDS(mcmc.out.summary, file = paste0(path_to_output, "mcmc.out_HomoBHMFA_summary", "_simu", task_id, ".Rdata"))
saveRDS(TrueVals, file = paste0(path_to_output, "TrueVals", "_simu", task_id, ".Rdata"))
saveRDS(NIMBLE_Setup, file = paste0(path_to_output, "NIMBLE_Setup", "_simu", task_id, ".Rdata"))
