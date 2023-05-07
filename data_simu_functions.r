#########################################################################################################
######## Project 3 Hierarchical Demographic Groups with CRC estimator for population size estimation #############
########################### Utility Functions ###########################################################

# library("copula")
if (!require("mvtnorm")) install.packages("mvtnorm") else (require("mvtnorm", quietly = TRUE)) 
if (!require("matrixsampling")) install.packages("matrixsampling") else (require("matrixsampling", quietly = TRUE)) 
if (!require("LaplacesDemon")) install.packages("LaplacesDemon") else (require("LaplacesDemon", quietly = TRUE)) 
if (!require("spam")) install.packages("spam") else (require("spam", quietly = TRUE)) 

######################################
### Model 2: detection probability ###
######################################
### Generate list effects mu_j
genListEffects <- function(p_j){
  # p_j = wanted detection probabilities
  logit_p_j <- log(p_j/(1-p_j))
  return(logit_p_j)
}

### Generate alpha - List interaction with Individual Heterogeneity
genAlpha <- function(J, alpha_sd = 100, no_interact = TRUE){ 
  # J = number of Datasets
  # alpha_sd = sd of loading factors
  # no_interact = no behavioral response
  
  if(no_interact){
    mat <- diag(J)
  }
  else{
    mat <- diag(J)
    n <- J*(J-1)/2
    alpha_jrs <- rnorm(n, 0, alpha_sd)
    mat[lower.tri(mat, diag = FALSE)] <- alpha_jrs 
  }
  return(mat) # loading matrix of alpha
}

### Generate individual heterogeneity vector 
genTheta <- function(D = n_latent_factor, sample_size, theta_mean = theta_val, theta_sd = theta_sd_val, theta_corr = theta_corr_val, no_interact = FALSE){
  # p = number of latent factors 
  # sample_size = number of individuals need to draw values
  # theta_mean = mean of individual heterogeneity latent factors, dim = 1*D
  # theta_sd = sd of individual heterogeneity latent factors, dim = 1*D
  # theta_corr = correlation of latent factors, dim = 1*D
  # no_interact = logical operation to have interaction between list and individual heterogeneity
  if(D == 1){
    if(!no_interact){
      theta_vec <- rnorm(sample_size, theta_mean, sd = theta_sd) # Generate theta for all individuals
      theta_mat <- matrix(rep(theta_vec, each = D), ncol = D, byrow = TRUE) # Repeat for J lists
    }
  }
  if(D > 1){
    if(no_interact){ 
      # Model Mt+h
      theta_vec <- rnorm(sample_size, theta_mean, sd = theta_sd) # Generate theta for all individuals
      theta_mat <- matrix(rep(theta_vec, each = D), ncol = D, byrow = TRUE) # Repeat for J lists
    }
    else{
      # Generalized Mth
      # theta_Sigma = Covariance matrix of the latent factors
      # option 1: generate theta Sigma and use inverse wishart distribution
      if(is.null(theta_sd) & is.null(theta_corr)){
        theta_Sigma_array <- rinvwishart(n = 1, nu = D+1, Omega = diag(D))
        theta_Sigma <- apply(theta_Sigma_array, 2, c)  
      }
      if(!is.null(theta_sd) & !is.null(theta_corr)){
        # option 2: import fixed theta Sigma and use multivariate normal distribution
        theta_cov <- theta_sd*theta_sd*theta_corr # compute covariance between the latent factors
        theta_vcov <- diag(D)
        theta_vcov <- ifelse(theta_vcov == 0, theta_cov, theta_vcov) # covariance between latent factors 
        diag(theta_vcov) <-  theta_sd^2 # variance of each latent factor
        theta_Sigma <- theta_vcov
      }
      theta_mean <- theta_mean
      theta_mat <- rmvnorm(sample_size, theta_mean, sigma = theta_Sigma)
    } 
  }  
  return(theta_mat)
} # theta_val matrix for each G

### Get individual-based detection probability
getDetectProb <- function(n_target, mu_jg, alpha, theta){ 
  # n_target = N_targets[g]
  # mu_jg = mu_jg as a list, each row represents a J
  # alpha = Alpha matrix for factor loading, dim = J*D (n lists * n latent factors)
  # theta = Individual heterogeneity for each k, dim = P_val_k * D  
  
  if(!is.null(alpha)){
    ## With individual heterogeneity, whether interaction or not is controlled inside of getAlpha
    theta_target <- theta[sample(nrow(theta), n_target),] # only keep the target population, dim = n_target*D
    theta_target <- matrix(theta_target, nrow = n_target, ncol = ncol(theta)) # make sure the dim of the theta for target population = n_target * D
    theta_target <- t(theta_target) # dim = D*n_target
    alpha <- as.matrix(alpha)
    loading_mat <- t(alpha%*%theta_target) # dim = n_target*J
    mu_jg_mat <- matrix(rep(mu_jg, each = n_target), nrow = n_target, byrow = FALSE) # dim = n_target * J
    p <- plogis(mu_jg_mat + loading_mat) # dim = number of OUD at a single G (n_target) * number of lists (J)
  }
  if(is.null(alpha) & is.null(theta)){
    ## Only lists effect
    mu_jg_mat <- matrix(rep(mu_jg, each = n_target), nrow = n_target, byrow = FALSE) # dim = n_target * J
    p <- plogis(mu_jg_mat) # dim = number of OUD at single K (n_target) * number of lists (J)  
  }
  
  return(p) # individual detection probability across lists at single group g
} # p_detect for each G


##################################################################
### Bernoulli draw full of all targets across lists and groups ###
##################################################################
genYfull <- function(p){
  # p = p_detect, for one group g, given p matrix sample 0/1 across J lists
  
  p_vec <- as.vector(p) # append by columns
  yfull_vec <- rbinom(n = length(p_vec), size = 1, prob = p_vec) # Across all detection prob, draw 0/1 for each
  yfull_mat <- matrix(yfull_vec, nrow = nrow(p), ncol = ncol(p), byrow = FALSE) # convert back to p matrix
  return(yfull_mat) # detection history in one location across lists
}

### Manipulate full of all target data shape 
makeYfullTable <- function(yfull){ 
  # yfull = yfull list of yfull_mat, list of all target people, some observed, some not
  
  datalist <- list()
  for(i in 1:length(yfull)) {
    yfull_i <- yfull[[i]]
    obs_label <- apply(yfull_i, 1, max) # Mark if observed or not
    dt_i <- data.frame(yfull_i, target_label = 1, obs_label = obs_label, grp_label = i) 
    datalist[[i]] <- dt_i  
  }
  YfullTable <- do.call(rbind,datalist)
  return(YfullTable)
}
