############################################X
#--Elk Density-Dependent Habitat Selection--X
#---------------Brian J. Smith--------------X
#---Model Fitting Function for parLapply()--X
#===========================================X
#-----------Last update 2021-11-29----------X
############################################X

#Load NIMBLE if necessary
library(nimble)

#Custom NIMBLE functions----

# ... lambda function ----
calc_lambda <- nimbleFunction(
  run = function(beta = double(1),
                 eta_out = double(0),
                 log_dens = double(0),
                 outside = double(1),
                 swe = double(1),
                 elev = double(1),
                 cos_asp = double(1),
                 sin_asp = double(1),
                 biomass = double(1),
                 biomass_logdens = double(1),
                 open = double(1),
                 rough = double(1),
                 open2 = double(1),
                 rough2 = double(1),
                 open_wolf = double(1),
                 rough_wolf = double(1),
                 open2_wolf = double(1),
                 rough2_wolf = double(1),
                 open_cougar = double(1),
                 rough_cougar = double(1),
                 open2_cougar = double(1),
                 rough2_cougar = double(1),
                 open_logdens = double(1),
                 rough_logdens = double(1),
                 open2_logdens = double(1),
                 rough2_logdens = double(1),
                 s = double(1),
                 Omega = double(0)) {
    
    returnType(double(1))
    
    lambda <- numeric(length = Omega, init = FALSE)
    
    for (i in 1:Omega) {
      lambda[i] <- exp(
        # Offset of log_dens
        log_dens +
          # Random slope for outside YNP
          eta_out * outside[i] +
          # Landscape conditions
          beta[1] * swe[i] +
          beta[2] * elev[i] +
          beta[3] * cos_asp[i] +
          beta[4] * sin_asp[i] +
          # Resources
          beta[5] * biomass[i] +
          beta[6] * biomass_logdens[i] +
          # Risk main effect
          beta[7] * open[i] +
          beta[8] * open2[i] +
          beta[9] * rough[i] +
          beta[10] * rough2[i] +
          # Risk from wolves
          beta[11] * open_wolf[i] +
          beta[12] * open2_wolf[i] +
          beta[13] * rough_wolf[i] +
          beta[14] * rough2_wolf[i] +
          # Risk from cougars
          beta[15] * open_cougar[i] +
          beta[16] * open2_cougar[i] +
          beta[17] * rough_cougar[i] +
          beta[18] * rough2_cougar[i] +
          # Risk density-dependence
          beta[19] * open_logdens[i] +
          beta[20] * open2_logdens[i] +
          beta[21] * rough_logdens[i] +
          beta[22] * rough2_logdens[i] +
          # Spatial autocorrelation
          s[i]
        )
    }
    
    return(lambda)
  }
)

# ... null model lambda function ----
calc_lambda_null <- nimbleFunction(
  run = function(log_dens = double(0),
                 Omega = double(0)) {
    
    returnType(double(1))
    
    lambda <- numeric(length = Omega, init = FALSE)
    
    for (i in 1:Omega) {
      lambda[i] <- exp(
        # Offset of log_dens
        log_dens)
    }
    
    return(lambda)
  }
)

# ... Exponential covariance function ----
# Copied from here:
# https://r-nimble.org/nimbleExamples/gaussian_process.html
expcov <- nimbleFunction(     
  run = function(dists = double(2), rho = double(0), sigma = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    sigma2 <- sigma*sigma
    for(i in 1:n)
      for(j in 1:n)
        result[i, j] <- sigma2*exp(-dists[i,j]/rho)
    return(result)
  })

# cExpcov <- compileNimble(expcov)
# cExpcov(elkConst$dists, rho = 0.2, sigma = 1)

# Custom distributions ----

# ... vectorized NB ----
dnegbin_vec <- nimbleFunction(
  run = function(x = double(1), 
                 lambda = double(1), 
                 r = double(0),
                 Omega = double(0),
                 log = logical(0, default = 0)) {
    
    returnType(double(0))
    
    # Calculate p parameter of negative binomial
    # Useful reference: https://georgederpa.github.io/teaching/countModels.html
    # p[i] = r/(r + lambda[i])
    p <- nimNumeric(length = Omega, init = FALSE)
    p[1:Omega] <- r/(r + lambda[1:Omega])
    
    # Initialize logProb
    logProb <- 0
    
    for (i in 1:Omega) {
      logProb <- logProb + dnbinom(x[i], prob = p[i], size = r, log = TRUE)
    }
    
    if (log) {
      return(logProb)
    } else {
      return(exp(logProb))
    }
  })

# # Compile
# C_dnegbin_vec <- compileNimble(dnegbin_vec)
# 
# # Test
# xx = c(0, 25, 11)
# ll = c(12, 12, 12)
# rr = 0.5
# C_dnegbin_vec(x = xx,
#        lambda = ll,
#        r = rr,
#        Omega = length(xx),
#        log = FALSE)

rnegbin_vec <- nimbleFunction(
  run = function(n = integer(0), 
                 lambda = double(1), 
                 r = double(0),
                 Omega = double(0)) {
    
    returnType(double(1))
    
    if(n != 1) {
      print("rnegbin_vec only allows n = 1; using n = 1.")
      n = 1
    }
    
    res <- numeric(Omega)
    p <- numeric(Omega)
    
    for (i in 1:Omega) {
      p[i] <- r/(r + lambda[i])
      res[i] <- rnbinom(1, prob = p[i], size = r)
    }
    
    return(res)
  })


# Model specification----
elk_mod <- nimbleCode({
  
  ########## Priors ########## 
  
  # beta (covariates for HS)
  # Laplace (= double exponential) Prior
  # scale = 1/sqrt(2) = 0.7071068
  for (i in 1:22) {
    beta[i] ~ dlaplace(0, scale = 0.7071068)
  }
  
  # r parameter of Negative Binomial
  nb_r ~ dunif(0, 1)
  
  # Random effect of outside park
  mu_out ~ dlaplace(0, scale = 0.7071068)
  sigma_out ~ dunif(0, 10)
  
  # Gaussian process (spatial autocorrelation)
  # Values taken from https://r-nimble.org/nimbleExamples/gaussian_process.html
  sigma ~ dunif(0, 100)  # prior for variance components based on Gelman (2006)
  rho ~ dunif(0, 100) # Consider (0, 10)
  
  ########## Likelihood ##########
  
  # # Gaussian process
  cov[1:Omega, 1:Omega] <- expcov(dists[1:Omega, 1:Omega], rho, sigma)
  s[1:Omega] ~ dmnorm(zeroes[1:Omega], cov = cov[1:Omega, 1:Omega])
  
  for (t in 1:length_t){
    
    # Random slope for outside YNP
    eta_out[t] ~ dnorm(mu_out, sd = sigma_out)
    
    # Expected abundance
    lambda[1:Omega, t] <- calc_lambda(beta = beta[1:22],
                                      eta_out = eta_out[t],
                                      log_dens = log_dens[t],
                                      outside = outside[1:Omega, t],
                                      swe = swe[1:Omega, t],
                                      elev = elev[1:Omega],
                                      cos_asp = cos_asp[1:Omega],
                                      sin_asp = sin_asp[1:Omega],
                                      biomass = biomass[1:Omega, t],
                                      biomass_logdens = biomass_logdens[1:Omega, t],
                                      open = open[1:Omega, t],
                                      rough = rough[1:Omega, t],
                                      open2 = open2[1:Omega, t],
                                      rough2 = rough2[1:Omega, t],
                                      open_wolf = open_wolf[1:Omega, t],
                                      rough_wolf = rough_wolf[1:Omega, t],
                                      open2_wolf = open2_wolf[1:Omega, t],
                                      rough2_wolf = rough2_wolf[1:Omega, t],
                                      open_cougar = open_cougar[1:Omega, t],
                                      rough_cougar = rough_cougar[1:Omega, t],
                                      open2_cougar = open2_cougar[1:Omega, t],
                                      rough2_cougar = rough2_cougar[1:Omega, t],
                                      open_logdens = open_logdens[1:Omega, t],
                                      rough_logdens = rough_logdens[1:Omega, t],
                                      open2_logdens = open2_logdens[1:Omega, t],
                                      rough2_logdens = rough2_logdens[1:Omega, t],
                                      s = s[1:Omega],
                                      Omega = Omega)
    
    n[1:Omega, t] ~ dnegbin_vec(lambda = lambda[1:Omega, t], 
                                r = nb_r,
                                Omega = Omega)
  }
})

# Null model specification----
elk_mod_null <- nimbleCode({
  
  ########## Priors ########## 
  
  # r parameter of Negative Binomial
  nb_r ~ dunif(0, 1)
  
  ########## Likelihood ##########
  
  for (t in 1:length_t){
    
    # Expected abundance
    lambda[1:Omega, t] <- calc_lambda_null(log_dens = log_dens[t],
                                           Omega = Omega)
    
    n[1:Omega, t] ~ dnegbin_vec(lambda = lambda[1:Omega, t], 
                                r = nb_r,
                                Omega = Omega)
  }
})

# Full MCMC in function ----
# Chain list has two elements: chain number and random seed
# Example
# chain_list <- list(ch1 = list(chain = 1,
#                               seed = 123),
#                    ch2 = list(chain = 2, 
#                               seed = 456))
mod_mcmc <- function(chain_list, data, niter, nburn, nthin, 
                     save_rate = 10000, dir = "models") {
  
  # Make sure dir exists
  dir.create(dir, showWarnings = FALSE)
  dir.create(paste0(dir, "/temp"), showWarnings = FALSE)
  
  # Save rate determines how often intermediate posterior samples are saved to
  # disk during model fitting
  
  # Check that save_rate >= niter
  if(save_rate > niter) {
    stop("'save_rate' must be <= 'niter'.")
  }
  
  # When to save progress?
  save_iter <- seq(from = nburn, to = niter, by = save_rate)
  
  # If last iteration is not niter, then append niter
  if (save_iter[length(save_iter)] != niter) {
    save_iter <- c(save_iter, niter)
  }
  
  # Construct vector of iterations needed
  niter_vec <- c(nburn, diff(save_iter))
  
  if(niter_vec[1] == 0) {
    niter_vec <- niter_vec[-1]
  }
  
  cat("Saving at iterations", paste(cumsum(niter_vec)), "\n")
  
  # ... source helper functions ----
  source("99_fun.R")
  
  # ... source model functions ----
  source("99_nimble_NB_model_function.R")
  
  # ... get chain and seed from list ----
  # (and open cutoff, if applicable)
  chain <- chain_list$chain
  rseed <- chain_list$seed
  
  if (is.null(chain_list$open_cutoff)) {
    oc <- 0
  } else {
    oc <- chain_list$open_cutoff
  }
  
  # ... Prep data ----
  nim_dat <- prep_nimble(dat = data, n_beta = 22, rand_seed = rseed,
                         spatial_effect = TRUE, open_cutoff = oc,
                         open_direction = chain_list$open_direction)
  
  # ... Create model definition ----
  elk <- nimbleModel(code = elk_mod, 
                     data = nim_dat$data,
                     constants = nim_dat$const, 
                     inits = nim_dat$init,
                     check = FALSE, calculate = TRUE)
  
  # ... Configure MCMC ----
  elkConf <- configureMCMC(elk, 
                           monitors = c("beta", 
                                        "sigma", "rho",
                                        "nb_r", 
                                        "mu_out", "sigma_out", 
                                        "eta_out"),
                           thin = 1,
                           monitors2 = "s",
                           thin2 = 1,
                           useConjugacy = FALSE)
  
  # ... ... block sample betas----
  elkConf$removeSampler("beta")
  elkConf$removeSampler("mu_out")
  # Sample conditions together (SWE, Elev, cos(Asp), sin(Asp))
  # As well as mean effect for outside park
  elkConf$addSampler(c(paste0("beta[", 1:4, "]"), "mu_out"), 
                     type = "AF_slice",
                     control = list(sliceWidths = rep(0.01, 
                                                      5),
                                    sliceAdaptFactorInterval = 200,
                                    sliceMaxSteps = 100,
                                    maxContractions = 100,
                                    maxContractionsWarning = FALSE))
  # Sample resources together (biomass, biomass x log(Dens))
  elkConf$addSampler(paste0("beta[", 5:6, "]"), 
                     type = "AF_slice",
                     control = list(sliceWidths = rep(0.01, 
                                                      2),
                                    sliceAdaptFactorInterval = 200,
                                    sliceMaxSteps = 100,
                                    maxContractions = 100,
                                    maxContractionsWarning = FALSE))
  # Sample risks together (open, rough, and all intxns)
  elkConf$addSampler(paste0("beta[", 7:22, "]"), 
                     type = "AF_slice",
                     control = list(sliceWidths = rep(0.01, 
                                                      16),
                                    sliceAdaptFactorInterval = 200,
                                    sliceMaxSteps = 100,
                                    maxContractions = 100,
                                    maxContractionsWarning = FALSE))
  
  
  # ... block sample etas ----
  elkConf$removeSampler("eta_out")
  elkConf$addSampler("eta_out", 
                     type = "AF_slice",
                     control = list(sliceWidths = rep(0.01, 
                                                      nim_dat$const$length_t),
                                    sliceAdaptFactorInterval = 200,
                                    sliceMaxSteps = 100,
                                    maxContractions = 100,
                                    maxContractionsWarning = FALSE))
  
  # ... ... block sample GP parameters ----
  elkConf$removeSampler("sigma")
  elkConf$removeSampler("rho")
  elkConf$addSampler(c("sigma", "rho"), type = "AF_slice",
                     control = list(sliceWidths = rep(0.01, 2),
                                    sliceAdaptFactorInterval = 200,
                                    sliceMaxSteps = 100,
                                    maxContractions = 100,
                                    maxContractionsWarning = FALSE))
  
  # ... ... change sampler on s ----
  elkConf$removeSamplers("s")
  ## reduce the initial proposal covariance scale for better mixing
  elkConf$addSampler("s", 'RW_block', control = list(scale = 0.1))
  
  # ... ... slice sample nb_r ----
  elkConf$removeSampler("nb_r")
  elkConf$addSampler("nb_r", type = "slice",
                     control = list(sliceWidth = 0.01,
                                    adaptInterval = 100,
                                    sliceMaxSteps = 100,
                                    maxContractions = 100,
                                    maxContractionsWarning = FALSE))
  
  # ... ... slice sample sigma_out ----
  elkConf$removeSampler("sigma_out")
  elkConf$addSampler("sigma_out", type = "slice",
                     control = list(sliceWidth = 0.01,
                                    adaptInterval = 100,
                                    sliceMaxSteps = 100,
                                    maxContractions = 100,
                                    maxContractionsWarning = FALSE))
  
  # ... Build and compile ----
  
  #Build the MCMC
  elkMCMC <- buildMCMC(elkConf)
  
  #Compile model
  Celk <- compileNimble(elk)
  
  #Compile MCMC
  CelkMCMC <- compileNimble(elkMCMC)
  
  # Run
  for (i in 1:length(niter_vec)) {
    if (i == 1) {
      reset <- TRUE
    } else {
      reset <- FALSE
    }
    
    CelkMCMC$run(niter = niter_vec[i], 
                 nburnin = 0,
                 thin = nthin, 
                 thin2 = nthin,
                 reset = reset)
    
    # Get samples
    samples1 <- as.matrix(CelkMCMC$mvSamples)
    samples2 <- as.matrix(CelkMCMC$mvSamples2)
    
    # If i == 1, we are still in burnin
    # If i == length(niter_vec), we are saving final chain
    if (i != length(niter_vec)) {
      
      cum_iter <- sum(niter_vec[1:i])
      
      # Save
      saveRDS(samples1, paste0(dir, "/temp/samples1_iter_", cum_iter, 
                               "_ch", chain, "_", 
                               Sys.Date(), ".rds"))
      saveRDS(samples2, paste0(dir, "/temp/samples2_iter_", cum_iter,
                               "_ch", chain, "_",
                               Sys.Date(), ".rds"))
      
    } else {
      if (i == length(niter_vec)) {
        # Save
        saveRDS(samples1, paste0(dir, "/samples1", 
                                 "_ch", chain, "_", 
                                 Sys.Date(), ".rds"))
        saveRDS(samples2, paste0(dir, "/samples2",
                                 "_ch", chain, "_",
                                 Sys.Date(), ".rds"))
      }
    }
  }
}

# Null MCMC in function ----
# Chain list has two elements: chain number and random seed
# Example
# chain_list <- list(ch1 = list(chain = 1,
#                               seed = 123),
#                    ch2 = list(chain = 2, 
#                               seed = 456))
null_mcmc <- function(chain_list, data, niter, nburn, nthin, 
                     save_rate = 10000, dir = "models") {
  
  # Make sure dir exists
  dir.create(dir, showWarnings = FALSE)
  dir.create(paste0(dir, "/temp"), showWarnings = FALSE)
  
  # Save rate determines how often intermediate posterior samples are saved to
  # disk during model fitting
  
  # Check that save_rate >= niter
  if(save_rate > niter) {
    stop("'save_rate' must be <= 'niter'.")
  }
  
  # When to save progress?
  save_iter <- seq(from = nburn, to = niter, by = save_rate)
  
  # If last iteration is not niter, then append niter
  if (save_iter[length(save_iter)] != niter) {
    save_iter <- c(save_iter, niter)
  }
  
  # Construct vector of iterations needed
  niter_vec <- c(nburn, diff(save_iter))
  
  if(niter_vec[1] == 0) {
    niter_vec <- niter_vec[-1]
  }
  
  cat("Saving at iterations", paste(cumsum(niter_vec)), "\n")
  
  # ... source helper functions ----
  source("99_fun.R")
  
  # ... source model functions ----
  source("99_nimble_NB_model_function.R")
  
  # ... get chain and seed from list ----
  # (and open cutoff, if applicable)
  chain <- chain_list$chain
  rseed <- chain_list$seed
  
  if (is.null(chain_list$open_cutoff)) {
    oc <- 0
  } else {
    oc <- chain_list$open_cutoff
  }
  
  # ... Prep data ----
  nim_dat <- prep_nimble(dat = data, n_beta = 1, rand_seed = rseed,
                         spatial_effect = FALSE, open_cutoff = oc,
                         null = TRUE)
  
  # ... Create model definition ----
  elk <- nimbleModel(code = elk_mod_null, 
                     data = nim_dat$data,
                     constants = nim_dat$const, 
                     inits = nim_dat$init,
                     check = FALSE, calculate = TRUE)
  
  # ... Configure MCMC ----
  elkConf <- configureMCMC(elk, 
                           monitors = c("nb_r"),
                           thin = 1,
                           # monitors2 = "lambda",
                           # thin2 = 1,
                           useConjugacy = FALSE)
  
  # ... ... slice sample nb_r ----
  elkConf$removeSampler("nb_r")
  elkConf$addSampler("nb_r", type = "slice",
                     control = list(sliceWidth = 0.01,
                                    adaptInterval = 100,
                                    sliceMaxSteps = 100,
                                    maxContractions = 100,
                                    maxContractionsWarning = FALSE))
  # ... Build and compile ----
  
  #Build the MCMC
  elkMCMC <- buildMCMC(elkConf)
  
  #Compile model
  Celk <- compileNimble(elk)
  
  #Compile MCMC
  CelkMCMC <- compileNimble(elkMCMC)
  
  # Run
  for (i in 1:length(niter_vec)) {
    if (i == 1) {
      reset <- TRUE
    } else {
      reset <- FALSE
    }
    
    CelkMCMC$run(niter = niter_vec[i], 
                 nburnin = 0,
                 thin = nthin, 
                 thin2 = nthin,
                 reset = reset)
    
    # Get samples
    samples1 <- as.matrix(CelkMCMC$mvSamples)
    # samples2 <- as.matrix(CelkMCMC$mvSamples2)
    
    # If i == 1, we are still in burnin
    # If i == length(niter_vec), we are saving final chain
    if (i != length(niter_vec)) {
      
      cum_iter <- sum(niter_vec[1:i])
      
      # Save
      saveRDS(samples1, paste0(dir, "/temp/samples1_null_iter_", cum_iter, 
                               "_ch", chain, "_", 
                               Sys.Date(), ".rds"))
      # saveRDS(samples2, paste0(dir, "/temp/samples2_null_iter_", cum_iter,
      #                          "_ch", chain, "_",
      #                          Sys.Date(), ".rds"))
      
    } else {
      if (i == length(niter_vec)) {
        # Save
        saveRDS(samples1, paste0(dir, "/samples1_null", 
                                 "_ch", chain, "_", 
                                 Sys.Date(), ".rds"))
        # saveRDS(samples2, paste0(dir, "/samples2_null",
        #                          "_ch", chain, "_",
        #                          Sys.Date(), ".rds"))
      }
    }
  }
  # Garbage cleanup in case we're low on memory.
  gc()
}
