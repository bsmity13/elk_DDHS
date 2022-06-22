############################################X
#--Elk Density-Dependent Habitat Selection--X
#---------------Brian J. Smith--------------X
#----------Bayesian Model Checking----------X
#===========================================X
#-----------Last update 2022-02-22----------X
############################################X

#Load packages----
library(tidyverse)

# Custom functions ----
source("99_fun.R")

#Load data----
# Original data with predictions
# Note, this was created in script "02_figures.R", but the folder "pred"
# is added to .gitignore. This object will not exist in cloned repos unless
# script 02 has already been run.
dat <- readRDS("pred/model_data_predictions.rds")

#Times
years <- sort(unique(dat$year))

dat$t <- as.numeric(factor(dat$year))

# MCMC samples ----
# Note: the model runs are added to .gitignore, so are not available
# in the GitHub repo (samples from final run are > 1.2 GB)

# For demonstration of the code (without asking the user to re-fit the model),
# I have provided a very small subset of the original chains. See script 
# '99_demo_chains.R' for more details.

# This switch determines whether to load the real data or the demo data
demo_data <- FALSE

if (demo_data) {
  # Demo data
  # 33 samples from each chain (99 total samples per node)
  
  samples1 <- readRDS("demo_models/example_samples1.rds")
  samples2 <- readRDS("demo_models/example_samples2.rds")
  
} else {
  # Real data
  
  samples1 <- list(
    ch1 = readRDS("models/samples1_ch1_2022-02-22.rds"),
    ch2 = readRDS("models/samples1_ch2_2022-02-22.rds"),
    ch3 = readRDS("models/samples1_ch3_2022-02-22.rds")
  )
  
  samples2 <- list(
    ch1 = readRDS("models/samples2_ch1_2022-02-22.rds"),
    ch2 = readRDS("models/samples2_ch2_2022-02-22.rds"),
    ch3 = readRDS("models/samples2_ch3_2022-02-22.rds")
  )
  
  # Get rid of burn-in and thin
  # Burn-in of 20k (samplers adapt for 15k)
  # Thinning by 20 leaves 4k samples for inference
  samples1 <- lapply(samples1, function(x) {
    return(x[seq(20020, nrow(x), by = 20),])
  })
  
  samples2 <- lapply(samples2, function(x) {
    return(x[seq(20020, nrow(x), by = 20),])
  })
  
  # Garbage cleanup
  gc()
}

# Null model samples ----

samples_null <- list(
  ch1 = readRDS("models/samples1_null_ch1_2022-01-12.rds"),
  ch2 = readRDS("models/samples1_null_ch2_2022-01-12.rds"),
  ch3 = readRDS("models/samples1_null_ch3_2022-01-12.rds")
)

# Get rid of burn-in and thin
# Burn-in of 20k (samplers adapt for 15k)
# Thinning by 20 leaves 4k samples for inference
samples_null <- lapply(samples_null, function(x) {
  return(x[seq(20020, nrow(x), by = 20), , drop = FALSE])
})

# Get coefficients ----
b <- lapply(samples1, function(x) {
  return(as.data.frame(x[1:nrow(x), 
                         grep("beta", colnames(samples1[[1]]))]))
}) %>% 
  bind_rows() %>% 
  as.data.frame()
names(b) <- paste0("b", 1:ncol(b))

e <- lapply(samples1, function(x) {
  return(as.data.frame(x[1:nrow(x), 
                         grep("eta_out", colnames(samples1[[1]]))]))
}) %>% 
  bind_rows() %>% 
  as.data.frame()
names(e) <- paste0("e", 1:ncol(e))

mu_eta <- unname(do.call(c, lapply(samples1, function(x) {
  return(x[, "mu_out"])
})))

s <- lapply(samples2, function(x) {
  return(as.data.frame(x[1:nrow(x), ]))
})  %>% 
  bind_rows() %>% 
  as.data.frame()
names(s) <- paste0("s", 1:ncol(s))

r <- unname(do.call(c, lapply(samples1, function(x) {
  return(x[, "nb_r"])
})))

# Null
r_null <- unname(do.call(c, lapply(samples_null, function(x) {
  return(x[, "nb_r"])
})))

# Number of posterior samples ----
npost <- length(r)

# Predict lambda for each iteration
lams <- lapply(1:npost, function(ii) {
  cat("\r", ii, "of", npost, "    ")
  
  lams <- predict_lambda_iter(dat = dat, i = ii,
                              b = b,
                              e = as.matrix(e),
                              s = as.matrix(s))
  
  return(lams)
})

# Calculate log-likelihood for each iteration ----
ll <- lapply(1:npost, function(ii) {
  
  cat("\r", ii, "of", npost, "    ")

  
  # ... negative binomial 'p' ----
  p_all <- r[ii]/(r[ii] + lams[[ii]])
  p_null <- r_null[ii]/(r_null[ii] + exp(dat$log_dens))
  
  # ... log-likelihood ----
  ll_all <- sum(dnbinom(dat$n, size = r[ii], prob = p_all, log = TRUE))
  ll_null <- sum(dnbinom(dat$n, size = r_null[ii], prob = p_null, log = TRUE))
  
  # ... iteration log-likelihood ----
  res <- data.frame(iter = ii, 
                    all = ll_all,
                    null = ll_null)
  return(res)
}) %>% 
  bind_rows()


# Pseudo-R2 ----
# McFadden's
PR2_McF <- 1 - (ll$all/ll$null)

mean(PR2_McF)
plot(density(PR2_McF), main = "McFadden's R2")

# Cox-Snell method
# See Eqn 1b in Nagelkerke 1991 (Biometrica)
PR2_CS <- 1 - exp((-2/nrow(dat)) * (ll$all - ll$null))

mean(PR2_CS)
plot(density(PR2_CS), main = "Cox-Snell R2")

# Nagelkerke method
# See Eqn 3 in Nagelkerke 1991 (Biometrica)
PR2_Nag <- PR2_CS/(1 - (exp(2/nrow(dat) * ll$null)))

mean(PR2_Nag)
plot(density(PR2_Nag), main = "Nagelkerke R2")
quantile(PR2_Nag, c(0.05, 0.95))
