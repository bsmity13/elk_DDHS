#######################################X
#----Spatially Explicit Elk Density----X
#----------Started 2019-09-11----------X
#------------Brian J. Smith------------X
#======================================X
#--------Last update 2021-11-29--------X
#######################################X

# The full model runs are stored in a directory called "models", which
# I am adding to .gitignore because samples from final run are > 1.2 GB.

# For the purposes of understanding the code, I am providing a very small
# subset of the final run to the GitHub repo for demonstration purposes. This
# subset will likely have different properties than the final version 
# presented in the manuscript, but can be used to illustrate the code.

# MCMC samples ----
# Load the real final samples (not on GitHub)
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

# Randomly sample 33 rows from each chain
set.seed(123456)

samples1 <- lapply(samples1, function(x) {
  return(x[sample.int(n = 4000, size = 33), ])
})

samples2 <- lapply(samples2, function(x) {
  return(x[sample.int(n = 4000, size = 33), ])
})

# Garbage cleanup
gc()

# Save
dir.create("demo_models", showWarnings = FALSE)
saveRDS(samples1, "demo_models/example_samples1.rds")
saveRDS(samples2, "demo_models/example_samples2.rds")
