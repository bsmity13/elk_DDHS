############################################X
#--Elk Density-Dependent Habitat Selection--X
#---------------Brian J. Smith--------------X
#---------------Model Fitting---------------X
#===========================================X
#-----------Last update 2022-01-10----------X
############################################X

# Load packages ----
library(nimble)
library(parallel)

# Create model directories ----
# Adding models to .gitignore, so directory will not exist in cloned repos
dir.create("models", showWarnings = FALSE)
dir.create("models/temp", showWarnings = FALSE)

# Load data ----
# GLM data
dat <- read.csv("data/all_data.csv")

# Source model functions ----
source("99_nimble_NB_model_function.R")

# List of models to fit ----
chain_list <- list(list(chain = 1,
                        seed = 123),
                   list(chain = 2,
                        seed = 456),
                   list(chain = 3, 
                        seed = 789))

# Start cluster ----

time_df <- data.frame(event = c("Start", "End"),
                      time = Sys.time())
time_df[1, 2] <- Sys.time()

# ... start cluster ----
this_cluster <- makeCluster(length(chain_list))

# ... run chains ----
# Seems like ~ 15 min for 25k iterations
niter <- 100000 
nburn <- 0
nthin <- 1
save_rate <- 50000

lapply(chain_list,
       null_mcmc, 
          data = dat,
          niter = niter,
          nburn = nburn,
          nthin = nthin,
          save_rate = save_rate)

# ... stop cluster ----
# stopCluster(this_cluster)


##Timing##
time_df[2, 2] <- Sys.time()


total.time <- round(difftime(time_df$time[2], time_df$time[1], 
                             units = "hours"), 2)
cat(paste("Total time was", total.time, "hours. \n"))

time_df$diff_h <- c(NA, total.time)

#Store MCMC parameters for timing
time_df$param <- c("niter", "nchains")
time_df$val <- c(niter, 3)

write.csv(time_df, "timing.csv", row.names = FALSE)
