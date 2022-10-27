############################################X
#--Elk Density-Dependent Habitat Selection--X
#---------------Brian J. Smith--------------X
#------------Openness Sensitivity-----------X
#===========================================X
#-----------Last update 2022-09-21----------X
############################################X

# Load packages ----
library(nimble)
library(parallel)

# Create model directories ----
# Adding models to .gitignore, so directory will not exist in cloned repos
dir.create("models", showWarnings = FALSE)
dir.create("models/sens", showWarnings = FALSE)
dir.create("models/sens/temp", showWarnings = FALSE)

# Load data ----
# GLM data
dat <- read.csv("data/all_data.csv")

# Determine openness cutoffs ----
hist(dat$open_orig)
# (oc1 <- unname(quantile(dat$open_orig, 0.25)))
# (oc2 <- unname(quantile(dat$open_orig, 0.50)))
# (oc3 <- unname(quantile(dat$open_orig, 0.99)))

(oc1 <- 0.30)
(oc2 <- 0.50)
(oc3 <- 0.9999)

# Source model functions ----
source("99_nimble_NB_model_function.R")

# List of models to fit ----
# Seeds are different to each other AND different from original model fit
chain_list <- list(#list(chain = "_oc1_1",
#                         seed = 112233,
#                         open_cutoff = oc1,
#                         open_direction = "less"),
#                    list(chain = "_oc1_2",
#                         seed = 445566,
#                         open_cutoff = oc1,
#                         open_direction = "less"),
#                    list(chain = "_oc2_1", 
#                         seed = 101112,
#                         open_cutoff = oc2,
#                         open_direction = "less"),
#                    list(chain = "_oc2_2", 
#                         seed = 131415,
#                         open_cutoff = oc2,
#                         open_direction = "less"),
                   list(chain = "_oc3_1", 
                        seed = 161718,
                        open_cutoff = oc3,
                        open_direction = "greater"),
                   list(chain = "_oc3_2", 
                        seed = 192021,
                        open_cutoff = oc3,
                        open_direction = "greater")
                   )

# Start cluster ----

time_df <- data.frame(event = c("Start", "End"),
                      time = Sys.time())
time_df[1, 2] <- Sys.time()

# ... start cluster ----
this_cluster <- makeCluster(length(chain_list))

# ... run chains ----
# Seems like ~ 2.5 h for first 250 iter, 30 min/250 after that 
## 41.83h for 20k
## 245.57h for 100k
niter <- 50000 
nburn <- 0
nthin <- 1
save_rate <- 10000

parLapply(cl = this_cluster, 
          X = chain_list, 
          fun = mod_mcmc, 
          data = dat,
          niter = niter,
          nburn = nburn,
          nthin = nthin,
          save_rate = save_rate,
          dir = "models/sens")

# ... stop cluster ----
stopCluster(this_cluster)


##Timing##
time_df[2, 2] <- Sys.time()


total.time <- round(difftime(time_df$time[2], time_df$time[1], 
                             units = "hours"), 2)
cat(paste("Total time was", total.time, "hours. \n"))

time_df$diff_h <- c(NA, total.time)

#Store MCMC parameters for timing
time_df$param <- c("niter", "nchains")
time_df$val <- c(niter, 6)

write.csv(time_df, "open_timing.csv", row.names = FALSE)

