############################################X
#--Elk Density-Dependent Habitat Selection--X
#---------------Brian J. Smith--------------X
#--------------Helper Functions-------------X
#===========================================X
#-----------Last update 2021-11-29----------X
############################################X

# Create interaction terms ----
intxn <- function(data) {
  # If terms already exist, don't make them again
  # Using just open2 and rough2 as checks
  if (!is.null(data$open2) & !is.null(data$rough2)) {
    warning("Quadratic terms already exist in `data`. Returning input.")
    return(data)
  }
  # Create quadratic terms
  data$open2 <- data$open^2
  data$rough2 <- data$rough^2
  # Create interaction terms
  data$biomass_logdens <- data$biomass * data$log_dens
  data$open_wolf <- data$open * data$wolf
  data$rough_wolf <- data$rough * data$wolf # this was 'open'
  data$open2_wolf <- data$open2 * data$wolf
  data$rough2_wolf <- data$rough2 * data$wolf
  data$open_cougar <- data$open * data$cougar
  data$rough_cougar <- data$rough * data$cougar # this was 'open'
  data$open2_cougar <- data$open2 * data$cougar
  data$rough2_cougar <- data$rough2 * data$cougar
  data$open_logdens <- data$open * data$log_dens
  data$rough_logdens <- data$rough * data$log_dens
  data$open2_logdens <- data$open2 * data$log_dens
  data$rough2_logdens <- data$rough2 * data$log_dens
  # Return
  return(data)
}

# Scale and center data ----
scale_dat <- function(dat, scale_df) {
  # If terms with "_orig" already exist, the data have already been scaled
  if (any(grepl("_orig", names(dat)))) {
    warning(paste("Columns named '*_orig' already exist.",
                  "Assuming data have already been scaled. Returning input."))
    return(dat)
  }
  
  for (i in 1:nrow(scale_df)) {
    term <- scale_df$term[i]
    mu <- scale_df$mean[i]
    sig <- scale_df$sd[i]
    # Copy unscaled data as "<term>_orig"
    orig_name <- paste0(term, "_orig")
    dat[[orig_name]] <- dat[[term]]
    # Scale and center
    dat[[term]] <- (dat[[term]] - mu)/sig
  }
  return(dat)
}

# Prep data for Nimble ----
prep_nimble <- function(dat, n_beta, rand_seed = rpois(1, 99), 
                        spatial_effect = TRUE,
                        open_cutoff = 0, 
                        open_direction = "less", 
                        null = FALSE) {
  # Set random seed (for generating initial values)
  set.seed(rand_seed)
  
  # If this is for a null model
  if (null) {
    # Data
    elkData <- list()
    elkData$n <- tapply(dat$n, list(dat$cell, dat$year), identity)
    elkData$log_dens <- tapply(dat$log_dens, list(dat$year), unique)
    # Constants
    elkConst <- list()
    elkConst$length_t <- length(unique(dat$year))
    elkConst$Omega <- length(unique(dat$cell))
    # Inits
    elkInits <- list("nb_r" = runif(1, 0, 1))
    
  } else {
    
    
    # Flag any cell where openness is ever below openness cutoff
    # (done for consistency with indexing over time)
    
    # Check open_direction
    stopifnot(open_direction %in% c("greater", "less"))
    
    if (open_direction == "less") {
      flagged <- sort(unique(dat$cell[which(dat$open_orig < open_cutoff)]))
    } else {
      flagged <- sort(unique(dat$cell[which(dat$open_orig > open_cutoff)]))
    }
    
    # Remove any cell that is flagged
    dat <- dat[which(!(dat$cell %in% flagged)),]
    
    if (spatial_effect) {
      # Calculate distance matrix
      dists <- as.matrix(dist(unique(dat[, c("x", "y")])))
      dists <- dists / max(dists)
    }
    
    # List of annual data
    dat_yr <- split(dat, dat$year)
    
    # Convert data into matrices
    elkData <- list()
    elkData$n <- tapply(dat$n, list(dat$cell, dat$year), identity)
    elkData$log_dens <- tapply(dat$log_dens, list(dat$year), unique)
    elkData$outside <- tapply(dat$outside, list(dat$cell, dat$year), identity)
    elkData$swe <- tapply(dat$swe, list(dat$cell, dat$year), identity)
    elkData$elev <- tapply(dat$elev, list(dat$cell), unique)
    elkData$cos_asp <- tapply(dat$cos_asp, list(dat$cell), unique)
    elkData$sin_asp <- tapply(dat$sin_asp, list(dat$cell), unique)
    elkData$biomass <- tapply(dat$biomass, list(dat$cell, dat$year), identity)
    elkData$biomass_logdens <- tapply(dat$biomass_logdens, list(dat$cell, dat$year), identity)
    elkData$open <- tapply(dat$open, list(dat$cell, dat$year), identity)
    elkData$rough <- tapply(dat$rough, list(dat$cell, dat$year), identity)
    elkData$open2 <- tapply(dat$open2, list(dat$cell, dat$year), identity)
    elkData$rough2 <- tapply(dat$rough2, list(dat$cell, dat$year), identity)
    elkData$open_wolf <- tapply(dat$open_wolf, list(dat$cell, dat$year), identity)
    elkData$rough_wolf <- tapply(dat$rough_wolf, list(dat$cell, dat$year), identity)
    elkData$open2_wolf <- tapply(dat$open2_wolf, list(dat$cell, dat$year), identity)
    elkData$rough2_wolf <- tapply(dat$rough2_wolf, list(dat$cell, dat$year), identity)
    elkData$open_cougar <- tapply(dat$open_cougar, list(dat$cell, dat$year), identity)
    elkData$rough_cougar <- tapply(dat$rough_cougar, list(dat$cell, dat$year), identity)
    elkData$open2_cougar <- tapply(dat$open2_cougar, list(dat$cell, dat$year), identity)
    elkData$rough2_cougar <- tapply(dat$rough2_cougar, list(dat$cell, dat$year), identity)
    elkData$open_logdens <- tapply(dat$open_logdens, list(dat$cell, dat$year), identity)
    elkData$rough_logdens <- tapply(dat$rough_logdens, list(dat$cell, dat$year), identity)
    elkData$open2_logdens <- tapply(dat$open2_logdens, list(dat$cell, dat$year), identity)
    elkData$rough2_logdens <- tapply(dat$rough2_logdens, list(dat$cell, dat$year), identity)
    
    elkConst <- list()
    elkConst$length_t <- length(unique(dat$year))
    elkConst$Omega <- length(unique(dat$cell))
    
    if (spatial_effect) {
      elkConst$dists <- dists
      elkConst$zeroes <- rep(0, elkConst$Omega)
    }
    
    elkInits <- list("beta" = runif((n_beta), -1, 1),
                     "eta_out" = rnorm(elkConst$length_t, 0, 1),
                     "mu_out" = 0,
                     "sigma_out" = 1,
                     "nb_r" = runif(1, 0, 1))
    
    if (spatial_effect) {
      elkInits$sigma <- 2
      elkInits$rho <- 0.1  
      elkInits$cov <- expcov(dists, elkInits$rho, elkInits$sigma)
      elkInits$s <-  t(chol(elkInits$cov))
      elkInits$s <- elkInits$s[ , 1]  # give nimble a vector, not one-column matrix
    }
  }
  return(list(data = elkData,
              const = elkConst,
              init = elkInits))
}

# Create prediction data ----
# Function allows you to specify only the terms of interest, with others
# defaulting to the mean. Does not create interactions, which can be
# created with function 'intxn()', and does not scale and center, which
# can be done with 'scale_dat()'.
create_pred_dat <- function(dat,
                            swe = mean(dat$swe_orig),
                            elev = mean(dat$elev_orig),
                            cos_asp = 1,
                            sin_asp = 0,
                            biomass = mean(dat$biomass_orig),
                            open = mean(dat$open_orig),
                            rough = mean(dat$rough_orig),
                            wolf = mean(dat$wolf),
                            cougar = mean(dat$cougar),
                            log_dens = mean(dat$log_dens),
                            outside = 0
) {
  res <- expand.grid(swe = swe,
                     elev = elev,
                     cos_asp = cos_asp,
                     sin_asp = sin_asp,
                     biomass = biomass,
                     open = open,
                     rough = rough,
                     wolf = wolf,
                     cougar = cougar,
                     log_dens = log_dens,
                     outside = outside)
  
  return(res)
  
}

# Calculate model prediction ----
# Vectorized across posterior samples
predict_lambda <- function(dat, years, i, b, e = NA, s, RE = TRUE, mu_eta = NA) {
  
  if (RE) {
    # Index for temporal random effect
    j <- which(years == dat$year[i])
    
    l <- exp(dat$log_dens[i] +
               # Random effect of in/out
               e[[j]] * dat$outside[i] +
               # Conditions
               b[[1]] * dat$swe[i] +
               b[[2]] * dat$elev[i] +
               b[[3]] * dat$cos_asp[i] +
               b[[4]] * dat$sin_asp[i] +
               # Resources
               b[[5]] * dat$biomass[i] +
               b[[6]] * dat$biomass_logdens[i] +
               # Risk main effect
               b[[7]] * dat$open[i] +
               b[[8]] * dat$open2[i] +
               b[[9]] * dat$rough[i] +
               b[[10]] * dat$rough2[i] +
               # Risk wolf effect
               b[[11]] * dat$open_wolf[i] +
               b[[12]] * dat$open2_wolf[i] +
               b[[13]] * dat$rough_wolf[i] +
               b[[14]] * dat$rough2_wolf[i] +
               # Risk cougar effect
               b[[15]] * dat$open_cougar[i] +
               b[[16]] * dat$open2_cougar[i] +
               b[[17]] * dat$rough_cougar[i] +
               b[[18]] * dat$rough2_cougar[i] +
               # Risk density-dependence
               b[[19]] * dat$open_logdens[i] +
               b[[20]] * dat$open2_logdens[i] +
               b[[21]] * dat$rough_logdens[i] +
               b[[22]] * dat$rough2_logdens[i] +
               s[, dat$cell[i]])
  } else {
    
    l <- exp(dat$log_dens[i] +
               # Random effect of in/out
               mu_eta * dat$outside[i] +
               # Conditions
               b[[1]] * dat$swe[i] +
               b[[2]] * dat$elev[i] +
               b[[3]] * dat$cos_asp[i] +
               b[[4]] * dat$sin_asp[i] +
               # Resources
               b[[5]] * dat$biomass[i] +
               b[[6]] * dat$biomass_logdens[i] +
               # Risk main effect
               b[[7]] * dat$open[i] +
               b[[8]] * dat$open2[i] +
               b[[9]] * dat$rough[i] +
               b[[10]] * dat$rough2[i] +
               # Risk wolf effect
               b[[11]] * dat$open_wolf[i] +
               b[[12]] * dat$open2_wolf[i] +
               b[[13]] * dat$rough_wolf[i] +
               b[[14]] * dat$rough2_wolf[i] +
               # Risk cougar effect
               b[[15]] * dat$open_cougar[i] +
               b[[16]] * dat$open2_cougar[i] +
               b[[17]] * dat$rough_cougar[i] +
               b[[18]] * dat$rough2_cougar[i] +
               # Risk density-dependence
               b[[19]] * dat$open_logdens[i] +
               b[[20]] * dat$open2_logdens[i] +
               b[[21]] * dat$rough_logdens[i] +
               b[[22]] * dat$rough2_logdens[i])
  }
  
  return(l)
}

# Calculate model prediction per iteration ----
# Vectorized across data
predict_lambda_iter <- function(dat, i, b, e, s, RE = TRUE) {
  
  if (RE) {
    # Predict lambda
    l <- exp(dat$log_dens +
               # Random effect of in/out
               unname(e[i, dat$t]) * dat$outside +
               # Conditions
               b[[1]][i] * dat$swe +
               b[[2]][i] * dat$elev +
               b[[3]][i] * dat$cos_asp +
               b[[4]][i] * dat$sin_asp +
               # Resources
               b[[5]][i] * dat$biomass +
               b[[6]][i] * dat$biomass_logdens +
               # Risk main effect
               b[[7]][i] * dat$open +
               b[[8]][i] * dat$open2 +
               b[[9]][i] * dat$rough +
               b[[10]][i] * dat$rough2 +
               # Risk wolf effect
               b[[11]][i] * dat$open_wolf +
               b[[12]][i] * dat$open2_wolf +
               b[[13]][i] * dat$rough_wolf +
               b[[14]][i] * dat$rough2_wolf +
               # Risk cougar effect
               b[[15]][i] * dat$open_cougar +
               b[[16]][i] * dat$open2_cougar +
               b[[17]][i] * dat$rough_cougar +
               b[[18]][i] * dat$rough2_cougar +
               # Risk density-dependence
               b[[19]][i] * dat$open_logdens +
               b[[20]][i] * dat$open2_logdens +
               b[[21]][i] * dat$rough_logdens +
               b[[22]][i] * dat$rough2_logdens +
               unname(s[i, dat$cell]))
  } else {
    # Predict lambda
    l <- exp(dat$log_dens +
               # Conditions
               b[[1]][i] * dat$swe +
               b[[2]][i] * dat$elev +
               b[[3]][i] * dat$cos_asp +
               b[[4]][i] * dat$sin_asp +
               # Resources
               b[[5]][i] * dat$biomass +
               b[[6]][i] * dat$biomass_logdens +
               # Risk main effect
               b[[7]][i] * dat$open +
               b[[8]][i] * dat$open2 +
               b[[9]][i] * dat$rough +
               b[[10]][i] * dat$rough2 +
               # Risk wolf effect
               b[[11]][i] * dat$open_wolf +
               b[[12]][i] * dat$open2_wolf +
               b[[13]][i] * dat$rough_wolf +
               b[[14]][i] * dat$rough2_wolf +
               # Risk cougar effect
               b[[15]][i] * dat$open_cougar +
               b[[16]][i] * dat$open2_cougar +
               b[[17]][i] * dat$rough_cougar +
               b[[18]][i] * dat$rough2_cougar +
               # Risk density-dependence
               b[[19]][i] * dat$open_logdens +
               b[[20]][i] * dat$open2_logdens +
               b[[21]][i] * dat$rough_logdens +
               b[[22]][i] * dat$rough2_logdens)
  }

  return(unname(l))
}

# Plot prediction map ----
map_lambda <- function(pred, yr, ynp, nr){
  # Subset data.frame
  r_df <- pred %>% 
    filter(year == yr)
  
  pp <- ggplot() + 
    geom_raster(data = r_df, aes(x = x, y = y, fill = lambda)) +
    geom_sf(data = ynp, color = "green", fill = NA) +
    geom_sf(data = nr, color = "goldenrod", fill = NA, size = 1) +
    ylim(4955000, 5020000) +
    ggtitle(yr) +
    xlab(NULL) +
    ylab(NULL) +
    theme_bw()
  
  return(pp)
}

map_n <- function(pred, yr, ynp, nr){
  # Subset data.frame
  r_df <- pred %>% 
    filter(year == yr)
  
  pp <- ggplot() + 
    geom_raster(data = r_df, aes(x = x, y = y, fill = n)) +
    geom_sf(data = ynp, color = "green", fill = NA) +
    geom_sf(data = nr, color = "goldenrod", fill = NA, size = 1) +
    ylim(4955000, 5020000) +
    ggtitle(yr) +
    xlab(NULL) +
    ylab(NULL) +
    theme_bw()
  
  return(pp)
}

