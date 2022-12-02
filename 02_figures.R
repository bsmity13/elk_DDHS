############################################X
#--Elk Density-Dependent Habitat Selection--X
#---------------Brian J. Smith--------------X
#------------------Figures------------------X
#===========================================X
#-----------Last update 2022-12-02----------X
############################################X

# For Ecology Letters:
#   single column (82 mm)
#   two-thirds page width (110 mm)
#   full page width (173 mm)

#Load packages----
library(tidyverse)
library(sf)
library(raster)
library(extrafont)
library(patchwork)
library(RColorBrewer)
library(coda)
library(ragg)
library(ncf)

# Custom functions ----
source("99_fun.R")

# Colors ----
# ... grayscale ----
# 50% credible interval
color_50 <- "black"
# 80% credible interval
color_80 <- "gray40"
# 90% credible interval
color_90 <- "gray70"

# ... two color ----
## A (orange)
# main
color_a <- "#e66100"
# 50% credible interval
color_a_50 <- "#e66100"
# 80% credible interval
color_a_80 <- "#ff9433"
# 90% credible interval
color_a_90 <- "#ffc766"

## B (blue)
# main
color_b <- "#001576"
# 50% credible interval
color_b_50 <- "#001576"
# 80% credible interval
color_b_80 <- "#0048a9"
# 90% credible interval
color_b_90 <- "#0c7bdc"

#Load data----
# Model data
dat <- read.csv("data/all_data.csv")

#Times
years <- sort(unique(dat$year))

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

# Load other data ----
# Scaling data.frame
scale_df <- read.csv("data/scale_df.csv")

#Good cells
good <- readRDS("data/good_cells.rds")

# Template raster
trast <- readRDS("data/r_template.rds")

#Northern range shapefile
nr <- st_read("geo/NR_Bound_Revised/NR_bound_revised_2021.shp") %>% 
  st_transform(crs = 26912)

#Yellowstone boundary
ynp <- st_read("geo/YNP_boundary/YNP_boundary.shp") %>% 
  st_transform(crs = 26912)

# Estimated population size (from model in appendix of Tallian et al. 2017)
N_SSM <- read.csv("data/NR_elk_SSM_1988-2020.csv")

# Scaling factor (elk/pixel -> elk/km2)
pixel_scale <- 1/((res(trast)[1])^2/(1000^2))

# Create figure directories ----
# Adding figures to .gitignore, so directory will not exist in cloned repos
dir.create("fig", showWarnings = FALSE)
dir.create("fig/map_resid", showWarnings = FALSE)
dir.create("fig/map_scaled", showWarnings = FALSE)
dir.create("fig/map_unscaled", showWarnings = FALSE)
dir.create("fig/tweet", showWarnings = FALSE)

# Coda diagnostics ----
coda1 <- lapply(samples1, as.mcmc)

gelman.diag(coda1)

write.csv(gelman.diag(coda1)$psrf, "out/gelman_rubin.csv")

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
sigma_eta <- unname(do.call(c, lapply(samples1, function(x) {
  return(x[, "sigma_out"])
})))

s <- lapply(samples2, function(x) {
  return(as.data.frame(x[1:nrow(x), ]))
})  %>% 
  bind_rows() %>% 
  as.data.frame()
names(s) <- paste0("s", 1:ncol(s))

rho <- unname(do.call(c, lapply(samples1, function(x) {
  return(x[, "rho"])
})))
sig <- unname(do.call(c, lapply(samples1, function(x) {
  return(x[, "sigma"])
})))

# Food/Safety correlation ----
# Correlation between food and product open x rough
# Use Spearman to check ranks because of product

cor(dat$biomass, dat$open * dat$rough, method = "spearman")

# Result is similar with Pearson's correlation
cor(dat$biomass, dat$open * dat$rough, method = "pearson")

# RE summary stats ----
# ... temporal ----
# eta mean
mean(mu_eta)
quantile(mu_eta, c(0.05, 0.95))

# eta variance
mean(sigma_eta^2)
quantile(sigma_eta^2, c(0.05, 0.95))

# eta_t
apply(e, 2, mean)
range(apply(e, 2, mean))

# ... spatial ----
# Decay
mean(rho)
quantile(rho, c(0.05, 0.95))

# SD
mean(sig)
quantile(sig, c(0.05, 0.95))

# Coefficient plots ----

# ... traceplots ----

# ... ... betas ----
{
  pdf(file = "fig/trace_beta.pdf", onefile = TRUE, 
      width = 12, height = 6)
  for (i in 1:max(grep("beta", colnames(samples1[[1]])))) { # i indexes variable
    
    # Figure out y-limits for this column
    ylim_list <- lapply(samples1, function(x, i){
      min = min(x[, i])
      max = max(x[, i])
      return(list(min = min, max = max))
    }, i = i)
    ylim <- range(unlist(ylim_list))
    
    for (j in 1:length(samples1)) { # j indexes chains
      if (j == 1) {
        plot(samples1[[j]][, i], type = "l", col = j, lty = 2,
             xlab = "Iteration",
             ylab = expression(beta),
             ylim = ylim,
             main = paste("beta", i))
      } else {
        lines(samples1[[j]][, i], col = j, lty = 2)
      }
    }
  }
  dev.off()
}

# Write all the other parameters in samples1 to a single PDF
{
  pdf(file = "fig/trace_other.pdf", onefile = TRUE,
      width = 12, height = 6)
  for (i in which(!grepl("beta", colnames(samples1[[1]])))) { # i indexes variable
    
    # Figure out y-limits for this column
    ylim_list <- lapply(samples1, function(x, i){
      min = min(x[, i])
      max = max(x[, i])
      return(list(min = min, max = max))
    }, i = i)
    ylim <- range(unlist(ylim_list))
    
    for (j in 1:length(samples1)) { # j indexes chains
      if (j == 1) {
        plot(samples1[[j]][, i], type = "l", col = j, lty = 2,
             xlab = "Iteration",
             ylab = expression(beta),
             ylim = ylim,
             main = colnames(samples1[[1]])[i])
      } else {
        lines(samples1[[j]][, i], col = j, lty = 2)
      }
    }
  }
  dev.off()
}

# ... ... s-es ----
# Write every 50th s traceplot to a single PDF
{
  pdf(file = "fig/trace_s.pdf", onefile = TRUE,
      width = 12, height = 6)
  for (i in 1:ncol(samples2[[1]])) { # i indexes variable
    # Figure out y-limits for this column
    ylim_list <- lapply(samples2, function(x, i){
      min = min(x[, i])
      max = max(x[, i])
      return(list(min = min, max = max))
    }, i = i)
    ylim <- range(unlist(ylim_list))
    for (j in 1:length(samples2)) { # j indexes chains
      if (i %% 50 == 0) { # Only plot every 50th 's'
        if (j == 1) {
          plot(samples2[[j]][, i], type = "l", col = j, lty = 2,
               ylim = ylim,
               xlab = "Iteration",
               ylab = "Cell Value (s)",
               main = paste("s", i))
        } else {
          lines(samples2[[j]][, i], col = j, lty = 2)
        }
      }
    }
  }
  dev.off()
}

# ... count coefficients ----
beta_data <- function(b) {
  coefs <- b %>% 
    pivot_longer(everything()) %>% 
    group_by(name) %>% 
    summarize(l90 = quantile(value, 0.05),
              l80 = quantile(value, 0.10),
              l50 = quantile(value, 0.25),
              mean = mean(value),
              u50 = quantile(value, 0.75),
              u80 = quantile(value, 0.90),
              u90 = quantile(value, 0.95),
    )
  
  # Plotting order
  # Order in model code:
  #   1 "SWE"
  #   2 "Elevation"
  #   3 "cos(Aspect)"
  #   4 "sin(Aspect)"
  #   5 "Biomass"
  #   6 "Biomass:log(Dens)"
  #   7 "Open"
  #   8 "Rough"
  #   9 "Open^2"
  #   10 "Rough^2"
  #   11 "Open:Wolf"
  #   12 "Rough:Wolf"
  #   13 "Open^2:Wolf"
  #   14 "Rough^2:Wolf"
  #   15 "Open:Cougar"
  #   16 "Rough:Cougar"
  #   17 "Open^2:Cougar"
  #   18 "Rough^2:Cougar"
  #   19 "Open:log(Dens)"
  #   20 "Rough:log(Dens)"
  #   21 "Open^2:log(Dens)"
  #   22 "Rough^2:log(Dens)"
  coef_order <- data.frame(name = paste0("b", 1:22), 
                           order = factor(c(1:7, 9, 8, 10, 11, 13, 
                                            12, 14, 15, 17, 16, 18,
                                            19, 21, 20, 22)),
                           label = c("SWE",
                                     "Elevation",
                                     "cos(Aspect)",
                                     "sin(Aspect)",
                                     "Biomass",
                                     "Biomass:log(Dens)",
                                     "Open",
                                     "Rough",
                                     "Open^2",
                                     "Rough^2",
                                     "Open:Wolf",
                                     "Rough:Wolf",
                                     "Open^2:Wolf",
                                     "Rough^2:Wolf",
                                     "Open:Cougar",
                                     "Rough:Cougar",
                                     "Open^2:Cougar",
                                     "Rough^2:Cougar",
                                     "Open:log(Dens)",
                                     "Rough:log(Dens)",
                                     "Open^2:log(Dens)",
                                     "Rough^2:log(Dens)"))
  coef_plot <- left_join(coefs, coef_order, by = "name")  
  return(coef_plot)
}

# Plot
bd <- beta_data(b)
cp <- bd %>% 
  ggplot(aes(y = order, x = mean)) +
  geom_vline(xintercept = 0, color = "firebrick", linetype = "dashed", size = 0.5) +
  geom_errorbar(aes(xmin = l90, xmax = u90), color = "gray70", width = 0.2, size = 0.5) +
  geom_errorbar(aes(xmin = l80, xmax = u80), color = "gray40", width = 0, size = 0.75) +
  geom_errorbar(aes(xmin = l50, xmax = u50), color = color_50, width = 0, size = 1) +
  geom_point() +
  scale_color_ordinal(name = "Openness Cutoff") +
  scale_y_discrete(breaks = bd$order, labels = parse(text = bd$label)) +
  ylab(NULL) +
  xlab(expression(beta)) +
  # coord_cartesian(ylim = c(-1.5, 1.75)) +
  theme_bw() +
  # theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  NULL

ggsave("fig/betas.tiff", plot = cp, width = 4, height = 6, device = agg_tiff,
       units = "in", dpi = 200, compression = "lzw")

# Spatial effects ----
mean_s <- apply(s, 2, mean)
s_rast <- trast
s_vals <- rep(NA, ncell(trast))
s_vals[good] <- mean_s
values(s_rast) <- s_vals
# plot(s_rast)
s_dat <- as.data.frame(s_rast, xy = TRUE)

(s_fig <- ggplot() +
    geom_raster(data = s_dat, aes(x = x, y = y, fill = layer)) +
    scale_fill_gradient2(name = expression(s[i]),
                         high = "firebrick",
                         low = "navy", 
                         na.value = "gray80") +
    xlab("Easting") +
    ylab("Northing") +
    theme_bw())
ggsave("fig/spatial_effects.tiff", plot = s_fig, width = 10, height = 8, 
       units = "in", dpi = 300, compression = "lzw", device = agg_tiff)

# Mean effect inside vs outside

# Outside YNP
out_dat <- dat %>% 
  dplyr::select(x, y, outside) %>% 
  distinct()
s_dat_out <- left_join(s_dat, out_dat, by = c("x", "y"))

s_dat_out %>% 
  filter(!is.na(outside)) %>% 
  group_by(outside) %>% 
  summarize(mean = mean(layer),
            q05 = quantile(layer, 0.05),
            q95 = quantile(layer, 0.95))


# Predict expected density ----
pred <- dat
pred$lambda <- NA
pred$lwr <- NA
pred$upr <- NA
pred$var <- NA
lambda_samps <- list()

for (i in 1:nrow(pred)) {
  cat("\r", i, "       ")
  
  lambda_samps[[i]] <- predict_lambda(pred, years, i, b, e, s)
  
  pred$lambda[i] <- mean(lambda_samps[[i]])
  pred$lwr[i] <- quantile(lambda_samps[[i]], 0.05)
  pred$upr[i] <- quantile(lambda_samps[[i]], 0.95)
  pred$var[i] <- var(lambda_samps[[i]])
  
}

# Save predictions
# dir.create("pred", showWarnings = FALSE)
# saveRDS(lambda_samps, "pred/lambda_samps.rds")
# saveRDS(pred, "pred/model_data_predictions.rds")

# Expected population inside vs outside ----
io <- pred %>% 
  group_by(year, outside) %>% 
  summarize(E_N = sum(lambda)) %>% 
  ungroup() %>% 
  group_by(year) %>% 
  mutate(year_total = sum(E_N),
         prop = E_N/year_total)

io %>% 
  ggplot(aes(x = year, y = prop, color = as.logical(outside))) +
  geom_line() +
  geom_point() +
  scale_color_discrete(name = "Where",
                       breaks = c(FALSE, TRUE),
                       labels = c("Inside YNP", "Outside YNP")) +
  xlab("Winter") +
  ylab("Proportion of Population") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

# Plot density maps ----
hs_plots <- list()

for (i in 1:length(years)){
  
  #Plot
  hs_plots[[i]] <- map_lambda(pred, years[i], ynp, nr)
  
  #Save on same scale
  ggsave(paste0("fig/map_scaled/N_it_raster_", 
                years[i], ".tiff"), 
         plot = hs_plots[[i]] +
           scale_fill_viridis_c(option = "A", 
                                expression("Expected\nDensity" * (elk/km^2)),
                                breaks = pretty(0:550),
                                limits = c(0, 551),
                                alpha = 0.75), 
         width = 10, height = 6, units = "in", compression = "lzw", 
         device = agg_tiff)
  
  #Save on different scale
  ggsave(paste0("fig/map_unscaled/N_it_raster_", 
                years[i], ".tiff"), 
         plot = hs_plots[[i]] +
           scale_fill_viridis_c(option = "D", 
                                name = expression("Expected\nDensity" * (elk/km^2)),
                                alpha = 0.75), 
         width = 10, height = 6, units = "in", compression = "lzw", 
         device = agg_tiff)
  
}

# Population size ----
mod_pop <- pred %>% 
  group_by(year) %>% 
  summarize(n = sum(lambda))
count <- dat %>% 
  group_by(year) %>% 
  summarize(N = sum(n))
# Plot
{ 
  agg_tiff("fig/population.tiff", width = 12, height = 6, units = "in",
           res = 150, compression = "lzw")
  par(mar = c(2, 3, 0.5, 0.5))
  plot(years, mod_pop$n, pch = 16, ylim = c(0, 75000),
       xlab = NA, ylab = "Total Elk")
  lines(years, mod_pop$n)
  points(N_SSM$year, N_SSM$mean, pch = 16, col = 2)
  lines(N_SSM$year, N_SSM$mean, col = 2)
  points(count$year, count$N, pch = 16, col = 3)
  lines(count$year, count$N, col = 3)
  legend("topright", pch = 16, lty = 1, col = 1:3,
         legend = c("Model Pred",
                    "SSM Estimate",
                    "Raw Counts"))
  dev.off()
  }


# Coefficient interpretation ----

# ... biomass conversion ----
# 1 lb = 0.45359 kg
# 1 acre = 0.40469 ha
# ==> 1 lb/acre = 0.45359/0.40469 kg/ha
lbac_kgha <- 0.45359/0.40469

# ... standard devations ----
# (for RSS comparison)
# Biomass
bio_sd_kgha <- dat %>% 
  mutate(biomass_kgha = exp(biomass_orig) * lbac_kgha) %>% 
  pull(biomass_kgha) %>% 
  sd()

bio_sd <- scale_df %>% 
  filter(term == "biomass") %>% 
  pull(sd)

# Openness
open_sd <- scale_df %>% 
  filter(term == "open") %>% 
  pull(sd)

# Roughness
rough_sd <- scale_df %>% 
  filter(term == "rough") %>% 
  pull(sd)

# ... mean density ----
# ... ... biomass ----
# Setup data
bio_dens <- create_pred_dat(dat, biomass = log(seq(1, 1000, length.out = 50)),
                            log_dens = mean(dat$log_dens)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(biomass_nat = exp(biomass_orig) * lbac_kgha) %>% 
  mutate(dens = round(exp(log_dens), 1))

# Predict lambda
bio_dens$q05 <- bio_dens$q10 <- bio_dens$q25 <- 
  bio_dens$lambda <- 
  bio_dens$q75 <- bio_dens$q90 <- bio_dens$q95 <- NA

for (i in 1:nrow(bio_dens)) {
  cat("               \r", i)
  
  lambda <- predict_lambda(bio_dens, years, i, b, e = NA, s, 
                           RE = FALSE, mu_eta = mu_eta)
  
  # Summarize posterior
  bio_dens$q05[i] <- quantile(lambda, 0.05)
  bio_dens$q10[i] <- quantile(lambda, 0.10)
  bio_dens$q25[i] <- quantile(lambda, 0.25)
  bio_dens$lambda[i] <- mean(lambda)
  bio_dens$q75[i] <- quantile(lambda, 0.75)
  bio_dens$q90[i] <- quantile(lambda, 0.90)
  bio_dens$q95[i] <- quantile(lambda, 0.95)
}

# Plot
(bio_dens_plot <- ggplot(bio_dens, aes(x = biomass_nat, y = lambda)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
                alpha = 0.5) +
    geom_line() +
    scale_x_continuous(breaks = c(0, 250, 500, 750, 1000)) +
    xlab("Biomass (kg/ha)") +
    ylab(expression("Expected Density" ~ (elk/km^2))) +
    theme_bw())
ggsave("fig/density_biomass.tiff", plot = bio_dens_plot,
       width = 5, height = 3, units = "in", compression = "lzw",
       device = agg_tiff)

# ... ... openness ----
# Setup data
open_dens <- create_pred_dat(dat, open = seq(0.5, 1, length.out = 50),
                             log_dens = quantile(dat$log_dens, c(0.5))) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1))

# Predict lambda
open_dens$q05 <- open_dens$q10 <- open_dens$q25 <- 
  open_dens$lambda <- 
  open_dens$q75 <- open_dens$q90 <- open_dens$q95 <- NA

for (i in 1:nrow(open_dens)) {
  cat("               \r", i)
  
  lambda <- predict_lambda(open_dens, years, i, b, e = NA, s, 
                           RE = FALSE, mu_eta = mu_eta)
  
  # Summarize posterior
  open_dens$q05[i] <- quantile(lambda, 0.05)
  open_dens$q10[i] <- quantile(lambda, 0.10)
  open_dens$q25[i] <- quantile(lambda, 0.25)
  open_dens$lambda[i] <- mean(lambda)
  open_dens$q75[i] <- quantile(lambda, 0.75)
  open_dens$q90[i] <- quantile(lambda, 0.90)
  open_dens$q95[i] <- quantile(lambda, 0.95)
}

# Plot
(open_dens_plot <- ggplot(open_dens, aes(x = open_orig, y = lambda)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
                alpha = 0.5) +
    geom_line() +
    xlab("Openness (%)") +
    ylab(expression("Expected Density" ~ (elk/km^2))) +
    theme_bw())

ggsave("fig/density_open.tiff", plot = open_dens_plot,
       width = 5, height = 3, units = "in", compression = "lzw",
       device = agg_tiff)

# ... ... roughness ----
# Setup data
rough_dens <- create_pred_dat(dat, rough = seq(0, 50, length.out = 50),
                              log_dens = mean(dat$log_dens)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1))

# Predict lambda
rough_dens$q05 <- rough_dens$q10 <- rough_dens$q25 <- 
  rough_dens$lambda <- 
  rough_dens$q75 <- rough_dens$q90 <- rough_dens$q95 <- NA

for (i in 1:nrow(rough_dens)) {
  cat("               \r", i)
  
  lambda <- predict_lambda(rough_dens, years, i, b, e = NA, s, 
                           RE = FALSE, mu_eta = mu_eta)
  
  # Summarize posterior
  rough_dens$q05[i] <- quantile(lambda, 0.05)
  rough_dens$q10[i] <- quantile(lambda, 0.10)
  rough_dens$q25[i] <- quantile(lambda, 0.25)
  rough_dens$lambda[i] <- mean(lambda)
  rough_dens$q75[i] <- quantile(lambda, 0.75)
  rough_dens$q90[i] <- quantile(lambda, 0.90)
  rough_dens$q95[i] <- quantile(lambda, 0.95)
}

# Plot
(rough_dens_plot <- ggplot(rough_dens, aes(x = rough_orig, y = lambda)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
                alpha = 0.5) +
    geom_line() +
    xlab("Roughness (m)") +
    ylab(expression("Expected Density" ~ (elk/km^2))) +
    theme_bw())

ggsave("fig/density_rough.tiff", plot = rough_dens_plot,
       width = 5, height = 3, units = "in", compression = "lzw",
       device = agg_tiff)

# ... ... SWE ----
# Setup data
swe_dens <- create_pred_dat(dat, swe = seq(0, 400, length.out = 50),
                            log_dens = mean(dat$log_dens)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1))

# Predict lambda
swe_dens$q05 <- swe_dens$q10 <- swe_dens$q25 <- 
  swe_dens$lambda <- 
  swe_dens$q75 <- swe_dens$q90 <- swe_dens$q95 <- NA

for (i in 1:nrow(swe_dens)) {
  cat("               \r", i)
  
  lambda <- predict_lambda(swe_dens, years, i, b, e = NA, s, 
                           RE = FALSE, mu_eta = mu_eta)
  
  # Summarize posterior
  swe_dens$q05[i] <- quantile(lambda, 0.05)
  swe_dens$q10[i] <- quantile(lambda, 0.10)
  swe_dens$q25[i] <- quantile(lambda, 0.25)
  swe_dens$lambda[i] <- mean(lambda)
  swe_dens$q75[i] <- quantile(lambda, 0.75)
  swe_dens$q90[i] <- quantile(lambda, 0.90)
  swe_dens$q95[i] <- quantile(lambda, 0.95)
}

# Plot
(swe_dens_plot <- ggplot(swe_dens, aes(x = swe_orig, y = lambda)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
                alpha = 0.5) +
    geom_line() +
    xlab(expression("SWE" ~ (kg/m^2))) +
    ylab(expression("Expected Density" ~ (elk/km^2))) +
    theme_bw())

# ... ... elevation ----
# Setup data
elev_dens <- create_pred_dat(dat, elev = seq(1700, 2800, length.out = 50),
                             log_dens = mean(dat$log_dens)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1))

# Predict lambda
elev_dens$q05 <- elev_dens$q10 <- elev_dens$q25 <- 
  elev_dens$lambda <- 
  elev_dens$q75 <- elev_dens$q90 <- elev_dens$q95 <- NA

for (i in 1:nrow(elev_dens)) {
  cat("               \r", i)
  
  lambda <- predict_lambda(elev_dens, years, i, b, e = NA, s, 
                           RE = FALSE, mu_eta = mu_eta)
  
  # Summarize posterior
  elev_dens$q05[i] <- quantile(lambda, 0.05)
  elev_dens$q10[i] <- quantile(lambda, 0.10)
  elev_dens$q25[i] <- quantile(lambda, 0.25)
  elev_dens$lambda[i] <- mean(lambda)
  elev_dens$q75[i] <- quantile(lambda, 0.75)
  elev_dens$q90[i] <- quantile(lambda, 0.90)
  elev_dens$q95[i] <- quantile(lambda, 0.95)
}

# Plot
(elev_dens_plot <- ggplot(elev_dens, aes(x = elev_orig, y = lambda)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
                alpha = 0.5) +
    geom_line() +
    xlab("Elevation (m)") +
    ylab(expression("Expected Density" ~ (elk/km^2))) +
    theme_bw())

# ... ... aspect ----
# Setup data
asp_dens <- create_pred_dat(dat, 
                            cos_asp = cos(seq(0, 2*pi, length.out = 50)),
                            log_dens = mean(dat$log_dens)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1),
         asp = seq(0, 2*pi, length.out = 50),
         sin_asp = sin(asp))

# Predict lambda
asp_dens$q05 <- asp_dens$q10 <- asp_dens$q25 <- 
  asp_dens$lambda <- 
  asp_dens$q75 <- asp_dens$q90 <- asp_dens$q95 <- NA

for (i in 1:nrow(asp_dens)) {
  cat("               \r", i)
  
  lambda <- predict_lambda(asp_dens, years, i, b, e = NA, s, 
                           RE = FALSE, mu_eta = mu_eta)
  
  # Summarize posterior
  asp_dens$q05[i] <- quantile(lambda, 0.05)
  asp_dens$q10[i] <- quantile(lambda, 0.10)
  asp_dens$q25[i] <- quantile(lambda, 0.25)
  asp_dens$lambda[i] <- mean(lambda)
  asp_dens$q75[i] <- quantile(lambda, 0.75)
  asp_dens$q90[i] <- quantile(lambda, 0.90)
  asp_dens$q95[i] <- quantile(lambda, 0.95)
}

# Plot
(asp_dens_plot <- ggplot(asp_dens, aes(x = asp, y = lambda)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
                alpha = 0.5) +
    geom_line() +
    scale_x_continuous(name = "Aspect",
                       breaks = c(0, pi/2, pi, 3*pi/2, 2*pi),
                       labels = c("N", "E", "S", "W", "N")) +
    ylab(expression("Expected Density" ~ (elk/km^2))) +
    theme_bw())

# ... DDHS ----
# ... ... biomass ----
# Biomass x1 value
biomass_x1 <- mean(dat$biomass_orig) + 0.5 * bio_sd
biomass_x1_kgha <- exp(biomass_x1) * lbac_kgha

# Biomass x2 value
biomass_x2 <- mean(dat$biomass_orig) - 0.5* bio_sd
biomass_x2_kgha <- exp(biomass_x2) * lbac_kgha

bio_ddhs_x1 <- create_pred_dat(dat, biomass = biomass_x1,
                               log_dens = seq(quantile(dat$log_dens, 0.05),
                                              quantile(dat$log_dens, 0.95),
                                              length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(biomass_nat = exp(biomass_orig) * lbac_kgha,
         dens = round(exp(log_dens), 1))

bio_ddhs_x2 <- create_pred_dat(dat, biomass = biomass_x2,
                               log_dens = seq(quantile(dat$log_dens, 0.05),
                                              quantile(dat$log_dens, 0.95),
                                              length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(biomass_nat = exp(biomass_orig) * lbac_kgha,
         dens = round(exp(log_dens), 1))

# Predict lambda
bio_ddhs_x1$q05 <- bio_ddhs_x1$q10 <- bio_ddhs_x1$q25 <- 
  bio_ddhs_x1$rss <- 
  bio_ddhs_x1$q75 <- bio_ddhs_x1$q90 <- bio_ddhs_x1$q95 <- NA

for (i in 1:nrow(bio_ddhs_x1)) {
  cat("               \r", i)
  
  lambda_x1 <- predict_lambda(bio_ddhs_x1, years, i, b, e = NA, s, 
                              RE = FALSE, mu_eta = mu_eta)
  lambda_x2 <- predict_lambda(bio_ddhs_x2, years, i, b, e = NA, s, 
                              RE = FALSE, mu_eta = mu_eta)
  rss <- lambda_x1/lambda_x2
  
  # Summarize posterior
  bio_ddhs_x1$q05[i] <- quantile(rss, 0.05)
  bio_ddhs_x1$q10[i] <- quantile(rss, 0.10)
  bio_ddhs_x1$q25[i] <- quantile(rss, 0.25)
  bio_ddhs_x1$rss[i] <- mean(rss)
  bio_ddhs_x1$q75[i] <- quantile(rss, 0.75)
  bio_ddhs_x1$q90[i] <- quantile(rss, 0.90)
  bio_ddhs_x1$q95[i] <- quantile(rss, 0.95)
}

# Plot
(bio_ddhs_plot <- ggplot(bio_ddhs_x1, aes(x = dens, y = rss)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
                alpha = 0.5) +
    geom_line() +
    xlab(expression("Average Elk Density" ~ (elk/km^2))) +
    ylab("RSS for Biomass") +
    theme_bw())

ggsave("fig/ddhs_biomass.tiff", plot = bio_ddhs_plot,
       width = 6, height = 3, units = "in", compression = "lzw",
       device = agg_tiff)

# ... ... openness ----
# Openness x1 value
open_x1 <- 1

# Openness x2 value
open_x2 <- 1 - open_sd

open_ddhs_x1 <- create_pred_dat(dat, open = open_x1,
                                log_dens = seq(quantile(dat$log_dens, 0.05),
                                               quantile(dat$log_dens, 0.95),
                                               length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1))

open_ddhs_x2 <- create_pred_dat(dat, open = open_x2,
                                log_dens = seq(quantile(dat$log_dens, 0.05),
                                               quantile(dat$log_dens, 0.95),
                                               length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1))

# Predict lambda
open_ddhs_x1$q05 <- open_ddhs_x1$q10 <- open_ddhs_x1$q25 <- 
  open_ddhs_x1$rss <- 
  open_ddhs_x1$q75 <- open_ddhs_x1$q90 <- open_ddhs_x1$q95 <- NA

for (i in 1:nrow(open_ddhs_x1)) {
  cat("               \r", i)
  
  lambda_x1 <- predict_lambda(open_ddhs_x1, years, i, b, e = NA, s, 
                              RE = FALSE, mu_eta = mu_eta)
  lambda_x2 <- predict_lambda(open_ddhs_x2, years, i, b, e = NA, s, 
                              RE = FALSE, mu_eta = mu_eta)
  rss <- lambda_x1/lambda_x2
  
  # Summarize posterior
  open_ddhs_x1$q05[i] <- quantile(rss, 0.05)
  open_ddhs_x1$q10[i] <- quantile(rss, 0.10)
  open_ddhs_x1$q25[i] <- quantile(rss, 0.25)
  open_ddhs_x1$rss[i] <- mean(rss)
  open_ddhs_x1$q75[i] <- quantile(rss, 0.75)
  open_ddhs_x1$q90[i] <- quantile(rss, 0.90)
  open_ddhs_x1$q95[i] <- quantile(rss, 0.95)
}

# Plot
(open_ddhs_plot <- ggplot(open_ddhs_x1, aes(x = dens, y = rss)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
                alpha = 0.5) +
    geom_line() +
    xlab(expression("Average Elk Density" ~ (elk/km^2))) +
    ylab("RSS for Openness") +
    theme_bw())

ggsave("fig/ddhs_open.tiff", plot = open_ddhs_plot,
       width = 6, height = 3, units = "in", compression = "lzw",
       device = agg_tiff)

# ... ... roughness ----
# Roughness x1 value
rough_x1 <- 23

# Roughness x2 value
rough_x2 <- rough_x1 - rough_sd

rough_ddhs_x1 <- create_pred_dat(dat, rough = rough_x1,
                                 log_dens = seq(quantile(dat$log_dens, 0.05),
                                                quantile(dat$log_dens, 0.95),
                                                length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1))

rough_ddhs_x2 <- create_pred_dat(dat, rough = rough_x2,
                                 log_dens = seq(quantile(dat$log_dens, 0.05),
                                                quantile(dat$log_dens, 0.95),
                                                length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1))

# Predict lambda
rough_ddhs_x1$q05 <- rough_ddhs_x1$q10 <- rough_ddhs_x1$q25 <- 
  rough_ddhs_x1$rss <- 
  rough_ddhs_x1$q75 <- rough_ddhs_x1$q90 <- rough_ddhs_x1$q95 <- NA

for (i in 1:nrow(rough_ddhs_x1)) {
  cat("               \r", i)
  
  lambda_x1 <- predict_lambda(rough_ddhs_x1, years, i, b, e = NA, s, 
                              RE = FALSE, mu_eta = mu_eta)
  lambda_x2 <- predict_lambda(rough_ddhs_x2, years, i, b, e = NA, s, 
                              RE = FALSE, mu_eta = mu_eta)
  rss <- lambda_x1/lambda_x2
  
  # Summarize posterior
  rough_ddhs_x1$q05[i] <- quantile(rss, 0.05)
  rough_ddhs_x1$q10[i] <- quantile(rss, 0.10)
  rough_ddhs_x1$q25[i] <- quantile(rss, 0.25)
  rough_ddhs_x1$rss[i] <- mean(rss)
  rough_ddhs_x1$q75[i] <- quantile(rss, 0.75)
  rough_ddhs_x1$q90[i] <- quantile(rss, 0.90)
  rough_ddhs_x1$q95[i] <- quantile(rss, 0.95)
}

# Plot
(rough_ddhs_plot <- ggplot(rough_ddhs_x1, aes(x = dens, y = rss)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
                alpha = 0.5) +
    geom_line() +
    xlab(expression("Average Elk Density" ~ (elk/km^2))) +
    ylab("RSS for Roughness") +
    theme_bw())

ggsave("fig/ddhs_rough.tiff", plot = rough_ddhs_plot,
       width = 6, height = 3, units = "in", compression = "lzw",
       device = agg_tiff)

# ... ... safety ----
# 1 SD of openness and roughness (0.5 SD each)

# Openness x1 value
open_safe_x1 <- 1

# Openness x2 value
open_safe_x2 <- open_safe_x1 - (0.5 * open_sd)

# Roughness x1 value
rough_safe_x1 <- 23

# Roughness x2 value
rough_safe_x2 <- rough_safe_x1 - (0.5 * rough_sd)

safe_ddhs_x1 <- create_pred_dat(dat, 
                                open = open_safe_x1,
                                rough = rough_safe_x1,
                                log_dens = seq(quantile(dat$log_dens, 0.05),
                                               quantile(dat$log_dens, 0.95),
                                               length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1))

safe_ddhs_x2 <- create_pred_dat(dat, 
                                open = open_safe_x2,
                                rough = rough_safe_x2,
                                log_dens = seq(quantile(dat$log_dens, 0.05),
                                               quantile(dat$log_dens, 0.95),
                                               length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1))

# Predict lambda
safe_ddhs_x1$q05 <- safe_ddhs_x1$q10 <- safe_ddhs_x1$q25 <- 
  safe_ddhs_x1$rss <- 
  safe_ddhs_x1$q75 <- safe_ddhs_x1$q90 <- safe_ddhs_x1$q95 <- NA

for (i in 1:nrow(safe_ddhs_x1)) {
  cat("               \r", i)
  
  lambda_x1 <- predict_lambda(safe_ddhs_x1, years, i, b, e = NA, s, 
                              RE = FALSE, mu_eta = mu_eta)
  lambda_x2 <- predict_lambda(safe_ddhs_x2, years, i, b, e = NA, s, 
                              RE = FALSE, mu_eta = mu_eta)
  rss <- lambda_x1/lambda_x2
  
  # Summarize posterior
  safe_ddhs_x1$q05[i] <- quantile(rss, 0.05)
  safe_ddhs_x1$q10[i] <- quantile(rss, 0.10)
  safe_ddhs_x1$q25[i] <- quantile(rss, 0.25)
  safe_ddhs_x1$rss[i] <- mean(rss)
  safe_ddhs_x1$q75[i] <- quantile(rss, 0.75)
  safe_ddhs_x1$q90[i] <- quantile(rss, 0.90)
  safe_ddhs_x1$q95[i] <- quantile(rss, 0.95)
}

# Plot
(safe_ddhs_plot <- ggplot(safe_ddhs_x1, aes(x = dens, y = rss)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
                alpha = 0.5) +
    geom_line() +
    xlab(expression("Average Elk Density" ~ (elk/km^2))) +
    ylab("RSS for Safety") +
    theme_bw())

ggsave("fig/ddhs_safety.tiff", plot = safe_ddhs_plot,
       width = 6, height = 3, units = "in", compression = "lzw",
       device = agg_tiff)

# ... predator density ----
# Area of *original* NR polygon clipped to inside YNP only
nr_ynp_m2 <- 994881034 #m2
nr_ynp_km2 <- nr_ynp_m2/1e6

# Wolf abundance to wolf density (wolves/100 km2)
wa_wd <- 1/(nr_ynp_km2/100)

# ... openness ----
# ... ... wolves ----
open_wolf_x1 <- create_pred_dat(dat, open = open_x1,
                                log_dens = quantile(dat$log_dens, 0.5),
                                wolf = seq(quantile(dat$wolf, 0.05),
                                           quantile(dat$wolf, 0.95),
                                           length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1),
         wolf_dens = wolf * wa_wd)

open_wolf_x2 <- create_pred_dat(dat, open = open_x2,
                                log_dens = quantile(dat$log_dens, 0.5),
                                wolf = seq(quantile(dat$wolf, 0.05),
                                           quantile(dat$wolf, 0.95),
                                           length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1),
         wolf_dens = wolf * wa_wd)

# Predict lambda
open_wolf_x1$q05 <- open_wolf_x1$q10 <- open_wolf_x1$q25 <- 
  open_wolf_x1$rss <- 
  open_wolf_x1$q75 <- open_wolf_x1$q90 <- open_wolf_x1$q95 <- NA

for (i in 1:nrow(open_wolf_x1)) {
  cat("               \r", i)
  
  lambda_x1 <- predict_lambda(open_wolf_x1, years, i, b, e = NA, s, 
                              RE = FALSE, mu_eta = mu_eta)
  lambda_x2 <- predict_lambda(open_wolf_x2, years, i, b, e = NA, s, 
                              RE = FALSE, mu_eta = mu_eta)
  rss <- lambda_x1/lambda_x2
  
  # Summarize posterior
  open_wolf_x1$q05[i] <- quantile(rss, 0.05)
  open_wolf_x1$q10[i] <- quantile(rss, 0.10)
  open_wolf_x1$q25[i] <- quantile(rss, 0.25)
  open_wolf_x1$rss[i] <- mean(rss)
  open_wolf_x1$q75[i] <- quantile(rss, 0.75)
  open_wolf_x1$q90[i] <- quantile(rss, 0.90)
  open_wolf_x1$q95[i] <- quantile(rss, 0.95)
}

# Plot
(open_wolf_plot <- ggplot(open_wolf_x1, aes(x = wolf_dens, y = rss)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
                alpha = 0.5) +
    geom_line() +
    coord_cartesian(ylim = c(1.25, 3.25)) +
    xlab(expression("Wolf Density" ~ (wolves/100~km^2))) +
    ylab("RSS for Openness") +
    theme_bw())

# ... ... cougars ----
open_cougar_x1 <- create_pred_dat(dat, open = open_x1,
                                  log_dens = quantile(dat$log_dens, 0.5),
                                  cougar = seq(quantile(dat$cougar, 0.05),
                                               quantile(dat$cougar, 0.95),
                                               length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1))

open_cougar_x2 <- create_pred_dat(dat, open = open_x2,
                                  log_dens = quantile(dat$log_dens, 0.5),
                                  cougar = seq(quantile(dat$cougar, 0.05),
                                               quantile(dat$cougar, 0.95),
                                               length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1))

# Predict lambda
open_cougar_x1$q05 <- open_cougar_x1$q10 <- open_cougar_x1$q25 <- 
  open_cougar_x1$rss <- 
  open_cougar_x1$q75 <- open_cougar_x1$q90 <- open_cougar_x1$q95 <- NA

for (i in 1:nrow(open_cougar_x1)) {
  cat("               \r", i)
  
  lambda_x1 <- predict_lambda(open_cougar_x1, years, i, b, e = NA, s, 
                              RE = FALSE, mu_eta = mu_eta)
  lambda_x2 <- predict_lambda(open_cougar_x2, years, i, b, e = NA, s, 
                              RE = FALSE, mu_eta = mu_eta)
  rss <- lambda_x1/lambda_x2
  
  # Summarize posterior
  open_cougar_x1$q05[i] <- quantile(rss, 0.05)
  open_cougar_x1$q10[i] <- quantile(rss, 0.10)
  open_cougar_x1$q25[i] <- quantile(rss, 0.25)
  open_cougar_x1$rss[i] <- mean(rss)
  open_cougar_x1$q75[i] <- quantile(rss, 0.75)
  open_cougar_x1$q90[i] <- quantile(rss, 0.90)
  open_cougar_x1$q95[i] <- quantile(rss, 0.95)
}

# Plot
(open_cougar_plot <- ggplot(open_cougar_x1, aes(x = cougar, y = rss)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
                alpha = 0.5) +
    geom_line() +
    coord_cartesian(ylim = c(1.25, 3.25)) +
    xlab(expression("Cougar Density" ~ (cougars/100~km^2))) +
    # ylab("RSS for Openness") +
    ylab(NULL) +
    theme_bw())

# ... ... combine ----
open_pred_plot <- open_wolf_plot + open_cougar_plot

ggsave("fig/predator_open.tiff", plot = open_pred_plot,
       width = 7, height = 3, units = "in", compression = "lzw",
       device = agg_tiff)

# ... roughness ----
# ... ... wolves ----
rough_wolf_x1 <- create_pred_dat(dat, rough = rough_x1,
                                 log_dens = quantile(dat$log_dens, 0.5),
                                 wolf = seq(quantile(dat$wolf, 0.05),
                                            quantile(dat$wolf, 0.95),
                                            length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1),
         wolf_dens = wolf * wa_wd)

rough_wolf_x2 <- create_pred_dat(dat, rough = rough_x2,
                                 log_dens = quantile(dat$log_dens, 0.5),
                                 wolf = seq(quantile(dat$wolf, 0.05),
                                            quantile(dat$wolf, 0.95),
                                            length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1),
         wolf_dens = wolf * wa_wd)

# Predict RSS
rough_wolf_x1$q05 <- rough_wolf_x1$q10 <- rough_wolf_x1$q25 <- 
  rough_wolf_x1$rss <- 
  rough_wolf_x1$q75 <- rough_wolf_x1$q90 <- rough_wolf_x1$q95 <- NA

for (i in 1:nrow(rough_wolf_x1)) {
  cat("               \r", i)
  
  lambda_x1 <- predict_lambda(rough_wolf_x1, years, i, b, e = NA, s, 
                              RE = FALSE, mu_eta = mu_eta)
  lambda_x2 <- predict_lambda(rough_wolf_x2, years, i, b, e = NA, s, 
                              RE = FALSE, mu_eta = mu_eta)
  rss <- lambda_x1/lambda_x2
  
  # Summarize posterior
  rough_wolf_x1$q05[i] <- quantile(rss, 0.05)
  rough_wolf_x1$q10[i] <- quantile(rss, 0.10)
  rough_wolf_x1$q25[i] <- quantile(rss, 0.25)
  rough_wolf_x1$rss[i] <- mean(rss)
  rough_wolf_x1$q75[i] <- quantile(rss, 0.75)
  rough_wolf_x1$q90[i] <- quantile(rss, 0.90)
  rough_wolf_x1$q95[i] <- quantile(rss, 0.95)
}

# Plot
(rough_wolf_plot <- ggplot(rough_wolf_x1, aes(x = wolf_dens, y = rss)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
                alpha = 0.5) +
    geom_line() +
    coord_cartesian(ylim = c(1, 2.0)) +
    xlab(expression("Wolf Density" ~ (wolves/100~km^2))) +
    ylab("RSS for Roughness") +
    theme_bw())

# ... ... cougars ----
rough_cougar_x1 <- create_pred_dat(dat, rough = rough_x1,
                                   log_dens = quantile(dat$log_dens, 0.5),
                                   cougar = seq(quantile(dat$cougar, 0.05),
                                                quantile(dat$cougar, 0.95),
                                                length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1))

rough_cougar_x2 <- create_pred_dat(dat, rough = rough_x2,
                                   log_dens = quantile(dat$log_dens, 0.5),
                                   cougar = seq(quantile(dat$cougar, 0.05),
                                                quantile(dat$cougar, 0.95),
                                                length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1))

# Predict lambda
rough_cougar_x1$q05 <- rough_cougar_x1$q10 <- rough_cougar_x1$q25 <- 
  rough_cougar_x1$rss <- 
  rough_cougar_x1$q75 <- rough_cougar_x1$q90 <- rough_cougar_x1$q95 <- NA

for (i in 1:nrow(rough_cougar_x1)) {
  cat("               \r", i)
  
  lambda_x1 <- predict_lambda(rough_cougar_x1, years, i, b, e = NA, s, 
                              RE = FALSE, mu_eta = mu_eta)
  lambda_x2 <- predict_lambda(rough_cougar_x2, years, i, b, e = NA, s, 
                              RE = FALSE, mu_eta = mu_eta)
  rss <- lambda_x1/lambda_x2
  
  # Summarize posterior
  rough_cougar_x1$q05[i] <- quantile(rss, 0.05)
  rough_cougar_x1$q10[i] <- quantile(rss, 0.10)
  rough_cougar_x1$q25[i] <- quantile(rss, 0.25)
  rough_cougar_x1$rss[i] <- mean(rss)
  rough_cougar_x1$q75[i] <- quantile(rss, 0.75)
  rough_cougar_x1$q90[i] <- quantile(rss, 0.90)
  rough_cougar_x1$q95[i] <- quantile(rss, 0.95)
}

# Plot
(rough_cougar_plot <- ggplot(rough_cougar_x1, aes(x = cougar, y = rss)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
                alpha = 0.5) +
    geom_line() +
    coord_cartesian(ylim = c(1, 2.0)) +
    xlab(expression("Cougar Density" ~ (cougars/100~km^2))) +
    # ylab("RSS for Roughness") +
    ylab(NULL) +
    theme_bw())

# ... ... combine ----
rough_pred_plot <- rough_wolf_plot + rough_cougar_plot

ggsave("fig/predator_rough.tiff", plot = rough_pred_plot,
       width = 7, height = 3, units = "in", compression = "lzw",
       device = agg_tiff)

# ... alternative roughness plot ----
# Plot location of vertex

# ... ... wolves ----
rough_wolf_alt <- create_pred_dat(dat, rough = seq(18, 33, length.out = 100),
                                  log_dens = quantile(dat$log_dens, 0.5),
                                  wolf = seq(quantile(dat$wolf, 0.05),
                                             quantile(dat$wolf, 0.95),
                                             length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1),
         wolf_dens = wolf * wa_wd)

# Predict lambda
rough_by_wolf <- matrix(NA, nrow = nrow(b), 
                        ncol = length(unique(rough_wolf_alt$wolf_dens)))

for (i in 1:nrow(b)) {
  cat("               \r", i)
  
  xx <- rough_wolf_alt %>% 
    mutate(lambda = predict_lambda_iter(rough_wolf_alt, i, b, 
                                        e = NA, s = NA, RE = FALSE)) %>% 
    group_by(wolf_dens) %>% 
    filter(lambda == max(lambda)) %>% 
    pull(rough_orig)
  
  rough_by_wolf[i, ] <- xx
}

rough_wolf_vert <- data.frame(wolf_dens = sort(unique(rough_wolf_alt$wolf_dens))) %>% 
  mutate(q05 = apply(rough_by_wolf, 2, quantile, 0.05),
         q10 = apply(rough_by_wolf, 2, quantile, 0.10),
         q25 = apply(rough_by_wolf, 2, quantile, 0.25),
         mean = apply(rough_by_wolf, 2, mean),
         q75 = apply(rough_by_wolf, 2, quantile, 0.75),
         q90 = apply(rough_by_wolf, 2, quantile, 0.90),
         q95 = apply(rough_by_wolf, 2, quantile, 0.95))

# Plot
(rough_wolf_alt_plot <- ggplot(rough_wolf_vert, aes(x = wolf_dens, y = mean)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
                alpha = 0.5) +
    geom_line() +
    coord_cartesian(ylim = c(18.5, 31)) +
    xlab(expression("Wolf Density" ~ (wolves/100~km^2))) +
    ylab("Roughness vertex (m)") +
    theme_bw())

# ... ... cougars ----
rough_cougar_alt <- create_pred_dat(dat, rough = seq(18, 33, length.out = 100),
                                    log_dens = quantile(dat$log_dens, 0.5),
                                    cougar = seq(quantile(dat$cougar, 0.05),
                                                 quantile(dat$cougar, 0.95),
                                                 length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1),
         cougar_dens = cougar)

# Predict lambda
rough_by_cougar <- matrix(NA, nrow = nrow(b), 
                          ncol = length(unique(rough_cougar_alt$cougar_dens)))

for (i in 1:nrow(b)) {
  cat("               \r", i)
  
  xx <- rough_cougar_alt %>% 
    mutate(lambda = predict_lambda_iter(rough_cougar_alt, i, b, 
                                        e = NA, s = NA, RE = FALSE)) %>% 
    group_by(cougar_dens) %>% 
    filter(lambda == max(lambda)) %>% 
    pull(rough_orig)
  
  rough_by_cougar[i, ] <- xx
  
}

rough_cougar_vert <- data.frame(cougar_dens = sort(unique(rough_cougar_alt$cougar_dens))) %>% 
  mutate(q05 = apply(rough_by_cougar, 2, quantile, 0.05),
         q10 = apply(rough_by_cougar, 2, quantile, 0.10),
         q25 = apply(rough_by_cougar, 2, quantile, 0.25),
         mean = apply(rough_by_cougar, 2, mean),
         q75 = apply(rough_by_cougar, 2, quantile, 0.75),
         q90 = apply(rough_by_cougar, 2, quantile, 0.90),
         q95 = apply(rough_by_cougar, 2, quantile, 0.95))

# Plot
(rough_cougar_alt_plot <- ggplot(rough_cougar_vert, aes(x = cougar_dens, y = mean)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
                alpha = 0.5) +
    geom_line() +
    coord_cartesian(ylim = c(18.5, 31)) +
    xlab(expression("Cougar Density" ~ (cougars/100~km^2))) +
    ylab(NULL) +
    theme_bw() +
    theme(axis.text.y = element_blank()))

# ... ... combine ----
rough_pred_alt <- rough_wolf_alt_plot + rough_cougar_alt_plot

ggsave("fig/predator_rough_alt.tiff", plot = rough_pred_alt,
       width = 7, height = 3, units = "in", compression = "lzw",
       device = agg_tiff)

# Hypothetical high vs low RSS ----
# Hold all temporal factors at 2008 values (arbitrary)
# Map RSS vs mean conditions
# Map log_dens from 5th quantile and 95th quantile

# High density
hi_dens <- dat %>% 
  filter(year == 2008) %>% 
  mutate(log_dens = quantile(dat$log_dens, 0.95))

hi_dens_mean <- hi_dens %>% 
  summarize(across(everything(), mean))

x1_hi <- hi_dens  %>% 
  dplyr::select(x, y, cell, N_ssm,
                swe = swe_orig, rough = rough_orig, elev = elev_orig,
                cos_asp, sin_asp, open = open_orig, biomass = biomass_orig,
                hunt, outside, wolf, cougar, log_dens) %>% 
  intxn() %>% 
  scale_dat(scale_df = scale_df)

x2_hi <- hi_dens_mean  %>% 
  dplyr::select(x, y, cell, N_ssm,
                swe = swe_orig, rough = rough_orig, elev = elev_orig,
                cos_asp, sin_asp, open = open_orig, biomass = biomass_orig,
                hunt, outside, wolf, cougar, log_dens) %>% 
  intxn() %>% 
  scale_dat(scale_df = scale_df)

# Low density
lo_dens <- dat %>% 
  filter(year == 2008) %>% 
  mutate(log_dens = quantile(dat$log_dens, 0.05))

lo_dens_mean <- lo_dens %>% 
  summarize(across(everything(), mean))

x1_lo <- lo_dens  %>% 
  dplyr::select(x, y, cell, N_ssm,
                swe = swe_orig, rough = rough_orig, elev = elev_orig,
                cos_asp, sin_asp, open = open_orig, biomass = biomass_orig,
                hunt, outside, wolf, cougar, log_dens) %>% 
  intxn() %>% 
  scale_dat(scale_df = scale_df)

x2_lo <- lo_dens_mean  %>% 
  dplyr::select(x, y, cell, N_ssm,
                swe = swe_orig, rough = rough_orig, elev = elev_orig,
                cos_asp, sin_asp, open = open_orig, biomass = biomass_orig,
                hunt, outside, wolf, cougar, log_dens) %>% 
  intxn() %>% 
  scale_dat(scale_df = scale_df)

x1_hi <- x1_hi %>% 
  mutate(q05 = NA,
         q10 = NA,
         q25 = NA,
         rss = NA,
         q75 = NA,
         q90 = NA,
         q95 = NA)

x1_lo <- x1_lo %>% 
  mutate(q05 = NA,
         q10 = NA,
         q25 = NA,
         rss = NA,
         q75 = NA,
         q90 = NA,
         q95 = NA)

for (i in 1:nrow(x1_hi)) {
  cat("               \r", i)
  
  ## High density
  lambda_x1_hi <- predict_lambda(x1_hi, years, i, b, e = NA, s, 
                                 RE = FALSE, mu_eta = mu_eta)
  # Only need to calculate x2 once
  if (i == 1) {
    lambda_x2_hi <- predict_lambda(x2_hi, years, i, b, e = NA, s, 
                                   RE = FALSE, mu_eta = mu_eta)
  }
  
  rss_hi <- lambda_x1_hi/lambda_x2_hi
  
  # Summarize posterior
  x1_hi$q05[i] <- quantile(rss_hi, 0.05)
  x1_hi$q10[i] <- quantile(rss_hi, 0.10)
  x1_hi$q25[i] <- quantile(rss_hi, 0.25)
  x1_hi$rss[i] <- mean(rss_hi)
  x1_hi$q75[i] <- quantile(rss_hi, 0.75)
  x1_hi$q90[i] <- quantile(rss_hi, 0.90)
  x1_hi$q95[i] <- quantile(rss_hi, 0.95)
  
  ## Low density
  lambda_x1_lo <- predict_lambda(x1_lo, years, i, b, e = NA, s, 
                                 RE = FALSE, mu_eta = mu_eta)
  
  # Only need to calculate x2 once
  if (i == 1) {
    lambda_x2_lo <- predict_lambda(x2_lo, years, i, b, e = NA, s, 
                                   RE = FALSE, mu_eta = mu_eta)
  }
  
  rss_lo <- lambda_x1_lo/lambda_x2_lo
  
  # Summarize posterior
  x1_lo$q05[i] <- quantile(rss_lo, 0.05)
  x1_lo$q10[i] <- quantile(rss_lo, 0.10)
  x1_lo$q25[i] <- quantile(rss_lo, 0.25)
  x1_lo$rss[i] <- mean(rss_lo)
  x1_lo$q75[i] <- quantile(rss_lo, 0.75)
  x1_lo$q90[i] <- quantile(rss_lo, 0.90)
  x1_lo$q95[i] <- quantile(rss_lo, 0.95)
  
  
}

# Combine hi and lo
map_rss <- bind_rows("High" = x1_hi,
                     "Low" = x1_lo,
                     .id = "dens_label") %>% 
  mutate(log_rss = log(rss))

# Map
(hi_v_lo_rss <- map_rss %>% 
    ggplot(aes(x = x, y = y, fill = log_rss)) +
    facet_wrap(~ dens_label) +
    geom_raster() +
    coord_sf(crs = 32612) +
    xlab(NULL) +
    ylab(NULL) +
    scale_fill_gradient2(name = "log-RSS") +
    theme_bw())

ggsave("fig/high_v_low_RSS.png", plot = hi_v_lo_rss,
       width = 10, height = 4, units = "in",
       device = agg_png)

diff_pal <- rev(brewer.pal(3, "BrBG"))
diff_pal[1] <- colorspace::darken(diff_pal[1], 0.4)
diff_pal[3] <- colorspace::darken(diff_pal[3], 0.4)

(diff_rss <- map_rss %>% 
    arrange(cell, dens_label) %>% 
    group_by(cell) %>% 
    mutate(diff = log_rss[2] - log_rss[1]) %>% 
    # OMIT THE CELL WITH GREATEST 'diff'
    # This cell is dominated by Dailey Lake
    filter(cell != 35) %>% 
    ggplot(aes(x = x, y = y, fill = diff)) +
    geom_raster() +
    coord_sf(crs = 32612)+
    xlab(NULL) +
    ylab(NULL) +
    scale_fill_gradient2(name = "Δ log-RSS", 
                         high = diff_pal[1],
                         mid = diff_pal[2], 
                         low = diff_pal[3], midpoint = 0) +
    theme_bw())

ggsave("fig/high_v_low_RSS_diff.tif", plot = diff_rss,
       width = 5, height = 4, units = "in",
       device = agg_tiff, compression = "lzw")

map_rss %>% 
  arrange(cell, dens_label) %>% 
  group_by(cell) %>% 
  mutate(diff = log_rss[2] - log_rss[1]) %>% 
  ungroup() %>% 
  group_by(dens_label) %>% 
  summarize(mean = mean(diff),
            min(diff),
            max(diff))

# Biomass drives the very dark pixel in the north
# (very low value; pixel is mostly Dailey Lake)
map_rss %>% 
  ggplot(aes(x = x, y = y, fill = biomass_orig)) +
  geom_raster() +
  coord_sf(crs = 32612)+
  xlab(NULL) +
  ylab(NULL)

# Manuscript figures ----
# ... Figure 1 ----
# Created in separate project, "Density Prediction Figures"

# ... Figure 2 ----
fig2 <- cp

ggsave("fig/ms/fig2.tif", plot = fig2, 
       width = 82, height = 170, units = "mm", dpi = 500,
       compression = "lzw", device = agg_tiff)

# ... Figure 3 ----
fig3 <- 
  (# Top row
    (open_wolf_plot + open_cougar_plot) /
      # Bottom row
      (rough_wolf_alt_plot + rough_cougar_alt_plot)) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")

ggsave("fig/ms/fig3.tif", plot = fig3, device = agg_tiff,
       width = 173, height = 100, units = "mm", dpi = 500,
       scale = 1.2, compression = "lzw")

# ... Figure 4 ----
driver_sd1 <- bind_rows("bio" = bio_ddhs_x1,
                        "open" = open_ddhs_x1,
                        "rough" = rough_ddhs_x1,
                        "safe" = safe_ddhs_x1,
                        .id = "driver") %>% 
  dplyr::select(driver, dens, rss, q95:q05) %>% 
  pivot_longer(cols = c(q95:q75), 
               names_to = "quantile_upr", values_to = "upr") %>% 
  pivot_longer(cols = c(q25:q05), 
               names_to = "quantile_lwr", values_to = "lwr") %>% 
  mutate(level_upr = case_when(
    quantile_upr == "q95" ~ "90%",
    quantile_upr == "q90" ~ "80%",
    quantile_upr == "q75" ~ "50%",
    TRUE ~ NA_character_
  ),
  level_lwr = case_when(
    quantile_lwr == "q05" ~ "90%",
    quantile_lwr == "q10" ~ "80%",
    quantile_lwr == "q25" ~ "50%",
    TRUE ~ NA_character_
  )) %>% 
  filter(level_upr == level_lwr) %>% 
  mutate(driver_level = paste(driver, level_upr, sep = "_"))

# Plot
(fig4 <- driver_sd1 %>% 
    filter(driver != "rough",
           driver != "open") %>% 
    ggplot(aes(x = dens, y = rss, 
               color = driver, fill = driver_level)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr),
                alpha = 0.2, show.legend = FALSE, color = NA) +
    geom_line(size = 0.7) +
    scale_color_manual(name = "Habitat Variable",
                       breaks = c("bio", "safe"),
                       labels = c("Food", "Safety"),
                       values = c(color_a, color_b)) +
    scale_fill_manual(name = "Habitat Variable",
                      breaks = c("bio_90%", "bio_80%", "bio_50%",
                                 "open_90%", "open_80%", "open_50%"),
                      values = c(color_a_90, color_a_80, color_a_50,
                                 color_b_90, color_b_80, color_b_50)) +
    xlab(expression("Average Elk Density" ~ (elk/km^2))) +
    ylab("RSS") +
    coord_cartesian(ylim = c(1.1, 2.1)) +
    theme_bw())

ggsave("fig/ms/fig4.tif", plot = fig4, device = agg_tiff,
       width = 173, height = 80, units = "mm", dpi = 500,
       compression = "lzw")

# Numbers
driver_sd1 %>% 
  filter(driver != "rough",
         driver != "open") %>% 
  filter(dens == 2)

driver_sd1 %>% 
  filter(driver != "rough",
         driver != "open") %>% 
  filter(dens == 9.3)

# ... Figure 5 ----
fig5 <- diff_rss

# Full page
ggsave("fig/ms/fig5.tif", plot = fig5, device = agg_tiff,
       width = 173, height = 100, units = "mm", dpi = 500,
       compression = "lzw")

# Tweetable
# ggsave("fig/tweet/fig5.png", plot = fig5,
#        width = 5, height = 4, units = "in",
#        device = agg_png)

# ... Figure S1 ----
# Study area map
# Made in QGIS

# ... Figure S2 ----
# Survey dates and times
# See data preparation project

# ... Figure S3 ----
# Residuals by survey month
# See script 03_spatio-tempo_corr.R

# ... Figure S4 ----
# Scale of analysis
# Made in PowerPoint

# ... Figure S5 ----
# Openness sensitivity analysis
# See script 07_open_sens_fig.R

# ... Figure S6 ----
# Mean effect of SWE, elev, asp
figS6 <- (
  (swe_dens_plot) + 
    (elev_dens_plot +
       ylab(NULL) +
       # theme(axis.text.y = element_blank())+
       NULL) + 
    (asp_dens_plot +
       ylab(NULL) +
       # theme(axis.text.y = element_blank()) +
       NULL)
) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(nrow = 1) &
  coord_cartesian(ylim = c(0, 30))

ggsave("fig/ms/figS6.tif", plot = figS6, device = agg_tiff,
       width = 170, height = 70, units = "mm", dpi = 500,
       scale = 1.2, compression = "lzw")

# ... Figure S7 ----
# Mean effect of (A) biomass, (B) openness, (C) roughness
figS7 <- (
  (bio_dens_plot) + 
    (open_dens_plot +
       ylab(NULL) +
       theme(axis.text.y = element_blank())) + 
    (rough_dens_plot +
       ylab(NULL) +
       theme(axis.text.y = element_blank()))) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(nrow = 1) &
  coord_cartesian(ylim = c(0, 7))

ggsave("fig/ms/figS7.tif", plot = figS7, device = agg_tiff,
       width = 170, height = 70, units = "mm", dpi = 500,
       scale = 1.2,
       compression = "lzw")

# ... Figure S8 ----
# Roughness RSS
figS8 <- rough_pred_plot +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")

ggsave("fig/ms/figS8.tif", plot = figS8, device = agg_tiff,
       width = 170, height = 70, units = "mm", dpi = 500,
       scale = 1.2, compression = "lzw")

# ... Figure S9 ----
# DDHS for roughness compared to openness and biomass
figS9a <- rough_ddhs_plot +
  ylab("RSS") +
  coord_cartesian(ylim = c(1.0, 2.0))

figS9b <- driver_sd1 %>% 
  filter(driver != "safe") %>% 
  ggplot(aes(x = dens, y = rss, color = driver)) +
  geom_line(size = 1) +
  scale_color_manual(name = "Habitat Variable",
                     breaks = c("bio", "open", "rough"),
                     labels = c("Biomass", "Openness", "Roughness"),
                     values = c(color_a, color_b, "gray70")) +
  xlab(expression("Average Elk Density" ~ (elk/km^2))) +
  ylab(NULL) +
  coord_cartesian(ylim = c(1.2, 2.5)) +
  theme_bw() +
  theme(axis.text.y = element_blank())

figS9 <- figS9a + figS9b + 
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")

ggsave("fig/ms/figS9.tif", plot = figS9, device = agg_tiff,
       width = 170, height = 100, units = "mm", dpi = 500,
       scale = 1.2,
       compression = "lzw")

# ... Figure S10 ----
# See script 05_oos_valid.R

# ... Figure S11 ----
# Residual spatial autocorrelation
# See script 03_spatio-tempo_corr.R

# ... Figure S12 ----
# Trend in eta
# Temporal random effect
figS12 <- data.frame(year = years,
                    q05 = apply(e, 2, quantile, 0.05),
                    q10 = apply(e, 2, quantile, 0.10),
                    q25 = apply(e, 2, quantile, 0.25),
                    mean = apply(e, 2, mean),
                    q75 = apply(e, 2, quantile, 0.75),
                    q90 = apply(e, 2, quantile, 0.90),
                    q95 = apply(e, 2, quantile, 0.95)) %>% 
  ggplot(aes(x = year, y = mean)) +
  geom_errorbar(aes(ymin = q05, ymax = q95), color = color_90, width = 0.2, size = 0.5) +
  geom_errorbar(aes(ymin = q10, ymax = q90), color = color_80, width = 0, size = 0.75) +
  geom_errorbar(aes(ymin = q25, ymax = q75), color = color_50, width = 0, size = 1) +
  geom_point(size = 1) +
  xlab("Year") +
  ylab(expression(eta["i, t"])) +
  theme_bw() 

ggsave("fig/ms/figS12.tif", plot = figS12, device = agg_tiff,
       width = 170, height = 70, units = "mm", dpi = 500,
       scale = 1, compression = "lzw")

# ... Figure S13 ----
# Supplemental figure with spatial effect
(figS13 <- s_fig +
   geom_sf(data = ynp, aes(geometry = geometry), 
           color = "#DAA520FF", 
           # fill = "#DAA52025", 
           fill = NA,
           size = 1.5) +
   coord_sf(xlim = range(s_dat$x),
            ylim = range(s_dat$y),
            default_crs = 26912,
            datum = 4326) +
   scale_fill_viridis_c(name = expression(s[i]), option = "D") +
   theme_bw())
ggsave("fig/ms/figS13.tif", plot = figS13, 
       width = 170, height = 160, units = "mm", dpi = 500,
       compression = "lzw", device = agg_tiff)

# ... Figure S14 ----
# Spatial covariance decay
# Distances are scaled by max
max_dist <- max(dist(unique(dat[, c("x", "y")])))

# Exponential covariance
expcov <- function(dist, rho, sigma) {
  return(sigma^2 * exp(-1 * dist/rho))
}
s14_dat <- data.frame(dist = seq(0, 1.25, length.out = 100))

s14_dat <- lapply(split(s14_dat, 1:100), function(x) {
  d <- x$dist
  Sigma <- expcov(d, rho, sig)
  x$q05 <- quantile(Sigma, 0.05)
  x$q10 <- quantile(Sigma, 0.10)
  x$q25 <- quantile(Sigma, 0.25)
  x$mean <- mean(Sigma)
  x$q75 <- quantile(Sigma, 0.75)
  x$q90 <- quantile(Sigma, 0.90)
  x$q95 <- quantile(Sigma, 0.95)
  return(x)
}) %>% 
  bind_rows() %>% 
  mutate(dist_km = dist * max_dist / 1000)

(figS14 <- ggplot(s14_dat, aes(x = dist_km, y = mean)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
                alpha = 0.5) +
    geom_line() +
    xlab("Distance (km)") +
    ylab(expression(Sigma["i, j"])) +
    theme_bw())

ggsave("fig/ms/figS14.tif", plot = figS14, 
       width = 170, height = 160, units = "mm", dpi = 500,
       compression = "lzw", device = agg_tiff)

# Figure S15 ----
# Combined random effects

# Outside?
out <- distinct(dat, cell, outside)
s_dat$cell <- NA
s_dat$cell[which(!is.na(s_dat$layer))] <- 1:1978

# Mean eta by year
eta_year <- data.frame(year = years,
                       mean = apply(e, 2, mean))

s_dat_ex <- s_dat %>% 
  left_join(out) %>% 
  as_tibble() %>% 
  mutate(year_list = lapply(x, function(x) {
    return(tibble(year = c(1990, 2007, 2013, 2017)))
  })) %>% 
  unnest(cols = year_list) %>% 
  arrange(year, x, y) %>% 
  left_join(eta_year) %>% 
  mutate(comb = case_when(
    outside == 1 ~ layer + mean,
    outside == 0 ~ layer,
    TRUE ~ NA_real_
  ))

figS15 <- ggplot() +
  geom_raster(data = s_dat_ex, aes(x = x, y = y, fill = comb)) +
  facet_wrap(~ year) +
  scale_fill_viridis_c(name = "Total\nRandom\nEffect") +
  xlab("Easting") +
  ylab("Northing") +
  coord_sf(crs = 26912) +
  theme_bw()

ggsave("fig/ms/figS15.tif", plot = figS15, 
       width = 170, height = 160, units = "mm", dpi = 500,
       compression = "lzw", device = agg_tiff)

# Figure S16 ----
# Cougar densities
# Created in data preparation project
# See 06_predators.R in that project
