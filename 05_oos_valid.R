############################################X
#--Elk Density-Dependent Habitat Selection--X
#---------------Brian J. Smith--------------X
#----------Out-of-Sample Validation---------X
#===========================================X
#-----------Last update 2022-02-22----------X
############################################X

#Load packages----
library(tidyverse)
library(sf)
library(raster)
library(patchwork)
library(ragg)

# Custom functions ----
source("99_fun.R")

# Colors ----
# 50% credible interval
color_50 <- "black"
# 80% credible interval
color_80 <- "gray40"
# 90% credible interval
color_90 <- "gray70"

#Load data----
# Validation data
dat <- read.csv("data/all_data_test.csv")

#Times
years <- sort(unique(dat$year))

# Discard unreliable data ----
# For years 1994, 2002, and 2004, outside the park data are unreliable
dat <- dat %>% 
  mutate(n = case_when(
    year %in% c(1994, 2002, 2004) & outside == 1 ~ NA_integer_,
    TRUE ~ n
  )) %>% 
  filter(!is.na(n))

# Original data with predictions
# Note, this was created in script "02_figures.R", but the folder "pred"
# is added to .gitignore. This object will not exist in cloned repos unless
# script 02 has already been run.
orig <- readRDS("pred/model_data_predictions.rds")

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

# Scaling factor (elk/pixel -> elk/km2)
pixel_scale <- 1/((res(trast)[1])^2/(1000^2))

# Create figure directories ----
# Adding figures to .gitignore, so directory will not exist in cloned repos
dir.create("fig/validation", showWarnings = FALSE)
dir.create("fig/validation/map_scaled", showWarnings = FALSE)
dir.create("fig/validation/map_unscaled", showWarnings = FALSE)
dir.create("fig/validation/map_comp", showWarnings = FALSE)

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

# Predict expected density ----
pred <- dat
pred$lambda <- NA
pred$lwr <- NA
pred$upr <- NA
pred$var <- NA
lambda_samps <- list()

for (i in 1:nrow(pred)) {
  cat("\r", i, "of", nrow(pred), "       ")
  
  lambda_samps[[i]] <- predict_lambda(pred, years, i, b, e, s)
  
  pred$lambda[i] <- mean(lambda_samps[[i]])
  pred$lwr[i] <- quantile(lambda_samps[[i]], 0.05)
  pred$upr[i] <- quantile(lambda_samps[[i]], 0.95)
  pred$var[i] <- var(lambda_samps[[i]])
  
}

# Plot maps
lambda_maps <- list()
n_maps <- list()
comp_maps <- list()

for (i in 1:length(years)){
  
  #Plot
  lambda_maps[[i]] <- map_lambda(pred, years[i], ynp, nr)
  n_maps[[i]] <- map_n(pred, years[i], ynp, nr)
  
  # Compare
  comp_maps[[i]] <- (lambda_maps[[i]] +
                       scale_fill_viridis_c(option = "A", 
                                            expression("Expected\nDensity" * (elk/km^2)),
                                            breaks = pretty(0:550),
                                            limits = c(0, 551),
                                            alpha = 0.75) +
                       ggtitle("Expected Density")) /
    (n_maps[[i]] +
       scale_fill_viridis_c(option = "D", 
                            expression("Actual\nCount" * (elk/km^2)),
                            breaks = pretty(0:550),
                            limits = c(0, 551),
                            alpha = 0.75) +
       ggtitle("Actual Count")) +
    plot_annotation(title = years[i])
  
  #Save on same scale
  ggsave(paste0("fig/validation/map_scaled/N_it_raster_", 
                years[i], ".tiff"), 
         plot = lambda_maps[[i]] +
           scale_fill_viridis_c(option = "A", 
                                expression("Expected\nDensity" * (elk/km^2)),
                                breaks = pretty(0:550),
                                limits = c(0, 551),
                                alpha = 0.75), 
         width = 10, height = 6, units = "in", compression = "lzw", 
         device = agg_tiff)
  
  #Save on different scale
  ggsave(paste0("fig/validation/map_unscaled/N_it_raster_", 
                years[i], ".tiff"), 
         plot = lambda_maps[[i]] +
           scale_fill_viridis_c(option = "D", 
                                name = expression("Expected\nDensity" * (elk/km^2)),
                                alpha = 0.75), 
         width = 10, height = 6, units = "in", compression = "lzw", 
         device = agg_tiff)
  
  # Save together
  ggsave(paste0("fig/validation/map_comp/comparison_map_", 
                years[i], ".tiff"), 
         plot = comp_maps[[i]], 
         width = 10, height = 12, units = "in", compression = "lzw", 
         device = agg_tiff)
  
}

comb <- bind_rows("In" = orig,
                  "Out" = pred,
                  .id = "sample")

# Residuals ----
# Ordinary residuals
comb$resid <- comb$n - comb$lambda
# Pearson residuals
comb$pearson <- comb$resid/sqrt(comb$var)

comb %>% 
  group_by(year, sample) %>% 
  summarize(resid = mean(resid),
            pearson = mean(pearson))

# Plot ordinary boxplot
(resid_boxplot <- comb %>% 
  ggplot(aes(x = year, y = resid, group = year, color = factor(sample))) +
  geom_boxplot() +
  xlab("Year") +
  ylab(expression("Residual" ~ (n - lambda))) +
  scale_color_discrete(name = "Sample") +
  theme_bw())
ggsave("fig/validation/resid_boxplot.tif", plot = resid_boxplot, device = agg_tiff,
       width = 8, height = 6, units = "in", compression = "lzw")

comb %>% 
  group_by(sample) %>% 
  summarize(mean(resid), quantile(resid, 0.25), quantile(resid, 0.75))

# Plot pearson boxplot
comb %>% 
  ggplot(aes(x = year, y = pearson, group = year, color = factor(sample))) +
  geom_boxplot() +
  xlab("Year") +
  ylab(expression("Pearson Residual" ~ ((n - lambda)/sigma))) +
  theme_bw()

# Spearman correlation ----
# By year
(sp_r <- comb %>% 
   group_by(year, sample) %>% 
   summarize(r = cor(n, lambda, method = "spearman")) %>% 
   ggplot(aes(x = year, y = r, color = factor(sample))) +
   geom_point() +
   coord_cartesian(ylim = c(0, 0.5)) +
   xlab("Year") +
   ylab("Spearman's Correlation") +
   scale_color_discrete(name = "Sample") +
   theme_bw())
ggsave("fig/validation/spearman.tif", plot = sp_r, device = agg_tiff,
       width = 8, height = 6, units = "in", compression = "lzw")

comb %>% 
  group_by(year, sample) %>% 
  summarize(r = cor(n, lambda, method = "spearman")) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  summarize(mean(r))

# Figure S10 ----
# Colors
cols <- c("navy", "firebrick")
(figS10 <- (resid_boxplot +
            scale_color_manual(name = "Sample",
                               breaks = c("In", "Out"),
                               values = cols)) / 
   (sp_r +
      scale_color_manual(name = "Sample",
                           breaks = c("In", "Out"),
                           values = cols)) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")"))

ggsave("fig/ms/figS10.tif", plot = figS10, 
       width = 173, height = 140, units = "mm", dpi = 500,
       compression = "lzw", device = agg_tiff)
