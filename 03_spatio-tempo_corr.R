############################################X
#--Elk Density-Dependent Habitat Selection--X
#---------------Brian J. Smith--------------X
#--Residual Spatio-temporal Autocorrelation-X
#===========================================X
#-----------Last update 2022-09-21----------X
############################################X

#Load packages----
library(tidyverse)
library(sf)
library(coda)
library(ncf)
library(ragg)

# Custom functions ----
source("99_fun.R")

#Load data----
# Model data
dat <- read.csv("data/all_data.csv")
#Times
years <- sort(unique(dat$year))
# Date summary
dates <- readRDS("../Elk_Dens_Data/out/dates_summary.rds")

# MCMC samples ----
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
samples1 <- lapply(samples1, function(x) {
  return(x[seq(20020, nrow(x), by = 20),])
})

samples2 <- lapply(samples2, function(x) {
  return(x[seq(20020, nrow(x), by = 20),])
})

# Garbage cleanup
gc()

# Load other data ----
#Northern range shapefile
nr <- st_read(paste0("../../../Data/GIS/YNP/Northern Range/NR_Bound_Revised/", 
                     "NR_bound_revised_2021.shp")) %>% 
  st_transform(crs = 26912)

#Yellowstone boundary
ynp <- st_read("../../../Data/GIS/YNP/YNP_boundary/YNP_boundary.shp") %>% 
  st_transform(crs = 26912)

# Extract coefficients ----
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
names(e) <- paste0("e", years)

s <- lapply(samples2, function(x) {
  return(as.data.frame(x[1:nrow(x), ]))
}) %>% 
  bind_rows() %>% 
  as.data.frame()
names(s) <- paste0("s", 1:ncol(s))

# Calculate residuals ----
pred <- dat
pred$lambda <- NA
pred$lwr <- NA
pred$upr <- NA
pred$var <- NA
lambda_samps <- list()

for (i in 1:nrow(pred)) {
  cat("\r", i, "       ")
  
  lambda_samps[[i]] <- predict_lambda(dat = pred, years = years, 
                                      i = i, b = b, e = e, s = s)
  
  pred$lambda[i] <- mean(lambda_samps[[i]])
  pred$lwr[i] <- quantile(lambda_samps[[i]], 0.05)
  pred$upr[i] <- quantile(lambda_samps[[i]], 0.95)
  pred$var[i] <- var(lambda_samps[[i]])
  
}

# Ordinary residual
pred$resid <- pred$n - pred$lambda
# Pearson residual
pred$pearson <- pred$resid/sqrt(pred$var)

# Temporal correlation ----
# Does date of survey affect residuals?

# Join dates to data
pred <- left_join(pred, dates, by = c("year" = "winter"))

# ... Figure S3 ----
(resid_temp <- pred %>% 
  ggplot(aes(x = mean_date, y = resid)) +
  geom_point(size = 1, alpha = 0.3) +
  geom_smooth(method = "gam", se = TRUE) +
  xlab("Survey Date") +
  ylab("Ordinary Residuals") +
  theme_bw())

ggsave("fig/ms/figS3.tiff", plot = resid_temp, device = agg_tiff,
       width = 170, height = 100, units = "mm", dpi = 500,
       scale = 1.2,
       compression = "lzw")

# Model
pred <- pred %>% 
  mutate(yday = lubridate::yday(mean_date)) %>% 
  mutate(yday_adj = case_when(
    yday > 200 ~ yday - 365,
    TRUE ~ yday
  ))

date_mod <- aov(pearson ~ factor(yday), data = pred)
summary(date_mod)
hsd <- TukeyHSD(date_mod)
range(hsd$`factor(yday)`[,"p adj"])

# Write output to table
hsd_df <- as.data.frame(hsd$`factor(yday)`)
hsd_df$comp <- paste0("'", row.names(hsd_df))
write.csv(hsd_df, "out/TableS2.csv", row.names = FALSE)

# Plot residual maps ----

plot_resid <- function(pred, yr){
  # Subset data.frame
  r_df <- pred %>% 
    filter(year == yr)
  
  pp <- ggplot() + 
    geom_raster(data = r_df, aes(x = x, y = y, fill = pearson)) +
    geom_sf(data = ynp, color = "green", fill = NA) +
    geom_sf(data = nr, color = "goldenrod", fill = NA, size = 1) +
    coord_sf(ylim = c(4955000, 5020000)) +
    ggtitle(paste(yr, "Pearson residual")) +
    xlab(NULL) +
    ylab(NULL) +
    theme_bw()
  
  return(pp)
}

resid_plots <- list()

for (i in 1:length(years)){
  
  #Plot regular
  resid_plots[[i]] <- plot_resid(pred, years[i]) +
    scale_fill_viridis_c(option = "A", 
                         name = "Pearson Residual\n(n - lambda)/sd", 
                         breaks = pretty(-50:1500, 7),
                         limits = c(-200, 1500),
                         alpha = 0.75) +
    NULL
  
  ggsave(paste0("fig/map_resid/resid_map", years[i], "_pearson.tiff"), 
         plot = resid_plots[[i]], device = agg_tiff, 
         width = 8, height = 6, compression = "lzw")
  
}

# Spline correlogram for each year ----
spl_corr <- list()
for (i in 1:length(years)){
  
  cat("\n", years[i], "\n")
  
  temp <- pred %>% 
    filter(year == years[i])
  
  spl_corr[[i]] <- spline.correlog(x = temp$x, y = temp$y,
                                   z = temp$pearson, resamp = 1000,
                                   xmax = 50000)
}

# Save

{
  pdf("fig/pearson_resid_spatial_autocorrelation.pdf",
      width = 8, height = 6, onefile = TRUE)
  
  for (i in 1:length(years)){
    plot(spl_corr[[i]], main = years[i], xlim = c(0, 7500), 
         xlab = "Distance (m)")
  }
  
  dev.off()
}

{
  tiff("fig/ms/figS11.tiff",
      width = 173, height = 200, units = "mm", res = 300,
      compression = "lzw")
  par(mfrow = c(4, 4), cex = 0.5)
  
  for (i in 1:length(years)){
    plot(spl_corr[[i]], main = years[i], xlim = c(0, 7500), 
         xlab = "Distance (m)")
  }
  
  dev.off()
}
