############################################X
#--Elk Density-Dependent Habitat Selection--X
#---------------Brian J. Smith--------------X
#--------Openness Sensitivity Figures-------X
#===========================================X
#-----------Last update 2022-10-08----------X
############################################X

#Load packages----
library(tidyverse)
library(coda)
library(ragg)
library(patchwork)

# Custom functions ----
source("99_fun.R")

# Colors ----
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


# Original data ----
# Model data
dat <- read.csv("data/all_data.csv")

#Times
years <- sort(unique(dat$year))

# Scaling data.frame
scale_df <- read.csv("data/scale_df.csv")

# SDs
# 0.5 SD of openness and roughness
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

# Openness cutoffs ----
(oc1 <- unname(quantile(dat$open_orig, 0.05)))
(oc2 <- unname(quantile(dat$open_orig, 0.10)))

(oc1 <- 0.30)
(oc2 <- 0.50)
(oc3 <- 0.9999)

# ... histogram ----
open_hist <- dat %>% 
  ggplot() +
  geom_histogram(aes(x = open_orig), bins = 50,
                 color = "black") +
  # OC1
  geom_vline(xintercept = oc1, size = 1,
             color = "firebrick", linetype = "dashed") +
  # geom_segment(aes(x = oc1, xend = oc1 - 0.1, y = 3000, yend = 3000),
  #              color = "firebrick", 
  #              arrow = arrow(length = unit(0.1, "inches"),
  #                            type = "closed")) +
  # OC2
  geom_vline(xintercept = oc2, size = 1,
             color = "blue", linetype = "dashed") +
  # geom_segment(aes(x = oc2, xend = oc2 - 0.1, y = 3000, yend = 3000),
  #              color = "blue", 
  #              arrow = arrow(length = unit(0.1, "inches"),
  #                            type = "closed")) +
  # OC3
  geom_vline(xintercept = oc3, size = 1,
             color = "#00BB00", linetype = "dashed") +
  # geom_segment(aes(x = oc3, xend = oc3 + 0.05, y = 3000, yend = 3000),
  #              color = "#00BB00", 
  #              arrow = arrow(length = unit(0.1, "inches"),
  #                            type = "closed")) +
  scale_x_continuous(name = "Openness (%)",
                     breaks = seq(0, 1, length.out = 5),
                     labels = seq(0, 100, length.out = 5)) +
  ylab("Number of Pixels") +
  coord_cartesian(xlim = c(0, 1.05)) +
  theme_bw() +
  theme(axis.title = element_text(size = 10)) +
  NULL

ggsave("fig/open_sens/open_hist.tiff", plot = open_hist, width = 8, height = 6, 
       device = agg_tiff, units = "in", dpi = 200, compression = "lzw")

# Actual data proportions ----
ncell <- length(unique(dat$cell))
flagged_oc1 <- sort(unique(dat$cell[which(dat$open_orig < oc1)]))
flagged_oc2 <- sort(unique(dat$cell[which(dat$open_orig < oc2)]))
flagged_oc3 <- sort(unique(dat$cell[which(dat$open_orig > oc3)]))

# Proportion remaining
round(prop_oc1 <- (ncell - length(flagged_oc1))/ncell, 3)
round(prop_oc2 <- (ncell - length(flagged_oc2))/ncell, 3)
round(prop_oc3 <- (ncell - length(flagged_oc3))/ncell, 3)

# MCMC samples ----
# Note: the model runs are added to .gitignore, so are not available
# in the GitHub repo (samples from final run are > 1.2 GB)

# (No demo data in GitHub repo)

# ... original data ----
samples1_orig <- list(
  ch1 = readRDS("models/samples1_ch1_2022-02-22.rds"),
  ch2 = readRDS("models/samples1_ch2_2022-02-22.rds")#,
  # ch3 = readRDS("models/samples1_ch3_2022-02-22.rds")
)

samples2_orig <- list(
  ch1 = readRDS("models/samples2_ch1_2022-02-22.rds"),
  ch2 = readRDS("models/samples2_ch2_2022-02-22.rds")#,
  # ch3 = readRDS("models/samples2_ch3_2022-02-22.rds")
)

# Get rid of burn-in and thin
# Burn-in of 20k (samplers adapt for 15k)
# Thinning by 20 leaves 4k samples for inference
samples1_orig <- lapply(samples1_orig, function(x) {
  return(x[seq(20020, nrow(x), by = 20),])
})

samples2_orig <- lapply(samples2_orig, function(x) {
  return(x[seq(20020, nrow(x), by = 20),])
})

# Garbage cleanup
gc()

# ... openness cutoff 1 ----
samples1_oc1 <- list(
  ch1 = readRDS("models/sens/samples1_ch_oc1_1_2022-04-30.rds"),
  ch2 = readRDS("models/sens/samples1_ch_oc1_2_2022-04-29.rds")
)

samples2_oc1 <- list(
  ch1 = readRDS("models/sens/samples2_ch_oc1_1_2022-04-30.rds"),
  ch2 = readRDS("models/sens/samples2_ch_oc1_2_2022-04-29.rds")
)

# Get rid of burn-in and thin
# Want to leave 4k samples for inference (same as original)
# Burn-in of 22k (samplers adapt for 15k)
# Thinning by 7 leaves 4k samples for inference
samples1_oc1 <- lapply(samples1_oc1, function(x) {
  return(x[seq(22007, nrow(x), by = 7),])
})

samples2_oc1 <- lapply(samples2_oc1, function(x) {
  return(x[seq(22007, nrow(x), by = 7),])
})

# ... openness cutoff 2 ----
samples1_oc2 <- list(
  ch1 = readRDS("models/sens/samples1_ch_oc2_1_2022-04-27.rds"),
  ch2 = readRDS("models/sens/samples1_ch_oc2_2_2022-04-27.rds")
)

samples2_oc2 <- list(
  ch1 = readRDS("models/sens/samples2_ch_oc2_1_2022-04-27.rds"),
  ch2 = readRDS("models/sens/samples2_ch_oc2_2_2022-04-27.rds")
)

# Get rid of burn-in and thin
# Want to leave 4k samples for inference (same as original)
# Burn-in of 22k (samplers adapt for 15k)
# Thinning by 7 leaves 4k samples for inference
samples1_oc2 <- lapply(samples1_oc2, function(x) {
  return(x[seq(22007, nrow(x), by = 7),])
})

samples2_oc2 <- lapply(samples2_oc2, function(x) {
  return(x[seq(22007, nrow(x), by = 7),])
})

# Garbage cleanup
gc()

# ... openness cutoff 3 ----
samples1_oc3 <- list(
  ch1 = readRDS("models/sens/samples1_ch_oc3_1_2022-09-23.rds"),
  ch2 = readRDS("models/sens/samples1_ch_oc3_2_2022-09-23.rds")
)

samples2_oc3 <- list(
  ch1 = readRDS("models/sens/samples2_ch_oc3_1_2022-09-23.rds"),
  ch2 = readRDS("models/sens/samples2_ch_oc3_2_2022-09-23.rds")
)

# Get rid of burn-in and thin
# Want to leave 4k samples for inference (same as original)
# Burn-in of 22k (samplers adapt for 15k)
# Thinning by 7 leaves 4k samples for inference
samples1_oc3 <- lapply(samples1_oc3, function(x) {
  return(x[seq(22007, nrow(x), by = 7),])
})

samples2_oc3 <- lapply(samples2_oc3, function(x) {
  return(x[seq(22007, nrow(x), by = 7),])
})

# Create figure directories ----
# Adding figures to .gitignore, so directory will not exist in cloned repos
dir.create("fig", showWarnings = FALSE)
dir.create("fig/open_sens", showWarnings = FALSE)

# Coda diagnostics ----
coda1_oc1 <- lapply(samples1_oc1, as.mcmc)

gelman.diag(coda1_oc1)

coda1_oc2 <- lapply(samples1_oc2, as.mcmc)

gelman.diag(coda1_oc2)

coda1_oc3 <- lapply(samples1_oc3, as.mcmc)

gelman.diag(coda1_oc3)

# Plot coefficients ----
# ... original data ----
b_orig <- lapply(samples1_orig, function(x) {
  return(as.data.frame(x[1:nrow(x), 
                         grep("beta", colnames(samples1_orig[[1]]))]))
}) %>% 
  bind_rows() %>% 
  as.data.frame()
names(b_orig) <- paste0("b", 1:ncol(b_orig))

# ... open cutoff 1 ----
b_oc1 <- lapply(samples1_oc1, function(x) {
  return(as.data.frame(x[1:nrow(x), 
                         grep("beta", colnames(samples1_oc1[[1]]))]))
}) %>% 
  bind_rows() %>% 
  as.data.frame()
names(b_oc1) <- paste0("b", 1:ncol(b_oc1))

# ... open cutoff 2 ----
b_oc2 <- lapply(samples1_oc2, function(x) {
  return(as.data.frame(x[1:nrow(x), 
                         grep("beta", colnames(samples1_oc2[[1]]))]))
}) %>% 
  bind_rows() %>% 
  as.data.frame()
names(b_oc2) <- paste0("b", 1:ncol(b_oc2))

# ... open cutoff 3 ----
b_oc3 <- lapply(samples1_oc3, function(x) {
  return(as.data.frame(x[1:nrow(x), 
                         grep("beta", colnames(samples1_oc3[[1]]))]))
}) %>% 
  bind_rows() %>% 
  as.data.frame()
names(b_oc3) <- paste0("b", 1:ncol(b_oc3))

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
  coef_order <- data.frame(name = paste0("b", 1:22), 
                           order = factor(1:22),
                           label = c("SWE",
                                     "Elevation",
                                     "cos(Aspect)",
                                     "sin(Aspect)",
                                     "Biomass",
                                     "Biomass:log(Dens)",
                                     "Openness",
                                     "Roughness",
                                     "Open^2",
                                     "Rough^2",
                                     "Open:Wolf",
                                     "Rough:Wolf",
                                     "(Open^2):Wolf",
                                     "(Rough^2):Wolf",
                                     "Open:Cougar",
                                     "Rough:Cougar",
                                     "(Open^2):Cougar",
                                     "(Rough^2):Cougar",
                                     "Open:log(Dens)",
                                     "Rough:log(Dens)",
                                     "(Open^2):log(Dens)",
                                     "(Rough^2):log(Dens)"))
  coef_plot <- left_join(coefs, coef_order, by = "name")  
  return(coef_plot)
}

# Plot
bd_orig <- beta_data(b_orig)
bd_oc1 <- beta_data(b_oc1)
bd_oc2 <- beta_data(b_oc2)
bd_oc3 <- beta_data(b_oc3)

bd <- bind_rows("Original" = bd_orig, 
                "> 30% Open" = bd_oc1,
                "> 50% Open" = bd_oc2, 
                "< 99.9% Open" = bd_oc3, 
                .id = "cutoff") %>% 
  mutate(cutoff = factor(cutoff, levels = c("Original", "> 30% Open", 
                                            "> 50% Open", "< 99.9% Open")))

cp <- bd %>% 
  ggplot(aes(y = order, x = mean)) +
  facet_wrap(~ cutoff) +
  geom_vline(xintercept = 0, color = "firebrick", linetype = "dashed", size = 0.5) +
  geom_errorbar(aes(xmin = l90, xmax = u90), color = "gray70", width = 0.2, size = 0.5) +
  geom_errorbar(aes(xmin = l80, xmax = u80), color = "gray40", width = 0, size = 0.75) +
  geom_errorbar(aes(xmin = l50, xmax = u50), color = color_50, width = 0, size = 1) +
  geom_point() +
  scale_color_ordinal(name = "Openness Cutoff") +
  scale_y_discrete(breaks = bd$order, labels = bd$label) +
  ylab(NULL) +
  xlab(expression(beta)) +
  coord_cartesian(xlim = c(-2, 2)) +
  theme_bw() +
  # theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  NULL
ggsave("fig/open_sens/betas.tiff", plot = cp, width = 8, height = 6, 
       device = agg_tiff, units = "in", dpi = 200, compression = "lzw")

# ... difference from orig ----
diff_oc1 <- beta_data(b_orig - b_oc1)
diff_oc2 <- beta_data(b_orig - b_oc2)
diff_oc3 <- beta_data(b_orig - b_oc3)

b_diff <- bind_rows("> 30% Open" = diff_oc1,
                    "> 50% Open" = diff_oc2,
                    "< 99.9% Open" = diff_oc3, 
                    .id = "cutoff") %>% 
  mutate(cutoff = factor(cutoff, levels = c("> 30% Open", 
                                            "> 50% Open",
                                            "< 99.9% Open")))

p_diff <- b_diff %>% 
  ggplot(aes(y = order, x = mean)) +
  facet_wrap(~ cutoff) +
  geom_vline(xintercept = 0, color = "firebrick", linetype = "dashed", size = 0.5) +
  geom_errorbar(aes(xmin = l90, xmax = u90), color = "gray70", width = 0.2, size = 0.5) +
  geom_errorbar(aes(xmin = l80, xmax = u80), color = "gray40", width = 0, size = 0.75) +
  geom_errorbar(aes(xmin = l50, xmax = u50), color = color_50, width = 0, size = 1) +
  geom_point() +
  scale_color_ordinal(name = "Openness Cutoff") +
  scale_y_discrete(breaks = bd$order, labels = bd$label) +
  ylab(NULL) +
  xlab(expression(beta)) +
  coord_cartesian(xlim = c(-3, 3)) +
  theme_bw() +
  # theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  NULL
ggsave("fig/open_sens/beta_diff.tiff", plot = p_diff, width = 8, height = 6, 
       device = agg_tiff, units = "in", dpi = 200, compression = "lzw")


# Coefficient interpretation ----

# ... biomass conversion ----
# 1 lb = 0.45359 kg
# 1 acre = 0.40469 ha
# ==> 1 lb/acre = 0.45359/0.40469 kg/ha
lbac_kgha <- 0.45359/0.40469

# ... mean density ----
# ... ... biomass ----
# Setup data
bio_dens <- create_pred_dat(dat, biomass = log(seq(1, 2000, length.out = 50)),
                            log_dens = quantile(dat$log_dens, c(0.5))) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(biomass_nat = exp(biomass_orig) * lbac_kgha) %>% 
  mutate(dens = round(exp(log_dens), 1))

# Predict lambda
bio_dens$q05 <- bio_dens$q10 <- bio_dens$q25 <- 
  bio_dens$lambda <- 
  bio_dens$q75 <- bio_dens$q90 <- bio_dens$q95 <- NA

# Replicate rows for each model
bio_dens_oc <- bind_rows("orig" = bio_dens,
                         "oc1" = bio_dens,
                         "oc2" = bio_dens,
                         "oc3" = bio_dens,
                         .id = "oc")

for (i in 1:nrow(bio_dens_oc)) {
  cat("               \r", i)
  
  lambda <- predict_lambda(dat = bio_dens_oc, years = years,
                           i = i, 
                           b = get(paste0("b_", bio_dens_oc$oc[i])), 
                           e = NA, s = 0, RE = FALSE, mu_eta = 0)
  
  # Summarize posterior
  bio_dens_oc$q05[i] <- quantile(lambda, 0.05)
  bio_dens_oc$q10[i] <- quantile(lambda, 0.10)
  bio_dens_oc$q25[i] <- quantile(lambda, 0.25)
  bio_dens_oc$lambda[i] <- mean(lambda)
  bio_dens_oc$q75[i] <- quantile(lambda, 0.75)
  bio_dens_oc$q90[i] <- quantile(lambda, 0.90)
  bio_dens_oc$q95[i] <- quantile(lambda, 0.95)
}

# Plot
(bio_dens_plot <- ggplot(bio_dens_oc, 
                         aes(x = biomass_nat, y = lambda, color = oc)) +
    # geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
    #             alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    # geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
    #             alpha = 0.5) +
    geom_line() +
    xlab("Biomass (kg/ha)") +
    ylab(expression("Expected Density" ~ (elk/km^2))) +
    scale_color_discrete(name = "Openness Cutoff",
                         breaks = c("orig", "oc1", "oc2", "oc3"),
                         labels = c("Original", "> 30% Open", 
                                    "> 50% Open", "< 99.9% Open")) +
    theme_bw())
ggsave("fig/open_sens/density_biomass.tiff", plot = bio_dens_plot,
       width = 7, height = 4, units = "in", compression = "lzw",
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

# Replicate rows for each model
open_dens_oc <- bind_rows("orig" = open_dens,
                          "oc1" = open_dens,
                          "oc2" = open_dens,
                          "oc3" = open_dens,
                          .id = "oc")

for (i in 1:nrow(open_dens_oc)) {
  cat("               \r", i)
  
  lambda <- predict_lambda(dat = open_dens_oc, 
                           years = years,
                           i = i, 
                           b = get(paste0("b_", open_dens_oc$oc[i])), 
                           e = NA, s = 0, RE = FALSE, mu_eta = 0)
  
  # Summarize posterior
  open_dens_oc$q05[i] <- quantile(lambda, 0.05)
  open_dens_oc$q10[i] <- quantile(lambda, 0.10)
  open_dens_oc$q25[i] <- quantile(lambda, 0.25)
  open_dens_oc$lambda[i] <- mean(lambda)
  open_dens_oc$q75[i] <- quantile(lambda, 0.75)
  open_dens_oc$q90[i] <- quantile(lambda, 0.90)
  open_dens_oc$q95[i] <- quantile(lambda, 0.95)
}

# Plot
(open_dens_plot <- ggplot(open_dens_oc, 
                          aes(x = open_orig, y = lambda,
                              color = oc)) +
    # geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
    #             alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    # geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
    #             alpha = 0.5) +
    geom_line() +
    xlab("Openness (%)") +
    ylab(expression("Expected Density" ~ (elk/km^2))) +
    coord_cartesian(ylim = c(0, 20)) +
    scale_color_discrete(name = "Openness Cutoff",
                         breaks = c("orig", "oc1", "oc2", "oc3"),
                         labels = c("Original", "> 30% Open", 
                                    "> 50% Open", "< 99.9% Open")) +
    theme_bw())

ggsave("fig/open_sens/density_open.tiff", plot = open_dens_plot,
       width = 7, height = 4, units = "in", compression = "lzw",
       device = agg_tiff)

# ... ... roughness ----
# Setup data
rough_dens <- create_pred_dat(dat, 
                              rough = seq(0, 50, length.out = 50),
                              log_dens = quantile(dat$log_dens, c(0.5))) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1))

# Predict lambda
rough_dens$q05 <- rough_dens$q10 <- rough_dens$q25 <- 
  rough_dens$lambda <- 
  rough_dens$q75 <- rough_dens$q90 <- rough_dens$q95 <- NA

# Replicate rows for each model
rough_dens_oc <- bind_rows("orig" = rough_dens,
                           "oc1" = rough_dens,
                           "oc2" = rough_dens,
                           "oc3" = rough_dens,
                           .id = "oc")

for (i in 1:nrow(rough_dens_oc)) {
  cat("               \r", i)
  
  lambda <- predict_lambda(dat = rough_dens_oc, 
                           years = years,
                           i = i, 
                           b = get(paste0("b_", rough_dens_oc$oc[i])), 
                           e = NA, s = 0, RE = FALSE, mu_eta = 0)
  
  # Summarize posterior
  rough_dens_oc$q05[i] <- quantile(lambda, 0.05)
  rough_dens_oc$q10[i] <- quantile(lambda, 0.10)
  rough_dens_oc$q25[i] <- quantile(lambda, 0.25)
  rough_dens_oc$lambda[i] <- mean(lambda)
  rough_dens_oc$q75[i] <- quantile(lambda, 0.75)
  rough_dens_oc$q90[i] <- quantile(lambda, 0.90)
  rough_dens_oc$q95[i] <- quantile(lambda, 0.95)
}

# Plot
(rough_dens_plot <- ggplot(rough_dens_oc, 
                           aes(x = rough_orig, y = lambda,
                               color = oc)) +
    # geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
    #             alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    # geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
    #             alpha = 0.5) +
    geom_line() +
    xlab("Roughness (m)") +
    ylab(expression("Expected Density" ~ (elk/km^2))) +
    scale_color_discrete(name = "Openness Cutoff",
                         breaks = c("orig", "oc1", "oc2", "oc3"),
                         labels = c("Original", "> 30% Open", 
                                    "> 50% Open", "< 99.9% Open")) +
    theme_bw())

ggsave("fig/open_sens/density_rough.tiff", plot = rough_dens_plot,
       width = 7, height = 4, units = "in", compression = "lzw",
       device = agg_tiff)

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

# Replicate rows for each model
bio_ddhs_x1_oc <- bind_rows("orig" = bio_ddhs_x1,
                            "oc1" = bio_ddhs_x1,
                            "oc2" = bio_ddhs_x1,
                            "oc3" = bio_ddhs_x1,
                            .id = "oc")

bio_ddhs_x2_oc <- bind_rows("orig" = bio_ddhs_x2,
                            "oc1" = bio_ddhs_x2,
                            "oc2" = bio_ddhs_x2,
                            "oc3" = bio_ddhs_x2,
                            .id = "oc")

for (i in 1:nrow(bio_ddhs_x1_oc)) {
  cat("               \r", i)
  
  lambda_x1 <- predict_lambda(dat = bio_ddhs_x1_oc, 
                              years = years,
                              i = i, 
                              b = get(paste0("b_", bio_ddhs_x1_oc$oc[i])), 
                              e = NA, s = 0, RE = FALSE, mu_eta = 0)
  lambda_x2 <- predict_lambda(dat = bio_ddhs_x2_oc, 
                              years = years,
                              i = i, 
                              b = get(paste0("b_", bio_ddhs_x2_oc$oc[i])), 
                              e = NA, s = 0, RE = FALSE, mu_eta = 0)
  rss <- lambda_x1/lambda_x2
  
  # Summarize posterior
  bio_ddhs_x1_oc$q05[i] <- quantile(rss, 0.05)
  bio_ddhs_x1_oc$q10[i] <- quantile(rss, 0.10)
  bio_ddhs_x1_oc$q25[i] <- quantile(rss, 0.25)
  bio_ddhs_x1_oc$rss[i] <- mean(rss)
  bio_ddhs_x1_oc$q75[i] <- quantile(rss, 0.75)
  bio_ddhs_x1_oc$q90[i] <- quantile(rss, 0.90)
  bio_ddhs_x1_oc$q95[i] <- quantile(rss, 0.95)
}

# Plot
(bio_ddhs_plot <- ggplot(bio_ddhs_x1_oc, aes(x = dens, y = rss,
                                             color = oc)) +
    # geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
    #             alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    # geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
    #             alpha = 0.5) +
    geom_line() +
    xlab(expression("Average Elk Density" ~ (elk/km^2))) +
    ylab("RSS for Food (Biomass)") +
    coord_cartesian(ylim = c(1, 3)) +
    scale_color_discrete(name = "Openness Cutoff",
                         breaks = c("orig", "oc1", "oc2", "oc3"),
                         labels = c("Original", "> 30% Open", 
                                    "> 50% Open", "< 99.9% Open")) +
    theme_bw())

ggsave("fig/open_sens/ddhs_biomass.tiff", plot = bio_ddhs_plot,
       width = 7, height = 4, units = "in", compression = "lzw",
       device = agg_tiff)

# ... ... openness ----
open_ddhs_x1 <- create_pred_dat(dat, open = 1,
                                log_dens = seq(quantile(dat$log_dens, 0.05),
                                               quantile(dat$log_dens, 0.95),
                                               length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1))

open_ddhs_x2 <- create_pred_dat(dat, open = 0.7,
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

# Replicate rows for each model
open_ddhs_x1_oc <- bind_rows("orig" = open_ddhs_x1,
                             "oc1" = open_ddhs_x1,
                             "oc2" = open_ddhs_x1,
                             "oc3" = open_ddhs_x1,
                             .id = "oc")

open_ddhs_x2_oc <- bind_rows("orig" = open_ddhs_x2,
                             "oc1" = open_ddhs_x2,
                             "oc2" = open_ddhs_x2,
                             "oc3" = open_ddhs_x2,
                             .id = "oc")

for (i in 1:nrow(open_ddhs_x1_oc)) {
  cat("               \r", i)
  
  lambda_x1 <- predict_lambda(dat = open_ddhs_x1_oc, 
                              years = years,
                              i = i, 
                              b = get(paste0("b_", open_ddhs_x1_oc$oc[i])), 
                              e = NA, s = 0, RE = FALSE, mu_eta = 0)
  lambda_x2 <- predict_lambda(dat = open_ddhs_x2_oc, 
                              years = years,
                              i = i, 
                              b = get(paste0("b_", open_ddhs_x2_oc$oc[i])), 
                              e = NA, s = 0, RE = FALSE, mu_eta = 0)
  rss <- lambda_x1/lambda_x2
  
  # Summarize posterior
  open_ddhs_x1_oc$q05[i] <- quantile(rss, 0.05)
  open_ddhs_x1_oc$q10[i] <- quantile(rss, 0.10)
  open_ddhs_x1_oc$q25[i] <- quantile(rss, 0.25)
  open_ddhs_x1_oc$rss[i] <- mean(rss)
  open_ddhs_x1_oc$q75[i] <- quantile(rss, 0.75)
  open_ddhs_x1_oc$q90[i] <- quantile(rss, 0.90)
  open_ddhs_x1_oc$q95[i] <- quantile(rss, 0.95)
}

# Plot
(open_ddhs_plot <- ggplot(open_ddhs_x1_oc, 
                          aes(x = dens, y = rss, color = oc)) +
    # geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
    #             alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    # geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
    #             alpha = 0.5) +
    geom_line() +
    xlab(expression("Average Elk Density" ~ (elk/km^2))) +
    ylab("RSS for Openness") +
    scale_color_discrete(name = "Openness Cutoff",
                         breaks = c("orig", "oc1", "oc2", "oc3"),
                         labels = c("Original", "> 30% Open", 
                                    "> 50% Open", "< 99.9% Open")) +
    theme_bw())

ggsave("fig/open_sens/ddhs_open.tiff", plot = open_ddhs_plot,
       width = 7, height = 4, units = "in", compression = "lzw",
       device = agg_tiff)

# ... ... roughness ----
rough_ddhs_x1 <- create_pred_dat(dat, rough = 20,
                                 log_dens = seq(quantile(dat$log_dens, 0.05),
                                                quantile(dat$log_dens, 0.95),
                                                length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1))

rough_ddhs_x2 <- create_pred_dat(dat, rough = 0,
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

# Replicate rows for each model
rough_ddhs_x1_oc <- bind_rows("orig" = rough_ddhs_x1,
                              "oc1" = rough_ddhs_x1,
                              "oc2" = rough_ddhs_x1,
                              "oc3" = rough_ddhs_x1,
                              .id = "oc")

rough_ddhs_x2_oc <- bind_rows("orig" = rough_ddhs_x2,
                              "oc1" = rough_ddhs_x2,
                              "oc2" = rough_ddhs_x2,
                              "oc3" = rough_ddhs_x2,
                              .id = "oc")

for (i in 1:nrow(rough_ddhs_x1_oc)) {
  cat("               \r", i)
  
  lambda_x1 <- predict_lambda(dat = rough_ddhs_x1_oc, 
                              years = years,
                              i = i, 
                              b = get(paste0("b_", rough_ddhs_x1_oc$oc[i])), 
                              e = NA, s = 0, RE = FALSE, mu_eta = 0)
  lambda_x2 <- predict_lambda(dat = rough_ddhs_x2_oc, 
                              years = years,
                              i = i, 
                              b = get(paste0("b_", rough_ddhs_x2_oc$oc[i])), 
                              e = NA, s = 0, RE = FALSE, mu_eta = 0)
  rss <- lambda_x1/lambda_x2
  
  # Summarize posterior
  rough_ddhs_x1_oc$q05[i] <- quantile(rss, 0.05)
  rough_ddhs_x1_oc$q10[i] <- quantile(rss, 0.10)
  rough_ddhs_x1_oc$q25[i] <- quantile(rss, 0.25)
  rough_ddhs_x1_oc$rss[i] <- mean(rss)
  rough_ddhs_x1_oc$q75[i] <- quantile(rss, 0.75)
  rough_ddhs_x1_oc$q90[i] <- quantile(rss, 0.90)
  rough_ddhs_x1_oc$q95[i] <- quantile(rss, 0.95)
}

# Plot
(rough_ddhs_plot <- ggplot(rough_ddhs_x1_oc, 
                           aes(x = dens, y = rss, color = oc)) +
    # geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
    #             alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    # geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
    #             alpha = 0.5) +
    geom_line() +
    xlab(expression("Average Elk Density" ~ (elk/km^2))) +
    ylab("RSS for Roughness") +
    scale_color_discrete(name = "Openness Cutoff",
                         breaks = c("orig", "oc1", "oc2", "oc3"),
                         labels = c("Original", "> 30% Open", 
                                    "> 50% Open", "< 99.9% Open")) +
    theme_bw())

ggsave("fig/open_sens/ddhs_rough.tiff", plot = rough_ddhs_plot,
       width = 7, height = 4, units = "in", compression = "lzw",
       device = agg_tiff)

# ... predator density ----
# Area of *original* NR polygon clipped to inside YNP only
nr_ynp_m2 <- 994881034 #m2
nr_ynp_km2 <- nr_ynp_m2/1e6

# Wolf abundance to wolf density (wolves/100 km2)
wa_wd <- 1/(nr_ynp_km2/100)

# ... openness ----

# ... ... wolves ----
open_wolf_x1 <- create_pred_dat(dat, open = 1,
                                log_dens = quantile(dat$log_dens, 0.5),
                                wolf = seq(quantile(dat$wolf, 0.05),
                                           quantile(dat$wolf, 0.95),
                                           length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1),
         wolf_dens = wolf * wa_wd)

open_wolf_x2 <- create_pred_dat(dat, open = 0.7,
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

# Replicate rows for each model
open_wolf_x1_oc <- bind_rows("orig" = open_wolf_x1,
                             "oc1" = open_wolf_x1,
                             "oc2" = open_wolf_x1,
                             "oc3" = open_wolf_x1,
                             .id = "oc")

open_wolf_x2_oc <- bind_rows("orig" = open_wolf_x2,
                             "oc1" = open_wolf_x2,
                             "oc2" = open_wolf_x2,
                             "oc3" = open_wolf_x2,
                             .id = "oc")

for (i in 1:nrow(open_wolf_x1_oc)) {
  cat("               \r", i)
  
  lambda_x1 <- predict_lambda(dat = open_wolf_x1_oc, 
                              years = years,
                              i = i, 
                              b = get(paste0("b_", open_wolf_x1_oc$oc[i])), 
                              e = NA, s = 0, RE = FALSE, mu_eta = 0)
  lambda_x2 <- predict_lambda(dat = open_wolf_x2_oc, 
                              years = years,
                              i = i, 
                              b = get(paste0("b_", open_wolf_x2_oc$oc[i])), 
                              e = NA, s = 0, RE = FALSE, mu_eta = 0)
  rss <- lambda_x1/lambda_x2
  
  # Summarize posterior
  open_wolf_x1_oc$q05[i] <- quantile(rss, 0.05)
  open_wolf_x1_oc$q10[i] <- quantile(rss, 0.10)
  open_wolf_x1_oc$q25[i] <- quantile(rss, 0.25)
  open_wolf_x1_oc$rss[i] <- mean(rss)
  open_wolf_x1_oc$q75[i] <- quantile(rss, 0.75)
  open_wolf_x1_oc$q90[i] <- quantile(rss, 0.90)
  open_wolf_x1_oc$q95[i] <- quantile(rss, 0.95)
}


# Plot
(open_wolf_plot <- ggplot(open_wolf_x1_oc, 
                          aes(x = wolf_dens, y = rss, color = oc)) +
    # geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
    #             alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    # geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
    #             alpha = 0.5) +
    geom_line() +
    coord_cartesian(ylim = c(0, 7)) +
    xlab(expression("Wolf Density" ~ (wolves/100~km^2))) +
    ylab("RSS for Openness") +
    scale_color_discrete(name = "Openness Cutoff",
                         breaks = c("orig", "oc1", "oc2", "oc3"),
                         labels = c("Original", "> 30% Open", 
                                    "> 50% Open", "< 99.9% Open")) +
    theme_bw())

# ... ... cougars ----
open_cougar_x1 <- create_pred_dat(dat, open = 1,
                                  log_dens = quantile(dat$log_dens, 0.5),
                                  cougar = seq(quantile(dat$cougar, 0.05),
                                               quantile(dat$cougar, 0.95),
                                               length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1))

open_cougar_x2 <- create_pred_dat(dat, open = 0.7,
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

# Replicate rows for each model
open_cougar_x1_oc <- bind_rows("orig" = open_cougar_x1,
                               "oc1" = open_cougar_x1,
                               "oc2" = open_cougar_x1,
                               "oc3" = open_cougar_x1,
                               .id = "oc")

open_cougar_x2_oc <- bind_rows("orig" = open_cougar_x2,
                               "oc1" = open_cougar_x2,
                               "oc2" = open_cougar_x2,
                               "oc3" = open_cougar_x2,
                               .id = "oc")

for (i in 1:nrow(open_cougar_x1_oc)) {
  cat("               \r", i)
  
  lambda_x1 <- predict_lambda(dat = open_cougar_x1_oc, 
                              years = years,
                              i = i, 
                              b = get(paste0("b_", open_cougar_x1_oc$oc[i])), 
                              e = NA, s = 0, RE = FALSE, mu_eta = 0)
  lambda_x2 <- predict_lambda(dat = open_cougar_x2_oc, 
                              years = years,
                              i = i, 
                              b = get(paste0("b_", open_cougar_x2_oc$oc[i])), 
                              e = NA, s = 0, RE = FALSE, mu_eta = 0)
  rss <- lambda_x1/lambda_x2
  
  # Summarize posterior
  open_cougar_x1_oc$q05[i] <- quantile(rss, 0.05)
  open_cougar_x1_oc$q10[i] <- quantile(rss, 0.10)
  open_cougar_x1_oc$q25[i] <- quantile(rss, 0.25)
  open_cougar_x1_oc$rss[i] <- mean(rss)
  open_cougar_x1_oc$q75[i] <- quantile(rss, 0.75)
  open_cougar_x1_oc$q90[i] <- quantile(rss, 0.90)
  open_cougar_x1_oc$q95[i] <- quantile(rss, 0.95)
}

# Plot
(open_cougar_plot <- ggplot(open_cougar_x1_oc, 
                            aes(x = cougar, y = rss, color = oc)) +
    # geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
    #             alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    # geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
    #             alpha = 0.5) +
    geom_line() +
    coord_cartesian(ylim = c(0, 7)) +
    xlab(expression("Cougar Density" ~ (cougars/100~km^2))) +
    # ylab("RSS for Openness") +
    ylab(NULL) +
    scale_color_discrete(name = "Openness Cutoff",
                         breaks = c("orig", "oc1", "oc2", "oc3"),
                         labels = c("Original", "> 30% Open", 
                                    "> 50% Open", "< 99.9% Open")) +
    theme_bw())

# ... ... combine ----
open_pred_plot <- open_wolf_plot + open_cougar_plot + 
  plot_layout(guides = "collect")

ggsave("fig/open_sens/predator_open.tiff", plot = open_pred_plot,
       width = 8, height = 3, units = "in", compression = "lzw",
       device = agg_tiff)

# ... roughness ----
# ... ... wolves ----
rough_wolf_x1 <- create_pred_dat(dat, rough = 20,
                                 log_dens = quantile(dat$log_dens, 0.5),
                                 wolf = seq(quantile(dat$wolf, 0.05),
                                            quantile(dat$wolf, 0.95),
                                            length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1),
         wolf_dens = wolf * wa_wd)

rough_wolf_x2 <- create_pred_dat(dat, rough = 0,
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

# Replicate rows for each model
rough_wolf_x1_oc <- bind_rows("orig" = rough_wolf_x1,
                              "oc1" = rough_wolf_x1,
                              "oc2" = rough_wolf_x1,
                              "oc3" = rough_wolf_x1,
                              .id = "oc")

rough_wolf_x2_oc <- bind_rows("orig" = rough_wolf_x2,
                              "oc1" = rough_wolf_x2,
                              "oc2" = rough_wolf_x2,
                              "oc3" = rough_wolf_x2,
                              .id = "oc")

for (i in 1:nrow(rough_wolf_x1_oc)) {
  cat("               \r", i)
  
  lambda_x1 <- predict_lambda(dat = rough_wolf_x1_oc, 
                              years = years,
                              i = i, 
                              b = get(paste0("b_", rough_wolf_x1_oc$oc[i])), 
                              e = NA, s = 0, RE = FALSE, mu_eta = 0)
  lambda_x2 <- predict_lambda(dat = rough_wolf_x2_oc, 
                              years = years,
                              i = i, 
                              b = get(paste0("b_", rough_wolf_x2_oc$oc[i])), 
                              e = NA, s = 0, RE = FALSE, mu_eta = 0)
  rss <- lambda_x1/lambda_x2
  
  # Summarize posterior
  rough_wolf_x1_oc$q05[i] <- quantile(rss, 0.05)
  rough_wolf_x1_oc$q10[i] <- quantile(rss, 0.10)
  rough_wolf_x1_oc$q25[i] <- quantile(rss, 0.25)
  rough_wolf_x1_oc$rss[i] <- mean(rss)
  rough_wolf_x1_oc$q75[i] <- quantile(rss, 0.75)
  rough_wolf_x1_oc$q90[i] <- quantile(rss, 0.90)
  rough_wolf_x1_oc$q95[i] <- quantile(rss, 0.95)
}

# Plot
(rough_wolf_plot <- ggplot(rough_wolf_x1_oc, 
                           aes(x = wolf_dens, y = rss, color = oc)) +
    # geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
    #             alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    # geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
    #             alpha = 0.5) +
    geom_line() +
    coord_cartesian(ylim = c(0, 20)) +
    xlab(expression("Wolf Density" ~ (wolves/100~km^2))) +
    ylab("RSS for Roughness") +
    scale_color_discrete(name = "Openness Cutoff",
                         breaks = c("orig", "oc1", "oc2", "oc3"),
                         labels = c("Original", "> 30% Open", 
                                    "> 50% Open", "< 99.9% Open")) +
    theme_bw())

# ... ... cougars ----
rough_cougar_x1 <- create_pred_dat(dat, rough = 20,
                                   log_dens = quantile(dat$log_dens, 0.5),
                                   cougar = seq(quantile(dat$cougar, 0.05),
                                                quantile(dat$cougar, 0.95),
                                                length.out = 50)) %>%
  intxn() %>% 
  scale_dat(scale_df = scale_df) %>% 
  mutate(dens = round(exp(log_dens), 1))

rough_cougar_x2 <- create_pred_dat(dat, rough = 0,
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

# Replicate rows for each model
rough_cougar_x1_oc <- bind_rows("orig" = rough_cougar_x1,
                                "oc1" = rough_cougar_x1,
                                "oc2" = rough_cougar_x1,
                                "oc3" = rough_cougar_x1,
                                .id = "oc")

rough_cougar_x2_oc <- bind_rows("orig" = rough_cougar_x2,
                                "oc1" = rough_cougar_x2,
                                "oc2" = rough_cougar_x2,
                                "oc3" = rough_cougar_x2,
                                .id = "oc")

for (i in 1:nrow(rough_cougar_x1_oc)) {
  cat("               \r", i)
  
  lambda_x1 <- predict_lambda(dat = rough_cougar_x1_oc, 
                              years = years,
                              i = i, 
                              b = get(paste0("b_", rough_cougar_x1_oc$oc[i])), 
                              e = NA, s = 0, RE = FALSE, mu_eta = 0)
  lambda_x2 <- predict_lambda(dat = rough_cougar_x2_oc, 
                              years = years,
                              i = i, 
                              b = get(paste0("b_", rough_cougar_x2_oc$oc[i])), 
                              e = NA, s = 0, RE = FALSE, mu_eta = 0)
  rss <- lambda_x1/lambda_x2
  
  # Summarize posterior
  rough_cougar_x1_oc$q05[i] <- quantile(rss, 0.05)
  rough_cougar_x1_oc$q10[i] <- quantile(rss, 0.10)
  rough_cougar_x1_oc$q25[i] <- quantile(rss, 0.25)
  rough_cougar_x1_oc$rss[i] <- mean(rss)
  rough_cougar_x1_oc$q75[i] <- quantile(rss, 0.75)
  rough_cougar_x1_oc$q90[i] <- quantile(rss, 0.90)
  rough_cougar_x1_oc$q95[i] <- quantile(rss, 0.95)
}

# Plot
(rough_cougar_plot <- ggplot(rough_cougar_x1_oc, 
                             aes(x = cougar, y = rss, color = oc)) +
    # geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
    #             alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    # geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
    #             alpha = 0.5) +
    geom_line() +
    coord_cartesian(ylim = c(0, 20)) +
    xlab(expression("Cougar Density" ~ (cougars/100~km^2))) +
    # ylab("RSS for Roughness") +
    ylab(NULL) +
    scale_color_discrete(name = "Openness Cutoff",
                         breaks = c("orig", "oc1", "oc2", "oc3"),
                         labels = c("Original", "> 30% Open", 
                                    "> 50% Open", "< 99.9% Open")) +
    theme_bw())

# ... ... combine ----
rough_pred_plot <- rough_wolf_plot + rough_cougar_plot + 
  plot_layout(guides = "collect")

ggsave("fig/open_sens/predator_rough.tiff", plot = rough_pred_plot,
       width = 8, height = 3, units = "in", compression = "lzw",
       device = agg_tiff)

# Main result figure ----

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

# Replicate rows for each model
safe_ddhs_x1_oc <- bind_rows("orig" = safe_ddhs_x1,
                             "oc1" = safe_ddhs_x1,
                             "oc2" = safe_ddhs_x1,
                             "oc3" = safe_ddhs_x1,
                             .id = "oc")

safe_ddhs_x2_oc <- bind_rows("orig" = safe_ddhs_x2,
                             "oc1" = safe_ddhs_x2,
                             "oc2" = safe_ddhs_x2,
                             "oc3" = safe_ddhs_x2,
                             .id = "oc")

for (i in 1:nrow(safe_ddhs_x1_oc)) {
  cat("               \r", i)
  
  lambda_x1 <- predict_lambda(dat = safe_ddhs_x1_oc, 
                              years = years,
                              i = i, 
                              b = get(paste0("b_", safe_ddhs_x1_oc$oc[i])), 
                              e = NA, s = 0, RE = FALSE, mu_eta = 0)
  lambda_x2 <- predict_lambda(dat = safe_ddhs_x2_oc, 
                              years = years,
                              i = i, 
                              b = get(paste0("b_", safe_ddhs_x2_oc$oc[i])), 
                              e = NA, s = 0, RE = FALSE, mu_eta = 0)
  rss <- lambda_x1/lambda_x2
  
  # Summarize posterior
  safe_ddhs_x1_oc$q05[i] <- quantile(rss, 0.05)
  safe_ddhs_x1_oc$q10[i] <- quantile(rss, 0.10)
  safe_ddhs_x1_oc$q25[i] <- quantile(rss, 0.25)
  safe_ddhs_x1_oc$rss[i] <- mean(rss)
  safe_ddhs_x1_oc$q75[i] <- quantile(rss, 0.75)
  safe_ddhs_x1_oc$q90[i] <- quantile(rss, 0.90)
  safe_ddhs_x1_oc$q95[i] <- quantile(rss, 0.95)
}

# Plot
(safe_ddhs_plot <- ggplot(safe_ddhs_x1_oc, 
                          aes(x = dens, y = rss, 
                              color = oc)) +
    # geom_ribbon(aes(ymin = q05, ymax = q95), fill = color_90,
    #             alpha = 0.5) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = color_80,
                alpha = 0.5) +
    # geom_ribbon(aes(ymin = q25, ymax = q75), fill = color_50,
    #             alpha = 0.5) +
    geom_line() +
    coord_cartesian(ylim = c(1, 3)) +
    xlab(expression("Elk Density" ~ (elk/km^2))) +
    ylab("RSS for Safety") +
    scale_color_discrete(name = "Openness Cutoff",
                         breaks = c("orig", "oc1", "oc2", "oc3"),
                         labels = c("Original", "> 30% Open", 
                                    "> 50% Open", "< 99.9% Open")) +
    theme_bw())

ggsave("fig/open_sens/ddhs_safety.tiff", plot = safe_ddhs_plot,
       width = 6, height = 3, units = "in", compression = "lzw",
       device = agg_tiff)

# Figure S5 ----
# Split calculations with OCs together into 30% and 50%
food_safe_30 <- bind_rows(bio = bio_ddhs_x1_oc,
                          safe = safe_ddhs_x1_oc,
                          .id = "driver") %>% 
  filter(oc == "oc1") %>% 
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

food_safe_50 <- bind_rows(bio = bio_ddhs_x1_oc,
                          safe = safe_ddhs_x1_oc,
                          .id = "driver") %>% 
  filter(oc == "oc2") %>% 
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

food_safe_99 <- bind_rows(bio = bio_ddhs_x1_oc,
                          safe = safe_ddhs_x1_oc,
                          .id = "driver") %>% 
  filter(oc == "oc3") %>% 
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

# ... S5A ----
S5A <- open_hist

# ... S5B ----
(S5B <- food_safe_30 %>% 
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
    coord_cartesian(ylim = c(1.2, 2.8)) +
    ggtitle("> 30% Open") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          plot.title = element_text(size = 8)))

# ... S5C ----
(S5C <- food_safe_50 %>% 
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
    xlab("Average Elk Density") +
    ylab("RSS") +
    coord_cartesian(ylim = c(1.2, 2.8)) +
    ggtitle("> 50% Open") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 10, hjust = 1),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(size = 8)))

# ... S5D ----
(S5D <- food_safe_99 %>% 
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
    xlab(expression((elk/km^2))) +
    ylab("RSS") +
    coord_cartesian(ylim = c(1.2, 2.8)) +
    ggtitle("< 99.9% Open") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 10, hjust = -0.75),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(size = 8)))

# ... combine ----
figS5 <- S5A / (S5B + S5C + S5D) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect")

ggsave("fig/ms/figS5.tif", plot = figS5, device = agg_tiff,
       width = 173, height = 130, units = "mm", dpi = 500,
       compression = "lzw")

