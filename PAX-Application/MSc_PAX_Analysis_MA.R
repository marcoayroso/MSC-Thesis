## MSc Thesis - R Code for PAX Application
## Author: Marco Ayroso
## Email: ayrosom@myumanitoba.ca


#____________________________________________
# Load libraries
#____________________________________________


library(tidyverse)
library(readxl)
# library(lme4)
# library(lmerTest)
library(mnormt)
library(compiler)
library(simstudy)
library(faux)
library(afex)
library(broom.mixed)
library(writexl)
library(psych) #for description
# library(nlme) #for multielvel models
library(brms) #for bayesian multilevel models
library(jtools)
library(glmmTMB)
library(tictoc)
library(cmdstanr)
# library(gamlss)
# library(gamlss.mx)
library(future)
library(dplyr)
library(future)
library(progress)
library(progressr)
handlers(global = TRUE)
library(ggpubr)

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/MSc Thesis/Thesis/PAX_G5")
PAX_dest <- ("~/Library/Mobile Documents/com~apple~CloudDocs/MSc Thesis/Thesis/PAX_G5")

# set up data

master <- read_csv(paste0(PAX_dest,"/","PAX_G5_Marco.csv"))
head(master)
colnames(master)

control <- master %>%
  filter(case_sch == 0)

treatment <- master %>%
  filter(case_sch == 1)

master <- master %>%
  mutate(`0` = rowSums(across(c(emotT1, condT1, hyprT1, peerT1))),
         `1` = rowSums(across(c(emotT2, condT2, hyprT2, peerT2))),
         `2` = rowSums(across(c(emotT3, condT3, hyprT3, peerT3)))) %>%
  pivot_longer(cols = c(`0`, `1`, `2`),
               names_to = "time",
               values_to = "SDQ") %>%
  mutate(
    time = as.numeric(time)
  ) %>%
  select(schcode, case_sch, STID, Female, Abst, cmci, efal_new, S_sefi2, "time", SDQ)

master$case_sch <- as.factor(master$case_sch)
master$Female <- as.factor(master$Female)

# Parallel chains
future::plan(multisession, workers = 7)

# Analysis

# MLM (Model A)

tic("MLM")

# Priors
priors_melsm <- c(
  prior(normal(0, 3), class = "b"),
  prior(normal(20, 5), class = "Intercept"),
  prior(normal(0, 0.2), class = "sd")
  # prior(normal(0, 0.3), class = "b", dpar = "sigma"),
  # prior(normal(log(5), 0.3), class = "Intercept", dpar = "sigma"),
  # prior(normal(0, 0.2), class = "sd", dpar = "sigma")
)

# Fit empty MELSM with brms
fit <- function(data) {
  data <- data %>%
    mutate(
      schcode = factor(schcode),
      STID = factor(STID)
    )
  
  brm(
    bf(
      SDQ ~ time * case_sch + (1 + time | schcode/STID)
      # center = TRUE
    ),
    data = data,
    prior = priors_melsm,
    family = gaussian(),
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    threads = threading(2, grainsize = 500),
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    backend = "cmdstanr",
    stan_model_args = list(stanc_options = list(O1 = TRUE)),
    refresh = 1
  )
}

fit_melsm <- fit(master)
toc()

summary(fit_melsm)
save(fit_melsm, file = "UPDATE_brms_mlm_fit_4chains_2000iter.RData")

# Empty MELSM (Model B)

tic("Empty MELSM")

# Priors
priors_melsm <- c(
  prior(normal(0, 3), class = "b"),
  prior(normal(20, 5), class = "Intercept"),
  prior(normal(0, 0.2), class = "sd")
  # prior(normal(0, 0.3), class = "b", dpar = "sigma"),
  # prior(normal(log(5), 0.3), class = "Intercept", dpar = "sigma"),
  # prior(normal(0, 0.2), class = "sd", dpar = "sigma")
)

# Fit empty MELSM with brms
fit <- function(data) {
  data <- data %>%
    mutate(
      schcode = factor(schcode),
      STID = factor(STID)
    )
  
  brm(
    bf(
      SDQ ~ time * case_sch + (1 + time | schcode/STID),
      sigma ~ time * case_sch + (1 + time | schcode/STID)
      # center = TRUE
    ),
    data = data,
    prior = priors_melsm,
    family = gaussian(),
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    threads = threading(2, grainsize = 500),
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    backend = "cmdstanr",
    stan_model_args = list(stanc_options = list(O1 = TRUE)),
    refresh = 1
  )
}

fit_melsm <- fit(master)
toc()

summary(fit_melsm)
save(fit_melsm, file = "brms_e_melsm_fit_4chains_2000iter.RData")


# Full MELSM (Model C)

tic("Full MELSM")

# Priors
priors_melsm <- c(
  prior(normal(0, 3), class = "b"),
  prior(normal(20, 5), class = "Intercept"),
  prior(normal(0, 0.2), class = "sd"),
  prior(normal(0, 0.3), class = "b", dpar = "sigma"),
  prior(normal(log(5), 0.3), class = "Intercept", dpar = "sigma"),
  prior(normal(0, 0.2), class = "sd", dpar = "sigma")
)


# Fit MELSM with brms
fit <- function(data) {
  data <- data %>%
    mutate(
      schcode = factor(schcode),
      STID = factor(STID)
    )
  
  brm(
    bf(
      SDQ ~ time * case_sch + S_sefi2 + Female + (1 + time | schcode/STID),
      sigma ~ time * case_sch + S_sefi2 + Female + (1 + time | schcode/STID)
      # center = TRUE
    ),
    data = data,
    prior = priors_melsm,
    family = gaussian(),
    chains = 4,
    iter = 2000,      
    warmup = 1000, 
    cores = 4,  
    threads = threading(2, grainsize = 500),
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    backend = "cmdstanr",
    stan_model_args = list(stanc_options = list(O1 = TRUE)),
    refresh = 1
  )
}

fit_melsm <- fit(master)
toc()

summary(fit_melsm)
save(fit_melsm, file = "UPDATE_brms_melsm_fit_4chains_2000iter_v2.RData")

# code for moderator effect below

# tic("Full MELSM")
# # Fit MELSM with brms
# fit <- function(data) {
#   data <- data %>%
#     mutate(
#       schcode = factor(schcode),
#       STID = factor(STID)
#     )
#   
#   brm(
#     bf(
#       SDQ ~ time * case_sch + S_sefi2 + case_sch:S_sefi2 + time:case_sch:S_sefi2 + (1 + time | schcode/STID),
#       sigma ~ time * case_sch + S_sefi2 + case_sch:S_sefi2 + time:case_sch:S_sefi2 + (1 + time | schcode/STID)
#       # center = TRUE
#     ),
#     data = data,
#     prior = priors_melsm,
#     family = gaussian(),
#     chains = 4,
#     iter = 2000,      
#     warmup = 1000, 
#     cores = 4,  
#     threads = threading(2, grainsize = 500),
#     control = list(adapt_delta = 0.99, max_treedepth = 15),
#     backend = "cmdstanr",
#     stan_model_args = list(stanc_options = list(O1 = TRUE)),
#     refresh = 1
#   )
# }
# 
# fit_melsm <- fit(master)
# toc()
# 
# summary(fit_melsm)
# save(fit_melsm, file = "UPDATE_brms_melsm_fit_4chains_2000iter_SSefi2.RData")
# 
# 
# tic("Full MELSM")
# 
# # Priors
# priors_melsm <- c(
#   prior(normal(0, 3), class = "b"),
#   prior(normal(20, 5), class = "Intercept"),
#   prior(normal(0, 0.2), class = "sd"),
#   prior(normal(0, 0.3), class = "b", dpar = "sigma"),
#   prior(normal(log(5), 0.3), class = "Intercept", dpar = "sigma"),
#   prior(normal(0, 0.2), class = "sd", dpar = "sigma")
# )
# 
# 
# # Fit MELSM with brms
# fit <- function(data) {
#   data <- data %>%
#     mutate(
#       schcode = factor(schcode),
#       STID = factor(STID)
#     )
#   
#   brm(
#     bf(
#       SDQ ~ time * case_sch + Female + case_sch:Female + time:case_sch:Female + (1 + time | schcode/STID),
#       sigma ~ time * case_sch + Female + case_sch:Female + time:case_sch:Female + (1 + time | schcode/STID)
#       # center = TRUE
#     ),
#     data = data,
#     prior = priors_melsm,
#     family = gaussian(),
#     chains = 4,
#     iter = 2000,      
#     warmup = 1000, 
#     cores = 4,  
#     threads = threading(2, grainsize = 500),
#     control = list(adapt_delta = 0.99, max_treedepth = 15),
#     backend = "cmdstanr",
#     stan_model_args = list(stanc_options = list(O1 = TRUE)),
#     refresh = 1
#   )
# }
# 
# fit_melsm_Female <- fit(master)
# toc()
# 
# summary(fit_melsm_Female)
# save(fit_melsm_Female, file = "UPDATE_brms_melsm_fit_4chains_2000iter_Female.RData")

# ========================================
# Model A diagnostics
# ========================================

load("UPDATE_brms_melsm_fit_4chains_2000iter_v2.RData")
fit_C <- fit_melsm
load("brms_e_melsm_fit_4chains_2000iter.RData")
fit_B <- fit_melsm
load("UPDATE_brms_mlm_fit_4chains_2000iter.RData")
fit_A <- fit_melsm

models <- list(A=fit_A, B=fit_B, C=fit_C)

# Largest RHat
round(max(summary(fit_A)$fixed[,"Rhat"], na.rm=T), 3)

# 1. PPC
ppA <- pp_check(fit_A, nsamples = 15L, type = "dens_overlay") +
  xlim(-40, 40) +
  ylim(0, 0.10) +
  labs(title = "Model A", x = "", y = "Density") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "none")
ppA

# 2. LOO
loo_A <- loo(fit_A, save_psis=F)
print(loo_A)
readline()

waic_A <- waic(fit_A)

# 3. Core: Residuals vs Fitted
mu_A <- colMeans(posterior_epred(fit_A))
resids_A <- fit_A$data$SDQ - mu_A
plot(mu_A, resids_A, pch=16, cex=0.7, las=1, bty="l", main="A: Resids vs Fitted",
     xlab="Fitted (μ)", ylab="Resids")
abline(h=0,col="blue",lwd=2)
grid()

mu_A     <- colMeans(posterior_epred(fit_A))
resids_A <- fit_A$data$SDQ - mu_A
df_A     <- data.frame(fitted = mu_A, resid = resids_A)

rpA <- ggplot(df_A, aes(x = fitted, y = resid)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "blue", linewidth = 1) +
  xlim(-5, 30) +
  ylim(-20, 30) +
  labs(title = "Model A",
       x = "",
       y = "Pearson Residuals") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold"))

# 4. QQ Plot
qq_th <- qnorm(ppoints(length(resids_A)))
qq_res <- quantile(resids_A, ppoints(length(resids_A)))
plot(qq_th, qq_res, pch=16, cex=0.7, las=1, bty="l", main="A: QQ Resids",
     xlab="", ylab="")
abline(0,1,col="blue",lwd=2)
grid()

qq_th  <- qnorm(ppoints(length(resids_A)))
qq_res <- quantile(resids_A, ppoints(length(resids_A)))
qq_df  <- data.frame(theoretical = qq_th,
                     sample      = qq_res)

qqA <- ggplot(qq_df, aes(x = theoretical, y = sample)) +
  geom_point(size = 1.2) +
  geom_abline(intercept = 0, slope = 1, colour = "blue", linewidth = 1) +
  ylim(-20, 30) +
  labs(title = "Model A",
       x = "",
       y = "Sample Quantiles") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold"))


# 5. Index Plot (resids vs observation index)
plot(1:length(resids_A), resids_A, pch=16, cex=0.6, las=1, bty="l", main="A: Index Plot",
     xlab="Observation Index", ylab="Resids")
abline(h=0,col="blue",lwd=2)
grid()

df_index_A <- data.frame(index = 1:length(resids_A), resid = resids_A)

ipA <- ggplot(df_index_A, aes(x = index, y = resid)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "blue", linewidth = 1) +
  scale_x_continuous(breaks = seq(0, 6000, 2000)) +
  ylim(-20, 30) +
  labs(title = "Model A",
       x = "",
       y = "Pearson Residuals") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold"))


# 6. Residuals Histogram
hist(resids_A, breaks=30, main="A: Resids Histogram", xlab="Resids")
abline(v=0,col="blue",lwd=2)

hgA <- ggplot(data.frame(resids = resids_A), aes(x = resids)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, 
                 fill = "grey", color = "black", alpha = 0.7) +
  geom_density(color = "blue", linewidth = 1) +
  ylim(0, 0.25) +
  labs(title = "Model A",
       x = "",
       y = "Density") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold"))

# 7. Linearity (resids vs time)
plot(fit_A$data$time, resids_A, pch=16, cex=0.7, las=1, bty="l", main="A: Resids vs Time",
     xlab="Time", ylab="Resids")
abline(h=0,col="blue",lwd=2)
grid()

# ========================================
# Model B diagnostics
# ========================================

# Largest RHat
round(max(summary(fit_B)$fixed[,"Rhat"], na.rm=T), 3)

# 1. PPC
ppB <- pp_check(fit_B, nsamples = 15L, type = "dens_overlay") + 
  xlim(-40, 40) +
  ylim(0, 0.10) +
  labs(title = "Model B", x = "SDQ Score", y = "") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "none")
ppB

# 2. LOO
loo_B <- loo(fit_B, save_psis=F)
print(loo_B)
readline()

waic_B <- waic(fit_B)


# 3. Core: Residuals vs Fitted
mu_B <- colMeans(posterior_epred(fit_B))
resids_B <- fit_B$data$SDQ - mu_B
plot(mu_B, resids_B, pch=16, cex=0.7, las=1, bty="l", main="B: Resids vs Fitted",
     xlab="Fitted (μ)", ylab="Resids")
abline(h=0,col="blue",lwd=2)
grid()

mu_B     <- colMeans(posterior_epred(fit_B))
resids_B <- fit_A$data$SDQ - mu_B
df_B     <- data.frame(fitted = mu_B, resid = resids_B)

rpB <- ggplot(df_B, aes(x = fitted, y = resid)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "blue", linewidth = 1) +
  xlim(-5, 30) +
  ylim(-20, 30) +
  labs(title = "Model B",
       x = "Predicted Values",
       y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold"))

# 4. QQ Plot
qq_th <- qnorm(ppoints(length(resids_B)))
qq_res <- quantile(resids_B, ppoints(length(resids_B)))
plot(qq_th, qq_res, pch=16, cex=0.7, las=1, bty="l", main="B: QQ Resids",
     xlab="Theoretical", ylab="Sample")
abline(0,1,col="blue",lwd=2)
grid()

qq_th  <- qnorm(ppoints(length(resids_B)))
qq_res <- quantile(resids_B, ppoints(length(resids_B)))
qq_df  <- data.frame(theoretical = qq_th,
                     sample      = qq_res)

qqB <- ggplot(qq_df, aes(x = theoretical, y = sample)) +
  geom_point(size = 1.2) +
  geom_abline(intercept = 0, slope = 1, colour = "blue", linewidth = 1) +
  ylim(-15, 30) +
  labs(title = "Model B",
       x = "Theoretical Quantiles",
       y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold"))


# 5. Index Plot (resids vs observation index)
plot(1:length(resids_B), resids_B, pch=16, cex=0.6, las=1, bty="l", main="B: Index Plot",
     xlab="Observation Index", ylab="Resids")
abline(h=0,col="blue",lwd=2)
grid()

df_index_B <- data.frame(index = 1:length(resids_B), resid = resids_B)

ipB <- ggplot(df_index_B, aes(x = index, y = resid)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "blue", linewidth = 1) +
  scale_x_continuous(breaks = seq(0, 6000, 2000)) +
  ylim(-20, 30) +
  labs(title = "Model B",
       x = "Observation Number",
       y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold"))


# 6. Residuals Histogram
hist(resids_B, breaks=30, main="B: Resids Histogram", xlab="Resids")
abline(v=0,col="blue",lwd=2)

hgB <- ggplot(data.frame(resids = resids_B), aes(x = resids)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, 
                 fill = "grey", color = "black", alpha = 0.7) +
  geom_density(color = "blue", linewidth = 1) +
  ylim(0, 0.25) +
  labs(title = "Model B",
       x = "Pearson Residuals",
       y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold"))


# 7. Linearity (resids vs time)
plot(fit_B$data$time, resids_B, pch=16, cex=0.7, las=1, bty="l", main="B: Resids vs Time",
     xlab="Time", ylab="Resids")
abline(h=0,col="blue",lwd=2)
grid()

# ========================================
# Model C diagnostics and comparison
# ========================================

# Largest RHat
round(max(summary(fit_C)$fixed[,"Rhat"], na.rm=T), 3)

# 1. PPC
ppC <- pp_check(fit_C, nsamples = 15L, type = "dens_overlay") + 
  xlim(-40, 40) +
  ylim(0, 0.10) +
  labs(title = "Model C", x = "", y = "") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "none")
ppC

# 2. LOO
loo_C <- loo(fit_C, save_psis=F)
print(loo_C)
readline()

waic_C <- waic(fit_C)


# 3. Core: Residuals vs Fitted
mu_C <- colMeans(posterior_epred(fit_C))
resids_C <- fit_C$data$SDQ - mu_C
plot(mu_C, resids_C, pch=16, cex=0.7, las=1, bty="l", main="C: Resids vs Fitted",
     xlab="Fitted (μ)", ylab="Resids")
abline(h=0,col="blue",lwd=2)
grid()

mu_C     <- colMeans(posterior_epred(fit_C))
resids_C <- fit_C$data$SDQ - mu_C
df_C     <- data.frame(fitted = mu_C, resid = resids_C)

rpC <- ggplot(df_C, aes(x = fitted, y = resid)) +
  geom_point(size = 1.2, alpha = 0.6) +
  xlim(-5, 30) +
  ylim(-20, 30) +
  geom_hline(yintercept = 0, colour = "blue", linewidth = 1) +
  labs(title = "Model C",
       x = "",
       y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold"))

# 4. QQ Plot
qq_th <- qnorm(ppoints(length(resids_C)))
qq_res <- quantile(resids_C, ppoints(length(resids_C)))
plot(qq_th, qq_res, pch=16, cex=0.7, las=1, bty="l", main="C: QQ Resids",
     xlab="Theoretical", ylab="Sample")
abline(0,1,col="blue",lwd=2)
grid()

qq_th  <- qnorm(ppoints(length(resids_C)))
qq_res <- quantile(resids_C, ppoints(length(resids_C)))
qq_df  <- data.frame(theoretical = qq_th,
                     sample      = qq_res)

qqC <- ggplot(qq_df, aes(x = theoretical, y = sample)) +
  geom_point(size = 1.2) +
  geom_abline(intercept = 0, slope = 1, colour = "blue", linewidth = 1) +
  ylim(-15, 30) +
  labs(title = "Model C",
       x = "",
       y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold"))


# 5. Index Plot (resids vs observation index)
plot(1:length(resids_C), resids_C, pch=16, cex=0.6, las=1, bty="l", main="C: Index Plot",
     xlab="Observation Index", ylab="Resids")
abline(h=0,col="blue",lwd=2)
grid()

df_index_C <- data.frame(index = 1:length(resids_C), resid = resids_C)

ipC <- ggplot(df_index_C, aes(x = index, y = resid)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "blue", linewidth = 1) +
  scale_x_continuous(breaks = seq(0, 6000, 2000)) +
  ylim(-20, 30) +
  labs(title = "Model C",
       x = "",
       y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold"))



# 6. Residuals Histogram
hist(resids_C, breaks=30, main="C: Resids Histogram", xlab="Resids")
abline(v=0,col="blue",lwd=2)

hgC <- ggplot(data.frame(resids = resids_C), aes(x = resids)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, 
                 fill = "grey", color = "black", alpha = 0.7) +
  geom_density(color = "blue", linewidth = 1) +
  ylim(0, 0.25) +
  labs(title = "Model C",
       x = "",
       y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold"))

# 7. Linearity (resids vs time)
plot(fit_C$data$time, resids_C, pch=16, cex=0.7, las=1, bty="l", main="C: Resids vs Time",
     xlab="Time", ylab="Resids")
abline(h=0,col="blue",lwd=2)
grid()

# 8. Linearity (resids vs SES)
plot(fit_C$data$S_sefi2, resids_C, pch=16, cex=0.7, las=1, bty="l", main="",
     xlab="SEFI2", ylab="Pearson Residuals")
abline(h=0,col="blue",lwd=2)
grid()

# Level-1: Wave residuals vs SES (your current)
df_lin1 <- data.frame(
  SEFI2 = fit_C$data$S_sefi2,
  resid = resids_C  # Level-1 Pearson residuals
)

p_level1 <- ggplot(df_lin1, aes(x = SEFI2, y = resid)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "blue", linewidth = 1) +
  labs(title = "Level-1",
       x = "SEFI2", y = "Wave Residuals") +
  ylim(-15, 30) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        plot.title = element_text(face = "bold"))

# Level-2: STUDENT residuals vs SCHOOL SES
student_ses <- df_resids %>%  # From extraction code
  group_by(STID) %>%
  summarise(
    SES = mean(S_sefi2, na.rm = TRUE),  # School SES
    level2_resid = mean(level2_student, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(level2_resid), !is.na(SES))

p_level2 <- ggplot(student_ses, aes(x = SES, y = level2_resid)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "blue", linewidth = 1) +
  labs(title = "Level-2",
       x = "SEFI2", y = "Student Residuals") +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        plot.title = element_text(face = "bold"))

SES_fig <- ggarrange(p_level1, p_level2,
                    ncol = 2, nrow = 1,
                    common.legend = TRUE, legend = "bottom")
SES_fig

# arrange

pp_fig <- ggarrange(ppA, ppB, ppC,
                    ncol = 3, nrow = 1,
                    common.legend = TRUE, legend = "bottom")
pp_fig

qq_fig <- ggarrange(qqA, qqB, qqC,
                    ncol = 3, nrow = 1,
                    common.legend = TRUE, legend = "bottom")
qq_fig

rp_fig <- ggarrange(rpA, rpB, rpC,
                    ncol = 3, nrow = 1,
                    common.legend = TRUE, legend = "bottom")
rp_fig

ip_fig <- ggarrange(ipA, ipB, ipC,
                    ncol = 3, nrow = 1,
                    common.legend = TRUE, legend = "bottom")
ip_fig

hg_fig <- ggarrange(hgA, hgB, hgC,
                    ncol = 3, nrow = 1,
                    common.legend = TRUE, legend = "bottom")
hg_fig


# posterior predictions (location and scale)
mu_post <- posterior_epred(fit_melsm, ndraws = 1)
sigma_post <- posterior_epred(fit_melsm, ndraws = 1, dpar = "sigma")

# add unique ID matching model data order
master$id <- 1:nrow(master)  # ensure matches posterior dims

# compute residuals with NA handling
df_resids <- master %>%
  mutate(
    # Level-1 raw (location)
    level1_raw = SDQ - mu_post[id],
    # Level-1 pearson (scale adjusted) 
    level1_pearson = (SDQ - mu_post[id]) / sigma_post[id]
  ) %>%
  filter(!is.na(SDQ), !is.na(id)) %>%
  group_by(STID) %>%
  mutate(
    # student mean deviation
    level2_student = mean(level1_pearson, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  group_by(schcode) %>%
  mutate(
    # school mean deviation  
    level3_school = mean(level1_pearson, na.rm = TRUE)
  ) %>%
  ungroup()

# extract unique residuals
level2_resids <- df_resids %>%
  filter(!is.na(STID), !is.na(level2_student)) %>%
  distinct(STID, level2_student) %>%
  pull(level2_student)

level3_resids <- df_resids %>%
  filter(!is.na(schcode), !is.na(level3_school)) %>%
  group_by(schcode) %>%
  summarise(mean_resid = mean(level3_school, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(mean_resid)) %>%
  pull(mean_resid)

# QQ Plots
qq_level2 <- data.frame(
  theoretical = qnorm(ppoints(length(level2_resids))),
  sample = quantile(level2_resids, ppoints(length(level2_resids)), na.rm = TRUE)
)

qq_level3 <- data.frame(
  theoretical = qnorm(ppoints(length(level3_resids))),
  sample = quantile(level3_resids, ppoints(length(level3_resids)), na.rm = TRUE)
)

# plots
p2 <- ggplot(qq_level2, aes(x = theoretical, y = sample)) +
  geom_point(size = 1.2, alpha = 0.7) +
  geom_abline(colour = "blue", linewidth = 1) +
  labs(title = "Level-2",
       x = "Sample Quantiles",
       y = "Theoretical Quantiles") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold"))


p3 <- ggplot(qq_level3, aes(x = theoretical, y = sample)) +
  geom_point(size = 1.2, alpha = 0.7) +
  geom_abline(colour = "blue", linewidth = 1) +
  ylim(-5, 5) +
  labs(title = "Level-3",
       x = "Sample Quantiles",
       y = "Theoretical Quantiles") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold"))


print(p2)
print(p3)

qq23_fig <- ggarrange(p2, p3,
                    ncol = 2, nrow = 1,
                    common.legend = TRUE, legend = "bottom")
qq23_fig


# student/school predicted means
student_fitted <- df_resids %>%
  group_by(STID) %>%
  summarise(fitted_student = mean(mu_post[id], na.rm = TRUE),
            level2_resid = first(level2_student),  # unique per student
            .groups = "drop") %>%
  filter(!is.na(level2_resid))

school_fitted <- df_resids %>%
  group_by(schcode) %>%
  summarise(fitted_school = mean(mu_post[id], na.rm = TRUE),
            level3_resid = first(level3_school),   # unique per school
            .groups = "drop") %>%
  filter(!is.na(level3_resid))

# level-2 residuals plot
rp_level2 <- ggplot(student_fitted, aes(x = fitted_student, y = level2_resid)) +
  geom_point(size = 1.2, alpha = 0.6) +
  xlim(-10, 30) + ylim(-30, 55) +
  geom_hline(yintercept = 0, colour = "blue", linewidth = 1) +
  labs(title = "Level-2",
       x = "Predicted Values", y = "Pearson Residuals") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold"))

# level-3 residuals plot  
rp_level3 <- ggplot(school_fitted, aes(x = fitted_school, y = level3_resid)) +
  geom_point(size = 1.2, alpha = 0.6) +
  xlim(-10, 30) + ylim(-30, 55) + 
  geom_hline(yintercept = 0, colour = "blue", linewidth = 1) +
  labs(title = "Level-3",
       x = "Predicted Values", y = "Pearson Residuals") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold"))

print(rp_level2)
print(rp_level3)

rp23_fig <- ggarrange(rp_level2, rp_level3,
                      ncol = 2, nrow = 1,
                      common.legend = TRUE, legend = "bottom")
rp23_fig

# Index plots (ordered by student/school ID)
df_index_level2 <- data.frame(
  index = 1:length(level2_resids), 
  resid = level2_resids
)

ip_level2 <- ggplot(df_index_level2, aes(x = index, y = resid)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "blue", linewidth = 1) +
  # scale_x_continuous(breaks = seq(0, length(level2_resids), length.out = 5)) +
  ylim(-30, 60) + 
  labs(title = "Level-2",
       x = "Observation Number", y = "Pearson Residuals") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold"))

df_index_level3 <- data.frame(
  index = 1:length(level3_resids), 
  resid = level3_resids
)

ip_level3 <- ggplot(df_index_level3, aes(x = index, y = resid)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "blue", linewidth = 1) +
  # scale_x_continuous(breaks = seq(0, length(level3_resids), length.out = 5)) +
  ylim(-30, 60) +
  labs(title = "Level-3",
       x = "Observation Number", y = "Pearson Residuals") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold"))

ip23_fig <- ggarrange(ip_level2, ip_level3,
                      ncol = 2, nrow = 1,
                      common.legend = TRUE, legend = "bottom")
ip23_fig

# Histograms
hg_level2 <- ggplot(data.frame(resids = level2_resids), aes(x = resids)) +
  geom_histogram(aes(y = after_stat(density)), bins = min(30, length(level2_resids)/5), 
                 fill = "grey", color = "black", alpha = 0.7) +
  geom_density(color = "blue", linewidth = 1) +
  xlim(-25, 25) +
  ylim(0, 0.5) +
  labs(title = "Level-2",
       x = "Pearson Residuals", y = "Density") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold"))

hg_level3 <- ggplot(data.frame(resids = level3_resids), aes(x = resids)) +
  geom_histogram(aes(y = after_stat(density)), bins = min(20, length(level3_resids)/3), 
                 fill = "grey", color = "black", alpha = 0.7) +
  geom_density(color = "blue", linewidth = 1) +
  xlim(-25, 25) +
  ylim(0, 0.5) +
  labs(title = "Level-3",
       x = "Pearson Residuals", y = "Density") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold"))

hg23_fig <- ggarrange(hg_level2, hg_level3,
                      ncol = 2, nrow = 1,
                      common.legend = TRUE, legend = "bottom")
hg23_fig

# display all
print(ip_level2)
print(ip_level3)
print(hg_level2)
print(hg_level3)

# determine proportion of variance

# True unconditional 3-level null (intercept only)
fit_null <- brm(
  bf(SDQ ~ 1 + (1 | schcode/STID)),
  data = master,
  family = gaussian(),
  prior = c(prior(normal(20, 5), class = "Intercept"),
            prior(normal(0, 0.2), class = "sd")),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, backend = "cmdstanr"
)

vc <- VarCorr(fit_null)

# extract posterior mean SDs
var_school  <- vc$schcode$sd[1,1]^2  
var_student <- vc$`schcode:STID`$sd[1,1]^2
var_wave    <- sigma(fit_null)^2

var_total <- var_school + var_student + var_wave
pct_school  <- 100 * var_school / var_total
pct_student <- 100 * var_student / var_total
pct_wave    <- 100 * var_wave / var_total

data.frame(
  Level = c("School", "Student", "Wave"),
  Variance = round(c(var_school, var_student, var_wave), 2),
  SD = round(sqrt(c(var_school, var_student, var_wave)), 2),
  `% Total` = round(c(pct_school, pct_student, pct_wave), 1)
)

# fixed intercept posterior summary (mean + 95% CI)
intercept_summary <- posterior_summary(fit_null, pars = "^b_Intercept")
intercept_mean <- intercept_summary[1, "Estimate"]
intercept_sd   <- intercept_summary[1, "Est.Error"]
intercept_ci   <- posterior_summary(fit_null, pars = "^b_Intercept", robust = FALSE)[1, c("Q2.5", "Q97.5")]

# single line
intercept_post <- fixef(fit_null)["Intercept", ]  # Mean [lower, upper]

# display with variances
data.frame(
  Parameter = c("Fixed Intercept", "School Var", "Student Var", "Wave Var"),
  Estimate = c(
    round(intercept_mean, 2),
    round(var_school, 2),
    round(var_student, 2),
    round(var_wave, 2)
  ),
  SD_CI = c(
    paste0(round(intercept_sd, 2), " [", round(intercept_ci[1], 2), ",", round(intercept_ci[2], 2), "]"),
    paste0(round(sqrt(var_school), 2)),
    paste0(round(sqrt(var_student), 2)),
    paste0(round(sigma(fit_null), 2))
  )
)

# extract intercept + SE (posterior SD)
intercept_se <- posterior_summary(fit_null, pars = "^b_Intercept")[1, "Est.Error"]

# full table with SE
data.frame(
  Parameter = c("Fixed Intercept", "School Var", "Student Var", "Wave Var"),
  Estimate = c(
    round(fixef(fit_null)["Intercept", 1], 2),     # mean
    round(var_school, 2),
    round(var_student, 2),
    round(var_wave, 2)
  ),
  SE = c(
    round(intercept_se, 2),                         # posterior SD = SE
    round(sqrt(var_school), 2),                     # SD for variances
    round(sqrt(var_student), 2),
    round(sigma(fit_null), 2)
  )
)

# MAR test

# corrected complete-case filters
vars_core <- c("SDQ", "time", "case_sch", "schcode", "STID")
master_cc_core <- master %>% 
  filter(complete.cases(across(all_of(vars_core))))

vars_full <- c("SDQ", "time", "case_sch", "S_sefi2", "Female", "schcode", "STID")
master_cc_full <- master %>% 
  filter(complete.cases(across(all_of(vars_full))))

vars_time <- c("time", "STID", "Female", "S_sefi2")
master_time <- master %>% 
  filter(complete.cases(across(all_of(vars_time))))

fit <- function(primary_data, cc_data = NULL) {
  data <- if (!is.null(cc_data)) cc_data else primary_data
  data <- data %>%
    mutate(schcode = factor(schcode), STID = factor(STID))
  
  brm(
    bf(SDQ ~ time * case_sch + (1 + time | schcode/STID)),
    data = data,
    prior = priors_melsm,
    family = gaussian(),
    chains = 4, iter = 2000, warmup = 1000,
    cores = 4, threads = threading(2, grainsize = 500),
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    backend = "cmdstanr",
    stan_model_args = list(stanc_options = list(O1 = TRUE)),
    refresh = 1
  )
}

fitC <- function(primary_data, cc_data = NULL) {
  data <- if (!is.null(cc_data)) cc_data else primary_data
  data <- data %>%
    mutate(schcode = factor(schcode), STID = factor(STID))
  
  brm(
    bf(SDQ ~ time * case_sch + S_sefi2 + Female + (1 + time | schcode/STID),
       sigma ~ time * case_sch +  S_sefi2 + Female + (1 + time | schcode/STID)),
    data = data,
    prior = priors_melsm,  # your existing priors
    family = gaussian(),
    chains = 4, iter = 2000, warmup = 1000,
    cores = 4, threads = threading(2, grainsize = 500),
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    backend = "cmdstanr",
    stan_model_args = list(stanc_options = list(O1 = TRUE)),
    refresh = 1
  )
}

fitC_full <- fitC(master_cc_full)
fitA_full <- fit(master_cc_full)

sumC_full <- summary(fitC_full)
sumC_cc   <- summary(fitC_cc)
tabC_full <- as.data.frame(sumC_full$fixed)
tabC_cc   <- as.data.frame(sumC_cc$fixed)

save(fitA_full, file = "fitA_full.RData")
save(fitC_full, file = "fitC_full_v2_noMC.RData")

load("fitC_full.RData")
fitC <- fitC_full

load("UPDATE_brms_melsm_fit_4chains_2000iter.RData")
fitC_cc <- fit_melsm

load("fitA_full.RData")
fitA <- fitA_full

load("fitA_cc.RData")

tab_compare_C <- cbind(
  Parameter = rownames(tabC),
  Est_C  = round(tabC$Estimate, 2),
  SE_C   = round(tabC$Est.Error, 2),
  L95_C  = round(tabC$`l-95% CI`, 2),
  U95_C  = round(tabC$`u-95% CI`, 2),
  Est_C_ss    = round(tabC_cc$Estimate, 2),
  SE_C_ss     = round(tabC_cc$Est.Error, 2),
  L95_C_ss    = round(tabC_cc$`l-95% CI`, 2),
  U95_C_ss    = round(tabC_cc$`u-95% CI`, 2)
)

tab_compare_A <- cbind(
  Parameter = rownames(tabA),
  Est_A  = round(tabA$Estimate, 2),
  SE_A   = round(tabA$Est.Error, 2),
  L95_A  = round(tabA$`l-95% CI`, 2),
  U95_A  = round(tabA$`u-95% CI`, 2),
  Est_A_cc    = round(tabA_cc$Estimate, 2),
  SE_A_cc     = round(tabA_cc$Est.Error, 2),
  L95_A_cc    = round(tabA_cc$`l-95% CI`, 2),
  U95_A_cc    = round(tabA_cc$`u-95% CI`, 2)
)

load("UPDATE_brms_melsm_fit_4chains_2000iter.RData")
fit_C_orig <- fit_melsm
load("brms_e_melsm_fit_4chains_2000iter.RData")
fit_B_orig <- fit_melsm
load("UPDATE_brms_mlm_fit_4chains_2000iter.RData")
fit_A_orig <- fit_melsm

loo(fit_A_orig)
waic(fit_A_orig)

loo(fit_B_orig)
waic(fit_B_orig)

loo(fit_C_orig)
waic(fit_C_orig)

loo(fitA)
waic(fitA)

loo(fitC)
waic(fitC)
