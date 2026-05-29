## MSc Thesis - R Code for Simulation Studies
## Author: Marco Ayroso
## Email: ayrosom@myumanitoba.ca


#____________________________________________
# Load libraries
#____________________________________________

library(tidyverse)
library(readxl)
library(compiler)
library(simstudy)
library(faux)
library(afex)
library(broom.mixed)
library(writexl)
library(nlme)
library(brms)
library(jtools)
library(remotes)
library(data.table)
library(knitr)
library(glmmTMB)
library(dplyr)
library(tidyr)
library(furrr)
library(future)
library(progressr)
library(ggpubr)
library(tictoc)
library(ggpubr)

# Parallelize
available_cores <- parallel::detectCores()
plan(multisession, workers = max(1, available_cores - 1))

future::plan(multisession, workers = 7)

# Enable progress bar
handlers("txtprogressbar")

setwd("/Users/marcoayroso/Library/Mobile Documents/com~apple~CloudDocs/MSc Thesis/Thesis/Simulation Studies")

#____________________________________________
# Simulation
#____________________________________________

# MLM-generated (study 2)
# *******For study 1, see comments for modification********

tic("Time elapsed:")

simulate <- function(theta_101, rep_id, p = NULL) {
  if (!is.null(p)) p()

  n_classrooms <- 200
  n_students_total <- 4000
  n_time <- 3
  n_students_per_class <- n_students_total / n_classrooms
  gamma_01j <- 0.5
  gamma_11j <- 0.1
  pi_000 <- 10
  pi_001 <- 0.001
  pi_100 <- 0.001
  pi_101 <- -1.0
  sigma_class_loc <- 0.3
  sigma_student_loc <- sqrt(0.3)
  sigma_resid_loc <- sqrt(0.5)
  theta_101_true <- 0
  missing_prob <- 0.10
  
  classroom_df <- data.frame(
    classroom = 1:n_classrooms,
    RCTGrp = rbinom(n_classrooms, 1, 0.5),
    r_00j_loc = rnorm(n_classrooms, 0, sigma_class_loc)
  )
  
  student_df <- do.call(rbind, lapply(1:n_classrooms, function(j) {
    data.frame(
      classroom = j,
      student = 1:n_students_per_class,
      gender = rbinom(n_students_per_class, 1, 0.505), # COMMENT OUT GENDER FOR STUDY 1
      u_0ij_loc = rnorm(n_students_per_class, 0, sigma_student_loc)
    )
  }))
  student_df <- merge(student_df, classroom_df, by = "classroom")
  
  data_df <- do.call(rbind, lapply(1:nrow(student_df), function(idx) {
    s_row <- student_df[idx, ]
    
    gamma_00j <- pi_000 + pi_001 * s_row$RCTGrp + s_row$r_00j_loc
    gamma_10j <- pi_100 + pi_101 * s_row$RCTGrp
    b_0ij <- gamma_00j + gamma_01j * s_row$gender + s_row$u_0ij_loc # COMMENT OUT GENDER FOR STUDY 1
    b_1ij <- gamma_10j + gamma_11j * s_row$gender # COMMENT OUT GENDER FOR STUDY 1
    
    time_points <- 0:(n_time - 1)
    sigma_e_tij <- sigma_resid_loc
    
    data.frame(
      classroom = s_row$classroom,
      student = s_row$student,
      measure_tij = time_points,
      RCTGrp = s_row$RCTGrp,
      gender = s_row$gender, # COMMENT OUT GENDER FOR STUDY 1
      b_0ij = b_0ij,
      b_1ij = b_1ij,
      sigma_e_tij = sigma_e_tij
    )
  }))
  
  data_df <- data_df %>%
    mutate(
      e_tij = rnorm(nrow(.), 0, sigma_e_tij),
      y_tij = b_0ij + b_1ij * measure_tij + e_tij
    )
  
  data_df <- data_df %>%
    group_by(measure_tij) %>%
    mutate(
      is_missing = rbinom(n(), 1, missing_prob),
      y_tij = ifelse(is_missing == 1, NA, y_tij),
      measure_tij = as.numeric(measure_tij)
    ) %>%
    ungroup()
  
  fit_data <- filter(data_df, !is.na(y_tij))
  

  # MLM fit with glmmTMB
  fit_mlm <- try(
    glmmTMB(
      y_tij ~ measure_tij * RCTGrp + gender * measure_tij + (1 + measure_tij | classroom/student), # COMMENT OUT GENDER FOR STUDY 1
      data = fit_data,
      family = gaussian()
    ), silent = TRUE
  )
  
  # set priors for brm
  priors <- c(
    prior(normal(0, 3), class = "b"),
    prior(normal(20, 5), class = "Intercept"),
    prior(normal(0, 0.2), class = "sd"),
    prior(normal(0, 0.3), class = "b", dpar = "sigma"),
    prior(normal(log(5), 0.3), class = "Intercept", dpar = "sigma"),
    prior(normal(0, 0.2), class = "sd", dpar = "sigma"),
    prior(normal(-1.0, 0.1), class = "b", coef = "measure_tij:RCTGrp"),
    prior(normal(0, 0.1), class = "b", dpar = "sigma", coef = "measure_tij:RCTGrp")

  )
  
  # MELSM fit with brms
  fit_melsm <- try(
    brm(
      bf(
        y_tij ~ measure_tij * RCTGrp + gender * measure_tij + (1 + measure_tij | classroom/student), # COMMENT OUT GENDER FOR STUDY 1
        sigma ~ measure_tij * RCTGrp + gender * measure_tij + (1 + measure_tij | classroom/student) # COMMENT OUT GENDER FOR STUDY 1
      ),
      data = fit_data,
      prior = priors,
      family = gaussian(),
      iter = 1000,
      warmup = 500,
      chains = 1,
      cores = 4,
      backend = "cmdstanr",
      threads = threading(2, grainsize = 500),
      control = list(adapt_delta = 0.9, max_treedepth = 10),
      stan_model_args = list(stanc_options = list(O1 = TRUE)),
      refresh = 1
    ), silent = TRUE
  )
  
  # Parameter extraction for glmmTMB
  extract_params_glmmTMB <- function(fit, analysis_model) {
    if (inherits(fit, "try-error") || is.null(fit)) {
      return(tibble(
        replication = rep_id, theta_101 = theta_101, analysis_model = analysis_model,
        param = c("pi101", "theta101", "gamma_01j", "gamma_11j", "tau_01j", "tau_11j"),
        estimate = rep(NA_real_, 6),
        se = rep(NA_real_, 6),
        ci_lower = rep(NA_real_, 6),
        ci_upper = rep(NA_real_, 6),
        p_value = rep(NA_real_, 6)
      ))
    }
    s <- summary(fit)
    cond <- s$coefficients$cond
    disp <- s$coefficients$disp
    
    get_val <- function(mat, row, col) if(row %in% rownames(mat)) mat[row, col] else NA_real_
    
    est_pi101 <- get_val(cond, "measure_tij:RCTGrp", "Estimate")
    se_pi101 <- get_val(cond, "measure_tij:RCTGrp", "Std. Error")
    p_pi101 <- get_val(cond, "measure_tij:RCTGrp", "Pr(>|z|)")
    ci_pi101_l <- est_pi101 - 1.96 * se_pi101
    ci_pi101_u <- est_pi101 + 1.96 * se_pi101
    
    est_gamma01j <- get_val(cond, "gender", "Estimate")
    se_gamma01j <- get_val(cond, "gender", "Std. Error")
    p_gamma01j <- get_val(cond, "gender", "Pr(>|z|)")
    ci_gamma01j_l <- est_gamma01j - 1.96 * se_gamma01j
    ci_gamma01j_u <- est_gamma01j + 1.96 * se_gamma01j
    
    gamma11_name <- intersect(c("gender:measure_tij", "measure_tij:gender"), rownames(cond))
    est_gamma11j <- se_gamma11j <- p_gamma11j <- NA_real_
    if(length(gamma11_name) > 0) {
      est_gamma11j <- cond[gamma11_name[1], "Estimate"]
      se_gamma11j <- cond[gamma11_name[1], "Std. Error"]
      p_gamma11j <- cond[gamma11_name[1], "Pr(>|z|)"]
    }
    ci_gamma11j_l <- est_gamma11j - 1.96 * se_gamma11j
    ci_gamma11j_u <- est_gamma11j + 1.96 * se_gamma11j
    
    est_theta101 <- se_theta101 <- p_theta101 <- ci_theta101_l <- ci_theta101_u <- NA_real_
    est_tau01j <- se_tau01j <- p_tau01j <- ci_tau01j_l <- ci_tau01j_u <- NA_real_
    est_tau11j <- se_tau11j <- p_tau11j <- ci_tau11j_l <- ci_tau11j_u <- NA_real_
    
    if(!is.null(disp)) {
      est_theta101 <- get_val(disp, "measure_tij:RCTGrp", "Estimate")
      se_theta101 <- get_val(disp, "measure_tij:RCTGrp", "Std. Error")
      p_theta101 <- get_val(disp, "measure_tij:RCTGrp", "Pr(>|z|)")
      ci_theta101_l <- est_theta101 - 1.96 * se_theta101
      ci_theta101_u <- est_theta101 + 1.96 * se_theta101
      
      est_tau01j <- get_val(disp, "gender", "Estimate")
      se_tau01j <- get_val(disp, "gender", "Std. Error")
      p_tau01j <- get_val(disp, "gender", "Pr(>|z|)")
      ci_tau01j_l <- est_tau01j - 1.96 * se_tau01j
      ci_tau01j_u <- est_tau01j + 1.96 * se_tau01j
      
      tau11_name <- intersect(rownames(disp), c("gender:measure_tij", "measure_tij:gender"))
      if(length(tau11_name) > 0) {
        est_tau11j <- disp[tau11_name[1], "Estimate"]
        se_tau11j <- disp[tau11_name[1], "Std. Error"]
        p_tau11j <- disp[tau11_name[1], "Pr(>|z|)"]
        ci_tau11j_l <- est_tau11j - 1.96 * se_tau11j
        ci_tau11j_u <- est_tau11j + 1.96 * se_tau11j
      }
    }
    
    tibble(
      replication = rep_id, theta_101 = theta_101, analysis_model = analysis_model,
      param = c("pi101", "theta101", "gamma_01j", "gamma_11j", "tau_01j", "tau_11j"),
      estimate = c(est_pi101, est_theta101, est_gamma01j, est_gamma11j, est_tau01j, est_tau11j),
      se = c(se_pi101, se_theta101, se_gamma01j, se_gamma11j, se_tau01j, se_tau11j),
      ci_lower = c(ci_pi101_l, ci_theta101_l, ci_gamma01j_l, ci_gamma11j_l, ci_tau01j_l, ci_tau11j_l),
      ci_upper = c(ci_pi101_u, ci_theta101_u, ci_gamma01j_u, ci_gamma11j_u, ci_tau01j_u, ci_tau11j_u),
      p_value = c(p_pi101, p_theta101, p_gamma01j, p_gamma11j, p_tau01j, p_tau11j)
    )
  }
  
  # Parameter extraction for brm
  extract_params_brm <- function(fit, analysis_model) {
    if (inherits(fit, "try-error") || is.null(fit)) {
      return(tibble(
        replication = rep_id, theta_101 = theta_101, analysis_model = analysis_model,
        param = c("pi101", "theta101", "gamma_01j", "gamma_11j"),
        estimate = rep(NA_real_, 4),
        se = rep(NA_real_, 4),
        ci_lower = rep(NA_real_, 4),
        ci_upper = rep(NA_real_, 4),
        p_value = rep(NA_real_, 4)
      ))
    }
    post <- posterior_summary(fit)
    
    get_post_param <- function(pattern) {
      idx <- grep(pattern, rownames(post), ignore.case = TRUE)
      if(length(idx) == 0) return(c(NA, NA, NA, NA))
      est <- post[idx[1], "Estimate"]
      se <- post[idx[1], "Est.Error"]
      lower <- post[idx[1], "Q2.5"]
      upper <- post[idx[1], "Q97.5"]
      c(est, se, lower, upper)
    }
    
    pi101 <- get_post_param("b_measure_tij:RCTGrp")
    theta101 <- get_post_param("b_sigma_measure_tij:RCTGrp")
    gamma01j <- get_post_param("b_gender$")
    gamma11j <- get_post_param("b_gender:measure_tij|b_measure_tij:gender")
    
    tibble(
      replication = rep_id, theta_101 = theta_101, analysis_model = analysis_model,
      param = c("pi101", "theta101", "gamma_01j", "gamma_11j"),
      estimate = c(pi101[1], theta101[1], gamma01j[1], gamma11j[1]),
      se = c(pi101[2], theta101[2], gamma01j[2], gamma11j[2]),
      ci_lower = c(pi101[3], theta101[3], gamma01j[3], gamma11j[3]),
      ci_upper = c(pi101[4], theta101[4], gamma01j[4], gamma11j[4]),
      p_value = rep(NA_real_, 4)
    )
  }
  
  res_mlm <- extract_params_glmmTMB(fit_mlm, "MLM")
  res_melsm <- extract_params_brm(fit_melsm, "MELSM")
  
  list(selected = bind_rows(res_mlm, res_melsm), sim_data = data_df)
}

# Simulate
n_sets <- 1
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/MSc Thesis/Thesis/Simulation Studies")
options(scipen = 999)

for (set_id in 1:n_sets) {
  set.seed(071714 + set_id)  # Set a different seed for each set to ensure different outputs
  
  params_grid <- expand.grid(theta_101 = 0, rep_id = 1:10)
  
  with_progress({
    p <- progressor(along = seq_len(nrow(params_grid)))
    results <- future_map2(
      params_grid$theta_101, params_grid$rep_id,
      ~ simulate(.x, .y, p),
      .options = furrr_options(seed = TRUE)
    )
  })
  
  results_df <- bind_rows(lapply(results, `[[`, "selected"))
  toc()

  write.csv(results_df, paste0("UPDATE_study_2_MLM_bayes_1000iter_500warmup_set_", set_id, ".csv"), row.names = FALSE)
}

# MELSM-generated (study 2)
# *******For study 1, see comments for modification********

tic("Time elapsed:")

simulate <- function(theta_101, rep_id, p = NULL) {
  if (!is.null(p)) p()
  
  n_classrooms <- 200
  n_students_total <- 4000
  n_time <- 3
  n_students_per_class <- n_students_total / n_classrooms
  pi_000 <- 10
  pi_001 <- 0.001
  pi_100 <- 0.001
  pi_101 <- -1.0
  sigma_class_loc <- 0.3
  sigma_student_loc <- sqrt(0.3)
  sigma_resid_loc <- sqrt(0.5)
  theta_000 <- 1.0
  theta_001 <- 0.001
  theta_100 <- 0.001
  sigma_class_scale <- 0.4
  sigma_student_scale <- 0.3
  gamma_01j <- 0.5
  gamma_11j <- 0.1
  tau_01j <- 0.8
  tau_11j <- 0.2
  missing_prob <- 0.10
  
  classroom_df <- data.frame(
    classroom = 1:n_classrooms,
    RCTGrp = rbinom(n_classrooms, 1, 0.5),
    r_00j_loc = rnorm(n_classrooms, 0, sigma_class_loc),
    omega_00j_scale = rnorm(n_classrooms, 0, sigma_class_scale)
  )
  student_df <- do.call(rbind, lapply(1:n_classrooms, function(j) {
    data.frame(
      classroom = j,
      student = 1:n_students_per_class,
      gender = rbinom(n_students_per_class, 1, 0.505), # COMMENT OUT GENDER FOR STUDY 1
      u_0ij_loc = rnorm(n_students_per_class, 0, sigma_student_loc),
      v_0ij_scale = rnorm(n_students_per_class, 0, sigma_student_scale)
    )
  }))
  student_df <- merge(student_df, classroom_df, by = "classroom")
  data_df <- do.call(rbind, lapply(1:nrow(student_df), function(idx) {
    s_row <- student_df[idx, ]
    gamma_00j <- pi_000 + pi_001 * s_row$RCTGrp + s_row$r_00j_loc
    gamma_10j <- pi_100 + pi_101 * s_row$RCTGrp
    b_0ij <- gamma_00j + gamma_01j * s_row$gender + s_row$u_0ij_loc # COMMENT OUT GENDER FOR STUDY 1
    b_1ij <- gamma_10j + gamma_11j * s_row$gender # COMMENT OUT GENDER FOR STUDY 1
    tau_00j <- theta_000 + theta_001 * s_row$RCTGrp + s_row$omega_00j_scale
    tau_10j <- theta_100 + theta_101 * s_row$RCTGrp
    a_0ij <- tau_00j + tau_01j * s_row$gender + s_row$v_0ij_scale # COMMENT OUT GENDER FOR STUDY 1
    a_1ij <- tau_10j + tau_11j * s_row$gender # COMMENT OUT GENDER FOR STUDY 1
    time_points <- 0:(n_time - 1)
    log_var_e <- a_0ij + a_1ij * time_points
    sigma_e_tij <- sqrt(exp(log_var_e))
    data.frame(
      classroom = s_row$classroom,
      student = s_row$student,
      measure_tij = time_points,
      RCTGrp = s_row$RCTGrp,
      gender = s_row$gender, # COMMENT OUT GENDER FOR STUDY 1
      b_0ij = b_0ij,
      b_1ij = b_1ij,
      a_0ij = a_0ij,
      a_1ij = a_1ij,
      sigma_e_tij = sigma_e_tij
    )
  }))
  
  data_df$e_tij <- rnorm(nrow(data_df), 0, data_df$sigma_e_tij)
  data_df$y_tij <- data_df$b_0ij + data_df$b_1ij * data_df$measure_tij + data_df$e_tij
  
  data_df <- data_df %>%
    group_by(measure_tij) %>%
    mutate(
      is_missing = rbinom(n(), 1, missing_prob),
      y_tij = ifelse(is_missing == 1, NA, y_tij),
      measure_tij = as.numeric(measure_tij)
    ) %>%
    ungroup()
  
  fit_data <- filter(data_df, !is.na(y_tij))
  
  # MLM fit - glmmTMB
  fit_mlm <- try(
    glmmTMB(
      y_tij ~ measure_tij * RCTGrp + gender * measure_tij + (1 + measure_tij | classroom/student), # COMMENT OUT GENDER FOR STUDY 1
      data = fit_data,
      family = gaussian()
    ), silent = TRUE
  )
  
  # priors
  priors <- c(
    prior(normal(0, 3), class = "b"),
    prior(normal(20, 5), class = "Intercept"),
    prior(normal(0, 0.2), class = "sd"),
    prior(normal(0, 0.3), class = "b", dpar = "sigma"),
    prior(normal(log(5), 0.3), class = "Intercept", dpar = "sigma"),
    prior(normal(0, 0.2), class = "sd", dpar = "sigma"),
    prior(normal(-1.0, 0.1), class = "b", coef = "measure_tij:RCTGrp"),
    prior(normal(0, 0.1), class = "b", dpar = "sigma", coef = "measure_tij:RCTGrp")
    
  )
  
  # MELSM fit - brm
  fit_melsm <- try(
    brm(
      bf(
        y_tij ~ measure_tij * RCTGrp + gender * measure_tij + (1 + measure_tij | classroom/student), # COMMENT OUT GENDER FOR STUDY 1
        sigma ~ measure_tij * RCTGrp + gender * measure_tij + (1 + measure_tij | classroom/student) # COMMENT OUT GENDER FOR STUDY 1
      ),
      data = fit_data,
      prior = priors,
      family = gaussian(),
      iter = 1000,
      warmup = 500,
      chains = 1,
      cores = 4,
      backend = "cmdstanr",
      threads = threading(2, grainsize = 500),
      control = list(adapt_delta = 0.9, max_treedepth = 10),
      stan_model_args = list(stanc_options = list(O1 = TRUE)),
      refresh = 1
    ), silent = TRUE
  )
  
  # Parameter extraction for glmmTMB
  extract_params_glmmTMB <- function(fit, analysis_model) {
    if (inherits(fit, "try-error") || is.null(fit)) {
      return(tibble(
        replication = rep_id, theta_101 = theta_101, analysis_model = analysis_model,
        param = c("pi101", "theta101", "gamma_01j", "gamma_11j", "tau_01j", "tau_11j"),
        estimate = rep(NA_real_, 6),
        se = rep(NA_real_, 6),
        ci_lower = rep(NA_real_, 6),
        ci_upper = rep(NA_real_, 6),
        p_value = rep(NA_real_, 6)
      ))
    }
    s <- summary(fit)
    cond <- s$coefficients$cond
    disp <- s$coefficients$disp
    
    get_val <- function(mat, row, col) if(row %in% rownames(mat)) mat[row, col] else NA_real_
    
    est_pi101 <- get_val(cond, "measure_tij:RCTGrp", "Estimate")
    se_pi101 <- get_val(cond, "measure_tij:RCTGrp", "Std. Error")
    p_pi101 <- get_val(cond, "measure_tij:RCTGrp", "Pr(>|z|)")
    ci_pi101_l <- est_pi101 - 1.96 * se_pi101
    ci_pi101_u <- est_pi101 + 1.96 * se_pi101
    
    est_gamma01j <- get_val(cond, "gender", "Estimate")
    se_gamma01j <- get_val(cond, "gender", "Std. Error")
    p_gamma01j <- get_val(cond, "gender", "Pr(>|z|)")
    ci_gamma01j_l <- est_gamma01j - 1.96 * se_gamma01j
    ci_gamma01j_u <- est_gamma01j + 1.96 * se_gamma01j
    
    gamma11_name <- intersect(c("gender:measure_tij", "measure_tij:gender"), rownames(cond))
    est_gamma11j <- se_gamma11j <- p_gamma11j <- NA_real_
    if(length(gamma11_name) > 0) {
      est_gamma11j <- cond[gamma11_name[1], "Estimate"]
      se_gamma11j <- cond[gamma11_name[1], "Std. Error"]
      p_gamma11j <- cond[gamma11_name[1], "Pr(>|z|)"]
    }
    ci_gamma11j_l <- est_gamma11j - 1.96 * se_gamma11j
    ci_gamma11j_u <- est_gamma11j + 1.96 * se_gamma11j
    
    est_theta101 <- se_theta101 <- p_theta101 <- ci_theta101_l <- ci_theta101_u <- NA_real_
    est_tau01j <- se_tau01j <- p_tau01j <- ci_tau01j_l <- ci_tau01j_u <- NA_real_
    est_tau11j <- se_tau11j <- p_tau11j <- ci_tau11j_l <- ci_tau11j_u <- NA_real_
    
    if(!is.null(disp)) {
      est_theta101 <- get_val(disp, "measure_tij:RCTGrp", "Estimate")
      se_theta101 <- get_val(disp, "measure_tij:RCTGrp", "Std. Error")
      p_theta101 <- get_val(disp, "measure_tij:RCTGrp", "Pr(>|z|)")
      ci_theta101_l <- est_theta101 - 1.96 * se_theta101
      ci_theta101_u <- est_theta101 + 1.96 * se_theta101
      
      est_tau01j <- get_val(disp, "gender", "Estimate")
      se_tau01j <- get_val(disp, "gender", "Std. Error")
      p_tau01j <- get_val(disp, "gender", "Pr(>|z|)")
      ci_tau01j_l <- est_tau01j - 1.96 * se_tau01j
      ci_tau01j_u <- est_tau01j + 1.96 * se_tau01j
      
      tau11_name <- intersect(rownames(disp), c("gender:measure_tij", "measure_tij:gender"))
      if(length(tau11_name) > 0) {
        est_tau11j <- disp[tau11_name[1], "Estimate"]
        se_tau11j <- disp[tau11_name[1], "Std. Error"]
        p_tau11j <- disp[tau11_name[1], "Pr(>|z|)"]
        ci_tau11j_l <- est_tau11j - 1.96 * se_tau11j
        ci_tau11j_u <- est_tau11j + 1.96 * se_tau11j
      }
    }
    
    tibble(
      replication = rep_id, theta_101 = theta_101, analysis_model = analysis_model,
      param = c("pi101", "theta101", "gamma_01j", "gamma_11j", "tau_01j", "tau_11j"),
      estimate = c(est_pi101, est_theta101, est_gamma01j, est_gamma11j, est_tau01j, est_tau11j),
      se = c(se_pi101, se_theta101, se_gamma01j, se_gamma11j, se_tau01j, se_tau11j),
      ci_lower = c(ci_pi101_l, ci_theta101_l, ci_gamma01j_l, ci_gamma11j_l, ci_tau01j_l, ci_tau11j_l),
      ci_upper = c(ci_pi101_u, ci_theta101_u, ci_gamma01j_u, ci_gamma11j_u, ci_tau01j_u, ci_tau11j_u),
      p_value = c(p_pi101, p_theta101, p_gamma01j, p_gamma11j, p_tau01j, p_tau11j)
    )
  }
  
  extract_params_brm <- function(fit, analysis_model) {
    if (inherits(fit, "try-error") || is.null(fit)) {
      return(tibble(
        replication = rep_id, theta_101 = theta_101, analysis_model = analysis_model,
        param = c("pi101", "theta101", "gamma_01j", "gamma_11j"),
        estimate = rep(NA_real_, 4),
        se = rep(NA_real_, 4),
        ci_lower = rep(NA_real_, 4),
        ci_upper = rep(NA_real_, 4),
        p_value = rep(NA_real_, 4)
      ))
    }
    post <- posterior_summary(fit)
    
    get_post_param <- function(pattern) {
      idx <- grep(pattern, rownames(post), ignore.case = TRUE)
      if(length(idx) == 0) return(c(NA, NA, NA, NA))
      est <- post[idx[1], "Estimate"]
      se <- post[idx[1], "Est.Error"]
      lower <- post[idx[1], "Q2.5"]
      upper <- post[idx[1], "Q97.5"]
      c(est, se, lower, upper)
    }
    
    pi101 <- get_post_param("b_measure_tij:RCTGrp")
    theta101 <- get_post_param("b_sigma_measure_tij:RCTGrp")
    gamma01j <- get_post_param("b_gender$")
    gamma11j <- get_post_param("b_gender:measure_tij|b_measure_tij:gender")
    
    tibble(
      replication = rep_id, theta_101 = theta_101, analysis_model = analysis_model,
      param = c("pi101", "theta101", "gamma_01j", "gamma_11j"),
      estimate = c(pi101[1], theta101[1], gamma01j[1], gamma11j[1]),
      se = c(pi101[2], theta101[2], gamma01j[2], gamma11j[2]),
      ci_lower = c(pi101[3], theta101[3], gamma01j[3], gamma11j[3]),
      ci_upper = c(pi101[4], theta101[4], gamma01j[4], gamma11j[4]),
      p_value = rep(NA_real_, 4)
    )
  }
  
  res_mlm <- extract_params_glmmTMB(fit_mlm, "MLM")
  res_melsm <- extract_params_brm(fit_melsm, "MELSM")
  
  list(selected = bind_rows(res_mlm, res_melsm), sim_data = data_df)
}


# Simulate
n_sets <- 100
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/MSc Thesis/Thesis/Simulation Studies")
options(scipen = 999)

for (set_id in 1:n_sets) {
  set.seed(071714 + set_id)  # Set a different seed for each set to ensure different outputs
  
  params_grid <- expand.grid(theta_101 = c(-0.5, 0, 0.5), rep_id = 1:10)
  
  with_progress({
    p <- progressor(along = seq_len(nrow(params_grid)))
    results <- future_map2(
      params_grid$theta_101, params_grid$rep_id,
      ~ simulate(.x, .y, p),
      .options = furrr_options(seed = TRUE)
    )
  })
  
  results_df <- bind_rows(lapply(results, `[[`, "selected"))
  toc()
  
  write.csv(results_df, paste0("UPDATE_study_1_MELSM_bayes_1000iter_500warmup_set_", set_id, ".csv"), row.names = FALSE)
}


#____________________________________________
# Post simulation
#____________________________________________

###            ###
### LOAD FILES ###
###            ###

study_1_MLM <- sprintf("UPDATE_study_1_MLM_bayes_1000iter_500warmup_sets_%d.csv", 1:100)

study_1_MLM <- study_1_MLM %>%
  lapply(read_csv) %>%
  bind_rows()

study_2_MLM <- sprintf("UPDATE_study_2_MLM_bayes_1000iter_500warmup_set_%d.csv", 1:100)

study_2_MLM <- study_2_MLM %>%
  lapply(read_csv) %>%
  bind_rows()

study_1_MELSM <- sprintf("UPDATE_study_1_MELSM_bayes_1000iter_500warmup_set_%d.csv", 1:100)

study_1_MELSM <- study_1_MELSM %>%
  lapply(read_csv) %>%
  bind_rows()

study_2_MELSM <- sprintf("UPDATE_study_2_MELSM_bayes_1000iter_500warmup_set_%d.csv", 1:100)

study_2_MELSM <- study_2_MELSM %>%
  lapply(read_csv) %>%
  bind_rows()

# write.csv(study_1_MLM, "study_1_MLM_final.csv", row.names = TRUE)
# write.csv(study_2_MLM, "study_2_MLM_final.csv", row.names = TRUE)
# write.csv(study_1_MELSM, "study_1_MELSM_final.csv", row.names = TRUE)
# write.csv(study_2_MELSM, "study_2_MELSM_final.csv", row.names = TRUE)

study_1_MLM <- study_1_MLM %>%
  mutate(
         Study = rep("Study 1 - MLM"))

study_2_MLM <- study_2_MLM %>%
  mutate(
         Study = rep("Study 2 - MLM"))

study_1_MELSM <- study_1_MELSM %>%
  mutate(
         Study = rep("Study 1 - MELSM"))

study_2_MELSM <- study_2_MELSM %>%
  mutate(
         Study = rep("Study 2 - MELSM"))

MLM_vis <- list(study_1_MLM, study_2_MLM)
MLM_vis <- bind_rows(MLM_vis)

MLM_vis_pi101 <- MLM_vis %>%
  filter(param == "pi101")

MLM_vis_theta101 <- MLM_vis %>%
  filter(param == "theta101")

MELSM_vis <- list(study_1_MELSM, study_2_MELSM)
MELSM_vis <- bind_rows(MELSM_vis)

MELSM_vis_theta101 <- MELSM_vis %>%
  filter(param == "theta101")

MELSM_vis_pi101 <- MELSM_vis %>%
  filter(param == "pi101")

pi101 <- bind_rows(MLM_vis_pi101, MELSM_vis_pi101)
theta101 <- bind_rows(MELSM_vis_theta101)

theta101_improve <- theta101 %>%
  filter(theta_101 == -0.5)

theta101_equal <- theta101 %>%
  filter(theta_101 == 0)

theta101_dispar <- theta101 %>%
  filter(theta_101 == 0.5)

all_data <- list(study_1_MLM, study_1_MELSM, study_2_MLM, study_2_MELSM)
all_data <- bind_rows(all_data)

all_data$Study <- as.factor(all_data$Study)

# write.csv(all_data, "master.csv", row.names = TRUE)

###               ###
### VISUALIZATION ###
###               ###

# all pi101 (homoscedastic and heteroscedastic)
ggplot(pi101, aes(x = estimate, color = Study, fill = Study)) +
  geom_density(alpha = 0.3) +   # smooth density curves
  geom_vline(aes(xintercept = -1.0), 
             linetype = "dashed", color = "black") +
  ggtitle(expression(paste(pi, "_101"))) +
  scale_x_continuous(limits = c(-1.4, -0.6), breaks = seq(-1.4, -0.6, by = 0.1)) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, by = 5)) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "Estimate", y = "Density")

# theta101 = -0.5
t101_i <- ggplot(theta101_improve, aes(x = estimate, color = Study, fill = Study)) +
  geom_density(alpha = 0.3) +   # smooth density curves
  geom_vline(aes(xintercept = -0.5), 
             linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(-0.6, -0.1), breaks = seq(-0.6, -0.1, by = 0.1)) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, by = 5)) +
  labs(x = "Estimate", y = "Density", title = "Overlayed Smooth Histograms") +
  ggtitle(expression(paste(theta, "_101 = -0.5"))) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "", y = "Density")

# theta101 = 0
t101_e <- ggplot(theta101_equal, aes(x = estimate, color = Study, fill = Study)) +
  geom_density(alpha = 0.3) +   # smooth density curves
  geom_vline(aes(xintercept = 0), 
             linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(-0.2, 0.2)) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, by = 5)) +
  labs(x = "Estimate", y = "Density", title = "Overlayed Smooth Histograms") +
  ggtitle(expression(paste(theta, "_101 = 0"))) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "Estimate", y = "")

# theta101 = 0.5
t101_d <- ggplot(theta101_dispar, aes(x = estimate, color = Study, fill = Study)) +
  geom_density(alpha = 0.3) +   # smooth density curves
  geom_vline(aes(xintercept = 0.5), 
             linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0.1, 0.6), breaks = seq(0.1, 0.6, by = 0.1)) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, by = 5)) +
  labs(x = "Estimate", y = "Density", title = "Overlayed Smooth Histograms") +
  ggtitle(expression(paste(theta, "_101 = 0.5"))) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "", y = "")

t101 <- ggarrange(t101_i, t101_e, t101_d,
                  ncol = 3, nrow = 1,
                  common.legend = TRUE, legend = "bottom")

t101

# load files

performance_metrics <- read_xlsx("Performance_Metrics_v3.xlsx", sheet = "pi_101")

performance_metrics_theta <- read_xlsx("Performance_Metrics_v3.xlsx", sheet = "theta_101")

# 1. Pivot longer so all metrics are in one column
plot_df <- performance_metrics %>%
  pivot_longer(
    cols = c("Mean estimate", "Bias", "MSE", Coverage, Power),
    names_to = "Metric",
    values_to = "Value"
  )

# 2. Factor order for nicer display
plot_df$Metric <- factor(plot_df$Metric,
                         levels = c("Mean estimate", "Bias", 
                                    "MSE", "Coverage", "Power"),
                         labels = c("Mean estimate", "Bias", 
                                    "Mean Squared Error", "Coverage", "Power"))

# MLM generated data only

MLM_data <- plot_df %>%
  filter(DGP == "MLM")

# Baseline values for each metric
baseline_df <- tibble(
  Metric = c("Mean estimate", "Bias", "Mean Squared Error", "Coverage", "Power"),
  baseline = c(-1.0, 0.0, 0.0, 0.95, 1.0)
)

baseline_df$Metric <- factor(baseline_df$Metric,
                             levels = c("Mean estimate", "Bias", 
                                        "Mean Squared Error", "Coverage", "Power"),
                             labels = c("Mean estimate", "Bias", 
                                        "Mean Squared Error", "Coverage", "Power"))


# 2. Factor order for nicer display
plot_df$Metric <- factor(plot_df$Metric,
                         levels = c("Mean estimate", "Bias", 
                                    "Mean Squared Error", "Coverage", "Power"),
                         labels = c("Mean estimate", "Bias", 
                                    "MSE", "Coverage", "Power"))

my_colours <- c("#365D8DFF", "#47C16EFF")

MLM_data$Model <- factor(MLM_data$Model,
                         levels = c("MLM", "MELSM"),
                         labels = c("MLM", "MELSM"))

# theta_101 prep

# 1. Pivot longer so all metrics are in one column
plot_df_theta <- performance_metrics_theta %>%
  pivot_longer(
    cols = c(`Mean estimate`, "Bias", "MSE", `Bias-adjusted coverage`, Power_bayes),
    names_to = "Metric",
    values_to = "Value"
  )

# 2. Factor order for nicer display
plot_df_theta$Metric <- factor(plot_df_theta$Metric,
                               levels = c("Mean estimate", "Bias", 
                                          "MSE", "Bias-adjusted coverage", "Power_bayes"),
                               labels = c("Mean estimate", "Bias", 
                                          "Mean Squared Error", "Bias-adjusted coverage", "Power"))

plot_df_theta$`Scenario/true value` <- factor(plot_df_theta$`Scenario/true value`,
                                              levels = c("𝜃_101 = -0.5", "𝜃_101 = 0", 
                                                         "𝜃_101 = 0.5"),
                                              labels = c("-0.5", "0", "0.5"))


# Baseline values for each metric
baseline_df_theta <- tibble(
  Metric = c("ME_l", "ME_0", "ME_u", "Bias", "Mean Squared Error", "Bias-adjusted coverage", "Power"),
  baseline = c(-0.50, 0, 0.50, 0.0, 0.0, 0.95, 1.0)
)

MELSM_data <- plot_df_theta %>%
  filter(DGP == "MELSM")

MELSM_data$Study <- as.factor(MELSM_data$Study)

my_colours <- c("#365D8DFF", "#47C16EFF")


# mean estimate
pi_1 <- ggplot(MLM_data[which(MLM_data$Metric=="Mean estimate"),], aes(x = factor(Study), y = Value, 
                                                                       color = Model, shape = Model, 
                                                                       group = Model)) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(data = baseline_df[which(baseline_df$Metric=="Mean estimate"),], 
             aes(yintercept = baseline), 
             linetype = "dashed", color = "black") +
  scale_color_manual(values = my_colours) +
  scale_y_continuous(limits = c(-1.1, -0.995), breaks = seq(-1.1, -1, by = 0.02)) +
  ggtitle("Mean Estimate") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "", y = "", color = "Model", shape = "Model")

# bias
pi_2 <- ggplot(MLM_data[which(MLM_data$Metric=="Bias"),], aes(x = factor(Study), y = Value, 
                                                              color = Model, shape = Model, 
                                                              group = Model)) +
  geom_point(size = 3) +
  geom_line() + 
  geom_hline(data = baseline_df[which(baseline_df$Metric=="Bias"),], 
             aes(yintercept = baseline), 
             linetype = "dashed", color = "black") +
  scale_color_manual(values = my_colours) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, by = 0.2)) +
  ggtitle("Bias") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "", y = "", color = "Model", shape = "Model")

# mean squared error
pi_3 <- ggplot(MLM_data[which(MLM_data$Metric=="Mean Squared Error"),], aes(x = factor(Study), y = Value, 
                                                                            color = Model, shape = Model, 
                                                                            group = Model)) +
  geom_point(size = 3) +
  geom_line() + 
  geom_hline(data = baseline_df[which(baseline_df$Metric=="Mean Squared Error"),], 
             aes(yintercept = baseline), 
             linetype = "dashed", color = "black") +
  scale_color_manual(values = my_colours) +
  scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by = 0.02)) +
  ggtitle("Mean Squared Error") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "Study", y = "", color = "Model", shape = "Model")

# coverage
pi_4 <- ggplot(MLM_data[which(MLM_data$Metric=="Coverage"),], aes(x = factor(Study), y = Value, 
                                                                  color = Model, shape = Model, 
                                                                  group = Model)) +
  geom_point(size = 3) +
  geom_line() + 
  geom_hline(data = baseline_df[which(baseline_df$Metric=="Coverage"),], 
             aes(yintercept = baseline), 
             linetype = "dashed", color = "black") +
  scale_color_manual(values = my_colours) +
  scale_y_continuous(limits = c(0.80, 1), breaks = seq(0.80, 1, by = 0.04)) +
  ggtitle("Coverage") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "", y = "", color = "Model", shape = "Model")

# power
pi_5 <- ggplot(MLM_data[which(MLM_data$Metric=="Power"),], aes(x = factor(Study), y = Value, 
                                                               color = Model, shape = Model, 
                                                               group = Model)) +
  geom_point(size = 3) +
  geom_line() + 
  geom_hline(data = baseline_df[which(baseline_df$Metric=="Power"),], 
             aes(yintercept = baseline), 
             linetype = "dashed", color = "black") +
  scale_color_manual(values = my_colours) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  ggtitle("Power") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "", y = "", color = "Model", shape = "Model")

pi_fig <- ggarrange(pi_1, pi_2, pi_3, pi_4, pi_5,
                    ncol = 5, nrow = 1,
                    common.legend = TRUE, legend = "bottom")

pi_fig


# MELSM generated data pi_101

MELSM_pi_data <- performance_metrics[which(performance_metrics$DGP=="MELSM"),]

MELSM_pi_data <- MELSM_pi_data %>%
  pivot_longer(
    cols = c("Mean estimate", "Bias", "MSE", Coverage, Power),
    names_to = "Metric",
    values_to = "Value"
  )

MELSM_pi_data$Metric <- factor(MELSM_pi_data$Metric,
                               levels = c("Mean estimate", "Bias", 
                                          "MSE", "Coverage", "Power"),
                               labels = c("Mean estimate", "Bias", 
                                          "MSE", "Coverage", "Power"))

MELSM_pi_data$Scenario <- factor(MELSM_pi_data$Scenario,
                                 levels = c("𝜃_101 = -0.5", "𝜃_101 = 0", 
                                            "𝜃_101 = 0.5"),
                                 labels = c("-0.5", "0", "0.5"))

my_colours_MELSM_pi <- c("#440154FF", "#365D8DFF", "#47C16EFF", "#FDE725FF")

MELSM_pi_data$`Model & Study` <- interaction(MELSM_pi_data$Model, MELSM_pi_data$Study)

MELSM_pi_data$`Model & Study` <- factor(MELSM_pi_data$`Model & Study`, 
                                        levels = c("MLM.1", "MELSM.1", "MLM.2", "MELSM.2"),
                                        labels = c("MLM - Study 1", "MELSM - Study 1", "MLM - Study 2", "MELSM - Study 2"))


# mean estimate
s1_pi_1 <- ggplot(MELSM_pi_data[which(MELSM_pi_data$Metric=="Mean estimate"),], aes(x = Scenario, y = Value, 
                                                                                    color = `Model & Study`, shape = `Model & Study`, 
                                                                                    group = `Model & Study`)) +
  geom_point(size = 3) +
  geom_line() + 
  geom_hline(data = baseline_df[which(baseline_df$Metric=="Mean estimate"),],
             aes(yintercept = baseline), 
             linetype = "dashed", color = "black") +
  scale_color_manual(values = my_colours_MELSM_pi) +
  scale_y_continuous(limits = c(-1.1, -0.995), breaks = seq(-1.1, -1, by = 0.02)) +
  ggtitle("Mean Estimate") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "", y = "")

# bias
s1_pi_2 <- ggplot(MELSM_pi_data[which(MELSM_pi_data$Metric=="Bias"),], aes(x = Scenario, y = Value, 
                                                                           color = `Model & Study`, shape = `Model & Study`, 
                                                                           group = `Model & Study`)) +
  geom_point(size = 3) +
  geom_line()  +
  geom_hline(data = baseline_df[which(baseline_df$Metric=="Bias"),], 
             aes(yintercept = baseline), 
             linetype = "dashed", color = "black") +
  scale_color_manual(values = my_colours_MELSM_pi) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, by = 0.2)) +
  ggtitle("Bias") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "", y = "")

# mean squared error
s1_pi_3 <- ggplot(MELSM_pi_data[which(MELSM_pi_data$Metric=="MSE"),], aes(x = Scenario, y = Value, 
                                                                          color = `Model & Study`, shape = `Model & Study`, 
                                                                          group = `Model & Study`)) +
  geom_point(size = 3) +
  geom_line()  +
  geom_hline(data = baseline_df[which(baseline_df$Metric=="Mean Squared Error"),], 
             aes(yintercept = baseline), 
             linetype = "dashed", color = "black") +
  scale_color_manual(values = my_colours_MELSM_pi) +
  scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by = 0.02)) +
  ggtitle("Mean Squared Error") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = expression(paste(theta, "_101")), y = "")

# coverage
s1_pi_4 <- ggplot(MELSM_pi_data[which(MELSM_pi_data$Metric=="Coverage"),], aes(x = Scenario, y = Value, 
                                                                               color = `Model & Study`, shape = `Model & Study`, 
                                                                               group = `Model & Study`)) +
  geom_point(size = 3) +
  geom_line()  +
  geom_hline(data = baseline_df[which(baseline_df$Metric=="Coverage"),], 
             aes(yintercept = baseline), 
             linetype = "dashed", color = "black") +
  scale_color_manual(values = my_colours_MELSM_pi) +
  scale_y_continuous(limits = c(0.8, 1), breaks = seq(0.8, 1, by = 0.04)) +
  ggtitle("Coverage") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "", y = "")

# power
s1_pi_5 <- ggplot(MELSM_pi_data[which(MELSM_pi_data$Metric=="Power"),], aes(x = Scenario, y = Value, 
                                                                            color = `Model & Study`, shape = `Model & Study`, 
                                                                            group = `Model & Study`)) +
  geom_point(size = 3) +
  geom_line() + 
  geom_hline(data = baseline_df_theta[which(baseline_df_theta$Metric=="Power"),], 
             aes(yintercept = baseline), 
             linetype = "dashed", color = "black") +
  scale_color_manual(values = my_colours_MELSM_pi) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  ggtitle("Power") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "", y = "")

MELSM_pi_fig <- ggarrange(s1_pi_1, s1_pi_2, s1_pi_3, s1_pi_4, s1_pi_5,
                          ncol = 5, nrow = 1,
                          common.legend = TRUE, legend = "bottom")

MELSM_pi_fig


# mean estimate
s1_theta_1 <- ggplot(MELSM_data[which(MELSM_data$Metric=="Mean estimate"),], aes(x = `Scenario/true value`, y = Value, 
                                                                                 color = Study, shape = Study, 
                                                                                 group = Study)) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(data = baseline_df_theta[1:3,], 
             aes(yintercept = baseline), 
             linetype = "dashed", color = "black") +
  scale_color_manual(values = my_colours) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.4)) +
  ggtitle("Mean Estimate") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "", y = "")

# bias
s1_theta_2 <- ggplot(MELSM_data[which(MELSM_data$Metric=="Bias"),], aes(x = `Scenario/true value`, y = Value, 
                                                                        color = Study, shape = Study, 
                                                                        group = Study)) +
  geom_point(size = 3) +
  geom_line() + 
  geom_hline(data = baseline_df_theta[which(baseline_df_theta$Metric=="Bias"),], 
             aes(yintercept = baseline), 
             linetype = "dashed", color = "black") +
  scale_color_manual(values = my_colours) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, by = 0.2)) +
  ggtitle("Bias") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "", y = "")

# mean squared error
s1_theta_3 <- ggplot(MELSM_data[which(MELSM_data$Metric=="Mean Squared Error"),], aes(x = `Scenario/true value`, y = Value, 
                                                                                      color = Study, shape = Study, 
                                                                                      group = Study)) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(data = baseline_df_theta[which(baseline_df_theta$Metric=="Mean Squared Error"),], 
             aes(yintercept = baseline), 
             linetype = "dashed", color = "black") +
  scale_color_manual(values = my_colours) +
  scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by = 0.02)) +
  ggtitle("Mean Squared Error") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = expression(paste(theta, "_101")), y = "")

# bias-adjusted coverage
s1_theta_4 <- ggplot(MELSM_data[which(MELSM_data$Metric=="Bias-adjusted coverage"),], aes(x = `Scenario/true value`, y = Value, 
                                                                                          color = Study, shape = Study, 
                                                                                          group = Study)) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(data = baseline_df_theta[which(baseline_df_theta$Metric=="Bias-adjusted coverage"),], 
             aes(yintercept = baseline), 
             linetype = "dashed", color = "black") +
  scale_color_manual(values = my_colours) +
  scale_y_continuous(limits = c(0.8, 1), breaks = seq(0.8, 1, by = 0.04)) +
  ggtitle("Bias-adjusted coverage") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "", y = "")

# power
s1_theta_5 <- ggplot(MELSM_data[which(MELSM_data$Metric=="Power"),], aes(x = `Scenario/true value`, y = Value, 
                                                                         color = Study, shape = Study, 
                                                                         group = Study)) +
  geom_point(size = 3) +
  geom_line() + 
  geom_hline(data = baseline_df_theta[which(baseline_df_theta$Metric=="Power"),], 
             aes(yintercept = baseline), 
             linetype = "dashed", color = "black") +
  scale_color_manual(values = my_colours) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  ggtitle("Power") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "", y = "")

theta_fig <- ggarrange(s1_theta_1, s1_theta_2, s1_theta_3, s1_theta_4, s1_theta_5,
                       ncol = 5, nrow = 1,
                       common.legend = TRUE, legend = "bottom")

theta_fig

pi_est_fig <- ggarrange(pi_1, s1_pi_1, s1_theta_1,
                       ncol = 3, nrow = 1,
                       common.legend = TRUE, legend = "bottom")

pi_est_fig

## =====================================================================
## Figures
## =====================================================================

# set theme
theme_set(theme_bw())

theme_panel <- theme(
  plot.title         = element_text(hjust = 0.5, size = 11, face = "bold"),
  panel.grid.major.x = element_blank()
)

# colour palettes
cols_model <- c(
  "MLM"   = "#365D8DFF",
  "MELSM" = "#47C16EFF"
)

cols_modelStudy <- c(
  "MLM - Study 1"   = "#440154FF",
  "MELSM - Study 1" = "#365D8DFF",
  "MLM - Study 2"   = "#47C16EFF",
  "MELSM - Study 2" = "#FDE725FF"
)

cols_study <- c(
  "1" = "#365D8DFF",
  "2" = "#47C16EFF"
)

## =====================================================================
## π_101
## =====================================================================

performance_metrics_pi <- read_xlsx("Performance_Metrics_3.xlsx", sheet = "pi_101")

plot_df_pi <- performance_metrics_pi %>%
  pivot_longer(
    cols = c("Mean estimate", "Bias", "MSE", Coverage, Power),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("Mean estimate", "Bias", "MSE", "Coverage", "Power"),
      labels = c("Mean estimate", "Bias", "Mean Squared Error", "Coverage", "Power")
    ),
    DGP = factor(
      DGP,
      levels = c("MLM", "MELSM"),
      labels = c("Homoscedastic (MLM DGP)", "Heteroscedastic (MELSM DGP)")
    )
  )

# Split by DGP
MLM_data      <- plot_df_pi %>% filter(DGP == "Homoscedastic (MLM DGP)")
MELSM_pi_data <- plot_df_pi %>% filter(DGP == "Heteroscedastic (MELSM DGP)")

MLM_data$Model      <- factor(MLM_data$Model, levels = c("MLM", "MELSM"))
MELSM_pi_data$Model <- factor(MELSM_pi_data$Model, levels = c("MLM", "MELSM"))

# scenario factor
MELSM_pi_data$Scenario <- factor(
  MELSM_pi_data$Scenario,
  levels = c("𝜃_101 = -0.5", "𝜃_101 = 0", "𝜃_101 = 0.5"),
  labels = c("-0.5", "0", "0.5")
)

# model × study for heteroscedastic panels
MELSM_pi_data$ModelStudy <- interaction(MELSM_pi_data$Model, MELSM_pi_data$Study)
MELSM_pi_data$ModelStudy <- factor(
  MELSM_pi_data$ModelStudy,
  levels = c("MLM.1", "MELSM.1", "MLM.2", "MELSM.2"),
  labels = c("MLM - Study 1", "MELSM - Study 1", "MLM - Study 2", "MELSM - Study 2")
)

# baseline reference lines
baseline_pi <- tibble(
  Metric   = c("Mean estimate", "Bias", "Mean Squared Error", "Coverage", "Power"),
  baseline = c(-1.0, 0.0, 0.0, 0.95, 1.0)
)

baseline_pi$Metric <- factor(
  baseline_pi$Metric,
  levels = c("Mean estimate", "Bias", "Mean Squared Error", "Coverage", "Power")
)

# π_101 mean estimate

pi_me_mlm <- ggplot(
  MLM_data %>% filter(Metric == "Mean estimate"),
  aes(x = factor(Study), y = Value,
      colour = Model, shape = Model, group = Model)
) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(
    data = baseline_pi %>% filter(Metric == "Mean estimate"),
    aes(yintercept = baseline),
    linetype = "dashed", colour = "black"
  ) +
  scale_color_manual(values = cols_model) +
  scale_y_continuous(limits = c(-1.02, -0.98), breaks = seq(-1.02, -0.98, by = 0.01)) +
  labs(x = "Study", y = "", colour = NULL, shape = NULL) +
  ggtitle("Homoscedastic (MLM DGP)") +
  theme_panel
pi_me_melsm <- ggplot(
  MELSM_pi_data %>% filter(Metric == "Mean estimate"),
  aes(x = Scenario, y = Value,
      colour = ModelStudy, shape = ModelStudy, group = ModelStudy)
) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(
    data = baseline_pi %>% filter(Metric == "Mean estimate"),
    aes(yintercept = baseline),
    linetype = "dashed", colour = "black"
  ) +
  scale_color_manual(
    values = cols_modelStudy,
    labels = c("MLM\nStudy 1", "MELSM\nStudy 1", "MLM\nStudy 2", "MELSM\nStudy 2"),
    guide = guide_legend(
      title = NULL,
      override.aes = list(shape = c(16, 17, 15, 18))
    )) +
  scale_shape_manual(
    values = c(16, 17, 15, 18),
    guide = "none"
  ) +
  scale_y_continuous(limits = c(-1.02, -0.98), breaks = seq(-1.02, -0.98, by = 0.01)) +
  labs(x = expression(theta[101]), y = "") +
  ggtitle("Heteroscedastic (MELSM DGP)") +
  theme_panel

fig_pi_mean <- ggarrange(
  pi_me_mlm, pi_me_melsm,
  ncol = 2, nrow = 1,
  legend = "bottom",
  labels = c("A", "B")
)

# π_101 bias

pi_bias_mlm <- ggplot(
  MLM_data %>% filter(Metric == "Bias"),
  aes(x = factor(Study), y = Value,
      colour = Model, shape = Model, group = Model)
) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(
    data = baseline_pi %>% filter(Metric == "Bias"),
    aes(yintercept = baseline),
    linetype = "dashed", colour = "black"
  ) +
  scale_color_manual(values = cols_model) +
  scale_y_continuous(limits = c(-0.1, 0.1), breaks = seq(-0.1, 0.1, 0.05)) +
  labs(x = "Study", y = "", colour = NULL, shape = NULL) +
  ggtitle("Homoscedastic (MLM DGP)") +
  theme_panel

pi_bias_melsm <- ggplot(
  MELSM_pi_data %>% filter(Metric == "Bias"),
  aes(x = Scenario, y = Value,
      colour = ModelStudy, shape = ModelStudy, group = ModelStudy)
) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(
    data = baseline_pi %>% filter(Metric == "Bias"),
    aes(yintercept = baseline),
    linetype = "dashed", colour = "black"
  ) +
  scale_color_manual(
    values = cols_modelStudy,
    labels = c("MLM\nStudy 1", "MELSM\nStudy 1", "MLM\nStudy 2", "MELSM\nStudy 2"),
    guide = guide_legend(
      title = NULL,
      override.aes = list(shape = c(16, 17, 15, 18))
    )) +
  scale_shape_manual(
    values = c(16, 17, 15, 18),
    guide = "none"
  ) +
  scale_y_continuous(limits = c(-0.1, 0.1), breaks = seq(-0.1, 0.1, 0.05)) +
  labs(x = expression(theta[101]), y = "",
       colour = "Model & Study", shape = "Model & Study") +
  ggtitle("Heteroscedastic (MELSM DGP)") +
  theme_panel

fig_pi_bias <- ggarrange(
  pi_bias_mlm, pi_bias_melsm,
  ncol = 2, nrow = 1,
  legend = "bottom",
  labels = c("A", "B")
)

# π_101 MSE

pi_mse_mlm <- ggplot(
  MLM_data %>% filter(Metric == "Mean Squared Error"),
  aes(x = factor(Study), y = Value,
      colour = Model, shape = Model, group = Model)
) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(
    data = baseline_pi %>% filter(Metric == "Mean Squared Error"),
    aes(yintercept = baseline),
    linetype = "dashed", colour = "black"
  ) +
  scale_color_manual(values = cols_model) +
  scale_y_continuous(limits = c(0, 0.04), breaks = seq(0, 0.04, 0.01)) +
  labs(x = "Study", y = "", colour = NULL, shape = NULL) +
  ggtitle("Homoscedastic (MLM DGP)") +
  theme_panel

pi_mse_melsm <- ggplot(
  MELSM_pi_data %>% filter(Metric == "Mean Squared Error"),
  aes(x = Scenario, y = Value,
      colour = ModelStudy, shape = ModelStudy, group = ModelStudy)
) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(
    data = baseline_pi %>% filter(Metric == "Mean Squared Error"),
    aes(yintercept = baseline),
    linetype = "dashed", colour = "black"
  ) +
  scale_color_manual(
    values = cols_modelStudy,
    labels = c("MLM\nStudy 1", "MELSM\nStudy 1", "MLM\nStudy 2", "MELSM\nStudy 2"),
    guide = guide_legend(
      title = NULL,
      override.aes = list(shape = c(16, 17, 15, 18))
    )) +
  scale_shape_manual(
    values = c(16, 17, 15, 18),
    guide = "none"
  ) +
  scale_y_continuous(limits = c(0, 0.04), breaks = seq(0, 0.04, 0.01)) +
  labs(x = expression(theta[101]), y = "",
       colour = "Model & Study", shape = "Model & Study") +
  ggtitle("Heteroscedastic (MELSM DGP)") +
  theme_panel

fig_pi_mse <- ggarrange(
  pi_mse_mlm, pi_mse_melsm,
  ncol = 2, nrow = 1,
  legend = "bottom",
  labels = c("A", "B")
)

# π_101 coverage

pi_cov_mlm <- ggplot(
  MLM_data %>% filter(Metric == "Coverage"),
  aes(x = factor(Study), y = Value,
      colour = Model, shape = Model, group = Model)
) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(
    data = baseline_pi %>% filter(Metric == "Coverage"),
    aes(yintercept = baseline),
    linetype = "dashed", colour = "black"
  ) +
  scale_color_manual(values = cols_model) +
  scale_y_continuous(limits = c(0.88, 1.0), breaks = seq(0.88, 1.0, 0.02)) +
  labs(x = "Study", y = "", colour = NULL, shape = NULL) +
  ggtitle("Homoscedastic (MLM DGP)") +
  theme_panel

pi_cov_melsm <- ggplot(
  MELSM_pi_data %>% filter(Metric == "Coverage"),
  aes(x = Scenario, y = Value,
      colour = ModelStudy, shape = ModelStudy, group = ModelStudy)
) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(
    data = baseline_pi %>% filter(Metric == "Coverage"),
    aes(yintercept = baseline),
    linetype = "dashed", colour = "black"
  ) +
  scale_color_manual(
    values = cols_modelStudy,
    labels = c("MLM\nStudy 1", "MELSM\nStudy 1", "MLM\nStudy 2", "MELSM\nStudy 2"),
    guide = guide_legend(
      title = NULL,
      override.aes = list(shape = c(16, 17, 15, 18))
    )) +
  scale_shape_manual(
    values = c(16, 17, 15, 18),
    guide = "none"
  ) +
  scale_y_continuous(limits = c(0.88, 1.0), breaks = seq(0.88, 1.0, 0.02)) +
  labs(x = expression(theta[101]), y = "",
       colour = "Model & Study", shape = "Model & Study") +
  ggtitle("Heteroscedastic (MELSM DGP)") +
  theme_panel

fig_pi_cov <- ggarrange(
  pi_cov_mlm, pi_cov_melsm,
  ncol = 2, nrow = 1,
  legend = "bottom",
  labels = c("A", "B")
)

# π_101 power

pi_pow_mlm <- ggplot(
  MLM_data %>% filter(Metric == "Power"),
  aes(x = factor(Study), y = Value,
      colour = Model, shape = Model, group = Model)
) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(
    data = baseline_pi %>% filter(Metric == "Power"),
    aes(yintercept = baseline),
    linetype = "dashed", colour = "black"
  ) +
  scale_color_manual(values = cols_model) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, 0.2)) +
  labs(x = "Study", y = "", colour = NULL, shape = NULL) +
  ggtitle("Homoscedastic (MLM DGP)") +
  theme_panel

pi_pow_melsm <- ggplot(
  MELSM_pi_data %>% filter(Metric == "Power"),
  aes(x = Scenario, y = Value,
      colour = ModelStudy, shape = ModelStudy, group = ModelStudy)
) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(
    data = baseline_pi %>% filter(Metric == "Power"),
    aes(yintercept = baseline),
    linetype = "dashed", colour = "black"
  ) +
  scale_color_manual(
    values = cols_modelStudy,
    labels = c("MLM\nStudy 1", "MELSM\nStudy 1", "MLM\nStudy 2", "MELSM\nStudy 2"),
    guide = guide_legend(
      title = NULL,
      override.aes = list(shape = c(16, 17, 15, 18))
    )) +
  scale_shape_manual(
    values = c(16, 17, 15, 18),
    guide = "none"
  ) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, 0.2)) +
  labs(x = expression(theta[101]), y = "",
       colour = "Model & Study", shape = "Model & Study") +
  ggtitle("Heteroscedastic (MELSM DGP)") +
  theme_panel

fig_pi_power <- ggarrange(
  pi_pow_mlm, pi_pow_melsm,
  ncol = 2, nrow = 1,
  legend = "bottom",
  labels = c("A", "B")
)

## =====================================================================
## θ_101 
## =====================================================================

performance_metrics_theta <- read_xlsx("Performance_Metrics_3.xlsx", sheet = "theta_101")

plot_df_theta <- performance_metrics_theta %>%
  pivot_longer(
    cols = c(`Mean estimate`, Bias, MSE, `Bias-adjusted coverage`, Power_bayes),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("Mean estimate", "Bias", "MSE", "Bias-adjusted coverage", "Power_bayes"),
      labels = c("Mean estimate", "Bias", "Mean Squared Error",
                 "Bias-adjusted coverage", "Power")
    ),
    `Scenario/true value` = factor(
      `Scenario/true value`,
      levels = c("𝜃_101 = -0.5", "𝜃_101 = 0", "𝜃_101 = 0.5"),
      labels = c("-0.5", "0", "0.5")
    )
  )

# restrict to MELSM DGP
MELSM_theta <- plot_df_theta %>% filter(DGP == "MELSM")
MELSM_theta$Study <- factor(MELSM_theta$Study)

# baselines for θ_101
baseline_theta <- tibble(
  Metric   = c("Bias", "Mean Squared Error", "Bias-adjusted coverage", "Power"),
  baseline = c(0.0, 0.0, 0.95, 1.0)
)

# θ_101 performance figure (2×2)

theta_bias <- ggplot(
  MELSM_theta %>% filter(Metric == "Bias"),
  aes(x = `Scenario/true value`, y = Value,
      colour = Study, shape = Study, group = Study)
) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(
    data = baseline_theta %>% filter(Metric == "Bias"),
    aes(yintercept = baseline),
    linetype = "dashed", colour = "black"
  ) +
  scale_color_manual(values = cols_study) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, 0.2)) +
  labs(x = expression(theta[101]), y = "", colour = "Study", shape = "Study") +
  ggtitle("Bias") +
  theme_panel

theta_mse <- ggplot(
  MELSM_theta %>% filter(Metric == "Mean Squared Error"),
  aes(x = `Scenario/true value`, y = Value,
      colour = Study, shape = Study, group = Study)
) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(
    data = baseline_theta %>% filter(Metric == "Mean Squared Error"),
    aes(yintercept = baseline),
    linetype = "dashed", colour = "black"
  ) +
  scale_color_manual(values = cols_study) +
  scale_y_continuous(limits = c(0, 0.10), breaks = seq(0, 0.10, 0.02)) +
  labs(x = expression(theta[101]), y = "", colour = "Study", shape = "Study") +
  ggtitle("Mean Squared Error") +
  theme_panel

theta_cov <- ggplot(
  MELSM_theta %>% filter(Metric == "Bias-adjusted coverage"),
  aes(x = `Scenario/true value`, y = Value,
      colour = Study, shape = Study, group = Study)
) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(
    data = baseline_theta %>% filter(Metric == "Bias-adjusted coverage"),
    aes(yintercept = baseline),
    linetype = "dashed", colour = "black"
  ) +
  scale_color_manual(values = cols_study) +
  scale_y_continuous(limits = c(0.88, 1.0), breaks = seq(0.88, 1.0, 0.02)) +
  labs(x = expression(theta[101]), y = "", colour = "Study", shape = "Study") +
  ggtitle("Bias-adjusted coverage") +
  theme_panel

theta_power <- ggplot(
  MELSM_theta %>% filter(Metric == "Power"),
  aes(x = `Scenario/true value`, y = Value,
      colour = Study, shape = Study, group = Study)
) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(
    data = baseline_theta %>% filter(Metric == "Power"),
    aes(yintercept = baseline),
    linetype = "dashed", colour = "black"
  ) +
  scale_color_manual(values = cols_study) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, 0.2)) +
  labs(x = expression(theta[101]), y = "", colour = "Study", shape = "Study") +
  ggtitle("Power") +
  theme_panel

theta_fig <- ggarrange(
  theta_bias, theta_mse,
  theta_cov,  theta_power,
  ncol = 2, nrow = 2,
  common.legend = TRUE, legend = "bottom",
  labels = c("A", "B", "C", "D")
)

# θ_101 mean estimates

theta_mean_df <- MELSM_theta %>%
  filter(Metric == "Mean estimate") %>%
  mutate(theta_true = as.numeric(as.character(`Scenario/true value`)))

fig_theta_mean <- ggplot(
  theta_mean_df,
  aes(x = theta_true, y = Value,
      colour = Study, shape = Study, group = Study)
) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey40") +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(values = cols_study) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.5, 0, 0.5)) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks = c(-0.5, 0, 0.5)) +
  labs(
    x = expression(theta[101]),
    y = "",
    colour = "Study", shape = "Study"
  ) +
  theme_panel

# save all figs

fig_pi_mean
fig_pi_bias
fig_pi_mse
fig_pi_cov
fig_pi_power
theta_fig
fig_theta_mean

###         ###
### METRICS ###
###         ###

options(scipen = 999)
metrics_df <- study_1_MELSM %>%
  mutate(true_value = case_when(
    param == "pi101"    ~ -1.0,           
    param == "theta101" ~ theta_101, # comment out if assessing MLM-generated data only
    TRUE ~ NA_real_
  ))

compute_metrics <- function(df) {
  est <- df$estimate
  ci_l <- df$ci_lower
  ci_u <- df$ci_upper
  pvals <- df$p_value
  tv <- unique(df$true_value)
  se <- df$se
  
  estimate <- mean(est)
  
  # Bias
  bias <- mean(est - tv, na.rm = TRUE)
  rel_bias <- if (!is.na(tv) && tv != 0) mean((est - tv)/tv, na.rm = TRUE) else NA_real_
  
  # Mean Squared Error
  mse <- mean((est - tv)^2, na.rm = TRUE)
  
  # Classical coverage
  coverage_classical <- mean(ci_l <= tv & ci_u >= tv, na.rm = TRUE)
  
  # Bias-adjusted coverage: shift CI around mean estimate
  ci_l_adj <- est - 1.96 * se  # lower CI
  ci_u_adj <- est + 1.96 * se  # upper CI
  coverage_bias_adj <- mean(ci_l_adj <= mean(est, na.rm = TRUE) & ci_u_adj >= mean(est, na.rm = TRUE), na.rm = TRUE)
  
  # Power
  power <- mean(pvals < 0.05 & !is.na(pvals), na.rm = TRUE)
  power_bayes <- mean((ci_l > 0) | (ci_u < 0), na.rm = TRUE)
  
  # variance
  variance <- exp(mean(est))
  
  tibble(
    Estimate = estimate,
    Bias = bias,
    RelBias = rel_bias,
    MSE = mse,
    Coverage = coverage_classical,
    Coverage_BiasAdj = coverage_bias_adj,
    Power = power,
    Power_b = power_bayes,
    Variance = variance
  )
}

# for MELSM
metrics_summary <- metrics_df %>%
  group_by(theta_101, analysis_model, param) %>%
  summarise(compute_metrics(cur_data()), .groups = "drop") %>%
  mutate(across(c(Bias, RelBias, MSE, Coverage, Coverage_BiasAdj, Power, Variance), ~ round(.x, 4))) %>%
  arrange(analysis_model, param, theta_101)

# # for MLM
# metrics_summary <- results_params_df %>%
#   group_by(analysis_model, param) %>%
#   summarise(compute_metrics(cur_data()), .groups = "drop") %>%
#   mutate(across(c(Bias, RelBias, MSE, Coverage, Coverage_BiasAdj, Power, Variance), ~ round(.x, 4))) %>%
#   arrange(analysis_model, param)

metrics_summary
