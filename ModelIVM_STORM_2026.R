#setwd("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion")
setwd("/Users/crid/Desktop/D_Franck/Plots")
rm(list = ls())#remet à jour tout l'espace de travail

reset=function() {
  par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
  plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
}

library(xlsx)
library(readxl)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(stringr)
# DRC model to estimate dose response link ---- 
library(drc)
library(patchwork) # load after GAM, package incompatibility  #pour montage figure
library(openxlsx)

library(rstan)

library(future)
library(progressr)
library(future.apply)

handlers(global = TRUE)
plan(multisession, workers = 9)  # tous les cœurs sauf 2

# Stan configuration for performance
rstan_options(auto_write = TRUE)
options(mc.cores = 9)

trapz=function(x=NULL,y){
  if(is.null(x)){
    x=1:length(y)
  }
  if(length(y)==1){z=0}
  else{z=caTools::trapz(x,y)}
  return(z)
}

P_tau.fun = function(tau, tau50, alpha, Dmax) {
  p_tau = Dmax / (1 + (tau / tau50)^alpha)
  return(p_tau)
}

residual_fun = function(model, formulation_name) {
  pred = fitted(model)
  resid = residuals(model)  # drm calcule directement les résidus
  
  data.frame(
    fitted = pred,
    residuals = resid,
    Formulation = formulation_name
  )
}

# Fonction RMSE Root Mean Squared Error
rmse_fun = function(model, formulation_name) {
  pred = fitted(model)
  resid = residuals(model)
  rmse = sqrt(mean(resid^2))
  data.frame(Strategy = formulation_name, RMSE = rmse)
}

#Stan code
{
  # STAN CODE: Binomial model (for strategies 2, 3, 4)
  stan_code_binomial = "
data {
  int<lower=1> N;                    // Number of observations
  vector<lower=0>[N] DAI;            // Days after injection
  array[N] int<lower=0> n_dead;      // Number of dead

  array[N] int<lower=0> n_total;     // Total number
  real<lower=0, upper=1> Dmax;       // Upper limit (fixed)
  
  // For prediction
  int<lower=1> N_pred;               // Number of prediction points
  vector<lower=0>[N_pred] DAI_pred;  // Prediction points
}

parameters {
  real<lower=0> tau50;               // ED50 (time where mortality = Dmax/2)
  real<lower=0> alpha;               // Slope
}

transformed parameters {
  vector<lower=0, upper=1>[N] p_mort;
  
  for (i in 1:N) {
    p_mort[i] = Dmax / (1.0 + pow(DAI[i] / tau50, alpha));
  }
}

model {
  // Weakly informative priors
  tau50 ~ lognormal(3, 1.5);
  alpha ~ lognormal(0.5, 1);
  
  // Binomial likelihood
  n_dead ~ binomial(n_total, p_mort);
}

generated quantities {
  vector[N_pred] p_pred;
  
  for (i in 1:N_pred) {
    p_pred[i] = Dmax / (1.0 + pow(DAI_pred[i] / tau50, alpha));
  }
}
"

# STAN CODE: Beta-Binomial model (for strategies 0, 1, 5 with overdispersion)

stan_code_betabinomial = "
data {
  int<lower=1> N;                    // Number of observations
  vector<lower=0>[N] DAI;            // Days after injection
  array[N] int<lower=0> n_dead;      // Number of dead
  array[N] int<lower=0> n_total;     // Total number
  real<lower=0, upper=1> Dmax;       // Upper limit (fixed)
  
  // For prediction
  int<lower=1> N_pred;               // Number of prediction points
  vector<lower=0>[N_pred] DAI_pred;  // Prediction points
}

parameters {
  real<lower=0> tau50;               // ED50
  real<lower=0> alpha;               // Slope
  real<lower=1> kappa;               // Overdispersion parameter (concentration)
}

transformed parameters {
  vector<lower=0, upper=1>[N] p_mort;
  
  for (i in 1:N) {
    p_mort[i] = Dmax / (1.0 + pow(DAI[i] / tau50, alpha));
  }
}

model {
  // Priors
  tau50 ~ lognormal(3, 1.5);
  alpha ~ lognormal(0.5, 1);
  kappa ~ pareto(1.0, 1.5);          // Overdispersion prior
  
  // Beta-Binomial likelihood
  // Beta-Binomial(n | N, a, b) where a = kappa*p, b = kappa*(1-p)
  for (i in 1:N) {
    real a_param = kappa * p_mort[i];
    real b_param = kappa * (1.0 - p_mort[i]);
    n_dead[i] ~ beta_binomial(n_total[i], a_param, b_param);
  }
}

generated quantities {
  vector[N_pred] p_pred;
  // Also generate prediction intervals accounting for overdispersion
  array[N_pred] int n_pred_sim;
  
  for (i in 1:N_pred) {
    real p_i = Dmax / (1.0 + pow(DAI_pred[i] / tau50, alpha));
    p_pred[i] = p_i;
    // Simulate from beta-binomial for prediction intervals
    real p_draw = beta_rng(kappa * p_i + 0.001, kappa * (1.0 - p_i) + 0.001);
    n_pred_sim[i] = binomial_rng(1000, p_draw);
  }
}
"
}

Estimate_MCMC_Ptau = function(strategy, n_chains = 8, n_iter = 6000,n_warmup = 1000,
                              CI_level = 0.95, DAI_pred = seq(0.1, 166, 1)) {
  
  use_betabinomial = strategy %in% c(0, 1, 5)
  likelihood_type = ifelse(use_betabinomial, "Beta-Binomial", "Binomial")
  
  # STEP 1: Get drm results
  drm_res = ModelEtau50Alpha(strategy)
  
  # Extract data
  data = drm_res$sub_PK; Dmax = drm_res$Dmax; model_drm = drm_res$model
  
  # Identify DAI column
  if ("DAI" %in% names(data)) {DAI = as.numeric(data$DAI) } else if ("DAI2" %in% names(data)) {
    DAI = as.numeric(data$DAI2)}
  prop_dead = data$prop_dead
  
  # Number of observations
  if ("ntot" %in% names(data)) {
    n_total = as.integer(data$ntot)
  } else if ("n_mosquitoes" %in% names(data)) {
    n_total = as.integer(data$n_mosquitoes)
  } else {
    n_total = rep(50L, length(DAI))
    cat("   Note: n_total estimated at 50 per observation\n")
  }
  n_dead = as.integer(round(prop_dead * n_total))
  
  # STEP 2: drm prediction (target ribbon)
  
  newdata_drm = data.frame(DAI = DAI_pred)
  if ("DAI2" %in% names(data)) {
    newdata_drm = data.frame(DAI2 = DAI_pred)
  }
  
  pred_drm = predict(model_drm, newdata = newdata_drm, interval = "confidence")
  
  if (is.matrix(pred_drm)) {
    ribbon_drm = data.frame(DAI = DAI_pred, Prediction = pred_drm[, "Prediction"],
                            Lower = pred_drm[, "Lower"], Upper = pred_drm[, "Upper"] )
  } else {
    ribbon_drm = data.frame(DAI = DAI_pred, Prediction = pred_drm, Lower = NA, Upper = NA )
  }
  
  # STEP 3: Stan MCMC fitting
  
  # Filter DAI = 0 (division issue)
  valid_idx = DAI > 0
  
  stan_data = list(N = sum(valid_idx),DAI = DAI[valid_idx],n_dead = n_dead[valid_idx],
                   n_total = n_total[valid_idx], Dmax = Dmax,N_pred = length(DAI_pred),
                   DAI_pred = DAI_pred)
  
  if (use_betabinomial) {
    stan_model_obj = stan_model(model_code = stan_code_betabinomial)
  } else {
    stan_model_obj = stan_model(model_code = stan_code_binomial)
  }
  
  # MCMC Sampling
  fit = sampling(stan_model_obj, data = stan_data, chains = n_chains, iter = n_iter,
                 warmup = n_warmup, control = list(adapt_delta = 0.95, max_treedepth = 12),
                 refresh = 500)
  samples = extract(fit)
  n_samples = length(samples$tau50)
  
  # STEP 4: Compute MCMC ribbon quantiles
  alpha_ci = (1 - CI_level) / 2
  p_pred_samples = samples$p_pred
  
  ribbon_mcmc = data.frame( DAI = DAI_pred,Prediction = colMeans(p_pred_samples),
                            Lower = apply(p_pred_samples, 2, quantile, probs = alpha_ci),
                            Upper = apply(p_pred_samples, 2, quantile, probs = 1 - alpha_ci))
  
  # STEP 5: Parameter summary from MCMC
  
  param_summary = data.frame(
    parameter = c("tau50", "alpha"),
    mean = c(mean(samples$tau50), mean(samples$alpha)),
    median = c(median(samples$tau50), median(samples$alpha)),
    lower = c(quantile(samples$tau50, alpha_ci), quantile(samples$alpha, alpha_ci)),
    upper = c(quantile(samples$tau50, 1 - alpha_ci), quantile(samples$alpha, 1 - alpha_ci)),
    sd = c(sd(samples$tau50), sd(samples$alpha))
  )
  
  cat(sprintf("\n   MCMC Parameters (%.0f%% CI):\n", CI_level*100))
  print(param_summary, row.names = FALSE)
  
  # STEP 6: Find exact tau50/alpha that produce ribbon bounds
  # For each prediction point, find which (tau50, alpha) sample gives the lower/upper
  # This is more precise than just using parameter quantiles
  
  # Method: For each DAI point, identify the samples that produce min/max predictions
  lower_params_by_point = matrix(NA, nrow = length(DAI_pred), ncol = 2)
  upper_params_by_point = matrix(NA, nrow = length(DAI_pred), ncol = 2)
  
  for (j in 1:length(DAI_pred)) {
    pred_at_j = p_pred_samples[, j]
    
    # Find sample index closest to lower quantile
    lower_target = quantile(pred_at_j, alpha_ci)
    idx_lower = which.min(abs(pred_at_j - lower_target))
    lower_params_by_point[j, ] = c(samples$tau50[idx_lower], samples$alpha[idx_lower])
    
    # Find sample index closest to upper quantile
    upper_target = quantile(pred_at_j, 1 - alpha_ci)
    idx_upper = which.min(abs(pred_at_j - upper_target))
    upper_params_by_point[j, ] = c(samples$tau50[idx_upper], samples$alpha[idx_upper])
  }
  
  # Average parameters that produce lower/upper bounds
  avg_lower_params = colMeans(lower_params_by_point)
  avg_upper_params = colMeans(upper_params_by_point)
  
  # cat("\n   Parameters that produce MCMC ribbon bounds (averaged):\n")
  # cat(sprintf("   For Lower bound: tau50 = %.4f, alpha = %.4f\n", 
  #             avg_lower_params[1], avg_lower_params[2]))
  # cat(sprintf("   For Upper bound: tau50 = %.4f, alpha = %.4f\n", 
  #             avg_upper_params[1], avg_upper_params[2]))
  
  # Verify these parameters reproduce the ribbon
  P_tau_lower_check = P_tau.fun(DAI_pred, avg_lower_params[1], avg_lower_params[2], Dmax)
  P_tau_upper_check = P_tau.fun(DAI_pred, avg_upper_params[1], avg_upper_params[2], Dmax)
  
  rmse_lower_mcmc = sqrt(mean((P_tau_lower_check - ribbon_mcmc$Lower)^2))
  rmse_upper_mcmc = sqrt(mean((P_tau_upper_check - ribbon_mcmc$Upper)^2))
  
  # cat(sprintf("\n   Verification RMSE (MCMC ribbon):\n"))
  # cat(sprintf("   Lower: %.6f\n", rmse_lower_mcmc))
  # cat(sprintf("   Upper: %.6f\n", rmse_upper_mcmc))
  
  # STEP 7: Comparison drm vs MCMC
  
  cat("\n>>> Step 7: Comparison drm vs MCMC...\n")
  
  cat(sprintf("   tau50 - drm:  %.3f [%.3f, %.3f]\n", 
              drm_res$tau50, drm_res$tau50_low, drm_res$tau50_up))
  cat(sprintf("   tau50 - MCMC: %.3f [%.3f, %.3f]\n",
              param_summary$median[1], param_summary$lower[1], param_summary$upper[1]))
  cat(sprintf("   alpha - drm:  %.3f [%.3f, %.3f]\n",
              abs(drm_res$alpha), drm_res$alpha_low, drm_res$alpha_up))
  cat(sprintf("   alpha - MCMC: %.3f [%.3f, %.3f]\n",
              param_summary$median[2], param_summary$lower[2], param_summary$upper[2]))
  
  # Central estimates from MCMC
  tau50 = param_summary$median[1]
  alpha = param_summary$median[2]
  
  # CI bounds from MCMC
  tau50_low = param_summary$lower[1]
  tau50_up = param_summary$upper[1]
  alpha_low = param_summary$lower[2]
  alpha_up = param_summary$upper[2]
  
  predprop = ribbon_mcmc
  names(predprop) = c("DAI", "Prediction", "Lower", "Upper")
  
  # Exact parameters for ribbon bounds
  ribbon_params = list(
    lower = list(tau50 = avg_lower_params[1], alpha = avg_lower_params[2]),
    upper = list(tau50 = avg_upper_params[1], alpha = avg_upper_params[2])
  )
  
  FunList = list("tau50" = tau50, "tau50_low" = tau50_low, "tau50_up" = tau50_up,"alpha" = alpha,
                 "alpha_low" = alpha_low,"alpha_up" = alpha_up, "Dmax" = Dmax, "predprop" = predprop,
                 "sub_PK" = data, "DAI" = drm_res$DAI,"model" = model_drm,"mcmc_fit" = fit,
                 "mcmc_samples" = samples, "param_summary" = param_summary, "ribbon_drm" = ribbon_drm,
                 "ribbon_mcmc" = ribbon_mcmc, "ribbon_params" = ribbon_params, 
                 "likelihood" = likelihood_type, "CI_level" = CI_level)
  return(FunList)
}

Estimate_Ptau_Optim = function(strategy,CI_level = 0.95, DAI_pred = seq(0.1, 166, 1),
                               optim_method = "L-BFGS-B") {
  
  # STEP 1: Get drm results
  drm_res = ModelEtau50Alpha(strategy)
  data = drm_res$sub_PK; Dmax = drm_res$Dmax; model_drm = drm_res$model
  if ("DAI" %in% names(data)) {
    DAI = as.numeric(data$DAI)
  } else if ("DAI2" %in% names(data)) {
    DAI = as.numeric(data$DAI2)
  } 
  
  # STEP 2: drm prediction (target ribbon)
  
  newdata_drm = data.frame(DAI = DAI_pred)
  if ("DAI2" %in% names(data)) { newdata_drm = data.frame(DAI2 = DAI_pred)}
  pred_drm = predict(model_drm, newdata = newdata_drm, interval = "confidence")
  
  if (is.matrix(pred_drm)) {
    ribbon_drm = data.frame(DAI = DAI_pred, Prediction = pred_drm[, "Prediction"], 
                            Lower = pred_drm[, "Lower"], Upper = pred_drm[, "Upper"])
  }
  
  # STEP 3: Optimization to find exact parameters for Lower bound
  # Objective function: minimize SSE (Sum of squared errors) between P_tau and ribbon_drm$Lower
  obj_lower = function(params) {
    tau50_try = params[1];alpha_try = params[2]
    # Constraints
    if (tau50_try <= 0 || alpha_try <= 0) return(1e10)
    p_try = P_tau.fun(DAI_pred, tau50_try, alpha_try, Dmax)
    # Sum of squared errors
    sse = sum((p_try - ribbon_drm$Lower)^2)
    return(sse)
  }
  
  # Multiple starting points to avoid local minima
  start_points_lower = list( c(drm_res$tau50_low, drm_res$alpha_low),
                             c(drm_res$tau50_low, drm_res$alpha_up),
                             c(drm_res$tau50_up, drm_res$alpha_low),
                             c(drm_res$tau50, drm_res$alpha_low),
                             c(drm_res$tau50_low, drm_res$alpha))
  best_lower = list(par = c(NA, NA), value = Inf)
  
  for (start in start_points_lower) {
    if (any(start <= 0)) next
    opt_result = tryCatch(optim(par = start, fn = obj_lower, method = optim_method,
                                lower = c(0.01, 0.01),upper = c(500, 50) ),
                          error = function(e) list(par = c(NA, NA), value = Inf))
    if (opt_result$value < best_lower$value) {
      best_lower = opt_result
    }
  }
  
  tau50_for_lower = best_lower$par[1];alpha_for_lower = best_lower$par[2]
  
  # cat(sprintf("   Lower bound params: tau50 = %.6f, alpha = %.6f\n",
  #             tau50_for_lower, alpha_for_lower))
  # cat(sprintf("   SSE = %.8f\n", best_lower$value))
  
  # STEP 4: Optimization to find exact parameters for Upper bound
  
  # Objective function: minimize SSE between P_tau and ribbon_drm$Upper
  obj_upper = function(params) {
    tau50_try = params[1];alpha_try = params[2]
    if (tau50_try <= 0 || alpha_try <= 0) return(1e10)
    p_try = P_tau.fun(DAI_pred, tau50_try, alpha_try, Dmax)
    sse = sum((p_try - ribbon_drm$Upper)^2)
    return(sse)
  }
  
  start_points_upper = list( c(drm_res$tau50_up, drm_res$alpha_up),
                             c(drm_res$tau50_up, drm_res$alpha_low),
                             c(drm_res$tau50_low, drm_res$alpha_up),
                             c(drm_res$tau50, drm_res$alpha_up),
                             c(drm_res$tau50_up, drm_res$alpha))
  best_upper = list(par = c(NA, NA), value = Inf)
  
  for (start in start_points_upper) {
    if (any(start <= 0)) next
    
    opt_result = tryCatch( optim(par = start, fn = obj_upper,method = optim_method,
                                 lower = c(0.01, 0.01), upper = c(500, 50) ),
                           error = function(e) list(par = c(NA, NA), value = Inf))
    
    if (opt_result$value < best_upper$value) {
      best_upper = opt_result
    }
  }
  tau50_for_upper = best_upper$par[1];alpha_for_upper = best_upper$par[2]
  
  # cat(sprintf("   Upper bound params: tau50 = %.6f, alpha = %.6f\n",
  #             tau50_for_upper, alpha_for_upper))
  # cat(sprintf("   SSE = %.8f\n", best_upper$value))
  # STEP 5: Verification - check that calibrated params reproduce ribbon
  
  P_tau_lower_calib = P_tau.fun(DAI_pred, tau50_for_lower, alpha_for_lower, Dmax)
  P_tau_upper_calib = P_tau.fun(DAI_pred, tau50_for_upper, alpha_for_upper, Dmax)
  
  rmse_lower = sqrt(mean((P_tau_lower_calib - ribbon_drm$Lower)^2))
  rmse_upper = sqrt(mean((P_tau_upper_calib - ribbon_drm$Upper)^2))
  
  max_error_lower = max(abs(P_tau_lower_calib - ribbon_drm$Lower))
  max_error_upper = max(abs(P_tau_upper_calib - ribbon_drm$Upper))
  
  # cat(sprintf("   Lower bound: RMSE = %.8f, Max error = %.8f\n", 
  #             rmse_lower, max_error_lower))
  # cat(sprintf("   Upper bound: RMSE = %.8f, Max error = %.8f\n", 
  #             rmse_upper, max_error_upper))
  
  # STEP 6: Build prediction ribbon using calibrated params
  
  predprop = data.frame(
    DAI = DAI_pred,
    Prediction = P_tau.fun(DAI_pred, drm_res$tau50, abs(drm_res$alpha), Dmax),
    Lower = P_tau_lower_calib, Upper = P_tau_upper_calib
  )
  
  # Calibrated parameters
  ribbon_params = list(lower = list(tau50 = tau50_for_lower, alpha = alpha_for_lower),
                       upper = list(tau50 = tau50_for_upper, alpha = alpha_for_upper))
  
  # STEP 7: Comparison with drm SE-based bounds
  cat("   drm (using tau50 +/- SE, alpha +/- SE):\n")
  cat(sprintf("      tau50: %.3f [%.3f, %.3f]\n", 
              drm_res$tau50, drm_res$tau50_low, drm_res$tau50_up))
  cat(sprintf("      alpha: %.3f [%.3f, %.3f]\n",
              abs(drm_res$alpha), drm_res$alpha_low, drm_res$alpha_up))
  
  cat("\n   Calibrated (exact params for ribbon):\n")
  cat(sprintf("      For Lower: tau50 = %.4f, alpha = %.4f\n",
              tau50_for_lower, alpha_for_lower))
  cat(sprintf("      For Upper: tau50 = %.4f, alpha = %.4f\n",
              tau50_for_upper, alpha_for_upper))
  
  FunList = list("predprop" = predprop, "sub_PK" = data, "DAI" = drm_res$DAI, "model" = model_drm,
                 "ribbon_drm" = ribbon_drm, "ribbon_params" = ribbon_params,
                 "calibration_rmse" = list(lower = rmse_lower, upper = rmse_upper),
                 "drm_bounds" = list(tau50 = drm_res$tau50,tau50_low = drm_res$tau50_low,
                                     tau50_up = drm_res$tau50_up,alpha = abs(drm_res$alpha),
                                     alpha_low = drm_res$alpha_low,alpha_up = drm_res$alpha_up),
                 "method" = "Optimization", "CI_level" = CI_level)
  return(FunList)
}

ModelEtau50Alpha=function(strategy){
  if(strategy==0){#Long lasting 0.6
    #fitting parameter for LAIF 0.6
    # load data ---- 
    # P4D_PK = readRDS("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion/dataIVM/IDR_tab4J_PK.rds")
    P4D_PK = readRDS("/Users/crid/Desktop/D_Franck/dataIVM/IDR_tab4J_PK.rds")
    # mortality data, proportion of dead at 4 days after exposition 
    # gam_all = readRDS("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion/dataIVM/pred_gam_all.rds")
    gam_all = readRDS("/Users/crid/Desktop/D_Franck/dataIVM/pred_gam_all.rds")
    CritereIVMformulation = "F31-0.6"
    sub_PK_F31 = P4D_PK %>% filter(strain == "vk5" & IVM_formulation == CritereIVMformulation)
    sub_PK_F31["DAI"] = as.numeric(sub_PK_F31$DAI)
    #fitting the formula
    SimulDrmLongLasting06= drm(prop_dead ~ DAI, data = sub_PK_F31,
                               fct = LL.4(fixed = c(NA, 0, 1, NA),
                                          names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                               type = "binomial")
    modeldrm=SimulDrmLongLasting06
    summary(SimulDrmLongLasting06)
    
    coef_summary = summary(SimulDrmLongLasting06)$coefficients
    Dmax = 1; tau50 = coef_summary[2, "Estimate"]; alpha = coef_summary[1, "Estimate"] 
    tau50_se = coef_summary[2, "Std. Error"]; alpha_se = coef_summary[1, "Std. Error"]
    
    
    gam_F31 = gam_all %>%
      filter(IVM_formulation == CritereIVMformulation) %>%
      mutate(DAI = Hour / 24, prop_dead = NA)
    Day_seq=length(unique(gam_F31$DAI))
    pred_F31 = predict(SimulDrmLongLasting06, newdata = gam_F31, interval = "confidence")
    predprop_F31 = data.frame(pred_F31) %>% cbind(gam_F31) %>%  mutate(Formulation = "F31-0.6")
    predprop_F31 = predprop_F31[, !names(predprop_F31) %in% 
                                  c("IVM_concentration", "IVM_formulation", "Hour", "prop_dead")]
    predprop=predprop_F31
    
    sub_PK_F31 = sub_PK_F31 %>%mutate(DAI = as.numeric(DAI))
    sub_PK_F31 = sub_PK_F31 %>% mutate(Day = as.numeric(Day))
    sub_PK=sub_PK_F31
    DAI = seq(0,166, 1)
  }
  
  if(strategy==1){#Long lasting 1
    #fitting parameter for LAIF 1
    # Only for Lamidi study 1 mg/kg
    # load data ---- 
    # 1. Data from study of Lamidi (Bapla area?) with cattle injected at 1.0 mg/kg  
    #lamidi_data = read_xlsx("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion/dataIVM/Survie_anopheles_sauvage110923.xlsx")
    lamidi_data = read_xlsx("/Users/crid/Desktop/D_Franck/dataIVM/Survie_anopheles_sauvage110923.xlsx")
    surv_gamb = lamidi_data %>% dplyr::filter(Espece=="gambiae" & Traitement=="IVM") %>% 
      mutate(DAI2=as.numeric(str_extract(DAI, "\\d+")), Status4J=if_else(Temps<=4, 1,0))
    
    # vizu 
    # estimation of proportion of dead mosquitoes at 4 days by cattle and DAI
    # join with GAM_pk estimation of IVM plasma concentration 
    sg_sumdead = surv_gamb %>% group_by(Bovin, DAI2) %>% dplyr::filter(Status4J==1) %>% 
      summarise(ndead=n())
    sg_sum = surv_gamb %>% group_by(Bovin, DAI2) %>%  summarise(ntot=n())
    sg_tot = inner_join(sg_sumdead, sg_sum) 
    sg_tot = sg_tot %>% mutate(prop_dead=ndead/ntot, Day=DAI2) %>% droplevels()
    
    #fitting the formula
    SimulDrmLongLasting1=drm(prop_dead~DAI2, data=sg_tot, fct=LL.4(fixed=c(NA, 0, 1, NA),
                                                                   names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), type="binomial") 
    modeldrm=SimulDrmLongLasting1
    summary(SimulDrmLongLasting1)
    
    coef_summary = summary(SimulDrmLongLasting1)$coefficients
    Dmax = 1; tau50 = coef_summary[2, "Estimate"]; alpha = coef_summary[1, "Estimate"] 
    tau50_se = coef_summary[2, "Std. Error"]; alpha_se = coef_summary[1, "Std. Error"]
    
    #pour faire des figures après?
    DAI = seq(0,166, 1)
    ndat = data.frame(DAI2=DAI, prop_dead=NA)
    pred_Pdead_4J = predict(SimulDrmLongLasting1, newdata = ndat, interval="confidence")
    
    pred_Pdead_4J = cbind(pred_Pdead_4J,DAI)
    
    predprop_F1.0 = data.frame(pred_Pdead_4J) %>% mutate(Formulation = "Lam-1.0")  # predprop_F1.0, 1.0 pour 1mg/kg
    predprop=predprop_F1.0
    
    sg_tot = sg_tot %>%mutate(DAI2 = as.numeric(DAI2))
    sg_tot = sg_tot %>% mutate(Day = as.numeric(Day))
    sub_PK= sg_tot
  }
  
  #fitting parameter for Bohemia and Rimdamal II
  # Oral_Formulation = read.csv("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion/DataIVM/Oral_formulationF.csv")
  Oral_Formulation = read.csv("/Users/crid/Desktop/D_Franck/dataIVM/Oral formulation/Oral_formulationF.csv")
  # View(Oral_Formulation)
  #gam_all = readRDS("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion/DataIVM/pred_gam_all.rds")
  gam_all = readRDS("/Users/crid/Desktop/D_Franck/dataIVM/pred_gam_all.rds")
  
  if(strategy==2){#Pour BOHEMIA
    BOHEMIA= Oral_Formulation %>%filter(Article == "BOHEMIA")
    BOHEMIA["DAI"] = as.numeric(BOHEMIA$DAI)
    SimulDrmBohemia= drm(prop_dead ~ DAI, data = BOHEMIA,
                         fct = LL.4(fixed = c(NA, 0, 0.8626, NA),
                                    names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
    modeldrm=SimulDrmBohemia
    summary(SimulDrmBohemia)
    coef_summary = summary(SimulDrmBohemia)$coefficients
    Dmax = 0.8626; tau50 = coef_summary[2, "Estimate"]; alpha = coef_summary[1, "Estimate"] 
    tau50_se = coef_summary[2, "Std. Error"]; alpha_se = coef_summary[1, "Std. Error"]
    
    #pour faire des figures après?
    #plot(EDM_BOHEMIA)
    DAI = seq(0,166, 1)
    ndat = data.frame(DAI=DAI, prop_dead=NA)
    
    pred_BOHEMIA = predict(SimulDrmBohemia, newdata = ndat, interval = "confidence", level = 0.45)
    
    pred_BOHEMIA = data.frame(pred_BOHEMIA) %>% cbind(ndat) %>% mutate(Formulation = "BOHEMIA")
    predprop=pred_BOHEMIA
    
    sub_PK=BOHEMIA %>% mutate(Formulation = "BOHEMIA")
  }
  
  if(strategy==3){#Pour RIMDAMAL smitETAl 2018
    RIMDAMAL= Oral_Formulation %>%filter(Article == "RIMDAMAL")
    RIMDAMAL["DAI"] = as.numeric(RIMDAMAL$DAI)
    SimulDrmRimdamal= drm(prop_dead ~ DAI, data = RIMDAMAL,
                          fct = LL.4(fixed = c(NA, 0, 0.8632, NA),
                                     names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
    modeldrm=SimulDrmRimdamal
    summary(SimulDrmRimdamal)
    
    coef_summary = summary(SimulDrmRimdamal)$coefficients
    Dmax = 0.8632; tau50 = coef_summary[2, "Estimate"]; alpha = coef_summary[1, "Estimate"] 
    tau50_se = coef_summary[2, "Std. Error"]; alpha_se = coef_summary[1, "Std. Error"]
    
    #plot(EDM_RIMDAMAL)
    DAI = seq(0,166, 1)
    ndat = data.frame(DAI=DAI, prop_dead=NA)
    
    pred_RIMDAMAL=predict(SimulDrmRimdamal,newdata=ndat,interval="confidence",level = 0.7)
    
    pred_RIMDAMAL = data.frame(pred_RIMDAMAL) %>% cbind(ndat) %>%  mutate(Formulation = "RIMDAMAL")
    predprop=pred_RIMDAMAL
    
    sub_PK=RIMDAMAL %>% mutate(Formulation = "RIMDAMAL")
  }
  
  if(strategy==4){#Pour KamauRIMDAMAL  Kamau et Al 2024
    KamauRIMDAMAL= Oral_Formulation %>%filter(Article == "KamauRIMDAMAL")
    KamauRIMDAMAL["DAI"] = as.numeric(KamauRIMDAMAL$DAI)
    SimulDrmKamauRIMDAMAL= drm(prop_dead ~ DAI, data = KamauRIMDAMAL,
                               fct = LL.4(fixed = c(NA, 0, 0.6744, NA),
                                          names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
    modeldrm=SimulDrmKamauRIMDAMAL
    summary(SimulDrmKamauRIMDAMAL)
    
    coef_summary = summary(SimulDrmKamauRIMDAMAL)$coefficients
    Dmax = 0.6744; tau50 = coef_summary[2, "Estimate"]; alpha = coef_summary[1, "Estimate"] 
    tau50_se = coef_summary[2, "Std. Error"]; alpha_se = coef_summary[1, "Std. Error"]
    
    #plot(EDM_KamauRIMDAMAL)
    DAI = seq(0,166, 1)
    ndat = data.frame(DAI=DAI, prop_dead=NA)
    
    pred_KamauRIMDAMAL=predict(SimulDrmKamauRIMDAMAL,newdata=ndat,interval="confidence",level = 0.3)
    
    pred_KamauRIMDAMAL = data.frame(pred_KamauRIMDAMAL) %>% cbind(ndat) %>% mutate(Formulation = "KamauRIMDAMAL")
    predprop=pred_KamauRIMDAMAL
    
    sub_PK=KamauRIMDAMAL %>% mutate(Formulation = "KamauRIMDAMAL")
  }
  
  if(strategy==5){#Long lasting 0.6 train Kis
    #fitting parameter for LAIF 0.6
    # load data ---- View(gam_all)
    # P4D_PK_kis = read.csv("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion/dataIVM/data_mortality_rate_4j_LAIF06_KIS.csv")
    P4D_PK_kis = read.csv("/Users/crid/Desktop/D_Franck/dataIVM/data_mortality_rate_4j_LAIF06_KIS.csv")
    # mortality data, proportion of dead at 4 days after exposition 
    # gam_all = readRDS("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion/dataIVM/pred_gam_all.rds")
    gam_all = readRDS("/Users/crid/Desktop/D_Franck/dataIVM/pred_gam_all.rds")
    CritereIVMformulation = "mdc-STM-001"
    sub_PK_F31 = P4D_PK_kis %>% filter(strain == "KIS" & IVM_formulation == CritereIVMformulation)
    sub_PK_F31["DAI"] = as.numeric(sub_PK_F31$DAI)
    #fitting the formula
    SimulDrmLongLasting06_Kis= drm(prop_dead ~ DAI, data = sub_PK_F31,
                                   fct = LL.4(fixed = c(NA, 0, 1, NA),
                                              names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                                   type = "binomial")
    modeldrm=SimulDrmLongLasting06_Kis
    summary(SimulDrmLongLasting06_Kis)
    
    coef_summary = summary(SimulDrmLongLasting06_Kis)$coefficients
    Dmax = 1; tau50 = coef_summary[2, "Estimate"]; alpha = coef_summary[1, "Estimate"] 
    tau50_se = coef_summary[2, "Std. Error"]; alpha_se = coef_summary[1, "Std. Error"]
    
    gam_F31 = gam_all %>% filter(IVM_formulation == "F31-0.6") %>% mutate(DAI = Hour / 24, prop_dead = NA)
    Day_seq=length(unique(gam_F31$DAI))
    pred_F31 = predict(SimulDrmLongLasting06_Kis, newdata = gam_F31, interval = "confidence")
    predprop_F31 = data.frame(pred_F31) %>% cbind(gam_F31) %>% mutate(Formulation = "F31-0.6")
    predprop_F31 = predprop_F31[, !names(predprop_F31) %in% 
                                  c("IVM_concentration", "IVM_formulation", "Hour", "prop_dead")]
    predprop=predprop_F31
    
    sub_PK_F31 = sub_PK_F31 %>%mutate(DAI = as.numeric(DAI))
    # sub_PK_F31 = sub_PK_F31 %>% mutate(Day = as.numeric(Day))
    sub_PK=sub_PK_F31
    DAI = seq(0,166, 1)
  }
  
  alpha_low = alpha - alpha_se; alpha_up  = alpha + alpha_se 
  tau50_low = tau50 - tau50_se; tau50_up  = tau50 + tau50_se
  
  FunList = list("tau50" = tau50, "tau50_low" = tau50_low, "tau50_up" = tau50_up,  
                 "alpha_low" = alpha_low, "alpha_up" = alpha_up, "Dmax" = Dmax, 
                 "predprop" = predprop, "sub_PK" = sub_PK, "alpha" = alpha, "DAI" = DAI,
                 "model"=modeldrm)
  return(FunList)
}

gen_param_seq=function(estimate, lower, upper, n = 2) {
  seq_vals = seq(lower, upper, length.out = n)
  return(seq_vals)
}

ParamIvmFormulation = ModelEtau50Alpha(strategy = 0)

tau50_seq = gen_param_seq(ParamIvmFormulation$tau50, ParamIvmFormulation$tau50_low, 
                          ParamIvmFormulation$tau50_up, n = 5)
alpha_seq = gen_param_seq(ParamIvmFormulation$alpha, ParamIvmFormulation$alpha_low,
                          ParamIvmFormulation$alpha_up, n = 5)

show_CI=function(strategy) {
  P = ModelEtau50Alpha(strategy)
  
  cat("Strategy", strategy, "\n")
  cat(sprintf("tau50: %.3f [%.3f, %.3f]\n", P$tau50, P$tau50_low, P$tau50_up))
  cat(sprintf("alpha: %.3f [%.3f, %.3f]\n", P$alpha, P$alpha_low, P$alpha_up))
  cat(sprintf("Dmax:  %.4f\n\n", P$Dmax))
}

for(s in 0:5) show_CI(s)

calculate_CHR = function(tau_times, Vtau, ParamIvmFormulation, mu_m) {
  p_interp = approx(Vtau, P_tau.fun(Vtau, tau50=ParamIvmFormulation$tau50, 
                                    alpha=ParamIvmFormulation$alpha,
                                    Dmax=ParamIvmFormulation$Dmax), tau_times)$y
  
  mu_m_ivm = mu_m - log(1 - p_interp)
  CHR = mu_m_ivm / mu_m
  
  return(CHR)  
}

RhoFunction=function(Nah,Ntau,Vtau,Vah,IVM_field_dependancy,strategy, tau50, alpha){
  rho=matrix(0,nrow = Nah, ncol = Ntau)
  
  if (IVM_field_dependancy) {
    tau_eff = 0.5 * tau50
  } else {
    tau_eff = tau50
  }
  
  for (a in 1:Nah) {
    for (j in 1:Ntau) {
      if (Vtau[j] == 0) {
        rho[a, j] = 0  # Éviter 0^(alpha-1) quand alpha < 1
      } else {
        rho[a,j]=alpha * (Vtau[j] / tau_eff)^(alpha - 1) /
          (tau_eff * (1 + (Vtau[j] / tau_eff)^alpha)) #10*(Vtau[j]>tau_eff)
      }
    }
  }
  
  #On calcule les proba d'être encore sous IVM pour chaque group
  PropHuman_in_IVM_Group1=Vtau; IdAgeTemoin=3
  PropHuman_in_IVM_Group1[1]=1
  for (j in 2:Ntau) {
    Id=1:j; xId=Vtau[Id]; yId=rho[IdAgeTemoin,Id]
    PropHuman_in_IVM_Group1[j]=exp(-trapz(xId,yId))
  }
  
  FunList = list( "rho"=rho, "PropHuman_in_IVM_Group1"=PropHuman_in_IVM_Group1)
  return(FunList)
}

TheIvmStrategy=function(t_begin_Camp,Number_of_cycle,VectTime_between_cycles,
                        Dur_cycle,Gap,time,Ntime){
  
  phi=rep(0,Ntime)
  if (Number_of_cycle >= 1) {
    for (t in 1:Ntime) {
      for (k in 0:(Number_of_cycle - 1)) { 
        if (k == 0) { ta = t_begin_Camp}
        else {ta = t_begin_Camp + sum(Dur_cycle + VectTime_between_cycles[1:k]) }
        tb = ta + Dur_cycle - 1
        if (ta - Gap <= time[t] & time[t] < ta) { phi[t] = (time[t]-ta+Gap)/Gap } 
        else if (ta <= time[t] & time[t] <= tb) { phi[t] = 1 } 
        else if (tb < time[t] & time[t] <= tb + Gap){phi[t] = (tb + Gap - time[t])/Gap}
      }
    }
  }
  
  FunList = list("PhiIvm" = phi)
  return(FunList)
}

ModelIVM=function(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,mu_h,delta_h,
                  bar_beta_h,beta_h,gamma_h,theta,beta_m,Wedge_m,mu_m, ah_max,dah,Nah,tmax,
                  Gap,dt,time,Ntime,p_f, IVM_Pregnancy, dtau,tau_max,Vtau,Ntau,dsigma,sigma_max,Vsigma,
                  Nsigma,strategy,IVM_field_dependancy,tau50, alpha, Dmax){
  q=0.48 #proportion of human male
  k_h=5.5*10^(-4)#1.469e-2;
  Wedge_h=13250 # recrutment rate Bobo dioulasso in 2012
  
  p_vec = P_tau.fun(Vtau, tau50=tau50, alpha=alpha, Dmax=Dmax)
  mu_m.ivm = mu_m - log(1 - p_vec)
  
  ModelRho=RhoFunction(Nah,Ntau,Vtau,Vah,IVM_field_dependancy,strategy, tau50, alpha)
  rho=ModelRho$rho
  
  #phi(t,a)
  Modelphi=TheIvmStrategy(t_begin_Camp,Number_of_cycle,VectTime_between_cycles,
                          Dur_cycle,Gap,time,Ntime)
  phi=Modelphi$PhiIvm
  {
    mphi=matrix(0,nrow = Ntime,ncol = Nah);#IVM exposition rate for human male
    fphi=matrix(0,nrow = Ntime,ncol = Nah);#IVM exposition rate for human female
    
    phi=(PropIVM/Dur_cycle)*phi
    for (t in 1:Ntime) {
      for (a in 1:Nah) {
        if(Vah[a]<=5){mphi[t,a]=0; fphi[t,a]=0}
        else if (Vah[a] >= 15 & Vah[a] <= 45){
          if (IVM_Pregnancy) {mphi[t,a] = q * phi[t]; fphi[t,a] = (1 - q)*(1 - p_f)*phi[t]
          } else {
            mphi[t,a] = q * phi[t]; fphi[t,a] = 0
          }
        }
        else{mphi[t,a]=q*phi[t]; fphi[t,a]=(1-q)*phi[t]}
      }
    }
  }
  
  # mSh is Susceptible human male and fSh Susceptible human female
  mSh=fSh=matrix(0,nrow = Ntime,ncol = Nah); mAh=fAh=matrix(0,nrow = Ntime,ncol = Nah)
  mIh=fIh=matrix(0,nrow = Ntime,ncol = Nah); mRh=fRh=matrix(0,nrow = Ntime,ncol = Nah)
  Sm=rep(0,Ntime); Im=matrix(0,nrow = Ntime,ncol = Nsigma);
  Sm_ivm=matrix(0,nrow = Ntime,ncol = Ntau);Im_ivm=array(0, dim = c(Nsigma,Ntau,Ntime))
  mSh_ivm=fSh_ivm=array(0, dim = c(Nah,Ntau,Ntime)); mAh_ivm=fAh_ivm=array(0, dim = c(Nah,Ntau,Ntime))
  mIh_ivm=fIh_ivm=array(0, dim = c(Nah,Ntau,Ntime)); mRh_ivm=fRh_ivm=array(0, dim = c(Nah,Ntau,Ntime))
  Nh=rep(0,Ntime);Nm=rep(0,Ntime)
  
  
  # Initial values (t=0), Ref Quentin Richard et al Bobo Dioulasso 2012
  {
    Int_values = 8136.10 * (c(12.9, 12.5, 11.5, 13.1, 11.9, 9.3, 7.3, 5.5, 4.4,
                              3.2, 2.6, 1.8, 1.4, 0.9, 0.7, 0.3, 0.2, 0.2) / 0.997)
    # apparently the data only count for 99.7 % of the population
    # Age labels corresponding to each range in mSh_values
    age_label_3 = c("0 to 5 years", "5 to 10 years", "10 to 15 years", 
                    "15 to 20 years", "20 to 25 years", "25 to 30 years", 
                    "30 to 35 years", "35 to 40 years", "40 to 45 years", 
                    "45 to 50 years", "50 to 55 years", "55 to 60 years", 
                    "60 to 65 years", "65 to 70 years", "70 to 75 years", 
                    "75 to 80 years", "80 to 85 years", "85 to 90 years")
    
    for (i in seq_along(age_label_3)) {
      
      if (i == 1) {
        # For the first age range (0 to 5 years)
        age_indices=which(Vah <= 5)
      } else if (i == length(age_label_3)) {
        # For the last age range (> 85 years)
        age_indices=which(Vah > 85)
      } else {
        # For intermediate age ranges
        age_indices=which(Vah > (i - 1) * 5 & Vah <= i * 5)
      }
      
      # Calculate n_groups (number of individuals in this age group)
      if (length(age_indices) > 0) {
        n_groups=length(age_indices)
        
        mSh[1, age_indices] = 0.45 * q * Int_values[i] / n_groups
        fSh[1, age_indices] = 0.45 * (1-q) * Int_values[i] / n_groups
        mAh[1, age_indices] = 0.3 * q * Int_values[i] / n_groups
        fAh[1, age_indices] = 0.3 * (1-q) * Int_values[i] / n_groups
        mIh[1, age_indices] = 0.1 * q * Int_values[i] / n_groups
        fIh[1, age_indices] = 0.1 * (1-q) * Int_values[i] / n_groups
        mRh[1, age_indices] = 0.15 * q * Int_values[i] / n_groups
        fRh[1, age_indices] = 0.15 * (1-q) * Int_values[i] / n_groups
      }
    }
    
  }
  
  #nu_hm et nu_hf
  {
    nu_0=1/11; nu_1=1/170;nu_2=1/3
    nu_hm=rep(NA,Nah) # disease mortality for men
    nu_hf = rep(NA, Nah) # disease mortality for female
    for (i in 1:Nah) {
      if (Vah[i] <= 5) { nu_hm[i] = nu_2; nu_hf[i] = nu_2 } 
      else if (Vah[i] > 5 && Vah[i] <= 15) {nu_hm[i] = nu_0; nu_hf[i] = nu_0 }
      else {nu_hm[i] = nu_1 }
      
      if (Vah[i] > 15 && Vah[i] <= 45) { nu_hf[i] = p_f * nu_0 + (1 - p_f) * nu_1 } 
      else if (Vah[i] > 40) { nu_hf[i] = nu_1}
    }
  }
  
  Sm[1]=6.5*10^6# Wedge_m/mu_m;
  Im[1,1]=0.75*10^6#rep(0.75*10^6, Nsigma)
  
  Sm_ivm[1,]=rep(0,Ntau); Im_ivm[,,1]=matrix(0,nrow = Nsigma,ncol = Ntau)
  mSh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau);fSh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau)
  mAh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau);fAh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau)
  mIh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau);fIh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau)
  mRh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau);fRh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau)
  
  t=1
  Nh[t]= (sum(mSh[t,]+fSh[t,]+mAh[t,]+fAh[t,]+mIh[t,]+fAh[t,]+mRh[t,]+fRh[t,])
          +sum(mSh_ivm[,,t]+fSh_ivm[,,t]+mAh_ivm[,,t]+fAh_ivm[,,t]+mIh_ivm[,,t]
               +fIh_ivm[,,t]+mRh_ivm[,,t]+fRh_ivm[,,t]))
  Nm[t]= Sm[t] + sum(Im[t,])+ sum(Sm_ivm[t,]) +sum(Im_ivm[,,t])
  #Nbrun=0
  tau_eff_age= 0 #initialisation
  lam_hI=rep(0,Ntau);lam_hS=rep(0,Ntau)
  
  for (t in 1:(Ntime-1)) {
    lam_m = theta * ( sum(beta_m[]*Im[t,]) + sum(beta_m[]*Im_ivm[,,t]))/Nh[t];
    lam_h = theta * (sum( beta_h[]*(mAh[t,]+fAh[t,]) ) 
                     + sum( bar_beta_h[]*(mIh[t,]+fIh[t,])))/Nh[t];
    
    for (tau in 1:Ntau) {
      lam_hI[tau]=theta*sum(mAh_ivm[,tau,t]+fAh_ivm[,tau,t]+mIh_ivm[,tau,t]+fIh_ivm[,tau,t])/Nh[t];
      lam_hS[tau]=theta*sum(mSh_ivm[,tau,t]+fSh_ivm[,tau,t]+mRh_ivm[,tau,t]+fRh_ivm[,tau,t])/Nh[t]; 
    }
    
    Sm[t+1] = (Wedge_m + Sm[t]/dt )/(1/dt  + mu_m + sum(lam_hI[] + lam_hS[]) + lam_h );
    
    sigma=1
    Im[t+1,sigma]=(Im[t,sigma]/dt + lam_h*Sm[t]/dsigma)/
      (1/dt + 1/dsigma +  mu_m + sum(lam_hI[] + lam_hS[]))
    
    for (sigma in 2:Nsigma ) {
      Im[t+1,sigma]=(Im[t,sigma]/dt + Im[t+1,sigma-1]/dsigma)/
        (1/dt + 1/dsigma +  mu_m + sum(lam_hI[] + lam_hS[]))
    }
    
    ### nouvelle naissance
    a=1 
    mSh[t+1,a]=(mSh[t,a]/dt + q*Wedge_h +k_h*mRh[t,a])/
      (1/dt + 1/dah +  mu_h[a] + mphi[t,a] + lam_m); 
    fSh[t+1,a]=(fSh[t,a]/dt + (1-q)*Wedge_h +k_h*fRh[t,a])/
      (1/dt + 1/dah +  mu_h[a] + fphi[t,a] + lam_m); 
    mAh[t+1,a]=0
    fAh[t+1,a]=0
    mIh[t+1,a]=0
    fAh[t+1,a]=0
    mRh[t+1,a]=0
    fRh[t+1,a]=0
    mSh_ivm[a,,t+1]=rep(0,Ntau);fSh_ivm[a,,t+1]=rep(0,Ntau)
    mAh_ivm[a,,t+1]=rep(0,Ntau);fAh_ivm[a,,t+1]=rep(0,Ntau)
    mIh_ivm[a,,t+1]=rep(0,Ntau);fIh_ivm[a,,t+1]=rep(0,Ntau)
    mRh_ivm[a,,t+1]=rep(0,Ntau);fRh_ivm[a,,t+1]=rep(0,Ntau)
    
    #on integre en age
    for (a in 2:Nah) {
      
      mSh[t+1,a]=(mSh[t,a]/dt +mSh[t+1,a-1]/dah +k_h*mRh[t,a]+sum(rho[a,]*mSh_ivm[a,,t]))/
        (1/dt + 1/dah +  mu_h[a] + mphi[t,a] + lam_m); 
      fSh[t+1,a]=(fSh[t,a]/dt +fSh[t+1,a-1]/dah +k_h*fRh[t,a]+sum(rho[a,]*fSh_ivm[a,,t]))/
        (1/dt + 1/dah +  mu_h[a] + fphi[t,a] + lam_m);
      mAh[t+1,a]=(mAh[t,a]/dt +mAh[t+1,a-1]/dah +lam_m*mSh[t,a]+sum(rho[a,]*mAh_ivm[a,,t]))/
        (1/dt + 1/dah +  mu_h[a] + mphi[t,a] + nu_hm[a]);
      fAh[t+1,a]=(fAh[t,a]/dt +fAh[t+1,a-1]/dah +lam_m*fSh[t,a]+sum(rho[a,]*fAh_ivm[a,,t]))/
        (1/dt + 1/dah +  mu_h[a] + fphi[t,a] + nu_hf[a]);
      mIh[t+1,a]=(mIh[t,a]/dt +mIh[t+1,a-1]/dah +nu_hm[a]*mAh[t,a]+sum(rho[a,]*mIh_ivm[a,,t]))/
        (1/dt + 1/dah +  mu_h[a] + delta_h[a]  + gamma_h[a]);
      fIh[t+1,a]=(fIh[t,a]/dt +fIh[t+1,a-1]/dah +nu_hf[a]*fAh[t,a]+sum(rho[a,]*fIh_ivm[a,,t]))/
        (1/dt + 1/dah +  mu_h[a] + delta_h[a]  + gamma_h[a]);
      mRh[t+1,a]=(mRh[t,a]/dt +mRh[t+1,a-1]/dah +gamma_h[a]*mIh[t,a]+sum(rho[a,]*mRh_ivm[a,,t]))/
        (1/dt + 1/dah + mu_h[a] +mphi[t,a] + k_h);
      fRh[t+1,a] = (fRh[t,a]/dt + fRh[t+1,a-1]/dah + gamma_h[a]*fIh[t,a] + sum(rho[a,]*fRh_ivm[a,,t])) /
        (1/dt + 1/dah + mu_h[a] + fphi[t,a] + k_h)
      
      #nouvelle IVM d'age a
      tau=1
      mSh_ivm[a,tau,t+1]=(mSh_ivm[a,tau,t]/dt + mphi[t,a]*mSh[t,a]
                          +mSh_ivm[a-1,tau,t+1]/dah +k_h*mRh_ivm[a,tau,t])/
        (1/dt + 1/dtau +1/dah + mu_h[a] + lam_m + rho[a,tau]);
      fSh_ivm[a,tau,t+1]=(fSh_ivm[a,tau,t]/dt + fphi[t,a]*fSh[t,a] 
                          +fSh_ivm[a-1,tau,t+1]/dah +k_h*fRh_ivm[a,tau,t])/
        (1/dt + 1/dtau +1/dah + mu_h[a] + lam_m + rho[a,tau]);
      mAh_ivm[a,tau,t+1]=(mAh_ivm[a,tau,t]/dt + mphi[t,a]*mAh[t,a] 
                          +mAh_ivm[a-1,tau,t+1]/dah +lam_m*mSh_ivm[a,tau,t])/
        (1/dt + 1/dtau +1/dah + mu_h[a] + nu_hm[a] + rho[a,tau]);
      fAh_ivm[a,tau,t+1]=(fAh_ivm[a,tau,t]/dt + fphi[t,a]*fAh[t,a] 
                          +fAh_ivm[a-1,tau,t+1]/dah +lam_m*fSh_ivm[a,tau,t])/
        (1/dt + 1/dtau +1/dah + mu_h[a] + nu_hf[a] + rho[a,tau]);
      mIh_ivm[a,tau,t+1]=(mIh_ivm[a,tau,t]/dt + 0 +mIh_ivm[a-1,tau,t+1]/dah 
                          +nu_hm[a]*mAh_ivm[a,tau,t])/
        (1/dt + 1/dtau +1/dah + mu_h[a] + gamma_h[a]+delta_h[a] + rho[a,tau]);
      fIh_ivm[a,tau,t+1]=(fIh_ivm[a,tau,t]/dt + 0 +fIh_ivm[a-1,tau,t+1]/dah 
                          +nu_hf[a]*fAh_ivm[a,tau,t])/
        (1/dt + 1/dtau +1/dah + mu_h[a] + gamma_h[a]+delta_h[a] + rho[a,tau]);
      mRh_ivm[a,tau,t+1]=(mRh_ivm[a,tau,t]/dt + mphi[t,a]*mRh[t,a] 
                          +mRh_ivm[a-1,tau,t+1]/dah +gamma_h[a]*mIh_ivm[a,tau,t])/
        (1/dt + 1/dtau +1/dah + mu_h[a] + k_h + rho[a,tau]);
      fRh_ivm[a,tau,t+1]=(fRh_ivm[a,tau,t]/dt + fphi[t,a]*fRh[t,a] 
                          +fRh_ivm[a-1,tau,t+1]/dah +gamma_h[a]*fIh_ivm[a,tau,t])/
        (1/dt + 1/dtau +1/dah + mu_h[a] + k_h + rho[a,tau]);
      
      #on integre les humains sous IVM
      for (tau in 2:Ntau) {
        mSh_ivm[a,tau,t+1]=(mSh_ivm[a,tau,t]/dt + mSh_ivm[a,tau-1,t+1]/dtau 
                            +mSh_ivm[a-1,tau,t+1]/dah +k_h*mRh_ivm[a,tau,t])/
          (1/dt + 1/dtau +1/dah + mu_h[a] + lam_m + rho[a,tau]);
        fSh_ivm[a,tau,t+1]=(fSh_ivm[a,tau,t]/dt + fSh_ivm[a,tau-1,t+1]/dtau 
                            +fSh_ivm[a-1,tau,t+1]/dah +k_h*fRh_ivm[a,tau,t])/
          (1/dt + 1/dtau +1/dah + mu_h[a] + lam_m + rho[a,tau]);
        mAh_ivm[a,tau,t+1]=(mAh_ivm[a,tau,t]/dt + mAh_ivm[a,tau-1,t+1]/dtau 
                            +mAh_ivm[a-1,tau,t+1]/dah +lam_m*mSh_ivm[a,tau,t])/
          (1/dt + 1/dtau +1/dah + mu_h[a] + nu_hm[a] + rho[a,tau]);
        fAh_ivm[a,tau,t+1]=(fAh_ivm[a,tau,t]/dt + fAh_ivm[a,tau-1,t+1]/dtau 
                            +fAh_ivm[a-1,tau,t+1]/dah +lam_m*fSh_ivm[a,tau,t])/
          (1/dt + 1/dtau +1/dah + mu_h[a] + nu_hf[a] + rho[a,tau]);
        mIh_ivm[a,tau,t+1]=(mIh_ivm[a,tau,t]/dt + mIh_ivm[a,tau-1,t+1]/dtau 
                            +mIh_ivm[a-1,tau,t+1]/dah +nu_hm[a]*mAh_ivm[a,tau,t])/
          (1/dt + 1/dtau +1/dah + mu_h[a] + gamma_h[a]+delta_h[a] + rho[a,tau]);
        fIh_ivm[a,tau,t+1]=(fIh_ivm[a,tau,t]/dt + fIh_ivm[a,tau-1,t+1]/dtau 
                            +fIh_ivm[a-1,tau,t+1]/dah +nu_hf[a]*fAh_ivm[a,tau,t])/
          (1/dt + 1/dtau +1/dah + mu_h[a] + gamma_h[a]+delta_h[a] + rho[a,tau]);
        mRh_ivm[a,tau,t+1]=(mRh_ivm[a,tau,t]/dt + mRh_ivm[a,tau-1,t+1]/dtau 
                            +mRh_ivm[a-1,tau,t+1]/dah +gamma_h[a]*mIh_ivm[a,tau,t])/
          (1/dt + 1/dtau +1/dah + mu_h[a] + k_h + rho[a,tau]);
        fRh_ivm[a,tau,t+1]=(fRh_ivm[a,tau,t]/dt + fRh_ivm[a,tau-1,t+1]/dtau 
                            +fRh_ivm[a-1,tau,t+1]/dah +gamma_h[a]*fIh_ivm[a,tau,t])/
          (1/dt + 1/dtau +1/dah + mu_h[a] + k_h + rho[a,tau]);
        
      }
    }
    
    tau = 1
    Sm_ivm[t+1, tau] = 0
    
    # sigma > 1 à tau=1
    for(sigma in 2:Nsigma) {
      Im_ivm[sigma, tau, t+1] = 0
    }
    # tau > 1  
    for(tau in 2:Ntau) {
      Sm_ivm[t+1, tau] = (lam_hS[tau]*Sm[t]+Sm_ivm[t, tau]/dt + Sm_ivm[t+1, tau-1]/dtau) / 
        (1/dt + 1/dtau + mu_m.ivm[tau] + lam_hI[tau] + lam_h)
      
      sigma = 1
      Im_ivm[sigma, tau, t+1] = (Im_ivm[sigma, tau, t]/dt +  Im_ivm[sigma, tau-1, t+1]/dtau + 
                                   (lam_hI[tau] + lam_h)*Sm_ivm[t, tau]/dsigma + 
                                   lam_hI[tau]*Sm[t]/dsigma ) / 
        (1/dt + 1/dtau + 1/dsigma + mu_m.ivm[tau])
      
      #  sigma > 1
      for(sigma in 2:Nsigma) {
        Im_ivm[sigma, tau, t+1] = (Im_ivm[sigma, tau, t]/dt + (lam_hI[tau] + lam_hS[tau])*Im[t, sigma]+
                                     Im_ivm[sigma, tau-1, t+1]/dtau + 
                                     Im_ivm[sigma-1, tau, t+1]/dsigma) / 
          (1/dt + 1/dtau + 1/dsigma + mu_m.ivm[tau])
      }
    }
    
    
    Nh[t+1]= (sum(mSh[t+1,]+fSh[t+1,]+mAh[t+1,]+fAh[t+1,]+mIh[t+1,]+fAh[t+1,]
                  +mRh[t+1,]+fRh[t+1,]) +sum(mSh_ivm[,,t+1]+fSh_ivm[,,t+1] +mAh_ivm[,,t+1]
                                             +fAh_ivm[,,t+1]+mIh_ivm[,,t+1]+fIh_ivm[,,t+1]+mRh_ivm[,,t+1]+fRh_ivm[,,t+1]))
    
    Nm[t+1]= Sm[t+1] + sum(Im[t+1,])+ sum(Sm_ivm[t+1,]) +sum(Im_ivm[,,t+1])
    
    #Nbrun=Nbrun+1
    #print(Nbrun) #le max etant Ntime
  }  
  
  Sh_Tot=time;Ah_Tot=time;Ih_Tot=time;Rh_Tot=time;Ih_15=time;Nh_15=time
  Sm_Tot=time;Im_Tot=time; Im_TotEIP=time
  
  mIh_a=fIh_a=matrix(0,nrow = Ntime,ncol = Nah);mAh_a=fAh_a=matrix(0,nrow = Ntime,ncol = Nah)
  Shg=matrix(0,nrow= Ntime,ncol= Nah);Ahg=matrix(0,nrow= Ntime,ncol= Nah)
  Ihg=matrix(0,nrow= Ntime,ncol= Nah);Rhg=matrix(0,nrow= Ntime,ncol= Nah)
  Shg_ivm=matrix(0,nrow= Ntime,ncol= Nah);Ahg_ivm=matrix(0,nrow= Ntime,ncol= Nah)
  Ihg_ivm=matrix(0,nrow= Ntime,ncol= Nah);Rhg_ivm=matrix(0,nrow= Ntime,ncol= Nah)
  Nm_IVM=matrix(0,nrow = Ntime,ncol = Ntau)
  Smg_ivm=time;Img=time; Img_ivm=time
  for (t in 1:Ntime){
    Sh_Tot[t]=sum(mSh[t,] + fSh[t,])+sum(mSh_ivm[,,t] +fSh_ivm[,,t])
    Ah_Tot[t]=sum(mAh[t,] + fAh[t,])+sum(mAh_ivm[,,t] +fAh_ivm[,,t])
    Ih_Tot[t]=sum(mIh[t,] + fIh[t,])+sum(mIh_ivm[,,t] +fIh_ivm[,,t])
    Rh_Tot[t]=sum(mRh[t,] + fRh[t,])+sum(mRh_ivm[,,t] +fRh_ivm[,,t])
    Ih_15[t] = sum(mIh[t,6:15] + fIh[t,6:15]) + sum(mIh_ivm[6:15,,t] + fIh_ivm[6:15,,t])
    Nh_15[t]= (sum(mSh[t,6:15]+fSh[t,6:15]+mAh[t,6:15]+fAh[t,6:15]+mIh[t,6:15]+fAh[t,6:15]+mRh[t,6:15]
                   +fRh[t,6:15]) +sum(mSh_ivm[6:15,,t]+fSh_ivm[6:15,,t]+mAh_ivm[6:15,,t]+fAh_ivm[6:15,,t]
                                      +mIh_ivm[6:15,,t]+fIh_ivm[6:15,,t]+mRh_ivm[6:15,,t]+fRh_ivm[6:15,,t]))
    
    Sm_Tot[t]=sum(Sm_ivm[t,]) + Sm[t] 
    Smg_ivm[t]=sum(Sm_ivm[t,])
    Img[t]=sum(Im[t,])
    Img_ivm[t]= sum(Im_ivm[, , t]) 
    Im_Tot[t]=sum(Im_ivm[,,t]) + sum(Im[t,])
    Im_TotEIP[t]=sum(Im_ivm[8:Nsigma,,t]) + sum(Im[t,8:Nsigma])
    for (a in 1:Nah){
      mIh_a[t,a]=mIh[t,a] + sum(mIh_ivm[a,,t]); fIh_a[t,a]=fIh[t,a] + sum(fIh_ivm[a,,t])
      mAh_a[t,a]=mAh[t,a] + sum(mAh_ivm[a,,t]); fAh_a[t,a]=fAh[t,a] + sum(fAh_ivm[a,,t])
      Shg[t,a]=mSh[t,a] + fSh[t,a]; Ahg[t,a]=mAh[t,a] + fAh[t,a]
      Ihg[t,a]=mIh[t,a] + fIh[t,a]; Rhg[t,a]=mRh[t,a] + fRh[t,a]
      Shg_ivm[t,a]=sum(mSh_ivm[a,,t]+fSh_ivm[a,,t]); Ahg_ivm[t,a]=sum(mAh_ivm[a,,t]+fAh_ivm[a,,t])
      Ihg_ivm[t,a]=sum(mIh_ivm[a,,t]+fIh_ivm[a,,t]); Rhg_ivm[t,a]=sum(mRh_ivm[a,,t]+fRh_ivm[a,,t])
    }
    
    for(tau in 1:Ntau) {
      Nm_IVM[t,tau]= Sm_ivm[t, tau] + sum(Im_ivm[, tau, t]) 
    }
  }
  
  IdBegin_Camp= 1+floor(t_begin_Camp/dt)
  IndexTime=IdBegin_Camp:Ntime
  Delta=trapz(time[IndexTime],Ih_Tot[IndexTime])
  
  PropIhTot=Ih_Tot/Nh;PropShTot=Sh_Tot/Nh;PropAhTot=Ah_Tot/Nh
  PropRhTot=Rh_Tot/Nh;PropIh_15=Ih_15/Nh_15
  PropImTot=Im_Tot/Nm;PropSmTot=Sm_Tot/Nm; PropImEIP=Im_TotEIP/Nm;
  
  IdTop=which(PropIhTot[IndexTime]==min(PropIhTot[IndexTime]))
  DurIvmEffect=time[IdTop]
  
  MinPropIhTot=min(PropIhTot[IndexTime])
  PropIhTot0=PropIhTot[t_begin_Camp/dt]
  
  IndexTimeHR = (IdBegin_Camp ):Ntime   
  IdTopHR = which(PropIhTot[IndexTimeHR] >= PropIhTot0)[1] 
  
  if (!is.na(IdTopHR)) {
    IdTopHR = IdTopHR 
    DurIvmEffectHR = time[IdTopHR]
  } else {
    IdTopHR = 0
    DurIvmEffectHR = 0
  }
  
  #Optimisation Timing, et cumulative gain prev, avec le deuxieme gain prev
  # DeltaOptim est le nombre cummule' de Ih (with IVM)
  Nyr=360/dt
  IndexOptim=IdBegin_Camp:(IdBegin_Camp+Nyr)
  DeltaOptim=sum(Ih_Tot[IndexOptim])
  DeltaOptim_15=sum(Ih_15[IndexOptim])
  Delta_Tot=sum(Ih_Tot[IndexOptim]+Ah_Tot[IndexOptim])
  
  cut_age = seq(0, ah_max, by = 5)
  
  time_idx = IdBegin_Camp + IdTop[1]  # Prendre le premier si plusieurs
  time_idx = min(time_idx, nrow(mIh_a))  # S'assurer qu'on ne dépasse pas
  
  # Regrouper en classes de 5 ans (dah = 1)
  mIh_age_groups = sapply(1:(length(cut_age) - 1), function(i) {
    idx_debut = cut_age[i] + 1
    idx_fin = min(cut_age[i + 1], ncol(mIh_a))  
    
    if (idx_debut <= idx_fin) {
      sum(mIh_a[time_idx, idx_debut:idx_fin])
    } else {
      0
    }
  })
  
  fIh_age_groups = sapply(1:(length(cut_age) - 1), function(i) {
    idx_debut = cut_age[i] + 1
    idx_fin = min(cut_age[i + 1], ncol(mIh_a))  
    
    if (idx_debut <= idx_fin) {
      sum(fIh_a[time_idx, idx_debut:idx_fin])
    } else {
      0
    }
  })
  
  
  FunList = list("time" = time, "mSh"=mSh,"Sh_Tot"=Sh_Tot,"Ah_Tot"=Ah_Tot,"Ih_Tot"=Ih_Tot,"Rh_Tot"=Rh_Tot,
                 "Nh"=Nh,"Nm"=Nm, "Sm_Tot"=Sm_Tot,"Im_Tot"=Im_Tot,"Delta"=Delta,"lam_h"=lam_h,"lam_m"=lam_m, 
                 "lam_hI"=lam_hI, "lam_hS"=lam_hS, "phi"=phi,"mIh_a"=mIh_age_groups, "fIh_a"=fIh_age_groups,
                 "mAh_a"=mAh_a, "fAh_a"=fAh_a,"PropRhTot"=PropRhTot,"PropIhTot"=PropIhTot,"PropImTot"=PropImTot,
                 "DurIvmEffect"=DurIvmEffect,"PropShTot"=PropShTot,"PropSmTot"=PropSmTot,"PropAhTot"=PropAhTot,
                 "MinPropIhTot"=MinPropIhTot, "PropIhTot0"=PropIhTot0, "cut_age"=cut_age, "Nh_15"=Nh_15,
                 "DurIvmEffectHR"=DurIvmEffectHR, "DeltaOptim"=DeltaOptim,"IndexOptim"=IndexOptim, "Ih_15"=Ih_15,
                 "PropImEIP"=PropImEIP, "Nm_IVM"=Nm_IVM,"Shg"=Shg,"Ahg"=Ahg,"Ihg"=Ihg,"Rhg"=Rhg,"Shg_ivm"=Shg_ivm,
                 "Ahg_ivm"=Ahg_ivm,"Ihg_ivm"=Ihg_ivm,"Rhg_ivm"=Rhg_ivm,"PropIh_15"=PropIh_15,"DeltaOptim_15"=DeltaOptim_15,
                 "Img_ivm"=Img_ivm, "Smg_ivm"=Smg_ivm, "Img"=Img,"Delta_Tot"=Delta_Tot)
  return(FunList)
}

#parametres de base
{ 
  Gap=0.25;
  N_parm= 2
  dt=Gap;tmax = 611;time=seq(0,tmax,by=dt); Ntime=length(time)
  dtau=1;tau_max= 300;Vtau=seq(0,tau_max,by=dtau);Ntau=length(Vtau);
  dah=1;ah_max=90;Vah=seq(0,ah_max,by=dah); Nah=length(Vah)
  dsigma=dtau;sigma_max=60;Vsigma=seq(0,sigma_max,by=dsigma); Nsigma=length(Vsigma)
  mu_m = 0.13  
  IVM_Pregnancy=0
  EIP=8
  theta=0.5
  Dur_cycle=7 
  t_begin_Camp=250
  tmaxmonths=(tmax-t_begin_Camp)/30
  
  #mu_h(a) mortalité naturelle 
  {
    # Values corresponding to age ranges, Ref Quentin Richard et al Bobo Dioulasso
    mu_h_values = c(66.8, 7.6, 1.7, 0.9, 1.3, 1.9, 2.4, 2.8, 3.6, 4.7, 6.3, 8.9,
                    13.2, 19.8, 31.1, 47.7, 71.3, 110.5, 186.7) / 1000 
    # Age labels corresponding to age ranges
    age_labels0=c(0,1,seq(5,90,by=5))
    age_labels = c("0 to 1 year", "1 to 5 years", "5 to 10 years", "10 to 15 years", 
                   "15 to 20 years", "20 to 25 years", "25 to 30 years", 
                   "30 to 35 years", "35 to 40 years", "40 to 45 years", 
                   "45 to 50 years", "50 to 55 years", "55 to 60 years", 
                   "60 to 65 years", "65 to 70 years", "70 to 75 years", 
                   "75 to 80 years", "80 to 85 years", "85 to 90 years")
    
    mu_h = rep(NA, Nah) 
    # Assign mu_h according to age_labels
    for (i in 1:Nah) {
      for (IdAgeGroup in 1:(length(age_labels0)-1)) {
        if(age_labels0[IdAgeGroup]<=Vah[i]&
           Vah[i]<age_labels0[IdAgeGroup+1])
        {mu_h[i]=mu_h_values[IdAgeGroup]}
      }
      if (i>=Nah) {mu_h[i]=mu_h_values[length(mu_h_values)]}
    }
    plot(Vah,mu_h)
  }
  
  #delta_h(a), Ref Quentin Richard et al Bobo Dioulasso. Motalité induite maladie
  {
    # Values of delta_h corresponding to the age ranges
    delta_h_values = c(1.07 * 10^(-3), 7.02 * 10^(-4), 4.55 * 10^(-4), 5.73 * 10^(-5))
    # Age labels for the age ranges
    Age_label_2 = c("0-1 years", "1-5 years", "5-15 years", ">15 years")
    # Create a vector to store delta_h according to age
    delta_h = rep(NA, Nah)
    # Assign nu_h according to Age_label_2
    for (i in 1:Nah) {
      if(0<=Vah[i]&Vah[i]<=1)
      {delta_h[i]=delta_h_values[1]}
      
      if(1<Vah[i]&Vah[i]<=5)
      {delta_h[i]=delta_h_values[2]}
      
      if(5<Vah[i]&Vah[i]<=15)
      {delta_h[i]=delta_h_values[3]}
      
      if(15<Vah[i])
      {delta_h[i]=delta_h_values[4]}
    }
    plot(Vah,delta_h)
  }
  
  #beta_h(a) et Wedge_M: infectiosité humaine et recrutement des moustiques
  {
    params_Mos_Prev =function(Prev_mos10) { #Prevealence des moustiques infectieux 5, 10%
      if (Prev_mos10) {
        beta_m  <<-  0.0416*(Vsigma>EIP)
        Wedge_m <<- 1.33 * 10^6   
        alpha_1 <<- 0.122       # 0.071[0.023; 0.175]
        alpha_2 <<- 0.17         # 0.302[0.16; 0.475]
      } else { #Prev_mos5_percent
        beta_m <<-  0.063*(Vsigma>EIP)
        Wedge_m <<- 1.33 * 10^6    
        alpha_1 <<- 0.071        # 0.071[0.023; 0.175]
        alpha_2 <<- 0.15  
      }
      beta_h <<- rep(NA, Nah)
      
      # Define the function G(a)
      G <<- function(a) {
        return(22.7 * a * exp(-0.0934 * a))
      }
      # Compute beta_h for each age in Vah
      beta_h <<- alpha_1 * (G(Vah) ^ alpha_2)  
      bar_beta_h <<- 0.8 * beta_h
      
      #plot(Vah,G(Vah))
      #plot(Vah,bar_beta_h)
    }
    
    Prev_mos10=1
    params_Mos_Prev(Prev_mos10)
    
    #Nom pour l'enregistrement des figures pdf
    Name_Prev_mos=function(Prev_mos10) {
      ifelse(Prev_mos10 == 1, "10", "5")
    }
  }
  
  #nu_hm(a), nu_hf(a): progression Asympto vers cas clinique dans la main fonction
  {
    p_f = 0.15 # Fertility probability for women aged 15 to 40 in Bobo Dioulasso,
    
    # plot(Vah,nu_hm)
    # plot(Vah,nu_hf)
  }
  
  #gamma_h(a): guerison
  {
    # Recovery time ranges for each age group (days)
    recovery_ranges = list(c(163, 345), c(555, 714), c(344, 400), c(181, 204), 
                           c(82, 92), c(56, 61), c(48, 55))
    # Age ranges  corresponding to the recovery time ranges
    age_ranges =list(c(0, 1),c(1, 5),c(5, 8),c(8, 18),c(18, 28),c(28, 43),c(43, ah_max))
    gamma_h = rep(NA, Nah)
    
    #  value of gamma_h
    for (i in 1:Nah) {
      for (j in 1: length(age_ranges)) {
        age_min = age_ranges[[j]][1]
        age_max = age_ranges[[j]][2]
        if(Vah[i] >= age_min & Vah[i] <= age_max)
        {gamma_h[i] = 1 / mean(recovery_ranges[[j]])}
      }
    }
    plot(Vah,gamma_h)
  }
}

#Figure S3: Without IVM
{
  Prev_mos10=0
  VectPrev_mos10=c(1)
  for (Prev_mos10 in VectPrev_mos10) {
    params_Mos_Prev(Prev_mos10)
    IVM_field_dependancy=0; strategy =0; PropIVM=0;IVM_Pregnancy=0
    VectTime_between_cycles=c(0,0,0); Number_of_cycle=1
    ParamIvmFormulation = ModelEtau50Alpha(strategy)
    tau50 = ParamIvmFormulation$tau50; alpha = ParamIvmFormulation$alpha; Dmax = ParamIvmFormulation$Dmax
    
    OutPut = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                      mu_h,delta_h,bar_beta_h,beta_h,gamma_h,theta,beta_m,Wedge_m,mu_m,
                      ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy, dtau,tau_max,Vtau,Ntau,
                      dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,tau50, alpha, Dmax)
    
    maxY=max(OutPut$PropRhTot, max(OutPut$PropIhTot, max(OutPut$PropAhTot, max(OutPut$PropShTot))))
    maxYm=max(OutPut$PropSmTot,max(OutPut$PropImEIP))
    
    Time_plot0=110
    tmaxmonths=(tmax-Time_plot0)/30
    
    save(OutPut, file = paste0("Model_Without_IVM_Prev_mos", Name_Prev_mos(Prev_mos10), ".RData"))
    
    load(paste0("Model_Without_IVM_Prev_mos", Name_Prev_mos(Prev_mos10), ".RData"))
    ls()
    #plot
    {
      FigName=paste0("Dynamic_without-IVM_Prev_mos", Name_Prev_mos(Prev_mos10), ".pdf")
      pdf(FigName,width=7,height=3)
      par(oma=c(1,.1,1,.1),mar=c(2,3.2,1,.5))
      par(mfrow=c(1,2))
      
      GrNumber=ifelse(Prev_mos10, 0, 2)
      ColVect = c( "#037153", "#D98E04","#ff3355", "#205072","#8B008B","#8B4513")
      LineVect=c(1, 3, 2, 4)
      LC=2
      
      #human population
      {
        GrNumber=GrNumber+1
        plot(-1, 1, type = "l", xlab = "", xlim = c(Time_plot0, tmax), ylab = "",
             ylim = c(0, 0.77), cex.lab = 0.6, yaxt = "n", xaxt = "n")
        axis(1, at = seq(Time_plot0, tmax, by = 60), labels = seq(0, tmaxmonths, by = 2),cex.axis = 0.7, las = 1)
        y_vals = c(0, 0.15,  0.30, 0.45, 0.6, 0.75)
        axis(2, at = y_vals, labels = y_vals * 100,cex.axis = 0.7, las = 1)
        par(xpd = NA)
        text(Time_plot0-0.05, 0.75 * 1.15, paste("(", LETTERS[GrNumber], ")", sep = ""), cex = 1)
        par(xpd = FALSE)
        lines(time, OutPut$PropShTot, lwd = LC, cex.lab = 1.2, lty = LineVect[1], col = ColVect[1])
        lines(time, OutPut$PropAhTot, lwd = LC, cex.lab = 1.2,    lty = LineVect[2], col = ColVect[2])
        lines(time, OutPut$PropIhTot, lwd = LC, cex.lab = 1.2, lty = LineVect[3], col = ColVect[3])
        lines(time, OutPut$PropRhTot, lwd = LC, cex.lab = 1.2, lty = LineVect[4], col = ColVect[4])
        
        mtext("Human population (all ages; %)", side = 2, adj = 0.5, cex = .85, line = 2, font = .8) 
        
        mtext(ifelse(Prev_mos10, "Time (Months)", ""), side = 1, adj = 0.5, cex = .9, line = 2, font = .8) 
        
        LegendTex=c(expression(S[h]),expression(A[h]), expression(I[h]),expression(R[h]))
        
        if(Prev_mos10){
          legend("topright", legend=LegendTex, xpd = NA, horiz = TRUE,lwd=LC,
                 bty = "n",col=ColVect, lty=LineVect, cex = 0.75, text.col= 1)
        }
      }
      
      #Vector population
      {
        GrNumber=GrNumber+1
        plot(-1, 1, type = "l", xlab = "", xlim = c(Time_plot0, tmax), ylab = "",
             ylim = c(0, 1.02), cex.lab = 0.6, yaxt = "n", xaxt = "n")
        axis(1, at = seq(Time_plot0, tmax, by = 60), labels = seq(0, tmaxmonths, by = 2), cex.axis = 0.7, las = 1)
        y_vals = c(0, 0.20,  0.4,  0.6,  0.8, 1)
        axis(2, at = y_vals, labels = y_vals * 100, cex.axis = 0.7, las = 1)
        par(xpd = NA)
        text(Time_plot0-0.05, 1.02 * 1.15, paste("(", LETTERS[GrNumber], ")", sep = ""), cex = 1)
        par(xpd = FALSE)
        lines(time, OutPut$PropSmTot, lwd = LC, cex.lab = 1.2, lty = LineVect[1], col = ColVect[1])
        lines(time, OutPut$PropImEIP, lwd = LC, cex.lab = 1.2,    lty = LineVect[3], col = ColVect[3])
        lines(time, OutPut$PropImTot, lwd = LC, cex.lab = 1.2,    lty = LineVect[2], col = ColVect[2])
        mtext("Mosquito population", side = 2, adj = 0.5, cex = .85, line = 2, font = .8) 
        
        mtext(ifelse(Prev_mos10, "Time (Months)" , ""), side = 1, adj = 0.5, cex = .9, line = 2, font = .8)
        
        LegendTex=c( expression(S[m]), expression(I[m] * " Infected"), expression(I[m] * " Infectious"))
        
        eq = 1 + floor(t_begin_Camp/dt)
        ratio_eq = OutPut$Nm[eq] / OutPut$Nh[eq]
        
        text(x = tmax * 0.7, y = 0.55, labels = paste0("Mosquito/Human\nratio: ", round(ratio_eq, 2)),
             col = "#205072", cex = 0.6,font = 1)
        if(Prev_mos10){
          legend("topright", legend = LegendTex, xpd = NA, horiz = TRUE, lwd = LC, bty = "n", 
                 col = c("#037153", "#D98E04", "#ff3355"), lty = c(1, 3, 2), cex = 0.6, text.col = 1)
        }
      }
    }
    dev.off() 
  }
  
  #Renitialiser pour la suite
  tmaxmonths=(tmax-t_begin_Camp)/30
}

GainPrev_with_CI <- function(strategy, PropIVM, IVM_field_dependancy, scenario_type,
                             CoefDelay, coefTimingLAIF, DeltaOptimInit, DeltaOptimInit_15,
                             IVM_Pregnancy, NbCycle_override = NULL) {
  
  if (strategy %in% c(0, 1, 5)) NbCycle = 4
  else if (strategy == 2) NbCycle = ifelse(is.null(NbCycle_override), 3, NbCycle_override)
  else NbCycle = ifelse(is.null(NbCycle_override), 4, NbCycle_override)
  
  ParamIvmFormulation = Estimate_MCMC_Ptau(strategy)
  tau50_est = ParamIvmFormulation$tau50
  alpha_est = ParamIvmFormulation$alpha
  Dmax = ParamIvmFormulation$Dmax
  tau50_seq = gen_param_seq(tau50_est, ParamIvmFormulation$tau50_low, 
                            ParamIvmFormulation$tau50_up, n = N_parm)
  alpha_seq = gen_param_seq(alpha_est, ParamIvmFormulation$alpha_low,
                            ParamIvmFormulation$alpha_up, n = N_parm)
  
  #toutes les combinaisons de paramètres
  param_grid = expand.grid(tau50 = tau50_seq, alpha = alpha_seq)
  
  rm(ParamIvmFormulation)
  gc()
  
  run_model <- function(tau50, alpha, nc, return_model = FALSE) {
    if (scenario_type == "optimal") {
      VectTime_between_cycles = c(0,0,0)
      Number_of_cycle = 1
      ModelOutput = ModelIVM(PropIVM, VectTime_between_cycles, Number_of_cycle, t_begin_Camp, Dur_cycle,
                             mu_h, delta_h, bar_beta_h, beta_h, gamma_h, theta, beta_m, Wedge_m, mu_m,
                             ah_max, dah, Nah, tmax, Gap, dt, time, Ntime, p_f, IVM_Pregnancy,
                             dtau, tau_max, Vtau, Ntau, dsigma, sigma_max, Vsigma, Nsigma, strategy, 
                             IVM_field_dependancy, tau50, alpha, Dmax)
      if (nc > 1) {
        TempsNextCycle0 = 0
        for (k in 2:nc) {
          Number_of_cycle = k
          TempsNextCycle = ModelOutput$DurIvmEffect
          d0 = TempsNextCycle - TempsNextCycle0 - Dur_cycle
          TempsNextCycle0 = TempsNextCycle
          coef_timing = ifelse(strategy %in% c(0, 1, 5) && k > 2, coefTimingLAIF, 1)
          VectTime_between_cycles[k-1] = CoefDelay * d0 * coef_timing
          ModelOutput = ModelIVM(PropIVM, VectTime_between_cycles, Number_of_cycle, t_begin_Camp, Dur_cycle, 
                                 mu_h, delta_h, bar_beta_h, beta_h, gamma_h, theta, beta_m, Wedge_m, 
                                 mu_m, ah_max, dah, Nah, tmax, Gap, dt, time, Ntime, p_f, IVM_Pregnancy, 
                                 dtau, tau_max, Vtau, Ntau, dsigma, sigma_max, Vsigma, Nsigma, strategy, 
                                 IVM_field_dependancy, tau50, alpha, Dmax)
        }
      }
    } else {
      VectTime_between_cycles = if(strategy %in% c(0, 1, 5)) c(60, 60, 60) else c(30, 30, 30)
      Number_of_cycle = nc
      ModelOutput = ModelIVM(PropIVM, VectTime_between_cycles, Number_of_cycle, t_begin_Camp, Dur_cycle,
                             mu_h, delta_h, bar_beta_h, beta_h, gamma_h, theta, beta_m, Wedge_m,
                             mu_m, ah_max, dah, Nah, tmax, Gap, dt, time, Ntime, p_f, IVM_Pregnancy,
                             dtau, tau_max, Vtau, Ntau, dsigma, sigma_max, Vsigma, Nsigma, strategy,
                             IVM_field_dependancy, tau50, alpha, Dmax)
    }
    
    GainPrev = 1 - ModelOutput$DeltaOptim / DeltaOptimInit
    GainPrev_15 = 1 - ModelOutput$DeltaOptim_15 / DeltaOptimInit_15
    
    if (return_model) {
      return(list(GainPrev = GainPrev, GainPrev_15 = GainPrev_15, Model = ModelOutput))
    } else {
      return(c(GainPrev, GainPrev_15))
    }
  }
  
  if (strategy %in% c(0, 1, 5)) {
    GainPrev_list = GainPrev_15_list = GainPrev_low_list = GainPrev_up_list = list()
    GainPrev_15_low_list = GainPrev_15_up_list = Model_list = list()
    
    for (nc in 1:NbCycle) {
      # Run avec paramètres estimés
      res_est = run_model(tau50_est, alpha_est, nc, return_model = TRUE)
      GainPrev_list[[nc]] = res_est$GainPrev
      GainPrev_15_list[[nc]] = res_est$GainPrev_15
      Model_list[[nc]] = res_est$Model
      
      # Runs pour IC 
      all_gains_list = future_lapply(1:nrow(param_grid), function(i) {
        run_model(param_grid$tau50[i], param_grid$alpha[i], nc, return_model = FALSE)
      }, future.seed = TRUE)
      
      all_gains = do.call(rbind, all_gains_list)
      GainPrev_low_list[[nc]] = min(all_gains[,1])
      GainPrev_up_list[[nc]] = max(all_gains[,1])
      GainPrev_15_low_list[[nc]] = min(all_gains[,2])
      GainPrev_15_up_list[[nc]] = max(all_gains[,2])
    }
    
    return(list(GainPrev = GainPrev_list, GainPrev_15 = GainPrev_15_list, 
                GainPrev_low = GainPrev_low_list, GainPrev_up = GainPrev_up_list,
                GainPrev_15_low = GainPrev_15_low_list, GainPrev_15_up = GainPrev_15_up_list, 
                Model = Model_list))
  } else {
    # Run avec paramètres estimés
    res_est = run_model(tau50_est, alpha_est, NbCycle, return_model = TRUE)
    
    # Runs pour IC
    all_gains_list = future_lapply(1:nrow(param_grid), function(i) {
      run_model(param_grid$tau50[i], param_grid$alpha[i], NbCycle, return_model = FALSE)
    }, future.seed = TRUE)
    
    all_gains = do.call(rbind, all_gains_list)
    
    return(list(GainPrev = res_est$GainPrev, GainPrev_15 = res_est$GainPrev_15,
                GainPrev_low = min(all_gains[,1]), GainPrev_up = max(all_gains[,1]),
                GainPrev_15_low = min(all_gains[,2]), GainPrev_15_up = max(all_gains[,2]),
                Model = res_est$Model))
  }
}

run_all_strategies=function(PropIVM, IVM_field_dependancy, scenario_type, CoefDelay,
                            coefTimingLAIF, DeltaOptimInit, DeltaOptimInit_15, IVM_Pregnancy,
                            Prev_mos10, oral_cycle_scenario = "initial") {
  NbCycle_oral = switch(oral_cycle_scenario, "initial" = NULL, "all3" = 3, "all4" = 4)
  results = list()
  
  for (s in c(0, 5, 1)) {
    res = GainPrev_with_CI(s, PropIVM, IVM_field_dependancy, scenario_type, CoefDelay, 
                           coefTimingLAIF, DeltaOptimInit, DeltaOptimInit_15, IVM_Pregnancy)
    nm = switch(as.character(s), "0" = "LongLasting06", "5" = "LongLasting06_kis", "1" = "LongLasting1")
    results[[paste0("Model", nm)]] = res$Model
    results[[paste0("GainPrev_", nm)]] = res$GainPrev
    results[[paste0("GainPrev_", nm, "_15")]] = res$GainPrev_15
    results[[paste0("GainPrev_", nm, "_low")]] = res$GainPrev_low 
    results[[paste0("GainPrev_", nm, "_up")]] = res$GainPrev_up
    results[[paste0("GainPrev_", nm, "_15_low")]] = res$GainPrev_15_low
    results[[paste0("GainPrev_", nm, "_15_up")]] = res$GainPrev_15_up
  }
  
  for (s in c(2, 3, 4)) {
    res = GainPrev_with_CI(s, PropIVM, IVM_field_dependancy, scenario_type, CoefDelay, coefTimingLAIF,
                           DeltaOptimInit, DeltaOptimInit_15, IVM_Pregnancy, NbCycle_oral)
    nm = switch(as.character(s), "2" = "Bohemia", "3" = "Rimdamal", "4" = "KamauRimdamal")
    results[[paste0("Model", nm)]] = res$Model
    results[[paste0("GainPrev_", nm)]] = res$GainPrev
    results[[paste0("GainPrev_", nm, "_15")]] = res$GainPrev_15
    results[[paste0("GainPrev_", nm, "_low")]] = res$GainPrev_low
    results[[paste0("GainPrev_", nm, "_up")]] = res$GainPrev_up
    results[[paste0("GainPrev_", nm, "_15_low")]] = res$GainPrev_15_low
    results[[paste0("GainPrev_", nm, "_15_up")]] = res$GainPrev_15_up
  }
  return(results)
}

run_scenario_IVM=function(scenario_type = "optimal", CoefDelay = 0.85, IVM_Pregnancy = 0,
                          Prev_mos10 = 1, VectIVM_field_dependancy = c(0), 
                          VectPropIVM = c(0.5, 0.7, 0.9), coefTimingLAIF = 1.25) {
  
  params_Mos_Prev(Prev_mos10)
  
  strategy = 0; PropIVM = 0; VectTime_between_cycles = c(0,0,0); Number_of_cycle = 1
  ParamIvmFormulation = ModelEtau50Alpha(strategy)
  tau50 = ParamIvmFormulation$tau50; alpha = ParamIvmFormulation$alpha; Dmax = ParamIvmFormulation$Dmax
  
  ModelNoIVM = ModelIVM(PropIVM, VectTime_between_cycles, Number_of_cycle, t_begin_Camp, Dur_cycle,
                        mu_h, delta_h, bar_beta_h, beta_h, gamma_h, theta, beta_m, Wedge_m, mu_m,
                        ah_max, dah, Nah, tmax, Gap, dt, time, Ntime, p_f, IVM_Pregnancy, dtau, tau_max, Vtau, Ntau,
                        dsigma, sigma_max, Vsigma, Nsigma, strategy, 0, tau50, alpha, Dmax)
  
  DeltaOptimInit = sum(ModelNoIVM$Ih_Tot[ModelNoIVM$IndexOptim])
  DeltaOptimInit_15 = sum(ModelNoIVM$Ih_15[ModelNoIVM$IndexOptim])
  
  VectNumber_of_cycle = c(1, 2, 3, 4)
  nf = length(VectIVM_field_dependancy); np = length(VectPropIVM); nnc = length(VectNumber_of_cycle)
  
  # Arrays scenario initial
  GainPrev_array_oral = GainPrev_array_oral_low = GainPrev_array_oral_up = array(0, dim = c(3, np, nf))
  GainPrev_array_oral_15 = GainPrev_array_oral_15_low = GainPrev_array_oral_15_up = array(0, dim = c(3, np, nf))
  GainPrev_array_LAIF = GainPrev_array_LAIF_low = GainPrev_array_LAIF_up = array(0, dim = c(3, np, nnc, nf))
  GainPrev_array_LAIF_15 = GainPrev_array_LAIF_15_low = GainPrev_array_LAIF_15_up = array(0, dim = c(3, np, nnc, nf))
  
  # Arrays scenario all3
  GainPrev_array_oral_all3 = GainPrev_array_oral_all3_low = GainPrev_array_oral_all3_up = array(0, dim = c(3, np, nf))
  GainPrev_array_oral_all3_15 = GainPrev_array_oral_all3_15_low = GainPrev_array_oral_all3_15_up = array(0, dim = c(3, np, nf))
  
  # Arrays scenario all4
  GainPrev_array_oral_all4 = GainPrev_array_oral_all4_low = GainPrev_array_oral_all4_up = array(0, dim = c(3, np, nf))
  GainPrev_array_oral_all4_15 = GainPrev_array_oral_all4_15_low = GainPrev_array_oral_all4_15_up = array(0, dim = c(3, np, nf))
  
  all_results = list()
  
  for (field_idx in 1:nf) {
    IVM_field_dependancy = VectIVM_field_dependancy[field_idx]
    
    for (j in 1:np) {
      PropIVM = VectPropIVM[j]
      
      # Scenario initial (BOHEMIA=3, RIMDAMAL=4)
      results = run_all_strategies(PropIVM, IVM_field_dependancy, scenario_type, CoefDelay,
                                   coefTimingLAIF, DeltaOptimInit, DeltaOptimInit_15, IVM_Pregnancy, Prev_mos10, "initial")
      
      # Stocker oral - initial
      GainPrev_array_oral[1, j, field_idx] = results$GainPrev_Bohemia
      GainPrev_array_oral[2, j, field_idx] = results$GainPrev_Rimdamal
      GainPrev_array_oral[3, j, field_idx] = results$GainPrev_KamauRimdamal
      GainPrev_array_oral_low[1, j, field_idx] = results$GainPrev_Bohemia_low
      GainPrev_array_oral_low[2, j, field_idx] = results$GainPrev_Rimdamal_low
      GainPrev_array_oral_low[3, j, field_idx] = results$GainPrev_KamauRimdamal_low
      GainPrev_array_oral_up[1, j, field_idx] = results$GainPrev_Bohemia_up
      GainPrev_array_oral_up[2, j, field_idx] = results$GainPrev_Rimdamal_up
      GainPrev_array_oral_up[3, j, field_idx] = results$GainPrev_KamauRimdamal_up
      GainPrev_array_oral_15[1, j, field_idx] = results$GainPrev_Bohemia_15
      GainPrev_array_oral_15[2, j, field_idx] = results$GainPrev_Rimdamal_15
      GainPrev_array_oral_15[3, j, field_idx] = results$GainPrev_KamauRimdamal_15
      GainPrev_array_oral_15_low[1, j, field_idx] = results$GainPrev_Bohemia_15_low
      GainPrev_array_oral_15_low[2, j, field_idx] = results$GainPrev_Rimdamal_15_low
      GainPrev_array_oral_15_low[3, j, field_idx] = results$GainPrev_KamauRimdamal_15_low
      GainPrev_array_oral_15_up[1, j, field_idx] = results$GainPrev_Bohemia_15_up
      GainPrev_array_oral_15_up[2, j, field_idx] = results$GainPrev_Rimdamal_15_up
      GainPrev_array_oral_15_up[3, j, field_idx] = results$GainPrev_KamauRimdamal_15_up
      
      # Stocker LAIF
      for (nc in 1:4) {
        GainPrev_array_LAIF[1, j, nc, field_idx] = results$GainPrev_LongLasting06[[nc]]
        GainPrev_array_LAIF[2, j, nc, field_idx] = results$GainPrev_LongLasting1[[nc]]
        GainPrev_array_LAIF[3, j, nc, field_idx] = results$GainPrev_LongLasting06_kis[[nc]]
        GainPrev_array_LAIF_low[1, j, nc, field_idx] = results$GainPrev_LongLasting06_low[[nc]]
        GainPrev_array_LAIF_low[2, j, nc, field_idx] = results$GainPrev_LongLasting1_low[[nc]]
        GainPrev_array_LAIF_low[3, j, nc, field_idx] = results$GainPrev_LongLasting06_kis_low[[nc]]
        GainPrev_array_LAIF_up[1, j, nc, field_idx] = results$GainPrev_LongLasting06_up[[nc]]
        GainPrev_array_LAIF_up[2, j, nc, field_idx] = results$GainPrev_LongLasting1_up[[nc]]
        GainPrev_array_LAIF_up[3, j, nc, field_idx] = results$GainPrev_LongLasting06_kis_up[[nc]]
        GainPrev_array_LAIF_15[1, j, nc, field_idx] = results$GainPrev_LongLasting06_15[[nc]]
        GainPrev_array_LAIF_15[2, j, nc, field_idx] = results$GainPrev_LongLasting1_15[[nc]]
        GainPrev_array_LAIF_15[3, j, nc, field_idx] = results$GainPrev_LongLasting06_kis_15[[nc]]
        GainPrev_array_LAIF_15_low[1, j, nc, field_idx] = results$GainPrev_LongLasting06_15_low[[nc]]
        GainPrev_array_LAIF_15_low[2, j, nc, field_idx] = results$GainPrev_LongLasting1_15_low[[nc]]
        GainPrev_array_LAIF_15_low[3, j, nc, field_idx] = results$GainPrev_LongLasting06_kis_15_low[[nc]]
        GainPrev_array_LAIF_15_up[1, j, nc, field_idx] = results$GainPrev_LongLasting06_15_up[[nc]]
        GainPrev_array_LAIF_15_up[2, j, nc, field_idx] = results$GainPrev_LongLasting1_15_up[[nc]]
        GainPrev_array_LAIF_15_up[3, j, nc, field_idx] = results$GainPrev_LongLasting06_kis_15_up[[nc]]
      }
      
      # Sauvegarder resultats individuels avec ModelNoIVM
      if (scenario_type == "optimal") {
        RDataFileName = paste0("Results_Ihopt_IVM_field_dep_", IVM_field_dependancy, "_PropIVM_",
                               gsub("\\.", "_", PropIVM), "_Pregnant", IVM_Pregnancy, "Prev_mos",
                               Name_Prev_mos(Prev_mos10), "_Timing", gsub("\\.", "_", CoefDelay), ".RData")
      } else {
        RDataFileName = paste0("Results_IhBase_IVM_field_dep_", IVM_field_dependancy, "_PropIVM_",
                               gsub("\\.", "_", PropIVM), "_Pregnant", IVM_Pregnancy, "Prev_mos", Name_Prev_mos(Prev_mos10), ".RData")
      }
      
      simulation_results = list(PropIVM = PropIVM, IVM_field_dependancy = IVM_field_dependancy,
                                ModelBohemia = results$ModelBohemia, ModelRimdamal = results$ModelRimdamal,
                                ModelKamauRimdamal = results$ModelKamauRimdamal, ModelLongLasting06 = results$ModelLongLasting06,
                                ModelLongLasting1 = results$ModelLongLasting1,ModelLongLasting06_kis = results$ModelLongLasting06_kis,
                                ModelNoIVM = ModelNoIVM, time = time)
      
      save(simulation_results, file = RDataFileName)
      all_results[[paste0("field", field_idx, "_prop", j)]] = simulation_results
      
      # Scenario all3
      results_all3 = run_all_strategies(PropIVM, IVM_field_dependancy, scenario_type, CoefDelay,
                                        coefTimingLAIF, DeltaOptimInit, DeltaOptimInit_15, IVM_Pregnancy, Prev_mos10, "all3")
      GainPrev_array_oral_all3[1, j, field_idx] = results_all3$GainPrev_Bohemia
      GainPrev_array_oral_all3[2, j, field_idx] = results_all3$GainPrev_Rimdamal
      GainPrev_array_oral_all3[3, j, field_idx] = results_all3$GainPrev_KamauRimdamal
      GainPrev_array_oral_all3_low[1, j, field_idx] = results_all3$GainPrev_Bohemia_low
      GainPrev_array_oral_all3_low[2, j, field_idx] = results_all3$GainPrev_Rimdamal_low
      GainPrev_array_oral_all3_low[3, j, field_idx] = results_all3$GainPrev_KamauRimdamal_low
      GainPrev_array_oral_all3_up[1, j, field_idx] = results_all3$GainPrev_Bohemia_up
      GainPrev_array_oral_all3_up[2, j, field_idx] = results_all3$GainPrev_Rimdamal_up
      GainPrev_array_oral_all3_up[3, j, field_idx] = results_all3$GainPrev_KamauRimdamal_up
      GainPrev_array_oral_all3_15[1, j, field_idx] = results_all3$GainPrev_Bohemia_15
      GainPrev_array_oral_all3_15[2, j, field_idx] = results_all3$GainPrev_Rimdamal_15
      GainPrev_array_oral_all3_15[3, j, field_idx] = results_all3$GainPrev_KamauRimdamal_15
      GainPrev_array_oral_all3_15_low[1, j, field_idx] = results_all3$GainPrev_Bohemia_15_low
      GainPrev_array_oral_all3_15_low[2, j, field_idx] = results_all3$GainPrev_Rimdamal_15_low
      GainPrev_array_oral_all3_15_low[3, j, field_idx] = results_all3$GainPrev_KamauRimdamal_15_low
      GainPrev_array_oral_all3_15_up[1, j, field_idx] = results_all3$GainPrev_Bohemia_15_up
      GainPrev_array_oral_all3_15_up[2, j, field_idx] = results_all3$GainPrev_Rimdamal_15_up
      GainPrev_array_oral_all3_15_up[3, j, field_idx] = results_all3$GainPrev_KamauRimdamal_15_up
      
      # Scenario all4
      results_all4 = run_all_strategies(PropIVM, IVM_field_dependancy, scenario_type, CoefDelay,
                                        coefTimingLAIF, DeltaOptimInit, DeltaOptimInit_15, IVM_Pregnancy, Prev_mos10, "all4")
      GainPrev_array_oral_all4[1, j, field_idx] = results_all4$GainPrev_Bohemia
      GainPrev_array_oral_all4[2, j, field_idx] = results_all4$GainPrev_Rimdamal
      GainPrev_array_oral_all4[3, j, field_idx] = results_all4$GainPrev_KamauRimdamal
      GainPrev_array_oral_all4_low[1, j, field_idx] = results_all4$GainPrev_Bohemia_low
      GainPrev_array_oral_all4_low[2, j, field_idx] = results_all4$GainPrev_Rimdamal_low
      GainPrev_array_oral_all4_low[3, j, field_idx] = results_all4$GainPrev_KamauRimdamal_low
      GainPrev_array_oral_all4_up[1, j, field_idx] = results_all4$GainPrev_Bohemia_up
      GainPrev_array_oral_all4_up[2, j, field_idx] = results_all4$GainPrev_Rimdamal_up
      GainPrev_array_oral_all4_up[3, j, field_idx] = results_all4$GainPrev_KamauRimdamal_up
      GainPrev_array_oral_all4_15[1, j, field_idx] = results_all4$GainPrev_Bohemia_15
      GainPrev_array_oral_all4_15[2, j, field_idx] = results_all4$GainPrev_Rimdamal_15
      GainPrev_array_oral_all4_15[3, j, field_idx] = results_all4$GainPrev_KamauRimdamal_15
      GainPrev_array_oral_all4_15_low[1, j, field_idx] = results_all4$GainPrev_Bohemia_15_low
      GainPrev_array_oral_all4_15_low[2, j, field_idx] = results_all4$GainPrev_Rimdamal_15_low
      GainPrev_array_oral_all4_15_low[3, j, field_idx] = results_all4$GainPrev_KamauRimdamal_15_low
      GainPrev_array_oral_all4_15_up[1, j, field_idx] = results_all4$GainPrev_Bohemia_15_up
      GainPrev_array_oral_all4_15_up[2, j, field_idx] = results_all4$GainPrev_Rimdamal_15_up
      GainPrev_array_oral_all4_15_up[3, j, field_idx] = results_all4$GainPrev_KamauRimdamal_15_up
    }
  }
  
  # Sauvegarder tous les arrays
  if (scenario_type == "optimal") {
    save_file = paste0("ResultatsPrev_Ihopt_field_dep_", paste(VectIVM_field_dependancy, collapse="_"),
                       "_Pregnant", IVM_Pregnancy, "Prev_mos", Name_Prev_mos(Prev_mos10), ".RData")
  } else {
    save_file = paste0("ResultatsPrev_IhBase_field_dep_", paste(VectIVM_field_dependancy, collapse="_"),
                       "_Pregnant", IVM_Pregnancy, "Prev_mos", Name_Prev_mos(Prev_mos10), ".RData")
  }
  
  save(GainPrev_array_oral, GainPrev_array_oral_low, GainPrev_array_oral_up,
       GainPrev_array_oral_15, GainPrev_array_oral_15_low, GainPrev_array_oral_15_up,
       GainPrev_array_LAIF, GainPrev_array_LAIF_low, GainPrev_array_LAIF_up,
       GainPrev_array_LAIF_15, GainPrev_array_LAIF_15_low, GainPrev_array_LAIF_15_up,
       GainPrev_array_oral_all3, GainPrev_array_oral_all3_low, GainPrev_array_oral_all3_up,
       GainPrev_array_oral_all3_15, GainPrev_array_oral_all3_15_low, GainPrev_array_oral_all3_15_up,
       GainPrev_array_oral_all4, GainPrev_array_oral_all4_low, GainPrev_array_oral_all4_up,
       GainPrev_array_oral_all4_15, GainPrev_array_oral_all4_15_low, GainPrev_array_oral_all4_15_up,
       VectPropIVM, VectNumber_of_cycle, VectIVM_field_dependancy,
       DeltaOptimInit, DeltaOptimInit_15, file = save_file)
  
  # Export TXT pour scenario initial
  export_GainPrev_to_txt(scenario_type, IVM_Pregnancy, Prev_mos10, VectIVM_field_dependancy,
                         VectPropIVM, VectNumber_of_cycle,
                         GainPrev_array_oral, GainPrev_array_oral_low, GainPrev_array_oral_up,
                         GainPrev_array_oral_15, GainPrev_array_oral_15_low, GainPrev_array_oral_15_up,
                         GainPrev_array_LAIF, GainPrev_array_LAIF_low, GainPrev_array_LAIF_up,
                         GainPrev_array_LAIF_15, GainPrev_array_LAIF_15_low, GainPrev_array_LAIF_15_up,
                         oral_cycle_scenario = "initial")
  
  # Export TXT pour scenario all3
  export_GainPrev_to_txt(scenario_type, IVM_Pregnancy, Prev_mos10, VectIVM_field_dependancy,
                         VectPropIVM, VectNumber_of_cycle,
                         GainPrev_array_oral_all3, GainPrev_array_oral_all3_low, GainPrev_array_oral_all3_up,
                         GainPrev_array_oral_all3_15, GainPrev_array_oral_all3_15_low, GainPrev_array_oral_all3_15_up,
                         GainPrev_array_LAIF, GainPrev_array_LAIF_low, GainPrev_array_LAIF_up,
                         GainPrev_array_LAIF_15, GainPrev_array_LAIF_15_low, GainPrev_array_LAIF_15_up,
                         oral_cycle_scenario = "all3")
  
  # Export TXT pour scenario all4
  export_GainPrev_to_txt(scenario_type, IVM_Pregnancy, Prev_mos10, VectIVM_field_dependancy,
                         VectPropIVM, VectNumber_of_cycle,
                         GainPrev_array_oral_all4, GainPrev_array_oral_all4_low, GainPrev_array_oral_all4_up,
                         GainPrev_array_oral_all4_15, GainPrev_array_oral_all4_15_low, GainPrev_array_oral_all4_15_up,
                         GainPrev_array_LAIF, GainPrev_array_LAIF_low, GainPrev_array_LAIF_up,
                         GainPrev_array_LAIF_15, GainPrev_array_LAIF_15_low, GainPrev_array_LAIF_15_up,
                         oral_cycle_scenario = "all4")
  
  return(list(all_results = all_results,
              GainPrev_array_oral = GainPrev_array_oral, GainPrev_array_LAIF = GainPrev_array_LAIF,
              GainPrev_array_oral_15 = GainPrev_array_oral_15, GainPrev_array_LAIF_15 = GainPrev_array_LAIF_15))
}

plot_comparison_generic=function(comparison_type, IVM_Pregnancy = 0, Prev_mos10 = 1,
                                 IVM_field_dependancy = 0, VectNumber_of_cycle = c(1,2,3,4)) {
  
  if (comparison_type == "opt_vs_base") {
    load(paste0("ResultatsPrev_Ihopt_field_dep_", IVM_field_dependancy, "_Pregnant", IVM_Pregnancy,
                "Prev_mos", Name_Prev_mos(Prev_mos10), ".RData"))
    GainPrevopt_array_oral = GainPrev_array_oral; GainPrevopt_array_LAIF = GainPrev_array_LAIF; VectPropIVM_opt = VectPropIVM
    
    load(paste0("ResultatsPrev_IhBase_field_dep_", IVM_field_dependancy, "_Pregnant", IVM_Pregnancy,
                "Prev_mos", Name_Prev_mos(Prev_mos10), ".RData"))
    GainPrevBase_array_oral = GainPrev_array_oral; GainPrevBase_array_LAIF = GainPrev_array_LAIF; VectPropIVM = VectPropIVM_opt
    
    comparison_label = "Optimal_vs_Baseline"; y_label = expression("Performance (optimal - baseline) (%)"); GrNumber = 0
    
  } else if (comparison_type == "field_dep") {
    load(paste0("ResultatsPrev_Ihopt_field_dep_1_Pregnant", IVM_Pregnancy, "Prev_mos", Name_Prev_mos(Prev_mos10), ".RData"))
    GainPrev1_array_oral = GainPrev_array_oral; GainPrev1_array_LAIF = GainPrev_array_LAIF; VectPropIVM_1 = VectPropIVM
    
    load(paste0("ResultatsPrev_Ihopt_field_dep_0_Pregnant", IVM_Pregnancy, "Prev_mos", Name_Prev_mos(Prev_mos10), ".RData"))
    GainPrev2_array_oral = GainPrev_array_oral; GainPrev2_array_LAIF = GainPrev_array_LAIF; VectPropIVM = VectPropIVM_1
    
    comparison_label = "Field_dependency"; y_label = expression("Performance IVM field dependency (%)"); GrNumber = 4
    
  } else if (comparison_type == "prev_mos") {
    load(paste0("ResultatsPrev_Ihopt_field_dep_", IVM_field_dependancy, "_Pregnant", IVM_Pregnancy, "Prev_mos10.RData"))
    GainPrev1_array_oral = GainPrev_array_oral; GainPrev1_array_LAIF = GainPrev_array_LAIF; VectPropIVM_1 = VectPropIVM
    
    load(paste0("ResultatsPrev_Ihopt_field_dep_", IVM_field_dependancy, "_Pregnant", IVM_Pregnancy, "Prev_mos5.RData"))
    GainPrev2_array_oral = GainPrev_array_oral; GainPrev2_array_LAIF = GainPrev_array_LAIF; VectPropIVM = VectPropIVM_1
    
    comparison_label = "Mosquito_prevalence"; y_label = expression("Performance mosquito prevalence (%)"); GrNumber = 0
    
  } else if (comparison_type == "pregnancy") {
    load(paste0("ResultatsPrev_Ihopt_field_dep_", IVM_field_dependancy, "_Pregnant1Prev_mos", Name_Prev_mos(Prev_mos10), ".RData"))
    GainPrev1_array_oral = GainPrev_array_oral; GainPrev1_array_LAIF = GainPrev_array_LAIF; VectPropIVM_1 = VectPropIVM
    
    load(paste0("ResultatsPrev_Ihopt_field_dep_", IVM_field_dependancy, "_Pregnant0Prev_mos", Name_Prev_mos(Prev_mos10), ".RData"))
    GainPrev2_array_oral = GainPrev_array_oral; GainPrev2_array_LAIF = GainPrev_array_LAIF; VectPropIVM = VectPropIVM_1
    
    comparison_label = "Pregnancy"; y_label = expression("Performance with pregnancy (%)"); GrNumber = 0
    
  } else if (comparison_type == "deltoptim_vs_15") {
    load(paste0("ResultatsPrev_Ihopt_field_dep_", IVM_field_dependancy, "_Pregnant", IVM_Pregnancy,
                "Prev_mos", Name_Prev_mos(Prev_mos10), ".RData"))
    GainPrev2_array_oral = GainPrev_array_oral; GainPrev2_array_LAIF = GainPrev_array_LAIF
    GainPrev1_array_oral = GainPrev_array_oral_15; GainPrev1_array_LAIF = GainPrev_array_LAIF_15
    
    comparison_label = "AllAges_vs_Under5"; y_label = expression("Performance ([5-15] years - all ages) (%)"); GrNumber = 4
  }
  
  # Calcul ylim
  all_values = c()
  for (nc in VectNumber_of_cycle) {
    for (oral_strat in 1:3) {
      for (prop_idx in 1:length(VectPropIVM)) {
        if (comparison_type == "opt_vs_base") {
          if (length(dim(GainPrevopt_array_oral)) == 2) val = (GainPrevopt_array_oral[oral_strat, prop_idx] - GainPrevBase_array_oral[oral_strat, prop_idx]) * 100
          else val = (GainPrevopt_array_oral[oral_strat, prop_idx, 1] - GainPrevBase_array_oral[oral_strat, prop_idx, 1]) * 100
        } else {
          if (length(dim(GainPrev1_array_oral)) == 2) val = (GainPrev1_array_oral[oral_strat, prop_idx] - GainPrev2_array_oral[oral_strat, prop_idx]) * 100
          else val = (GainPrev1_array_oral[oral_strat, prop_idx, 1] - GainPrev2_array_oral[oral_strat, prop_idx, 1]) * 100
        }
        all_values = c(all_values, val)
      }
    }
    for (laif_strat in 1:3) {
      for (prop_idx in 1:length(VectPropIVM)) {
        if (comparison_type == "opt_vs_base") {
          if (length(dim(GainPrevopt_array_LAIF)) == 3) val = (GainPrevopt_array_LAIF[laif_strat, prop_idx, nc] - GainPrevBase_array_LAIF[laif_strat, prop_idx, nc]) * 100
          else val = (GainPrevopt_array_LAIF[laif_strat, prop_idx, nc, 1] - GainPrevBase_array_LAIF[laif_strat, prop_idx, nc, 1]) * 100
        } else {
          if (length(dim(GainPrev1_array_LAIF)) == 3) val = (GainPrev1_array_LAIF[laif_strat, prop_idx, nc] - GainPrev2_array_LAIF[laif_strat, prop_idx, nc]) * 100
          else val = (GainPrev1_array_LAIF[laif_strat, prop_idx, nc, 1] - GainPrev2_array_LAIF[laif_strat, prop_idx, nc, 1]) * 100
        }
        all_values = c(all_values, val)
      }
    }
  }
  
  all_values = all_values[!is.na(all_values)]
  y_min0 = min(all_values); y_min = y_min0 * ifelse(y_min0 < 0, 1.4, 0.75); y_max = max(all_values) * 1.4
  
  FigName = paste0("Comparison_", comparison_label, "_field_dep", IVM_field_dependancy,
                   "_Pregnant", IVM_Pregnancy, "Prev_mos", Name_Prev_mos(Prev_mos10), ".pdf")
  
  pdf(FigName, width = 14, height = 4.5)
  par(mfrow = c(1, 4), oma = c(9, 1, 2, 0.2))
  
  for (nc in VectNumber_of_cycle) {
    GrNumber = GrNumber + 1
    if(GrNumber %in% c(1, 5, 9, 13, 17)) mar = c(4.1, 3.35, 1.5, .1) else mar = c(4.1, 2.6, 1.5, .1)
    par(mar = mar)
    
    Summary_stats = data.frame()
    
    for (oral_strat in 1:3) {
      values_strat = c()
      for (prop_idx in 1:length(VectPropIVM)) {
        if (comparison_type == "opt_vs_base") {
          if (length(dim(GainPrevopt_array_oral)) == 2) val = (GainPrevopt_array_oral[oral_strat, prop_idx] - GainPrevBase_array_oral[oral_strat, prop_idx]) * 100
          else val = (GainPrevopt_array_oral[oral_strat, prop_idx, 1] - GainPrevBase_array_oral[oral_strat, prop_idx, 1]) * 100
        } else {
          if (length(dim(GainPrev1_array_oral)) == 2) val = (GainPrev1_array_oral[oral_strat, prop_idx] - GainPrev2_array_oral[oral_strat, prop_idx]) * 100
          else val = (GainPrev1_array_oral[oral_strat, prop_idx, 1] - GainPrev2_array_oral[oral_strat, prop_idx, 1]) * 100
        }
        values_strat = c(values_strat, val)
      }
      Summary_stats = rbind(Summary_stats, data.frame(Strategy = oral_strat, Min = min(values_strat, na.rm = TRUE),
                                                      Max = max(values_strat, na.rm = TRUE), Median = median(values_strat, na.rm = TRUE)))
    }
    
    for (laif_strat in 1:3) {
      values_strat = c()
      for (prop_idx in 1:length(VectPropIVM)) {
        if (comparison_type == "opt_vs_base") {
          if (length(dim(GainPrevopt_array_LAIF)) == 3) val = (GainPrevopt_array_LAIF[laif_strat, prop_idx, nc] - GainPrevBase_array_LAIF[laif_strat, prop_idx, nc]) * 100
          else val = (GainPrevopt_array_LAIF[laif_strat, prop_idx, nc, 1] - GainPrevBase_array_LAIF[laif_strat, prop_idx, nc, 1]) * 100
        } else {
          if (length(dim(GainPrev1_array_LAIF)) == 3) val = (GainPrev1_array_LAIF[laif_strat, prop_idx, nc] - GainPrev2_array_LAIF[laif_strat, prop_idx, nc]) * 100
          else val = (GainPrev1_array_LAIF[laif_strat, prop_idx, nc, 1] - GainPrev2_array_LAIF[laif_strat, prop_idx, nc, 1]) * 100
        }
        values_strat = c(values_strat, val)
      }
      Summary_stats = rbind(Summary_stats, data.frame(Strategy = laif_strat + 3, Min = min(values_strat, na.rm = TRUE),
                                                      Max = max(values_strat, na.rm = TRUE), Median = median(values_strat, na.rm = TRUE)))
    }
    
    x_labels = c("Delta(BOH,BOH)", "Delta(RII[S],RII[S])", "Delta(RII[K],RII[K])",
                 "Delta(0.6-vk5,0.6-vk5)", "Delta('1.0','1.0')", "Delta(0.6-kis,0.6-kis)")
    
    if (comparison_type == "opt_vs_base") { y_min = -5.5; y_max = 5.5 }
    y_min_lim = ifelse(y_min < 0, y_min, -6); y_max_lim = ifelse(y_max > 0, y_max, 6)
    ticks = sort(unique(c(pretty(c(y_min_lim, y_max_lim), n=6), 0)))
    
    plot(-1, 1, type="n", xlim=c(0.5, 6.5), ylim=c(y_min_lim, y_max_lim), xaxt="n", yaxt="n", xlab="", ylab="")
    abline(h=0, lty=2, col="gray"); box(); axis(2, at=ticks); axis(1, at = 1:6, labels = FALSE)
    
    for (i in 1:length(x_labels)) mtext(parse(text=x_labels[i]), side=1, line=1, at=i, cex=0.7, las=2)
    
    if(GrNumber %in% c(1, 5, 9, 13, 17)) mtext(y_label, side = 2, adj = 0.5, cex = .8, line = 2)
    
    strategy_colors = c("#D98E04", "#205072", "#8B008B", "#ff3355", "#037153", "#8B4513")
    
    for (i in 1:nrow(Summary_stats)) {
      strat = Summary_stats$Strategy[i]
      arrows(strat, Summary_stats$Min[i], strat, Summary_stats$Max[i], angle = 90, code = 3, length = 0.02, col = strategy_colors[strat], lwd = 2)
      points(strat, Summary_stats$Median[i], pch = 15, col = strategy_colors[strat], cex = 2)
    }
    
    par(xpd=NA)
    if (comparison_type %in% c("field_dep", "prev_mos")) text(0.2, ifelse(y_max > 0, y_max, 6)*1.5, paste0("(", LETTERS[GrNumber], ")"), cex=1.5, adj=0)
    else text(0.2, ifelse(y_max > 0, y_max, 6)*1.2, paste0("(", LETTERS[GrNumber], ")"), cex=1.5, adj=0)
    mtext(bquote("#cycle for mdc-STM-001: "~ n[c]==.(nc)), side=3, cex = 0.95, line=1)
    par(xpd=FALSE)
  }
  
  if (comparison_type %in% c("field_dep", "deltoptim_vs_15")) {
    reset()
    LEGEND = c("mdc-STM-001-0.6-vk5", "mdc-STM-001-0.6-kis", "mdc-STM-001-1.0",
               expression("RIMDAMAL II"[S]), expression("RIMDAMAL II"[K]), "BOHEMIA")
    strategy_colors = c("#ff3355", "#8B4513", "#037153", "#205072", "#8B008B", "#D98E04")
    legend("bottom", legend = LEGEND, col = strategy_colors, pch = 15, xpd = NA, horiz = FALSE, inset = c(0, 0.01), ncol = 3, cex = 1, bty = "n")
  }
  dev.off()
}

plot_formulations_comparison=function(IVM_field_dependancy = 0, IVM_Pregnancy = 0,
                                      Prev_mos10 = 1) {
  load_file = paste0("ResultatsPrev_Ihopt_field_dep_", IVM_field_dependancy,
                     "_Pregnant", IVM_Pregnancy, "Prev_mos", Name_Prev_mos(Prev_mos10), ".RData")
  load(load_file)
  
  if (length(dim(GainPrev_array_oral)) == 3) {
    GainPrev2_array_oral = GainPrev_array_oral[,,1]
    GainPrev2_array_LAIF = GainPrev_array_LAIF[,,,1]
  } else {
    GainPrev2_array_oral = GainPrev_array_oral
    GainPrev2_array_LAIF = GainPrev_array_LAIF
  }
  
  FigName = paste0("Comparisons_Strategies_points_field_dep", IVM_field_dependancy, "_Pregnant",
                   IVM_Pregnancy, "Prev_mos", Name_Prev_mos(Prev_mos10), ".pdf")
  
  pdf(FigName, width = 7.5, height = 6.2)
  par(mfrow = c(2, 2), oma = c(7, 1, 2, 0.2))
  
  colors = c("#ff3355", "#037153", "#D98E04")
  names(colors) = c("Oral vs Oral", "mdc-STM-001 vs Oral", "mdc-STM-001 vs mdc-STM-001")
  Vectpch = c(1, 15, 16)
  names(Vectpch) = c("Oral vs Oral", "mdc-STM-001 vs Oral", "mdc-STM-001 vs mdc-STM-001")
  
  VectPropIVM_target = c(0.5, 0.7, 0.9)
  strategies_LAIF_names = c("0.6-vk5", "1.0", "0.6-kis")
  GrNumber = 0
  
  for (nc in 1:4) {
    GrNumber = GrNumber + 1
    par(mar = c(2, 3.35, 1.5, .1))
    
    All_comparisons = data.frame()
    
    for (PropIVM_target in VectPropIVM_target) {
      IndexPropIVM = which.min(abs(VectPropIVM - PropIVM_target))
      
      # Oral vs Oral
      val1 = (GainPrev2_array_oral[2, IndexPropIVM] - GainPrev2_array_oral[1, IndexPropIVM]) * 100
      val2 = (GainPrev2_array_oral[3, IndexPropIVM] - GainPrev2_array_oral[1, IndexPropIVM]) * 100
      val3 = (GainPrev2_array_oral[2, IndexPropIVM] - GainPrev2_array_oral[3, IndexPropIVM]) * 100
      
      All_comparisons = rbind(All_comparisons,
                              data.frame(PropIVM=PropIVM_target, Comparison="Δ(RII[S],BOH)", Value=val1, Type="Oral vs Oral"),
                              data.frame(PropIVM=PropIVM_target, Comparison="Δ(RII[K],BOH)", Value=val2, Type="Oral vs Oral"),
                              data.frame(PropIVM=PropIVM_target, Comparison="Δ(RII[S],RII[K])", Value=val3, Type="Oral vs Oral"))
      
      # mdc-STM-001 vs Oral
      for (laif in 1:3) {
        for (oral in 1:3) {
          val = (GainPrev2_array_LAIF[laif, IndexPropIVM, nc] - GainPrev2_array_oral[oral, IndexPropIVM]) * 100
          oral_name = c("BOH", "RII[S]", "RII[K]")[oral]
          laif_name = strategies_LAIF_names[laif]
          comp_name = paste0("Δ(", laif_name, ",", oral_name, ")")
          All_comparisons = rbind(All_comparisons,
                                  data.frame(PropIVM=PropIVM_target, Comparison=comp_name, Value=val, Type="mdc-STM-001 vs Oral"))
        }
      }
      
      # mdc-STM-001 vs mdc-STM-001
      val4 = (GainPrev2_array_LAIF[2, IndexPropIVM, nc] - GainPrev2_array_LAIF[1, IndexPropIVM, nc]) * 100
      val5 = (GainPrev2_array_LAIF[2, IndexPropIVM, nc] - GainPrev2_array_LAIF[3, IndexPropIVM, nc]) * 100
      val6 = (GainPrev2_array_LAIF[3, IndexPropIVM, nc] - GainPrev2_array_LAIF[1, IndexPropIVM, nc]) * 100
      
      All_comparisons = rbind(All_comparisons,
                              data.frame(PropIVM=PropIVM_target, Comparison="Δ(1.0,0.6-vk5)", Value=val4, Type="mdc-STM-001 vs mdc-STM-001"),
                              data.frame(PropIVM=PropIVM_target, Comparison="Δ(1.0,0.6-kis)", Value=val5, Type="mdc-STM-001 vs mdc-STM-001"),
                              data.frame(PropIVM=PropIVM_target, Comparison="Δ(0.6-kis,0.6-vk5)", Value=val6, Type="mdc-STM-001 vs mdc-STM-001"))
    }
    
    Comparisons_list = unique(All_comparisons$Comparison)
    Summary_stats = data.frame()
    
    for (comp in Comparisons_list) {
      df_comp = subset(All_comparisons, Comparison == comp)
      Summary_stats = rbind(Summary_stats, data.frame(Comparison=comp, Min=min(df_comp$Value), Max=max(df_comp$Value),
                                                      Median=median(df_comp$Value), Type=df_comp$Type[1]))
    }
    
    ylim = c(-12, 68)
    plot(-1, 1, type="n", xlim=c(0.5, length(Comparisons_list)+0.5), ylim=ylim, xaxt="n", yaxt="n", xlab="", ylab="")
    abline(h=0, lty=2, col="gray"); box()
    ticks = c(-10, 0, 20, 40, 60) 
    axis(2, at=ticks, cex.axis = 0.7)
    
    if(GrNumber %in% c(1, 3)) mtext("Performance of formulations (all ages; %)", side = 2, adj = 0.5, cex = .78, line = 2)
    
    axis(1, at=1:length(Comparisons_list), labels=FALSE)
    if(GrNumber %in% c(3, 4)) {
      for (i in 1:length(Comparisons_list)) {
        label_text = Comparisons_list[i]
        label_text = gsub("Δ\\(", "Delta(", label_text)
        label_text = gsub("RII\\[S\\]", "RII[S]", label_text)
        label_text = gsub("RII\\[K\\]", "RII[K]", label_text)
        label_text = gsub("1\\.0", "'1.0'", label_text)
        mtext(parse(text=label_text), side=1, line=1, at=i, cex=0.7, las=2)
      }
    }
    
    for (i in 1:nrow(Summary_stats)) {
      row = Summary_stats[i,]
      arrows(i, row$Min, i, row$Max, angle = 90, code = 3, length = 0.015, col = colors[row$Type], lwd = 2)
      points(i, row$Median, pch = Vectpch[row$Type], bg = colors[row$Type], col = colors[row$Type], cex = 1.6)
    }
    
    par(xpd=NA)
    text(0.2, ylim[2]*1.15, paste0("(", LETTERS[GrNumber], ")"), cex=1.3, adj=0)
    mtext(bquote("#cycle for mdc-STM-001: "~ n[c]==.(nc)), side=3, cex = 0.95, line=1)
    par(xpd=FALSE)
  }
  
  reset()
  legend("bottom", legend=names(colors), col=colors, pch=Vectpch, xpd=NA, horiz=TRUE, inset=c(0,-0.01), cex=0.8, bty="n")
  dev.off()
}


export_GainPrev_to_txt <- function(scenario_type, IVM_Pregnancy, Prev_mos10, VectIVM_field_dependancy,
                                   VectPropIVM, VectNumber_of_cycle,GainPrev_array_oral, 
                                   GainPrev_array_oral_low, GainPrev_array_oral_up, GainPrev_array_oral_15,
                                   GainPrev_array_oral_15_low, GainPrev_array_oral_15_up,
                                   GainPrev_array_LAIF, GainPrev_array_LAIF_low, GainPrev_array_LAIF_up,
                                   GainPrev_array_LAIF_15, GainPrev_array_LAIF_15_low, GainPrev_array_LAIF_15_up,
                                   oral_cycle_scenario = "initial") {
  
  filename = paste0("GainPrev_", scenario_type, "_", oral_cycle_scenario, "_Pregnant", IVM_Pregnancy,
                    "_Prev_mos", Name_Prev_mos(Prev_mos10), ".csv")
  
  df = data.frame()
  strategy_names_oral = c("BOHEMIA", "RIMDAMAL_S", "RIMDAMAL_K")
  strategy_names_LAIF = c("LAIF_06vk5", "LAIF_10", "LAIF_06kis")
  
  # Nombre de cycles oral selon le scenario
  if (oral_cycle_scenario == "initial") {
    nc_oral = c(3, 4, 4)  # BOHEMIA=3, RIMDAMAL_S=4, RIMDAMAL_K=4
  } else if (oral_cycle_scenario == "all3") {
    nc_oral = c(3, 3, 3)
  } else if (oral_cycle_scenario == "all4") {
    nc_oral = c(4, 4, 4)
  }
  
  for (field_idx in 1:length(VectIVM_field_dependancy)) {
    for (j in 1:length(VectPropIVM)) {
      # Strategies orales
      for (s in 1:3) {
        df = rbind(df, data.frame(
          Strategy = strategy_names_oral[s],
          PropIVM = VectPropIVM[j],
          Number_of_cycle = nc_oral[s],
          IVM_field_dep = VectIVM_field_dependancy[field_idx],
          GainPrev = GainPrev_array_oral[s, j, field_idx],
          GainPrev_low = GainPrev_array_oral_low[s, j, field_idx],
          GainPrev_up = GainPrev_array_oral_up[s, j, field_idx],
          GainPrev_15 = GainPrev_array_oral_15[s, j, field_idx],
          GainPrev_15_low = GainPrev_array_oral_15_low[s, j, field_idx],
          GainPrev_15_up = GainPrev_array_oral_15_up[s, j, field_idx]
        ))
      }
      # Strategies LAIF
      for (s in 1:3) {
        for (nc in VectNumber_of_cycle) {
          df = rbind(df, data.frame(
            Strategy = strategy_names_LAIF[s],
            PropIVM = VectPropIVM[j],
            Number_of_cycle = nc,
            IVM_field_dep = VectIVM_field_dependancy[field_idx],
            GainPrev = GainPrev_array_LAIF[s, j, nc, field_idx],
            GainPrev_low = GainPrev_array_LAIF_low[s, j, nc, field_idx],
            GainPrev_up = GainPrev_array_LAIF_up[s, j, nc, field_idx],
            GainPrev_15 = GainPrev_array_LAIF_15[s, j, nc, field_idx],
            GainPrev_15_low = GainPrev_array_LAIF_15_low[s, j, nc, field_idx],
            GainPrev_15_up = GainPrev_array_LAIF_15_up[s, j, nc, field_idx]
          ))
        }
      }
    }
  }
  df[, 5:10] = round(df[, 5:10], 4)
  write.table(df, file = filename, sep = ",", row.names = FALSE, quote = FALSE)
}

load_simulation_data=function(scenario_type, IVM_field_dependancy, PropIVM, 
                              IVM_Pregnancy, Prev_mos10, CoefDelay) {
  if (scenario_type == "optimal") {
    RDataFileName = paste0("Results_Ihopt_IVM_field_dep_", IVM_field_dependancy, "_PropIVM_",
                           gsub("\\.", "_", PropIVM), "_Pregnant", IVM_Pregnancy, "Prev_mos",
                           Name_Prev_mos(Prev_mos10), "_Timing", gsub("\\.", "_", CoefDelay), ".RData")
  } else {
    RDataFileName = paste0("Results_IhBase_IVM_field_dep_", IVM_field_dependancy, "_PropIVM_",
                           gsub("\\.", "_", PropIVM), "_Pregnant", IVM_Pregnancy, "Prev_mos", Name_Prev_mos(Prev_mos10), ".RData")
  }
  load(RDataFileName)
  return(simulation_results)
}


plot_scenario_dynamics=function(scenario_type = "optimal", IVM_field_dependancy = 0, IVM_Pregnancy = 0,
                                Prev_mos10 = 1, VectPropIVM = c(0.5, 0.7, 0.9), plot_type = "Ih", 
                                CoefDelay = 0.85) {
  
  ColVect = c("#ff3355", "#037153", "#D98E04", "#205072", "#8B008B", "#8B4513"); LineVect = c(1, 3, 2, 4); LC = 3.5; GrNumber = 0
  
  FigName = paste0(ifelse(scenario_type == "optimal", "Scenario_", "Base_Scenario_"), 
                   ifelse(plot_type == "Ih", "Ih", "Im"), "opt_field_dep", IVM_field_dependancy, 
                   "_Pregnant", IVM_Pregnancy, "Prev_mos", Name_Prev_mos(Prev_mos10), ".pdf")
  
  pdf(FigName, width = 15, height = 8.7)
  par(oma = c(6, 1, 2, .1), mar = c(3, 3.5, 2, .1), mfrow = c(length(VectPropIVM), 4))
  
  for (PropIVM in VectPropIVM) {
    simulation_results = load_simulation_data(scenario_type, IVM_field_dependancy, 
                                              PropIVM, IVM_Pregnancy, Prev_mos10, CoefDelay)
    ModelNoIVM = simulation_results$ModelNoIVM
    time = simulation_results$time
    
    for (nc in 1:4) {
      GrNumber = GrNumber + 1; Max_y = ifelse(plot_type == "Ih", 0.103, 0.082)
      
      plot(-1, 1, type = "l", xlab = "", xlim = c(t_begin_Camp, tmax), ylab = "", ylim = c(0, Max_y), yaxt = "n", xaxt = "n")
      axis(1, at = seq(t_begin_Camp, tmax, by = 60), labels = seq(0, tmaxmonths, by = 2), las = 1)
      y_vals = if(plot_type == "Ih") c(0, 0.02, 0.04, 0.06, 0.08, 0.10) else c(0, 0.02, 0.04, 0.06, 0.08)
      axis(2, at = y_vals, labels = y_vals * 100, las = 1)
      
      if (plot_type == "Ih") { lines(time, ModelNoIVM$PropIhTot, lwd = 0.5, lty = 2, col = "black") }
      else { lines(time, ModelNoIVM$PropImEIP, lwd = 0.5, lty = 2, col = "black") }
      
      if (plot_type == "Ih") {
        lines(time, simulation_results$ModelLongLasting06[[nc]]$PropIhTot, lwd = LC, lty = LineVect[1], col = ColVect[1])
        lines(time, simulation_results$ModelLongLasting06_kis[[nc]]$PropIhTot, lwd = LC, lty = LineVect[4], col = ColVect[6])
        lines(time, simulation_results$ModelLongLasting1[[nc]]$PropIhTot, lwd = LC, lty = LineVect[2], col = ColVect[2])
        lines(time, simulation_results$ModelBohemia$PropIhTot, lwd = LC, lty = LineVect[3], col = ColVect[3])
        lines(time, simulation_results$ModelRimdamal$PropIhTot, lwd = LC, lty = LineVect[4], col = ColVect[4])
        lines(time, simulation_results$ModelKamauRimdamal$PropIhTot, lwd = LC, lty = LineVect[4], col = ColVect[5])
      } else {
        lines(time, simulation_results$ModelLongLasting06[[nc]]$PropImEIP, lwd = LC, lty = LineVect[1], col = ColVect[1])
        lines(time, simulation_results$ModelLongLasting06_kis[[nc]]$PropImEIP, lwd = LC, lty = LineVect[4], col = ColVect[6])
        lines(time, simulation_results$ModelLongLasting1[[nc]]$PropImEIP, lwd = LC, lty = LineVect[2], col = ColVect[2])
        lines(time, simulation_results$ModelBohemia$PropImEIP, lwd = LC, lty = LineVect[3], col = ColVect[3])
        lines(time, simulation_results$ModelRimdamal$PropImEIP, lwd = LC, lty = LineVect[4], col = ColVect[4])
        lines(time, simulation_results$ModelKamauRimdamal$PropImEIP, lwd = LC, lty = LineVect[4], col = ColVect[5])
      }
      
      par(xpd = NA); text(t_begin_Camp - 0.05, Max_y * 1.15, paste0("(", LETTERS[GrNumber], ")"), cex = 1.6); par(xpd = FALSE)
      if (GrNumber %in% 1:4) mtext(bquote("#cycle for mdc-STM-001: " ~ n[c] == .(nc)), side = 3, adj = 0.5, cex = 1.05, line = 1)
      if (GrNumber %in% c(1, 5, 9)) mtext(ifelse(plot_type == "Ih", "Clinical cases (all ages; %)", "Infectious mosquitoes (%)"), side = 2, adj = 0.5, cex = .95, line = 2)
      if (GrNumber %in% c(1, 5, 9)) mtext(bquote("Coverage: " ~ .(PropIVM * 100) ~ "%"), side = 2, adj = 0.5, cex = 1.05, line = 3)
      if (GrNumber %in% 9:12) mtext("Time (Months)", side = 1, adj = 0.5, cex = 1.05, line = 3)
    }
  }
  
  LEGEND = c("No IVM", "mdc-STM-001-0.6-vk5", "mdc-STM-001-0.6-kis", "mdc-STM-001-1.0", expression("RIMDAMAL II"[S]), expression("RIMDAMAL II"[K]), "BOHEMIA")
  reset(); legend("bottom", legend = LEGEND, xpd = NA, horiz = FALSE, bty = "n", lty = c(2, 1, 4, 3, 4, 4, 2), lwd = c(0.5, LC, LC, LC, LC, LC, LC),
                  col = c("black", "#ff3355", "#8B4513", "#037153", "#205072", "#8B008B", "#D98E04"), ncol = 4, cex = 1.2)
  dev.off()
}

plot_gainprev_barplot=function(IVM_field_dependancy = 0, IVM_Pregnancy = 0, Prev_mos10 = 1) {
  
  load_file = paste0("ResultatsPrev_Ihopt_field_dep_", IVM_field_dependancy, "_Pregnant", IVM_Pregnancy, "Prev_mos", Name_Prev_mos(Prev_mos10), ".RData")
  load(load_file)
  
  if (length(dim(GainPrev_array_oral)) == 3) {
    GP_oral = GainPrev_array_oral[,,1]; GP_oral_low = GainPrev_array_oral_low[,,1]; 
    GP_oral_up = GainPrev_array_oral_up[,,1]; GP_LAIF = GainPrev_array_LAIF[,,,1]; 
    GP_LAIF_low = GainPrev_array_LAIF_low[,,,1]; GP_LAIF_up = GainPrev_array_LAIF_up[,,,1]
  } else {
    GP_oral = GainPrev_array_oral; GP_oral_low = GainPrev_array_oral_low; 
    GP_oral_up = GainPrev_array_oral_up; GP_LAIF = GainPrev_array_LAIF; 
    GP_LAIF_low = GainPrev_array_LAIF_low; GP_LAIF_up = GainPrev_array_LAIF_up
  }
  
  VectPropIVM_sel = c(0.5, 0.7, 0.9); prop_idx_list = sapply(VectPropIVM_sel, function(v) which.min(abs(VectPropIVM - v)))
  ColVect = c("#D98E04", "#8B008B", "#205072", "#ff3355", "#8B4513", "#037153"); GrNumber = 0
  
  FigName = paste0("Scenarios_GainPrev_Barplot_field_dep", IVM_field_dependancy, "_Pregnant", IVM_Pregnancy, "Prev_mos", Name_Prev_mos(Prev_mos10), ".pdf")
  pdf(FigName, width = 7.5, height = 6.2); par(oma = c(3, 2, 2, 1), mar = c(2, 2, 1.7, 1), mfrow = c(2, 2))
  
  for (i in 1:4) {
    GrNumber = GrNumber + 1
    plot(1, type = "n", xlim = c(0.55, 3.45), ylim = c(0, 0.85), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    axis(2, at = seq(0, 0.8, 0.2), labels = seq(0, 80, 20), las = 1); axis(1, at = 1:3, labels = c("50", "70", "90"))
    
    par(xpd = NA); text(0.3, 0.85 * 1.15, paste0("(", LETTERS[GrNumber], ")"), cex = 1.3); par(xpd = FALSE)
    mtext(bquote("#cycle for mdc-STM-001: " ~ n[c] == .(i)), side = 3, adj = 0.5, cex = 0.95, line = 0.5)
    if (GrNumber %in% c(3, 4)) mtext("Coverage of the target population (%)", side = 1, adj = 0.5, cex = .85, line = 3)
    if (GrNumber %in% c(1, 3)) mtext("Relative gain (all ages, %)", side = 2, adj = 0.5, cex = 0.95, line = 2.2)
    
    bar_width = 0.12; bar_spacing = 0.02
    
    for (pos_idx in 1:3) {
      prop_idx = prop_idx_list[pos_idx]; if(is.na(prop_idx)) next
      base_x = pos_idx
      
      values = c(GP_oral[1, prop_idx], GP_oral[3, prop_idx], GP_oral[2, prop_idx],
                 GP_LAIF[1, prop_idx, i], GP_LAIF[3, prop_idx, i], GP_LAIF[2, prop_idx, i])
      values_low = c(GP_oral_low[1, prop_idx], GP_oral_low[3, prop_idx], GP_oral_low[2, prop_idx],
                     GP_LAIF_low[1, prop_idx, i], GP_LAIF_low[3, prop_idx, i], GP_LAIF_low[2, prop_idx, i])
      values_up = c(GP_oral_up[1, prop_idx], GP_oral_up[3, prop_idx], GP_oral_up[2, prop_idx], 
                    GP_LAIF_up[1, prop_idx, i], GP_LAIF_up[3, prop_idx, i], GP_LAIF_up[2, prop_idx, i])
      
      total_width = 6 * bar_width + 5 * bar_spacing; start_offset = -total_width / 2
      x_pos = sapply(0:5, function(k) base_x + start_offset + k*(bar_width + bar_spacing) + bar_width/2)
      
      for (b in 1:6) {
        if(!is.na(values[b]) && values[b] > 0) {
          # Arrows AVANT barplot (IC)
          segments(x_pos[b] - bar_width/2 + 0.01, values_low[b], x_pos[b] + bar_width/2 - 0.01, values_low[b], col = ColVect[b], lwd = 1.5)
          segments(x_pos[b] - bar_width/2 + 0.01, values_up[b], x_pos[b] + bar_width/2 - 0.01, values_up[b], col = ColVect[b], lwd = 1.5)
          segments(x_pos[b], values_low[b], x_pos[b], values_up[b], col = ColVect[b], lwd = 1)
          # Barplot PAR-DESSUS
          rect(x_pos[b] - bar_width/2, 0, x_pos[b] + bar_width/2, values[b], col = ColVect[b], border = "black", lwd = 0.5)
        }
      }
    }
    abline(h = 0.2, lty = 2, col = "gray", lwd = 1)
    
    if (GrNumber == 1) {
      legend("topleft",
             legend = c("BOHEMIA", expression("RIMDAMAL II"[K]), expression("RIMDAMAL II"[S]),
                        bquote("mdc-STM-001-0.6-vk5"), bquote("mdc-STM-001-0.6-kis"),
                        bquote("mdc-STM-001-1.0")),
             horiz = FALSE, fill = ColVect, bty = "n", cex = 1, border = "black")
    } 
  }
  dev.off()
}

export_to_excel=function(scenario_type = "optimal", IVM_field_dependancy = 0, IVM_Pregnancy = 0,
                         Prev_mos10 = 1, VectPropIVM = c(0.7), CoefDelay = 0.85) {
  wb_laif = createWorkbook(); wb_oral = createWorkbook()
  
  for (PropIVM in VectPropIVM) {
    sim = load_simulation_data(scenario_type, IVM_field_dependancy, PropIVM, IVM_Pregnancy, Prev_mos10, CoefDelay)
    
    # LAIF strategies
    for (nm in c("LongLasting06", "LongLasting1", "LongLasting06_kis")) {
      for (nc in 1:4) {
        sheet_name = paste0(nm, "_Prop", PropIVM*10, "_nc", nc); addWorksheet(wb_laif, sheet_name)
        df = data.frame(time = sim$time, Sm_Tot = sim[[paste0("Model", nm)]][[nc]]$Sm_Tot - sim[[paste0("Model", nm)]][[nc]]$Smg_ivm,
                        Im_Tot = sim[[paste0("Model", nm)]][[nc]]$Im_Tot - sim[[paste0("Model", nm)]][[nc]]$Img_ivm,
                        Smg_ivm = sim[[paste0("Model", nm)]][[nc]]$Smg_ivm, Img_ivm = sim[[paste0("Model", nm)]][[nc]]$Img_ivm)
        for (a in 1:Nah) {
          df[[paste0("Shg_age", a)]] = sim[[paste0("Model", nm)]][[nc]]$Shg[, a]; df[[paste0("Ahg_age", a)]] = sim[[paste0("Model", nm)]][[nc]]$Ahg[, a]
          df[[paste0("Ihg_age", a)]] = sim[[paste0("Model", nm)]][[nc]]$Ihg[, a]; df[[paste0("Rhg_age", a)]] = sim[[paste0("Model", nm)]][[nc]]$Rhg[, a]
          df[[paste0("Shg_ivm_age", a)]] = sim[[paste0("Model", nm)]][[nc]]$Shg_ivm[, a]; df[[paste0("Ahg_ivm_age", a)]] = sim[[paste0("Model", nm)]][[nc]]$Ahg_ivm[, a]
          df[[paste0("Ihg_ivm_age", a)]] = sim[[paste0("Model", nm)]][[nc]]$Ihg_ivm[, a]; df[[paste0("Rhg_ivm_age", a)]] = sim[[paste0("Model", nm)]][[nc]]$Rhg_ivm[, a]
        }
        writeData(wb_laif, sheet_name, df)
      }
    }
    
    # Oral strategies
    for (nm in c("Bohemia", "Rimdamal", "KamauRimdamal")) {
      sheet_name = paste0(nm, "_Prop", PropIVM*10); addWorksheet(wb_oral, sheet_name)
      df = data.frame(time = sim$time, Sm_Tot = sim[[paste0("Model", nm)]]$Sm_Tot - sim[[paste0("Model", nm)]]$Smg_ivm,
                      Im_Tot = sim[[paste0("Model", nm)]]$Im_Tot - sim[[paste0("Model", nm)]]$Img_ivm,
                      Smg_ivm = sim[[paste0("Model", nm)]]$Smg_ivm, Img_ivm = sim[[paste0("Model", nm)]]$Img_ivm)
      for (a in 1:Nah) {
        df[[paste0("Shg_age", a)]] = sim[[paste0("Model", nm)]]$Shg[, a]; df[[paste0("Ahg_age", a)]] = sim[[paste0("Model", nm)]]$Ahg[, a]
        df[[paste0("Ihg_age", a)]] = sim[[paste0("Model", nm)]]$Ihg[, a]; df[[paste0("Rhg_age", a)]] = sim[[paste0("Model", nm)]]$Rhg[, a]
        df[[paste0("Shg_ivm_age", a)]] = sim[[paste0("Model", nm)]]$Shg_ivm[, a]; df[[paste0("Ahg_ivm_age", a)]] = sim[[paste0("Model", nm)]]$Ahg_ivm[, a]
        df[[paste0("Ihg_ivm_age", a)]] = sim[[paste0("Model", nm)]]$Ihg_ivm[, a]; df[[paste0("Rhg_ivm_age", a)]] = sim[[paste0("Model", nm)]]$Rhg_ivm[, a]
      }
      writeData(wb_oral, sheet_name, df)
    }
  }
  
  scenario_label = ifelse(scenario_type == "optimal", "Opt", "Base")
  saveWorkbook(wb_laif, paste0("Data_LAIF_", scenario_label, "_field", IVM_field_dependancy, "_Preg", IVM_Pregnancy, "_Mos", Name_Prev_mos(Prev_mos10), ".xlsx"), overwrite = TRUE)
  saveWorkbook(wb_oral, paste0("Data_Oral_", scenario_label, "_field", IVM_field_dependancy, "_Preg", IVM_Pregnancy, "_Mos", Name_Prev_mos(Prev_mos10), ".xlsx"), overwrite = TRUE)
}

# ============================================================================
# RUNS ET PLOTS
# ============================================================================

# 1) Run et plot scénario optimal, dynamic Figure 2, Barplot avec IC Figure 3 et Comparaison Figure 4
results_opt = run_scenario_IVM(scenario_type = "optimal", IVM_Pregnancy = 0,
                               Prev_mos10 = 1, VectIVM_field_dependancy = c(0))
plot_scenario_dynamics(scenario_type = "optimal", plot_type = "Ih")
plot_scenario_dynamics(scenario_type = "optimal", plot_type = "Im")
export_to_excel(scenario_type = "optimal")
plot_gainprev_barplot()
plot_formulations_comparison()

# 2) Run scenario Base et plot scénario Figure S5  : Base Scenarios_Formulations for Ih and Im
results_base = run_scenario_IVM(scenario_type = "base", IVM_Pregnancy = 0,
                                Prev_mos10 = 1, VectIVM_field_dependancy = c(0))
plot_scenario_dynamics(scenario_type = "base", plot_type = "Ih")
plot_scenario_dynamics(scenario_type = "base", plot_type = "Im")
export_to_excel(scenario_type = "base")

# 3)Comparaison Optimal vs Base Figure S9,Pas besoin de Run pour les data: deja run en 1) et 2)
plot_comparison_generic("opt_vs_base")

# 4)Comparaison Field dependancy 0 vs 1 Figure S14, 
#data pour IVM_field_dependancy=0 deja, run juste en 1
# run pour IVM_field_dependancy=1 ensuite comparaison
results_opt = run_scenario_IVM(scenario_type = "optimal",IVM_Pregnancy = 0, 
                               Prev_mos10 = 1, VectIVM_field_dependancy = c(1))

plot_comparison_generic("field_dep")

# 5) Comparaison IVM_Pregnancy 0 vs 1 Figure S14, 
#data pour IVM_Pregnancy=0 deja,  run juste en 1
# run pour IVM_Pregnancy=1 ensuite comparaison
results_opt = run_scenario_IVM(scenario_type = "optimal",IVM_Pregnancy = 1, 
                               Prev_mos10 = 1, VectIVM_field_dependancy = c(0))

plot_comparison_generic("pregnancy")

# 6) Comparaison Prev_mos10 0 vs 1 , 
#data pour Prev_mos10=1 deja run en 1)
# run pour Prev_mos10=0 ensuite comparaison

results_opt = run_scenario_IVM(scenario_type = "optimal",IVM_Pregnancy = 0, 
                               Prev_mos10 = 0, VectIVM_field_dependancy = c(0))

plot_comparison_generic("prev_mos")

# 7) Comparaison gainPrev all ages vs 5-15ans Figure S9 tous les data deja run

plot_comparison_generic("deltoptim_vs_15")

#Figure S9 :   Scenarios_Formulations Barplots formulations Figure S13
###########################
CoefDelay=0.85 
IVM_field_dependancy=0;IVM_Pregnancy=0; Prev_mos10 = 1
sink(paste("DurEntreCycleBarplotsIVM_field_dep",IVM_field_dependancy,"_Timing",
           gsub("\\.", "_", CoefDelay),".txt",sep=""))

VectDataset=c("No IVM", "PropIVM 0.7")
LegendBarplot=c("Without IVM", "With IVM")
ColBarplot=c("No IVM" = "#148F77", "PropIVM 0.7" = "#F39C12")
VectPropIVM=c(0, 0.7)
{
  data_Ih_a_Bohemia = data.frame()
  data_Ih_a_Rimdamal = data.frame()
  data_Ih_a_KamauRimdamal = data.frame()
  data_Ih_a_LAIF0_6 = data.frame()
  data_Ih_a_LAIF0_6_kis = data.frame()
  data_Ih_a_LAIF1_0 = data.frame()
  
  Nbb=0
  for (PropIVM in VectPropIVM) {
    Nbb= Nbb+1
    
    print(paste0("CoefDelay=",CoefDelay))
    print(paste0("PropIVM=",PropIVM))
    
    #On tourne le modèle pour Bohemia
    strategy =2;NbCycle=3
    ParamIvmFormulation = ModelEtau50Alpha(strategy)
    tau50 = ParamIvmFormulation$tau50; alpha = ParamIvmFormulation$alpha; Dmax = ParamIvmFormulation$Dmax
    
    Number_of_cycle=1;VectTime_between_cycles=1;#pour le premier cycle, cette variable est inutile
    ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                           mu_h,delta_h,bar_beta_h,beta_h,gamma_h,theta,beta_m,Wedge_m,mu_m,
                           ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy, dtau,tau_max,Vtau,Ntau,
                           dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,tau50, alpha, Dmax)
    VectTime_between_cycles=rep(0,NbCycle-1)
    TempsNextCycle0=0
    for (Number_of_cycle in 2:NbCycle) {
      #On récupère le temps pour le cycle suivant
      TempsNextCycle=ModelOutput$DurIvmEffect #temps durée de décoirssance depuis begin camp
      d0=TempsNextCycle-TempsNextCycle0- Dur_cycle
      TempsNextCycle0=TempsNextCycle
      VectTime_between_cycles[Number_of_cycle-1]=CoefDelay*d0
      
      ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                             mu_h,delta_h,bar_beta_h,beta_h,gamma_h,theta,beta_m,Wedge_m,mu_m,
                             ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy, dtau,tau_max,Vtau,Ntau,
                             dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,tau50, alpha, Dmax)
    }
    #on récupère les valeurs ici
    print("Bohemia")
    print(VectTime_between_cycles)
    ModelBohemia=ModelOutput
    
    data_Ih_a_Bohemia = rbind(data_Ih_a_Bohemia, 
                              data.frame(Age = ModelOutput$cut_age[-length(ModelOutput$cut_age)], # Bornes inférieures
                                         mIh_a = ModelBohemia$mIh_a,  fIh_a = ModelBohemia$fIh_a, 
                                         Dataset = VectDataset[Nbb])) 
    
    
    
    #On tourne le modèle pour Rimdamal
    strategy =3;NbCycle=4
    ParamIvmFormulation = ModelEtau50Alpha(strategy)
    tau50 = ParamIvmFormulation$tau50; alpha = ParamIvmFormulation$alpha; Dmax = ParamIvmFormulation$Dmax
    
    Number_of_cycle=1;VectTime_between_cycles=1;#pour le premier cycle, cette variable est inutile
    ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                           mu_h,delta_h,bar_beta_h,beta_h,gamma_h,theta,beta_m,Wedge_m,mu_m,
                           ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy, dtau,tau_max,Vtau,Ntau,
                           dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,tau50, alpha, Dmax)
    VectTime_between_cycles=rep(0,NbCycle-1)
    TempsNextCycle0=0
    for (Number_of_cycle in 2:NbCycle) {
      #On récupère le temps pour le cycle suivant
      TempsNextCycle=ModelOutput$DurIvmEffect #temps durée de décoirssance depuis begin camp
      d0=TempsNextCycle-TempsNextCycle0- Dur_cycle
      TempsNextCycle0=TempsNextCycle
      VectTime_between_cycles[Number_of_cycle-1]=CoefDelay*d0
      
      ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                             mu_h,delta_h,bar_beta_h,beta_h,gamma_h,theta,beta_m,Wedge_m,mu_m,
                             ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy, dtau,tau_max,Vtau,Ntau,
                             dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,tau50, alpha, Dmax)
    }
    #on récupère les valeurs ici
    print("Rimdamal")
    print(VectTime_between_cycles)
    #VectTime_between_cycles=0.7*c(39.50, 79.25-39.50, 79.25-(79.25-39.50))
    ModelRimdamal=ModelOutput
    
    data_Ih_a_Rimdamal = rbind(data_Ih_a_Rimdamal, 
                               data.frame(Age = ModelOutput$cut_age[-length(ModelOutput$cut_age)], # Bornes inférieures
                                          mIh_a = ModelRimdamal$mIh_a,  fIh_a = ModelRimdamal$fIh_a, 
                                          Dataset = VectDataset[Nbb]))
    
    #On tourne le modèle pour KamauRimdamal
    strategy =4;NbCycle=4
    ParamIvmFormulation = ModelEtau50Alpha(strategy)
    tau50 = ParamIvmFormulation$tau50; alpha = ParamIvmFormulation$alpha; Dmax = ParamIvmFormulation$Dmax
    
    Number_of_cycle=1;VectTime_between_cycles=1;#pour le premier cycle, cette variable est inutile
    ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                           mu_h,delta_h,bar_beta_h,beta_h,gamma_h,theta,beta_m,Wedge_m,mu_m,
                           ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy, dtau,tau_max,Vtau,Ntau,
                           dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,tau50, alpha, Dmax)
    VectTime_between_cycles=rep(0,NbCycle-1)
    TempsNextCycle0=0
    for (Number_of_cycle in 2:NbCycle) {
      #On récupère le temps pour le cycle suivant
      TempsNextCycle=ModelOutput$DurIvmEffect #temps durée de décoirssance depuis begin camp
      d0=TempsNextCycle-TempsNextCycle0- Dur_cycle
      TempsNextCycle0=TempsNextCycle
      VectTime_between_cycles[Number_of_cycle-1]=CoefDelay*d0
      
      ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                             mu_h,delta_h,bar_beta_h,beta_h,gamma_h,theta,beta_m,Wedge_m,mu_m,
                             ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy, dtau,tau_max,Vtau,Ntau,
                             dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,tau50, alpha, Dmax)
    }
    #on récupère les valeurs ici
    print("KamauRimdamal")
    print(VectTime_between_cycles)
    #VectTime_between_cycles=0.7*c(39.50, 79.25-39.50, 79.25-(79.25-39.50))
    ModelKamauRimdamal=ModelOutput
    
    data_Ih_a_KamauRimdamal = rbind(data_Ih_a_KamauRimdamal, 
                                    data.frame(Age = ModelOutput$cut_age[-length(ModelOutput$cut_age)], # Bornes inférieures
                                               mIh_a = ModelKamauRimdamal$mIh_a,  fIh_a = ModelKamauRimdamal$fIh_a, 
                                               Dataset = VectDataset[Nbb]))
    
    
    
    #On tourne le modèle pour LongLasting 06
    strategy =0;NbCycle=4
    ParamIvmFormulation = ModelEtau50Alpha(strategy)
    tau50 = ParamIvmFormulation$tau50; alpha = ParamIvmFormulation$alpha; Dmax = ParamIvmFormulation$Dmax
    
    ModelLongLasting06=list()
    Number_of_cycle=1;VectTime_between_cycles=1;#pour le premier cycle, cette variable est inutile
    ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                           mu_h,delta_h,bar_beta_h,beta_h,gamma_h,theta,beta_m,Wedge_m,mu_m,
                           ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy, dtau,tau_max,Vtau,Ntau,
                           dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,tau50, alpha, Dmax)
    ModelLongLasting06[[Number_of_cycle]]=ModelOutput
    
    data_Ih_a_LAIF0_6 = rbind(data_Ih_a_LAIF0_6, 
                              data.frame(Age = ModelOutput$cut_age[-length(ModelOutput$cut_age)], # Bornes inférieures
                                         mIh_a = ModelLongLasting06[[Number_of_cycle]]$mIh_a,
                                         fIh_a = ModelLongLasting06[[Number_of_cycle]]$fIh_a, 
                                         Dataset = VectDataset[Nbb], Ncycle=Number_of_cycle))
    
    VectTime_between_cycles=rep(0,NbCycle-1)
    TempsNextCycle0=0
    for (Number_of_cycle in 2:NbCycle) {
      #On récupère le temps pour le cycle suivant
      TempsNextCycle=ModelOutput$DurIvmEffect #temps durée de décoirssance depuis begin camp
      d0=TempsNextCycle-TempsNextCycle0- Dur_cycle
      TempsNextCycle0=TempsNextCycle
      VectTime_between_cycles[Number_of_cycle-1]=CoefDelay*d0
      
      ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                             mu_h,delta_h,bar_beta_h,beta_h,gamma_h,theta,beta_m,Wedge_m,mu_m,
                             ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy, dtau,tau_max,Vtau,Ntau,
                             dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,tau50, alpha, Dmax)
      ModelLongLasting06[[Number_of_cycle]]=ModelOutput
      
      data_Ih_a_LAIF0_6 = rbind(data_Ih_a_LAIF0_6, 
                                data.frame(Age = ModelOutput$cut_age[-length(ModelOutput$cut_age)], # Bornes inférieures
                                           mIh_a = ModelLongLasting06[[Number_of_cycle]]$mIh_a,
                                           fIh_a = ModelLongLasting06[[Number_of_cycle]]$fIh_a, 
                                           Dataset = VectDataset[Nbb], Ncycle=Number_of_cycle))
      
    }
    #on récupère les valeurs ici
    print("LongLasting0.6")
    print(VectTime_between_cycles)
    
    #On tourne le modèle pour LongLasting 06 Kis
    strategy =5;NbCycle=4
    ParamIvmFormulation = ModelEtau50Alpha(strategy)
    tau50 = ParamIvmFormulation$tau50; alpha = ParamIvmFormulation$alpha; Dmax = ParamIvmFormulation$Dmax
    
    ModelLongLasting06_kis=list()
    Number_of_cycle=1;VectTime_between_cycles=1;#pour le premier cycle, cette variable est inutile
    ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                           mu_h,delta_h,bar_beta_h,beta_h,gamma_h,theta,beta_m,Wedge_m,mu_m,
                           ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy, dtau,tau_max,Vtau,Ntau,
                           dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,tau50, alpha, Dmax)
    ModelLongLasting06_kis[[Number_of_cycle]]=ModelOutput
    
    data_Ih_a_LAIF0_6_kis = rbind(data_Ih_a_LAIF0_6_kis, 
                                  data.frame(Age = ModelOutput$cut_age[-length(ModelOutput$cut_age)], # Bornes inférieures
                                             mIh_a = ModelLongLasting06_kis[[Number_of_cycle]]$mIh_a,
                                             fIh_a = ModelLongLasting06_kis[[Number_of_cycle]]$fIh_a, 
                                             Dataset = VectDataset[Nbb], Ncycle=Number_of_cycle))
    
    VectTime_between_cycles=rep(0,NbCycle-1)
    TempsNextCycle0=0
    for (Number_of_cycle in 2:NbCycle) {
      #On récupère le temps pour le cycle suivant
      TempsNextCycle=ModelOutput$DurIvmEffect #temps durée de décoirssance depuis begin camp
      d0=TempsNextCycle-TempsNextCycle0- Dur_cycle
      TempsNextCycle0=TempsNextCycle
      VectTime_between_cycles[Number_of_cycle-1]=CoefDelay*d0
      
      ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                             mu_h,delta_h,bar_beta_h,beta_h,gamma_h,theta,beta_m,Wedge_m,mu_m,
                             ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy, dtau,tau_max,Vtau,Ntau,
                             dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,tau50, alpha, Dmax)
      ModelLongLasting06_kis[[Number_of_cycle]]=ModelOutput
      
      data_Ih_a_LAIF0_6_kis = rbind(data_Ih_a_LAIF0_6_kis, 
                                    data.frame(Age = ModelOutput$cut_age[-length(ModelOutput$cut_age)], # Bornes inférieures
                                               mIh_a = ModelLongLasting06_kis[[Number_of_cycle]]$mIh_a,
                                               fIh_a = ModelLongLasting06_kis[[Number_of_cycle]]$fIh_a, 
                                               Dataset = VectDataset[Nbb], Ncycle=Number_of_cycle))
      
    }
    #on récupère les valeurs ici
    print("LongLasting0.6 Kis")
    print(VectTime_between_cycles)
    
    #On tourne le modèle pour LongLasting 1.0
    strategy =1;NbCycle=4
    ParamIvmFormulation = ModelEtau50Alpha(strategy)
    tau50 = ParamIvmFormulation$tau50; alpha = ParamIvmFormulation$alpha; Dmax = ParamIvmFormulation$Dmax
    
    ModelLongLasting1=list()
    Number_of_cycle=1;VectTime_between_cycles=1;#pour le premier cycle, cette variable est inutile
    ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                           mu_h,delta_h,bar_beta_h,beta_h,gamma_h,theta,beta_m,Wedge_m,mu_m,
                           ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy, dtau,tau_max,Vtau,Ntau,
                           dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,tau50, alpha, Dmax)
    ModelLongLasting1[[Number_of_cycle]]=ModelOutput
    
    data_Ih_a_LAIF1_0 = rbind(data_Ih_a_LAIF1_0, 
                              data.frame(Age = ModelOutput$cut_age[-length(ModelOutput$cut_age)], # Bornes inférieures
                                         mIh_a = ModelLongLasting1[[Number_of_cycle]]$mIh_a,
                                         fIh_a = ModelLongLasting1[[Number_of_cycle]]$fIh_a, 
                                         Dataset = VectDataset[Nbb], Ncycle=Number_of_cycle))
    
    VectTime_between_cycles=rep(0,NbCycle-1)
    TempsNextCycle0=0
    for (Number_of_cycle in 2:NbCycle) {
      #On récupère le temps pour le cycle suivant
      TempsNextCycle=ModelOutput$DurIvmEffect #temps durée de décoirssance depuis begin camp
      d0=TempsNextCycle-TempsNextCycle0- Dur_cycle
      TempsNextCycle0=TempsNextCycle
      VectTime_between_cycles[Number_of_cycle-1]=CoefDelay*d0
      
      ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                             mu_h,delta_h,bar_beta_h,beta_h,gamma_h,theta,beta_m,Wedge_m,mu_m,
                             ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy, dtau,tau_max,Vtau,Ntau,
                             dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,tau50, alpha, Dmax)
      ModelLongLasting1[[Number_of_cycle]]=ModelOutput
      
      data_Ih_a_LAIF1_0 = rbind(data_Ih_a_LAIF1_0, 
                                data.frame(Age = ModelOutput$cut_age[-length(ModelOutput$cut_age)], # Bornes inférieures
                                           mIh_a = ModelLongLasting1[[Number_of_cycle]]$mIh_a,
                                           fIh_a = ModelLongLasting1[[Number_of_cycle]]$fIh_a, 
                                           Dataset = VectDataset[Nbb], Ncycle=Number_of_cycle))
      
    }
    #on récupère les valeurs ici
    print("LongLasting1")
    print(VectTime_between_cycles)
    
  }
  
  save(data_Ih_a_Bohemia, data_Ih_a_Rimdamal, data_Ih_a_LAIF0_6,data_Ih_a_LAIF0_6_kis, 
       data_Ih_a_LAIF1_0, data_Ih_a_KamauRimdamal,
       VectDataset,  LegendBarplot,  ColBarplot,  VectPropIVM, CoefDelay,  IVM_field_dependancy,
       ah_max, file = paste0("ResultsBarplots","_Timing",gsub("\\.", "_", CoefDelay),"_field_dep",
                             IVM_field_dependancy ,"_Pregnant", IVM_Pregnancy,"Prev_mos",
                             Name_Prev_mos(Prev_mos10), ".RData"))
  
  Timing=CoefDelay
  load(paste0("ResultsBarplots_Timing",gsub("\\.", "_", Timing),"_field_dep",
              IVM_field_dependancy ,"_Pregnant", IVM_Pregnancy,"Prev_mos",
              Name_Prev_mos(Prev_mos10),".RData"))
  ls()
  
  Plot_mIh_LAIF06=list()
  Plot_fIh_LAIF06=list()
  Plot_mIh_LAIF06_kis=list()
  Plot_fIh_LAIF06_kis=list()
  Plot_mIh_LAIF1_0=list()
  Plot_fIh_LAIF1_0=list()
  #Oral and LAIF
  { 
    GrNumber=0
    #plot mIh_a
    GrNumber=GrNumber+1
    #Barplots Bohemia
    {
      maxPop=4200
      PlotmIh_a_Bohemia= ggplot(data = data_Ih_a_Bohemia, aes(x = Age + 2.5, y = mIh_a, fill = Dataset)) +
        geom_bar(
          stat = "identity", 
          color = "black", 
          position = position_dodge(width = 4), # Décalage ajusté
          width = 4                             # Largeur des barres ajustée
        ) +
        annotate(
          "rect", xmin = 5, xmax = 90, ymin = 0, ymax = Inf, fill =   "#76448A",  alpha = 0.15) +
        theme_light() +
        scale_x_continuous(breaks = seq(0, ah_max, by = 5), minor_breaks = NULL, expand = c(0.02, 0),
                           limits = c(0, 90))+
        scale_y_continuous(breaks = seq(0, maxPop, by = 1000), 
                           labels = c("0", expression(10^3), expression(2%*%10^3), expression(3%*%10^3),
                                      expression(4%*%10^3)),limits = c(0,maxPop)) +
        labs(
          x = "",
          y = "Male clinical cases",
          title = "",
          fill = ""
        ) + 
        ggtitle("BOHEMIA") + 
        labs(subtitle = paste0("(", LETTERS[GrNumber], ")")) +
        scale_fill_manual(
          values = ColBarplot,
          labels = LegendBarplot 
        ) +
        theme(
          legend.position = "none", 
          legend.title = element_text(size = 12, face = "bold"), 
          legend.text = element_text(size = 10), 
          legend.key.size = unit(0.3, "cm"), #dimension des carrés en legende
          plot.title = element_text(hjust = 0.5, size = 12), 
          plot.subtitle =element_text(hjust = 0, size = 14) , 
          panel.grid.major.x = element_line(color = "gray", linewidth = 0.2)
        ) +
        guides(fill = guide_legend(override.aes = list(alpha = 0.8)))
    }
    GrNumber=GrNumber+1
    #Barplots Rimdamal
    {
      maxPop=4200
      PlotmIh_a_Rimdamal=ggplot(data = data_Ih_a_Rimdamal, aes(x = Age + 2.5, y = mIh_a, fill = Dataset)) +
        geom_bar(
          stat = "identity", 
          color = "black", 
          position = position_dodge(width = 4), # Décalage ajusté
          width = 4                             # Largeur des barres ajustée
        ) +
        annotate(
          "rect", xmin = 5, xmax = 90, ymin = 0, ymax = Inf, fill =   "#76448A",  alpha = 0.15) +
        theme_light() +
        scale_x_continuous(breaks = seq(0, ah_max, by = 5), minor_breaks = NULL, expand = c(0.02, 0),
                           limits = c(0, 90))+
        scale_y_continuous(breaks = seq(0, maxPop, by = 1000), 
                           labels = c("0", expression(10^3), expression(2%*%10^3), expression(3%*%10^3),
                                      expression(4%*%10^3)),limits = c(0,maxPop)) +
        labs(
          x = "",
          y = "",
          title = "",
          fill = ""
        ) +
        ggtitle(expression("RIMDAMAL II"[S])) + 
        labs(subtitle = paste0("(", LETTERS[GrNumber], ")")) +
        scale_fill_manual(
          values = ColBarplot,
          labels = LegendBarplot 
        ) +
        theme(
          legend.position = "none", 
          legend.title = element_text(size = 12, face = "bold"), 
          legend.text = element_text(size = 10), 
          legend.key.size = unit(0.3, "cm"), #dimension des carrés en legende
          plot.title = element_text(hjust = 0.5, size = 12), 
          plot.subtitle =element_text(hjust = 0, size = 14) , 
          panel.grid.major.x = element_line(color = "gray", linewidth = 0.2)
        ) +
        guides(fill = guide_legend(override.aes = list(alpha = 0.8)))
      
    }
    
    GrNumber=GrNumber+2
    #plot fIh_a
    #Barplots Bohemia
    {
      maxPop=4200
      PlotfIh_a_Bohemia= ggplot(data = data_Ih_a_Bohemia, aes(x = Age + 2.5, y = fIh_a, fill = Dataset)) +
        geom_bar(
          stat = "identity", 
          color = "black", 
          position = position_dodge(width = 4), # Décalage ajusté
          width = 4                             # Largeur des barres ajustée
        ) +
        annotate(
          "rect", xmin = 5, xmax = 15, ymin = 0, ymax = Inf, fill =  "#76448A",  alpha = 0.15) +
        annotate(
          "rect", xmin = 45, xmax = 90, ymin = 0, ymax = Inf, fill = "#76448A",  alpha = 0.15) +
        theme_light() +
        scale_x_continuous(breaks = seq(0, ah_max, by = 5), minor_breaks = NULL, expand = c(0.02, 0),
                           limits = c(0, 90))+
        scale_y_continuous(breaks = seq(0, maxPop, by = 1000), 
                           labels = c("0", expression(10^3), expression(2%*%10^3), expression(3%*%10^3),
                                      expression(4%*%10^3)),limits = c(0,maxPop)) +
        labs(
          x = "Age (Years)",
          y = "Female clinical cases",
          title = "",
          fill = ""
        ) +
        labs(subtitle =paste0("(", LETTERS[GrNumber], ")") ) +
        scale_fill_manual(
          values = ColBarplot,
          labels = LegendBarplot 
        ) +
        theme(
          legend.position = "none", 
          legend.title = element_text(size = 12, face = "bold"), 
          legend.text = element_text(size = 10), 
          legend.key.size = unit(0.3, "cm"), #dimension des carrés en legende
          plot.title = element_text(hjust = 0.5, size = 16), 
          plot.subtitle = element_text(hjust = 0, size = 14), 
          panel.grid.major.x = element_line(color = "gray", linewidth = 0.2)
        ) +
        guides(fill = guide_legend(override.aes = list(alpha = 0.8)))
    }
    GrNumber=GrNumber+1
    #Barplots Rimdamal
    {
      maxPop=4200
      PlotfIh_a_Rimdamal= ggplot(data = data_Ih_a_Rimdamal, aes(x = Age + 2.5, y = fIh_a, fill = Dataset)) +
        geom_bar(
          stat = "identity", 
          color = "black", 
          position = position_dodge(width = 4), # Décalage ajusté
          width = 4                             # Largeur des barres ajustée
        ) +
        annotate(
          "rect", xmin = 5, xmax = 15, ymin = 0, ymax = Inf, fill =  "#76448A",  alpha = 0.15) +
        annotate(
          "rect", xmin = 45, xmax = 90, ymin = 0, ymax = Inf, fill = "#76448A",  alpha = 0.15) +
        theme_light() +
        scale_x_continuous(breaks = seq(0, ah_max, by = 5), minor_breaks = NULL, expand = c(0.02, 0),
                           limits = c(0, 90))+
        scale_y_continuous(breaks = seq(0, maxPop, by = 1000), 
                           labels = c("0", expression(10^3), expression(2%*%10^3), expression(3%*%10^3),
                                      expression(4%*%10^3)),limits = c(0,maxPop)) +
        labs(
          x = "Age (Years)",
          y = "",
          title = "",
          fill = ""
        ) +
        labs(subtitle =paste0("(", LETTERS[GrNumber], ")") ) +
        scale_fill_manual(
          values = ColBarplot,
          labels = LegendBarplot 
        ) +
        theme(
          legend.position = "none", 
          legend.title = element_text(size = 12, face = "bold"), 
          legend.text = element_text(size = 10), 
          legend.key.size = unit(0.3, "cm"), #dimension des carrés en legende
          plot.title = element_text(hjust = 0.5, size = 16), 
          plot.subtitle = element_text(hjust = 0, size = 14), 
          panel.grid.major.x = element_line(color = "gray", linewidth = 0.2)
        ) +
        guides(fill = guide_legend(override.aes = list(alpha = 0.8)))
    }
    
    # Barplot LAIF0.6 Nc=3
    GrNumber = 2
    for (Number_of_cycle in 3:3) {
      GrNumber = GrNumber + 1
      data_Ih_a_LAIF0_6_filt = data_Ih_a_LAIF0_6 %>% filter(Ncycle == Number_of_cycle)
      
      # plot mIh_a
      {
        maxPop = 4200
        Plot_mIh_LAIF06[[Number_of_cycle]] = ggplot(data = data_Ih_a_LAIF0_6_filt, aes(x = Age + 2.5, y = mIh_a, fill = Dataset)) +
          geom_bar(
            stat = "identity",
            color = "black",
            position = position_dodge(width = 4),
            width = 4
          ) +
          annotate("rect", xmin = 5, xmax = 90, ymin = 0, ymax = Inf, fill = "#76448A", alpha = 0.15) +
          theme_light() +
          scale_x_continuous(
            breaks = seq(0, ah_max, by = 5), minor_breaks = NULL, expand = c(0.02, 0),
            limits = c(0, 90)
          ) +
          scale_y_continuous(breaks = seq(0, maxPop, by = 1000), 
                             labels = c("0", expression(10^3), expression(2%*%10^3), expression(3%*%10^3),
                                        expression(4%*%10^3)),limits = c(0,maxPop)) +
          labs(
            x = "",
            y = ifelse(GrNumber == 1, "Male clinical cases", ""),
            title = "",
            fill = ""
          ) +
          ggtitle(bquote("#cycle for mdc-STM-001-0.6-vk5: " ~ n[c] == .(Number_of_cycle))) +
          labs(subtitle = paste0("(", LETTERS[GrNumber], ")")) +
          scale_fill_manual(
            values = ColBarplot,
            labels = LegendBarplot
          ) +
          theme(
            #  axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = "none",  
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10),
            legend.key.size = unit(0.3, "cm"),
            plot.title = element_text(hjust = 0.5, size = 12),
            plot.subtitle = element_text(hjust = 0, size = 14),
            panel.grid.major.x = element_line(color = "gray", linewidth = 0.2)
          ) +
          guides(fill = guide_legend(override.aes = list(alpha = 0.8)))
      }
    }
    GrNumber = 5
    for (Number_of_cycle in 3:3) {
      GrNumber = GrNumber + 1
      data_Ih_a_LAIF0_6_filt = data_Ih_a_LAIF0_6 %>% filter(Ncycle == Number_of_cycle)
      
      # plot fIh_a
      {maxPop = 4200
        Plot_fIh_LAIF06[[Number_of_cycle]] = ggplot(data = data_Ih_a_LAIF0_6_filt, aes(x = Age + 2.5, y = fIh_a, fill = Dataset)) +
          geom_bar(
            stat = "identity",
            color = "black",
            position = position_dodge(width = 4),
            width = 4
          ) +
          annotate("rect", xmin = 5, xmax = 15, ymin = 0, ymax = Inf, fill = "#76448A", alpha = 0.15) +
          annotate("rect", xmin = 45, xmax = 90, ymin = 0, ymax = Inf, fill = "#76448A", alpha = 0.15) +
          theme_light() +
          scale_x_continuous(
            breaks = seq(0, ah_max, by = 5), minor_breaks = NULL, expand = c(0.02, 0),
            limits = c(0, 90)
          ) +
          scale_y_continuous(breaks = seq(0, maxPop, by = 1000), 
                             labels = c("0", expression(10^3), expression(2%*%10^3), expression(3%*%10^3),
                                        expression(4%*%10^3)),limits = c(0,maxPop)) +
          labs(
            x = "Age (Years)",
            y = ifelse(GrNumber == 4, "Female clinical cases", ""),
            title = "",
            fill = ""
          ) +
          labs(subtitle = paste0("(", LETTERS[GrNumber], ")")) +
          scale_fill_manual(
            values = ColBarplot,
            labels = LegendBarplot
          ) +
          theme(
            #  axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = "none",  
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10),
            legend.key.size = unit(0.3, "cm"),
            plot.title = element_text(hjust = 0.5, size = 12),
            plot.subtitle = element_text(hjust = 0, size = 14),
            panel.grid.major.x = element_line(color = "gray", linewidth = 0.2)
          ) +
          guides(fill = guide_legend(override.aes = list(alpha = 0.8)))
      }
    }
    
    plot_list_LAIForal = list( PlotmIh_a_Bohemia, PlotmIh_a_Rimdamal, Plot_mIh_LAIF06[[3]],  
                               PlotfIh_a_Bohemia, PlotfIh_a_Rimdamal, Plot_fIh_LAIF06[[3]])
    
    combined_plot_LAIForal = wrap_plots(plot_list_LAIForal, ncol = 3, nrow = 2) +
      plot_layout(guides = "collect") +
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10)
      ) +
      plot_annotation(
        theme = theme(legend.position = "bottom"))
    
    ggsave(
      filename = paste0("Barplot_Ih_LAIF_oral","_field_dep",IVM_field_dependancy ,"_Pregnant",
                        IVM_Pregnancy,"Prev_mos", Name_Prev_mos(Prev_mos10), ".pdf"),
      plot = combined_plot_LAIForal, width = 13, height = 6.5, units = "in", dpi = 300)
  }
  
}

######################FIGURE S2 prop death mosquitoes, mum_IVM, RH pour LOngLasting, Bohemia, Rimdal
{
  strategy=0#Long lasting 0.6
  ParamLongLasting06=ModelEtau50Alpha(strategy)
  strategy=1#Long lasting 1
  ParamLongLasting1=ModelEtau50Alpha(strategy)
  strategy=2#Bohemia
  ParamBohemia=ModelEtau50Alpha(strategy)
  strategy=3#Rimdamal
  ParamRimdamal=ModelEtau50Alpha(strategy)
  strategy=4#KamauRIMDAMAL
  ParamKamauRIMDAMAL=ModelEtau50Alpha(strategy)
  strategy=5#Long lasting 0.6 train Kis
  ParamLongLasting06Kis=ModelEtau50Alpha(strategy)
  
  #Figure parameters LOng Lasting (propd death mosquitoes, mum_ivm, HR)
  # on combine 3 formulations LongLasting (0, 1, 5) 
  combined_dataLongLasting = bind_rows(
    ParamLongLasting06$predprop %>% mutate(Formulation = "F31-0.6"),
    ParamLongLasting1$predprop %>% mutate(Formulation = "Lam-1.0"), 
    ParamLongLasting06Kis$predprop %>% mutate(Formulation = "F31-0.6-Kis") )
  
  # Correction pour combined_PkLongLasting - standardisation des colonnes
  combined_PkLongLasting = bind_rows( ParamLongLasting06$sub_PK %>% 
                                        mutate(Formulation = "F31-0.6", N_bovin = as.integer(as.character(N_bovin))),
                                      ParamLongLasting1$sub_PK %>% 
                                        mutate(Formulation = "Lam-1.0", N_bovin = if("N_bovin" %in% names(.)) as.integer(as.character(N_bovin)) 
                                               else as.integer(as.character(Bovin))),
                                      ParamLongLasting06Kis$sub_PK %>% 
                                        mutate(Formulation = "F31-0.6-Kis", N_bovin = as.integer(as.character(N_bovin)))  )
  
  combined_PkLongLasting$DAI[is.na(combined_PkLongLasting$DAI)] = 
    combined_PkLongLasting$DAI2[is.na(combined_PkLongLasting$DAI)]
  
  combined_dataOral = bind_rows(  ParamBohemia$predprop %>% mutate(Formulation = "BOHEMIA"),
                                  ParamRimdamal$predprop %>% mutate(Formulation = "RIMDAMAL"),
                                  ParamKamauRIMDAMAL$predprop %>% mutate(Formulation = "KamauRIMDAMAL") )
  combined_PkOral = bind_rows( ParamBohemia$sub_PK %>% mutate(Formulation = "BOHEMIA"),
                               ParamRimdamal$sub_PK %>% mutate(Formulation = "RIMDAMAL"),
                               ParamKamauRIMDAMAL$sub_PK %>% mutate(Formulation = "KamauRIMDAMAL") )
  
  if("DAI2" %in% names(combined_PkOral)) {
    combined_PkOral$DAI[is.na(combined_PkOral$DAI)] =    combined_PkOral$DAI2[is.na(combined_PkOral$DAI)]
  }
  
  #plot Proportion of dead mosquitoes: LongLasting et Orale
  {
    Plot_Prop_dead_mosquitoesLongLasting = ggplot() +
      geom_ribbon( data = combined_dataLongLasting %>%
                     filter(!is.na(Lower) & !is.na(Upper) & Lower >= 0 & Upper <= 1),
                   aes(x = DAI, ymin = pmax(Lower, 0), ymax = pmin(Upper, 1), fill = Formulation),
                   alpha = 0.3 ) +
      geom_line( data = combined_dataLongLasting %>% filter(!is.na(Prediction)),
                 aes(x = DAI, y = Prediction, colour = Formulation), linewidth = 1.2 ) +
      geom_point( data = combined_PkLongLasting %>% filter(!is.na(prop_dead)),
                  aes(x = DAI, y = prop_dead, colour = Formulation, shape = Formulation),
                  size = 2.2, stroke = 0.8 ) +
      scale_y_continuous(expand = c(0.02, 0), limits = c(0, 1.05)) +
      scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150),
                         labels = c("0", "30", "60", "90", "120", "150"),
                         expand = c(0.02, 0), minor_breaks = NULL ) +
      scale_color_manual( values = c("F31-0.6" = "#ff3355", "Lam-1.0" = "#037153", "F31-0.6-Kis" = "#8B4513"),
                          labels = c("F31-0.6" = "mdc-STM-001-0.6-vk5", "Lam-1.0" = "mdc-STM-001-1.0", "F31-0.6-Kis" = "mdc-STM-001-0.6-Kis") ) +
      scale_fill_manual(  values = c("F31-0.6" = "#ff3355", "Lam-1.0" = "#037153", "F31-0.6-Kis" = "#8B4513"),
                          labels = c("F31-0.6" = "mdc-STM-001-0.6-vk5", "Lam-1.0" = "mdc-STM-001-1.0", "F31-0.6-Kis" = "mdc-STM-001-0.6-Kis")) +
      scale_shape_manual(  values = c("F31-0.6" = 16, "Lam-1.0" = 18, "F31-0.6-Kis" = 17),
                           labels = c("F31-0.6" = "mdc-STM-001-0.6-vk5", "Lam-1.0" = "mdc-STM-001-1.0", "F31-0.6-Kis" = "mdc-STM-001-0.6-Kis") ) +
      guides( fill = guide_legend(order = 1), #, override.aes = list(alpha = 0.3)
              colour = guide_legend(order = 1), shape = guide_legend(order = 1)  ) +
      theme_light() +
      labs( x = "Time post exposure to IVM (days)",
            y = "Proportion of dead mosquitoes"  ) +
      theme( legend.position = "bottom",  legend.title = element_blank(),
             legend.text = element_text(size = 8),  legend.key.width = unit(1.5, "cm"),
             plot.title = element_text(hjust = 0, size = 16), axis.text = element_text(size = 12),
             axis.title = element_text(size = 11) ) +
      ggtitle("(A)")
    
    
    Plot_Prop_dead_mosquitoesOral = ggplot() +
      geom_ribbon(data = combined_dataOral %>%
                    filter(!is.na(Lower) & !is.na(Upper) & Lower >= 0 & Upper <= 1),
                  aes(x = DAI, ymin = pmax(Lower, 0), ymax = pmin(Upper, 1), fill = Formulation),
                  alpha = 0.3) +
      geom_line( data = combined_dataOral %>% filter(!is.na(Prediction)),
                 aes(x = DAI, y = Prediction, colour = Formulation),linewidth = 1.2) +
      geom_point(data = combined_PkOral %>% filter(!is.na(prop_dead)),
                 aes(x = DAI, y = prop_dead, colour = Formulation, shape = Formulation),
                 size = 2.2, stroke = 0.8) +
      scale_y_continuous(expand = c(0.02, 0), limits = c(0, 1.05)) +
      scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150),
                         labels = c("0", "30", "60", "90", "120", "150"),
                         expand = c(0.02, 0), minor_breaks = NULL ) +
      scale_color_manual( values = c("BOHEMIA" = "#D98E04", "RIMDAMAL" = "#205072", "KamauRIMDAMAL" = "#8B008B"),
                          labels = c("BOHEMIA" = "BOHEMIA", "KamauRIMDAMAL" = expression("RIMDAMAL II" [K]), "RIMDAMAL" = expression("RIMDAMAL II" [S])) ) +
      scale_fill_manual( values = c("BOHEMIA" = "#D98E04", "RIMDAMAL" = "#205072", "KamauRIMDAMAL" = "#8B008B"),
                         labels = c("BOHEMIA" = "BOHEMIA", "KamauRIMDAMAL" = expression("RIMDAMAL II" [K]), "RIMDAMAL" = expression("RIMDAMAL II" [S])) ) +
      scale_shape_manual(values = c("BOHEMIA" = 16, "RIMDAMAL" = 18, "KamauRIMDAMAL" = 17),
                         labels = c("BOHEMIA" = "BOHEMIA", "KamauRIMDAMAL" = expression("RIMDAMAL II" [K]), "RIMDAMAL" = expression("RIMDAMAL II" [S]))) +
      guides(fill = guide_legend(order = 1), #, override.aes = list(alpha = 0.3, linewidth = 2)
             colour = guide_legend(order = 1),  shape = guide_legend(order = 1)) +
      theme_light() +
      labs( x = "Time post exposure to IVM (days)", y = "Proportion of dead mosquitoes") +
      theme(
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 8), legend.key.width = unit(1.5, "cm"),
        plot.title = element_text(hjust = 0, size = 16), axis.text = element_text(size = 12),
        axis.title = element_text(size = 11) ) +
      ggtitle("(D)")
    
  }
  
  #plot mu.m.ivm and HR: LongLasting and Oral
  Day_plot_LL = 155  
  Vtau_plot_LL = seq(0, Day_plot_LL, by = 0.5)  
  LegendTex_LL = c("F31-0.6", "Lam-1.0", "F31-0.6-Kis") 
  
  p_vecLongLasting06 = P_tau.fun(Vtau_plot_LL, tau50=ParamLongLasting06$tau50, alpha=ParamLongLasting06$alpha,
                                 Dmax=ParamLongLasting06$Dmax)
  mu_m_ivmLongLasting06 = mu_m - log(1 - p_vecLongLasting06)
  HRLongLasting06 = mu_m_ivmLongLasting06/mu_m
  
  p_vecLongLasting1 = P_tau.fun(Vtau_plot_LL, tau50=ParamLongLasting1$tau50, alpha=ParamLongLasting1$alpha,
                                Dmax=ParamLongLasting1$Dmax)
  mu_m_ivmLongLasting1 = mu_m - log(1 - p_vecLongLasting1)
  HRLongLasting1 = mu_m_ivmLongLasting1/mu_m
  
  p_vecLongLasting06Kis = P_tau.fun(Vtau_plot_LL, tau50=ParamLongLasting06Kis$tau50, alpha=ParamLongLasting06Kis$alpha,
                                    Dmax=ParamLongLasting06Kis$Dmax)
  mu_m_ivmLongLasting06Kis = mu_m - log(1 - p_vecLongLasting06Kis)
  HRLongLasting06Kis = mu_m_ivmLongLasting06Kis/mu_m
  
  # Construction des dataframes pour Long Lasting
  mu_m_ivm_dataLongLasting = rbind(  data.frame(tau = Vtau_plot_LL, mu_m_ivm = mu_m_ivmLongLasting06, Formulation = LegendTex_LL[1]),
                                     data.frame(tau = Vtau_plot_LL, mu_m_ivm = mu_m_ivmLongLasting1, Formulation = LegendTex_LL[2]),
                                     data.frame(tau = Vtau_plot_LL, mu_m_ivm = mu_m_ivmLongLasting06Kis, Formulation = LegendTex_LL[3])
  )
  
  HR_dataLongLasting = rbind( data.frame(tau = Vtau_plot_LL, HR = HRLongLasting06, Formulation = LegendTex_LL[1]),
                              data.frame(tau = Vtau_plot_LL, HR = HRLongLasting1, Formulation = LegendTex_LL[2]),
                              data.frame(tau = Vtau_plot_LL, HR = HRLongLasting06Kis, Formulation = LegendTex_LL[3])
  )
  
  # Calculs pour les formulations orales
  Day_plot_Oral = 155  
  Vtau_plot_Oral = seq(0, Day_plot_Oral, by = 0.1)
  LegendTex_Oral = c("BOHEMIA", "RIMDAMAL", "KamauRIMDAMAL") 
  
  p_vecBohemia = P_tau.fun(Vtau_plot_Oral, tau50=ParamBohemia$tau50, alpha=ParamBohemia$alpha,
                           Dmax=ParamBohemia$Dmax)
  mu_m_ivmBohemia = mu_m - log(1 - p_vecBohemia)
  HRBohemia = mu_m_ivmBohemia/mu_m
  
  p_vecRimdamal = P_tau.fun(Vtau_plot_Oral, tau50=ParamRimdamal$tau50, alpha=ParamRimdamal$alpha,
                            Dmax=ParamRimdamal$Dmax)
  mu_m_ivmRimdamal = mu_m - log(1 - p_vecRimdamal)
  HRRimdamal = mu_m_ivmRimdamal/mu_m
  
  p_vecKamauRIMDAMAL = P_tau.fun(Vtau_plot_Oral, tau50=ParamKamauRIMDAMAL$tau50, alpha=ParamKamauRIMDAMAL$alpha,
                                 Dmax=ParamKamauRIMDAMAL$Dmax)
  mu_m_ivmKamauRIMDAMAL = mu_m - log(1 - p_vecKamauRIMDAMAL)
  HRKamauRIMDAMAL = mu_m_ivmKamauRIMDAMAL/mu_m
  
  # Construction des dataframes pour formulations orales
  mu_m_ivm_dataOral = rbind( data.frame(tau = Vtau_plot_Oral, mu_m_ivm = mu_m_ivmBohemia, Formulation = LegendTex_Oral[1]),
                             data.frame(tau = Vtau_plot_Oral, mu_m_ivm = mu_m_ivmRimdamal, Formulation = LegendTex_Oral[2]),
                             data.frame(tau = Vtau_plot_Oral, mu_m_ivm = mu_m_ivmKamauRIMDAMAL, Formulation = LegendTex_Oral[3])
  )
  
  HR_Oral = rbind( data.frame(tau = Vtau_plot_Oral, HR = HRBohemia, Formulation = LegendTex_Oral[1]),
                   data.frame(tau = Vtau_plot_Oral, HR = HRRimdamal, Formulation = LegendTex_Oral[2]),
                   data.frame(tau = Vtau_plot_Oral, HR = HRKamauRIMDAMAL, Formulation = LegendTex_Oral[3])
  )
  
  # Plots mu_m_ivm 
  Plot_mu_m_ivmLongLasting = ggplot(mu_m_ivm_dataLongLasting, aes(x = tau, y = mu_m_ivm, color = Formulation)) +
    geom_line(linewidth = 1.2) +
    scale_y_continuous(breaks = seq(0, 8, 2), labels = as.character(seq(0, 8, 2)), 
                       expand = c(0.01, 0),  limits = c(0, 8.5)) +
    scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150), 
                       labels = c("0", "30","60", "90", "120","150"),
                       limits = c(0, 155), expand = c(0.01, 0), minor_breaks = NULL) +
    theme_light() +
    scale_color_manual(values = c("F31-0.6" = "#ff3355", "Lam-1.0" = "#037153", "F31-0.6-Kis" = "#8B4513"),
                       labels = c("F31-0.6" = "mdc-STM-001-0.6-vk5", "Lam-1.0" = "mdc-STM-001-1.0", "F31-0.6-Kis" = "mdc-STM-001-0.6-Kis")) +
    labs( x = "Time post exposure to IVM (days)", y = "Mosquitoes mortality rate",
          color = "IVM formulation" ) +
    theme( legend.position = "bottom",legend.title = element_blank(),
           legend.text = element_text(size = 9), plot.title = element_text(hjust = 0, size = 16),
           axis.text = element_text(size = 12), axis.title = element_text(size = 11)  ) +
    ggtitle("(B)")
  
  # Plot mu_m_ivm oral 
  Plot_mu_m_ivm_oral = ggplot(mu_m_ivm_dataOral, aes(x = tau, y = mu_m_ivm, color = Formulation)) +
    geom_line(linewidth = 1.2) +
    scale_y_continuous(breaks = seq(0, 4, 1),labels = as.character(seq(0, 4, 1)), 
                       expand = c(0.01, 0), limits = c(0, 4.25)) +
    scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150), 
                       labels = c("0", "30","60", "90", "120","150"),
                       limits = c(0, 155), expand = c(0.01, 0), minor_breaks = NULL) +
    theme_light() +
    scale_color_manual(values = c("BOHEMIA" = "#D98E04", "RIMDAMAL" = "#205072", "KamauRIMDAMAL" = "#8B008B"),
                       labels = c("BOHEMIA" = "BOHEMIA", "KamauRIMDAMAL" = expression("RIMDAMAL II" [K]), "RIMDAMAL" = expression("RIMDAMAL II" [S]))) +
    labs(  x = "Time post exposure to IVM (days)", y = "Mosquitoes mortality rate",
           color = "IVM formulation" ) +
    theme(  legend.position = "bottom", legend.title = element_blank(),
            legend.text = element_text(size = 9), plot.title = element_text(hjust = 0, size = 16),
            axis.text = element_text(size = 12), axis.title = element_text(size = 11)  ) +
    ggtitle("(E)")
  
  # Plots HR 
  Plot_HRLongLasting = ggplot(HR_dataLongLasting, aes(x = tau, y = HR, color = Formulation)) +
    geom_line(linewidth = 1.2) +
    scale_y_continuous(breaks = c(0, 4, 20, 40, 60),
                       labels = c("0", "4", "20", "40","60"), 
                       expand = c(0.01, 0), limits = c(0, 60)) +
    scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150), 
                       labels = c("0", "30","60", "90", "120","150"),
                       limits = c(0, 155), expand = c(0.01, 0), minor_breaks = NULL) +
    # geom_hline(yintercept = 4, linetype = "dashed", color = "black", linewidth = 0.6) +
    theme_light() +
    scale_color_manual(values = c("F31-0.6" = "#ff3355", "Lam-1.0" = "#037153", "F31-0.6-Kis" = "#8B4513"),
                       labels = c("F31-0.6" = "mdc-STM-001-0.6-vk5", "Lam-1.0" = "mdc-STM-001-1.0", "F31-0.6-Kis" = "mdc-STM-001-0.6-Kis")) +
    labs(
      x = "Time post exposure to IVM (days)", y = "Crude hazard ratio",
      color = "IVM formulation"
    ) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      plot.title = element_text(hjust = 0, size = 16),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 11)
    ) +
    ggtitle("(C)")
  
  # Plot HR oral 
  Plot_HR_oral = ggplot(HR_Oral, aes(x = tau, y = HR, color = Formulation)) +
    geom_line(linewidth = 1.2) +
    scale_y_continuous(breaks = c(0,  4,  10 ,  20, 30),
                       labels = c("0",  "4", "10", "20", "30"), 
                       expand = c(0.02, 0), limits = c(0, 30)) +
    scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150), 
                       labels = c("0", "30","60", "90", "120","150"),
                       limits = c(0, 155), expand = c(0.01, 0), minor_breaks = NULL) +
    # geom_hline(yintercept = 4, linetype = "dashed", color = "black", linewidth = 0.6) +
    theme_light() +
    scale_color_manual(values = c("BOHEMIA" = "#D98E04", "RIMDAMAL" = "#205072", "KamauRIMDAMAL" = "#8B008B"),
                       labels = c("BOHEMIA" = "BOHEMIA", "RIMDAMAL" = expression("RIMDAMAL II" [S]), "KamauRIMDAMAL" = expression("RIMDAMAL II" [K]))) +
    labs(
      x = "Time post exposure to IVM (days)", y = "Crude hazard ratio",
      color = "IVM formulation") +
    theme( legend.position = "bottom", legend.title = element_blank(),
           legend.text = element_text(size = 9), plot.title = element_text(hjust = 0, size = 16),
           axis.text = element_text(size = 12),  axis.title = element_text(size = 11)  ) +
    ggtitle("(F)")
  
  # Création de la liste des graphiques
  plot_list_data = list(Plot_Prop_dead_mosquitoesLongLasting, Plot_mu_m_ivmLongLasting, Plot_HRLongLasting,
                        Plot_Prop_dead_mosquitoesOral, Plot_mu_m_ivm_oral, Plot_HR_oral)
  
  # plot_list_data = list(Plot_Prop_dead_mosquitoesLongLasting, Plot_mu_m_ivmLongLasting,
  #                       Plot_Prop_dead_mosquitoesOral, Plot_mu_m_ivm_oral)
  
  # Combinaison finale avec patchwork
  combined_plot = plot_list_data[[1]] + plot_list_data[[2]] + plot_list_data[[3]] +
    plot_list_data[[4]] + plot_list_data[[5]]+ plot_list_data[[6]] +  plot_layout(ncol = 3, nrow = 2)
}
combined_plot
ggsave("plot_data_Estimaton.pdf", combined_plot, width = 16.5, height = 10.2, units = "in", dpi = 300)


##Residual plots Figure S3
{
  # Long Lasting (strategies 0, 1, 5)
  resid_LL06 = residual_fun(ParamLongLasting06$model, "F31-0.6")
  resid_LL1 = residual_fun(ParamLongLasting1$model, "Lam-1.0")
  resid_LL06Kis = residual_fun(ParamLongLasting06Kis$model, "F31-0.6-Kis")
  
  combined_resid_LL = bind_rows(resid_LL06, resid_LL1, resid_LL06Kis)
  
  # Oral (strategies 2, 3, 4)
  resid_Bohemia = residual_fun(ParamBohemia$model, "BOHEMIA")
  resid_Rimdamal = residual_fun(ParamRimdamal$model, "RIMDAMAL")
  resid_Kamau = residual_fun(ParamKamauRIMDAMAL$model, "KamauRIMDAMAL")
  
  combined_resid_Oral = bind_rows(resid_Bohemia, resid_Rimdamal, resid_Kamau)
  
  mean_resid = data.frame(
    Strategy = c("mdc-STM-001-0.6-vk5", "mdc-STM-001-1.0", "mdc-STM-001-0.6-Kis", 
                 "BOHEMIA", "RIMDAMAL II (S)", "RIMDAMAL II (K)"),
    Mean_residuals = c(
      mean(resid_LL06$residuals),
      mean(resid_LL1$residuals),
      mean(resid_LL06Kis$residuals),
      mean(resid_Bohemia$residuals),
      mean(resid_Rimdamal$residuals),
      mean(resid_Kamau$residuals)
    )
  )
  
  print(mean_resid) 
  # if mean_resid close to 0, then Residual analysis confirmed adequate model fit, no systematic patterns, and unbiased predictions.
  
  # Plot Residuals vs Fitted - Long Lasting
  Plot_Resid_LL = ggplot(combined_resid_LL, aes(x = fitted, y = residuals, color = Formulation, shape = Formulation)) +
    geom_point(size = 2, stroke = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
    scale_y_continuous(breaks = seq(-0.5, 0.5, 0.25), labels = as.character(seq(-0.5, 0.5, 0.25)), 
                       expand = c(0.01, 0), limits = c(-0.65, 0.55)) +
    scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1), labels = c("0.2", "0.4","0.6", "0.8", "1"),
                       limits = c(0.2, 1), expand = c(0.01, 0)) +
    scale_color_manual(values = c("F31-0.6" = "#ff3355", "Lam-1.0" = "#037153", "F31-0.6-Kis" = "#8B4513"),
                       labels = c("F31-0.6" = "mdc-STM-001-0.6-vk5", "Lam-1.0" = "mdc-STM-001-1.0", "F31-0.6-Kis" = "mdc-STM-001-0.6-Kis")) +
    scale_shape_manual(values = c("F31-0.6" = 16, "Lam-1.0" = 18, "F31-0.6-Kis" = 17),
                       labels = c("F31-0.6" = "mdc-STM-001-0.6-vk5", "Lam-1.0" = "mdc-STM-001-1.0", "F31-0.6-Kis" = "mdc-STM-001-0.6-Kis")) +
    theme_light() +
    labs(x = "Fitted values", y = "Residuals") +
    theme(legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = 9), plot.title = element_text(hjust = 0, size = 16),
          axis.text = element_text(size = 12), axis.title = element_text(size = 11)) +
    ggtitle("(A)")
  
  # Plot Residuals vs Fitted - Oral
  Plot_Resid_Oral = ggplot(combined_resid_Oral, aes(x = fitted, y = residuals, color = Formulation, shape = Formulation)) +
    geom_point(size = 2.5, stroke = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
    scale_y_continuous(breaks = seq(-0.05, 0.075, 0.025), labels = as.character(seq(-0.05, 0.075, 0.025)), 
                       expand = c(0.01, 0), limits = c(-0.055, 0.08)) +
    scale_x_continuous(breaks = c(0,0.2, 0.4, 0.6, 0.8, 1),labels = c("0", "0.2", "0.4","0.6", "0.8", "1"),
                       limits = c(0, 1), expand = c(0.01, 0)) +
    scale_color_manual(values = c("BOHEMIA" = "#D98E04", "RIMDAMAL" = "#205072", "KamauRIMDAMAL" = "#8B008B"),
                       labels = c("BOHEMIA" = "BOHEMIA", "RIMDAMAL" = expression("RIMDAMAL II"[S]), "KamauRIMDAMAL" = expression("RIMDAMAL II"[K]))) +
    scale_shape_manual(values = c("BOHEMIA" = 16, "RIMDAMAL" = 18, "KamauRIMDAMAL" = 17),
                       labels = c("BOHEMIA" = "BOHEMIA", "RIMDAMAL" = expression("RIMDAMAL II"[S]), "KamauRIMDAMAL" = expression("RIMDAMAL II"[K]))) +
    theme_light() +
    labs(x = "Fitted values", y = "Residuals") +
    theme(legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = 9), plot.title = element_text(hjust = 0, size = 16),
          axis.text = element_text(size = 12), axis.title = element_text(size = 11)) +
    ggtitle("(B)")
  
  combined_diagnostic = Plot_Resid_LL + Plot_Resid_Oral + plot_layout(ncol = 2, nrow = 1)
}
combined_diagnostic

ggsave("diagnostic_residuals.pdf", combined_diagnostic, width = 8.5, height = 4, units = "in", dpi = 300)

#RMSE Values
rmse_results = bind_rows(
  rmse_fun(ParamLongLasting06$model, "mdc-STM-001-0.6-vk5"),
  rmse_fun(ParamLongLasting1$model, "mdc-STM-001-1.0"),
  rmse_fun(ParamLongLasting06Kis$model, "mdc-STM-001-0.6-Kis"),
  rmse_fun(ParamBohemia$model, "BOHEMIA"),
  rmse_fun(ParamRimdamal$model, "RIMDAMAL II (S)"),
  rmse_fun(ParamKamauRIMDAMAL$model, "RIMDAMAL II (K)")
)

print(rmse_results)

################# FIGURE S4, the initial population Sh0, natural mortality and nu_h 
{
  plot_list_parms <- list()
  
  #plot Initial population Sh0
  {
    # Initial human population Sh0
    Int_Humanpop <- 8136.10 * (c(12.9, 12.5, 11.5, 13.1, 11.9, 9.3, 7.3, 5.5, 4.4,
                                 3.2, 2.6, 1.8, 1.4, 0.9, 0.7, 0.3, 0.2, 0.2) / 0.997)
    
    # Create a dataframe for Sh0
    Sh0 <- data.frame(
      Age = seq(2.5, 87.5, by = 5),  # Adjust sequence to match length of Int_Humanpop
      Population = Int_Humanpop)
    
    plot_Sh0=ggplot(Sh0, aes(x = Age, y = Population, fill = "Population S[h0]")) +
      geom_bar(stat = "identity") +
      theme_light() +
      scale_x_continuous(breaks = seq(0, 90, by = 5), minor_breaks = NULL,
                         expand = c(0.02, 0), limits = c(0, 90) ) +
      scale_y_continuous( breaks = c(0, 5000, 15000, 30000, 50000, 75000, 100000),
                          labels = c("0", expression(5 %*% 10^3), expression(1.5 %*% 10^4),
                                     expression(3 %*% 10^4), expression(5 %*% 10^4), 
                                     expression(7.5 %*% 10^4), expression(10^5)),
                          minor_breaks = NULL, limits = c(0, 107000) ) +
      labs( x = "Age (years)",  y = "Population", title = "(A)", fill = "") +
      scale_fill_manual( values = c("Population S[h0]" = "#148F77"),
                         labels = expression("Initial human susceptible population "*S[h0]))+
      theme( legend.position = "bottom", # Position de la légende
             legend.title = element_text(size = 12, face = "bold"), 
             legend.text = element_text(size = 12), 
             legend.key.size = unit(0.3, "cm"), # dimension des carrés
             plot.title = element_text(hjust = 0, size = 16), 
             panel.grid.major.x = element_line(color = "gray", linewidth = 0.1)
      ) +
      guides(fill = guide_legend(override.aes = list(alpha = 0.8)))
  }
  
  #plot  mu_h
  {
    # Natural mortality rate of human
    mu_h = c(66.8, 7.6, 1.7, 0.9, 1.3, 1.9, 2.4, 2.8, 3.6, 4.7, 6.3, 8.9, 13.2, 19.8,
             31.1, 47.7, 71.3, 110.5, 186.7) / 1000
    
    # Corresponding ages
    ages <- seq(2.5, 92.5, by = 5)
    
    # Create a dataframe for the mortality rate plot
    mu_h <- data.frame(
      Age = ages,
      MortalityRate = mu_h
    )
    
    
    Age_labels <- c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 
                    60, 65, 70, 75, 80, 85,90)
    
    plot_mu_h=ggplot(mu_h, aes(x = Age, y = MortalityRate, fill = "mu_h")) +
      geom_bar(stat = "identity") +
      theme_light() +
      labs( x = "Age (years)",  y = "Mortality rate", title = "(B)", fill = "") +
      scale_x_continuous(breaks =  seq(0, 95, by = 5), labels = Age_labels,
                         minor_breaks = NULL, expand = c(0.02, 0),limits = c(0, 95)) +
      scale_y_continuous(breaks = seq(0, 0.2, by = 0.05),
                         minor_breaks = NULL, limits = c(0, 0.21)) +
      scale_fill_manual( values = c("mu_h" = "#148F77"),
                         labels = expression("Natural mortality rate for humans "*mu[h])) +
      theme( legend.position = "bottom", # Position de la légende
             legend.title = element_text(size = 12, face = "bold"), 
             legend.text = element_text(size = 12), 
             legend.key.size = unit(0.3, "cm"), # dimension des carrés
             plot.title = element_text(hjust = 0, size = 16), 
             panel.grid.major.x = element_line(color = "gray", linewidth = 0.1)
      ) +
      guides(fill = guide_legend(override.aes = list(alpha = 0.8)))
  }
  
  #plot  nu_hm nu_hf et gamma_h
  {
    #nu_hm(a), nu_hf(a): progression Asympto vers cas clinique
    {
      p_f = 0.15 # Fertility probability for women aged 15 to 40 in Bobo Dioulasso,
      nu_0=1/10; nu_1=1/50;nu_2=1/5
      eps=0.0005
      nu_hm=rep(NA,Nah) # disease mortality for men
      nu_hf = rep(NA, Nah) # disease mortality for female
      for (i in 1:Nah) {
        if (Vah[i] <= 5) { nu_hm[i] = nu_2; nu_hf[i] = nu_2+eps } 
        else if (Vah[i] > 5 && Vah[i] <= 15) {nu_hm[i] = nu_0; nu_hf[i] = nu_0+eps }
        else {nu_hm[i] = nu_1 }
        
        if (Vah[i] > 15 && Vah[i] <= 40) { nu_hf[i] = p_f * nu_0 + (1 - p_f) * nu_1 } 
        else if (Vah[i] > 40) { nu_hf[i] = nu_1+eps}
      }
    }
    
    #gamma_h(a): guerison
    {
      # Recovery time ranges for each age group (days)
      recovery_ranges = list(c(163, 345), c(555, 714), c(344, 400), c(181, 204), 
                             c(82, 92), c(56, 61), c(48, 55))
      # Age ranges  corresponding to the recovery time ranges
      age_ranges =list(c(0, 1),c(1, 5),c(5, 8),c(8, 18),c(18, 28),c(28, 43),c(43, ah_max))
      gamma_h = rep(NA, Nah)
      
      #  value of gamma_h
      for (i in 1:Nah) {
        for (j in 1: length(age_ranges)) {
          age_min = age_ranges[[j]][1]
          age_max = age_ranges[[j]][2]
          if(Vah[i] >= age_min & Vah[i] <= age_max)
          {gamma_h[i] = 1 / mean(recovery_ranges[[j]])}
        }
      }
    }
    
    df_segments <- data.frame(xstart = Vah[-length(Vah)], xend = Vah[-1],
                              nu_hm = nu_hm[-length(nu_hm)], nu_hf= nu_hf[-length(nu_hf)],
                              gamma_h=gamma_h[-length(gamma_h)])
    
    cut_age = seq(0, ah_max, by = 5)
    IdSeq <- seq(1, length(Vah), by = 2)
    
    df_points <- data.frame(Age = Vah[IdSeq],nu_hm = nu_hm[IdSeq], 
                            nu_hf= nu_hf[IdSeq],gamma_h=gamma_h[IdSeq])
    
    coeff=10
    plot_nu_gamma=ggplot() +
      geom_segment(data = df_segments, aes(x = xstart, xend = xend, y = nu_hf, yend = nu_hf, color = "nu_hf"), size = 1) +
      geom_point(data = df_points, aes(x = Age, y = nu_hf, color = "nu_hf", shape = "nu_hf"), size = 1.5) +
      geom_segment(data = df_segments, aes(x = xstart, xend = xend, y = nu_hm, yend = nu_hm, color = "nu_hm"), size = 1) +
      geom_point(data = df_points, aes(x = Age, y = nu_hm, color = "nu_hm", shape = "nu_hm"), size = 1.5) +
      geom_segment(data =df_segments,aes(x=xstart,xend=xend,y=gamma_h*coeff,yend=gamma_h*coeff,color="gamma_h"),size = 1) +
      geom_point(data = df_points, aes(x = Age, y = gamma_h * coeff, color = "gamma_h", shape = "gamma_h"), size = 1.5) +
      theme_light() +
      scale_x_continuous(breaks = cut_age, minor_breaks = NULL, expand = c(0.02, 0),
                         limits = c(0, 90))+
      scale_y_continuous( name = "Asymptomatic progression rate",
                          sec.axis = sec_axis(~ . / coeff, name = "Recovery rate")  ) +
      scale_color_manual(
        values = c("nu_hf" = "#76448A", "nu_hm" = "#148F77", "gamma_h" = "#F39C12"),
        labels = c(expression(gamma[h]), expression(nu[h]^bold(italic("♀"))), expression(nu[h]^bold(italic("♂"))))
      ) +
      scale_shape_manual(values = c("nu_hf" = 8, "nu_hm" = 16, "gamma_h" = 18),
                         labels = "" ) +
      labs( x = "Age (years)",y = "", title = "(C)", fill = "" )+
      theme( 
        legend.position = "bottom",
        legend.title = element_blank(), 
        legend.text = element_text(size = 10), 
        legend.spacing.x = unit(2, 'cm'),
        legend.key.width = unit(1.1, "cm"),
        legend.key.size = unit(0.5, "cm"), #dimension des carrés en legende
        plot.title = element_text(hjust = 0, size = 16), 
        plot.subtitle = element_text(hjust = 0.5, size = 14), 
        panel.grid.major.x = element_line(color = "gray", linewidth = 0.2)
      ) +
      guides(color = guide_legend(override.aes = list(shape = c(18, 8, 16), size = 3)),
             shape = "none" )
  }
  
  plot_list_parms=list(plot_Sh0, plot_mu_h, plot_nu_gamma)
  
  combined_plot=plot_list_parms[[1]]+ plot_list_parms[[2]]+plot_list_parms[[3]]+plot_layout(ncol = 3)
}
combined_plot
ggsave("plot_parameters_&_Initial_pop.pdf", combined_plot, width = 12, height = 4, units = "in",dpi = 300)

#Figure S15: Time since first IVM administration
{
  VectTime_between_cycles_LongLasting06 = c(60.56, 56.31, 56.31) # Strategy 0
  VectTime_between_cycles_LongLasting06_kis = c(63.32, 59.76, 60.03) # Strategy 5
  VectTime_between_cycles_LongLasting1 = c(63.53, 59.76, 60.29)  # Strategy 1  
  VectTime_between_cycles_Bohemia = c(39.95, 24.22)             # Strategy 2 - toujours 3 cycles
  VectTime_between_cycles_Rimdamal =  c(46.53, 30.60, 21.67)     # Strategy 3 - toujours 4 cycles
  VectTime_between_cycles_KamauRimdamal =  c(39.31, 23.58, 17.00)    # Strategy 4 - toujours 4 cycles
  
  
  calculate_multi_cycle_probability= function(strategy,n_cycles_override=NULL,max_time = 250) {
    ParamStrategy=ModelEtau50Alpha(strategy)
    tau50=ParamStrategy$tau50 
    alpha=ParamStrategy$alpha
    ParamRho = RhoFunction(Nah,Ntau,Vtau,Vah,IVM_field_dependancy,strategy, tau50, alpha)
    single_cycle_prob = ParamRho$PropHuman_in_IVM_Group1
    if(strategy == 0) {
      cycle_intervals = VectTime_between_cycles_LongLasting06
      n_cycles = ifelse(is.null(n_cycles_override), 4, n_cycles_override)
    } else if(strategy == 1) {
      cycle_intervals = VectTime_between_cycles_LongLasting1
      n_cycles = ifelse(is.null(n_cycles_override), 4, n_cycles_override)
    } else if(strategy == 2) {
      cycle_intervals = VectTime_between_cycles_Bohemia
      n_cycles = 3  # Toujours 3 cycles pour BOHEMIA
    } else if(strategy == 3) {
      cycle_intervals = VectTime_between_cycles_Rimdamal
      n_cycles = 4  # Toujours 4 cycles pour RIMDAMAL II_S
    } else if(strategy == 4) {
      cycle_intervals = VectTime_between_cycles_KamauRimdamal
      n_cycles = 4  # Toujours 4 cycles pour RIMDAMAL II_K
    } else if(strategy == 5) {
      cycle_intervals = VectTime_between_cycles_LongLasting06_kis
      n_cycles = ifelse(is.null(n_cycles_override), 4, n_cycles_override)
    }
    
    if(n_cycles > 1) {
      cycle_intervals = cycle_intervals[1:(n_cycles-1)]
    } else {
      cycle_intervals = c()
    }
    
    time_full = seq(0, max_time, 1)
    prob_full = rep(0, length(time_full))
    
    cycle_starts = c(0)
    for(i in 1:length(cycle_intervals)) {
      cycle_starts = c(cycle_starts, cycle_starts[length(cycle_starts)] + cycle_intervals[i])
    }
    
    for(t_idx in 1:length(time_full)) {
      current_time = time_full[t_idx]
      
      current_cycle = 1
      time_since_last_cycle = current_time
      
      for(cycle in 1:n_cycles) {
        if(cycle <= length(cycle_starts) && current_time >= cycle_starts[cycle]) {
          current_cycle = cycle
          time_since_last_cycle = current_time - cycle_starts[cycle]
        }
      }
      
      if(current_cycle > n_cycles || (length(cycle_starts) >= n_cycles && current_time > cycle_starts[n_cycles])) {
        if(length(cycle_starts) >= n_cycles) {
          final_time_in_last_cycle = current_time - cycle_starts[n_cycles]
          if(final_time_in_last_cycle < length(single_cycle_prob)) {
            prob_full[t_idx] = single_cycle_prob[final_time_in_last_cycle + 1]
          } else {
            prob_full[t_idx] = single_cycle_prob[length(single_cycle_prob)]
          }
        } else {
          prob_full[t_idx] = 0
        }
      } else {
        # Calculer la probabilité pour le cycle actuel
        if(time_since_last_cycle < length(single_cycle_prob)) {
          prob_full[t_idx] = single_cycle_prob[time_since_last_cycle + 1]
        } else {
          prob_full[t_idx] = single_cycle_prob[length(single_cycle_prob)]
        }
      }
    }
    
    return(list(time = time_full, probability = prob_full, cycle_starts = cycle_starts[1:n_cycles]))
  }
  
  FigName = "Nc_PropHuman_in_IVM_withoutAge.pdf"
  pdf(FigName, width = 13.5, height = 3.8)
  
  par(oma = c(3.5, 4, 3, 1), mar = c(3, 2.5, 2.5, 1))
  par(mfrow = c(1, 4))
  
  ColVect = c("#ff3355",  "#8B4513", "#037153", "#205072", "#8B008B", "#D98E04")
  strategy_to_color = c(1, 3, 6, 4, 5, 2)  # indices dans ColVect pour strategies 0,1,2,3,4,5
  strategy_to_lty = c(1, 1, 2, 4, 3, 1)   # lty pour strategies 0,1,2,3,4,5
  
  FormulationNames = c("mdc-STM-001-0.6-vk5", "mdc-STM-001-0.6-kis", "mdc-STM-001-1.0", 
                       "RIMDAMAL II_S", "RIMDAMAL II_K", "BOHEMIA")
  LC = 2
  
  for(panel in 1:4) {
    n_cycles_variable = panel  # 1, 2, 3, ou 4 cycles
    
    # Créer le graphique
    plot(-1, 1, type = "l", xlab = "", xlim = c(0, 250), ylab = "", ylim = c(0, 1.05), 
         cex.lab = 1, xaxt = "n", yaxt = "n")
    
    # Axes
    axis(2, at = seq(0, 1, by = 0.2), labels = seq(0, 1, by = 0.2), cex.axis = 0.7, las = 1)
    axis(1, at = seq(0, 250, by = 50), labels = seq(0, 250, by = 50), cex.axis = 0.7, las = 1)
    
    # Tracer chaque stratégie (0, 1, 2, 3, 4, 5)
    for(strategy in c(0, 1, 2, 3, 4, 5)) {
      if(strategy %in% c(0, 1, 5)) {
        # Pour strategies 0, 1, et 5, utiliser le nombre de cycles variable
        result = calculate_multi_cycle_probability(strategy, n_cycles_override = n_cycles_variable, max_time = 250)
      } else {
        # Pour strategies 2, 3, et 4, garder leurs cycles fixes
        result = calculate_multi_cycle_probability(strategy, max_time = 250)
      }
      
      color_idx = strategy_to_color[strategy + 1]
      lty_val = strategy_to_lty[strategy + 1]
      
      lines(result$time, result$probability, col = ColVect[color_idx], lty = lty_val, lwd = LC)
      
      # # Ajouter des lignes verticales pour marquer le début de chaque cycle
      # for(start_time in result$cycle_starts[-1]) { # Exclure le temps 0
      #   if(start_time <= 250) {
      #     abline(v = start_time, col = ColVect[color_idx], lty = 3, lwd = 0.8, alpha = 0.7)
      #   }
      # }
    }
    
    # Titre spécifique pour chaque panel
    mtext(bquote("#cycles for mdc-STM-001: " ~ n[c] == .(n_cycles_variable)),
          side = 3, adj = 0.5, cex = 0.85, line = 1, font = 1)
    
    # Label des panels
    par(xpd = NA)
    text(-0.05 * 250, 1.05 * 1.15, paste("(", LETTERS[panel], ")", sep = ""), cex = 1.6)
    par(xpd = FALSE)
    
    mtext("Time post first exposure to IVM (Days)", side = 1, adj = 0.5, cex = 0.8, line = 2, font = 1)
    # Labels Y seulement sur le premier panel
    if(panel == 1) {
      mtext("Probability of remaining exposed to IVM", side = 2, adj = 0.5, cex = 0.8, line = 2.5, font = 1)
    }
  }
  
  # Plot legend
  LEGEND = c("mdc-STM-001-0.6-vk5", "mdc-STM-001-0.6-kis", "mdc-STM-001-1.0", 
             expression("RIMDAMAL II"[S]), expression("RIMDAMAL II"[K]), "BOHEMIA") 
  reset()
  legend("bottom", legend=LEGEND,
         xpd = NA, horiz = FALSE, #inset = c(-3,-0.6),
         bty = "n", lty = c(1, 1, 1, 4, 3, 2), lw=LC,
         col = c("#ff3355","#8B4513", "#037153", "#205072","#8B008B", "#D98E04"),ncol =3,  cex = 0.85)
}
dev.off()

#Global sensitivity
{
  GainPrev_sensitivity <- function(strategy, PropIVM, IVM_field_dependancy, 
                                   Number_of_cycle, IVM_Pregnancy, p_f,
                                   DeltaOptimInit) {
    
    ParamIvmFormulation = ModelEtau50Alpha(strategy)
    tau50 = ParamIvmFormulation$tau50
    alpha = ParamIvmFormulation$alpha
    Dmax = ParamIvmFormulation$Dmax
    
    if (strategy %in% c(0, 1, 5)) {
      VectTime_between_cycles = rep(60, Number_of_cycle - 1)
    } else {
      VectTime_between_cycles = rep(30, Number_of_cycle - 1)
    }
    if (Number_of_cycle == 1) VectTime_between_cycles = 0
    
    ModelOutput = ModelIVM(PropIVM, VectTime_between_cycles, Number_of_cycle, t_begin_Camp, Dur_cycle,
                           mu_h, delta_h, bar_beta_h, beta_h, gamma_h, theta, beta_m, Wedge_m, mu_m,
                           ah_max, dah, Nah, tmax, Gap, dt, time, Ntime, p_f, IVM_Pregnancy,
                           dtau, tau_max, Vtau, Ntau, dsigma, sigma_max, Vsigma, Nsigma, strategy,
                           IVM_field_dependancy, tau50, alpha, Dmax)
    
    GainPrev = 1 - ModelOutput$DeltaOptim / DeltaOptimInit
    return(GainPrev)
  }
  
  run_sensitivity_analysis <- function(Prev_mos10 = 1) {
    
    VectStrategy = c(0, 1, 2, 3, 4, 5)
    VectIVM_field_dependancy = c(0, 1)
    VectPropIVM = c(0.5, 0.7, 0.9)
    VectNumber_of_cycle = c(1, 2, 3, 4)
    VectIVM_Pregnancy = c(0, 1)
    Vectp_f = c(0.1, 0.15, 0.2)
    
    strategy_names = c("LAIF_06vk5", "LAIF_10", "BOHEMIA", "RIMDAMAL_S", "RIMDAMAL_K", "LAIF_06kis")
    names(strategy_names) = as.character(0:5)
    
    params_Mos_Prev(Prev_mos10)
    
    SizePerPf = length(VectStrategy) * length(VectIVM_field_dependancy) * 
      length(VectPropIVM) * length(VectNumber_of_cycle) * length(VectIVM_Pregnancy)
    
    print(paste("Simulations par p_f:", SizePerPf))
    print(paste("Total simulations:", SizePerPf * length(Vectp_f)))
    
    for (p_f_val in Vectp_f) {
      
      print(paste("=== Traitement p_f =", p_f_val, "==="))
      
      # Calculer DeltaOptimInit pour ce p_f (run sans IVM)
      strategy = 0; PropIVM = 0; IVM_field_dependancy = 0
      ParamRef = ModelEtau50Alpha(strategy)
      Number_of_cycle=1
      IVM_Pregnancy=0
      VectTime_between_cycles= c(0,0,0)
      ModelNoIVM = ModelIVM(PropIVM, VectTime_between_cycles,Number_of_cycle, t_begin_Camp, Dur_cycle,mu_h, delta_h, bar_beta_h,
                            beta_h, gamma_h, theta, beta_m, Wedge_m, mu_m,ah_max, dah, Nah, tmax,
                            Gap, dt, time, Ntime, p_f_val, IVM_Pregnancy,dtau, tau_max,Vtau, Ntau, dsigma, 
                            sigma_max, Vsigma, Nsigma, strategy, IVM_field_dependancy, 
                            ParamRef$tau50, ParamRef$alpha, ParamRef$Dmax)
      DeltaOptimInit = sum(ModelNoIVM$Ih_Tot[ModelNoIVM$IndexOptim])
      
      # Initialiser colonnes
      ColStrategy = rep(NA, SizePerPf)
      ColIVM_field_dep = rep(NA, SizePerPf)
      ColPropIVM = rep(NA, SizePerPf)
      ColNumber_of_cycle = rep(NA, SizePerPf)
      ColIVM_Pregnancy = rep(NA, SizePerPf)
      ColGainPrev = rep(NA, SizePerPf)
      
      RunNumber = 0
      
      for (strategy in VectStrategy) {
        for (IVM_field_dep in VectIVM_field_dependancy) {
          for (PropIVM in VectPropIVM) {
            for (nc in VectNumber_of_cycle) {
              for (IVM_Preg in VectIVM_Pregnancy) {
                
                RunNumber = RunNumber + 1
                
                if (RunNumber %% 5 == 0) {
                  print(paste("p_f =", p_f_val, "- Run", RunNumber, "/", SizePerPf))
                }
                
                tryCatch({
                  GainPrev = GainPrev_sensitivity(strategy, PropIVM, IVM_field_dep,
                                                  nc, IVM_Preg, p_f_val, DeltaOptimInit)
                  
                  ColStrategy[RunNumber] = strategy_names[as.character(strategy)]
                  ColIVM_field_dep[RunNumber] = IVM_field_dep
                  ColPropIVM[RunNumber] = PropIVM
                  ColNumber_of_cycle[RunNumber] = nc
                  ColIVM_Pregnancy[RunNumber] = IVM_Preg
                  ColGainPrev[RunNumber] = GainPrev
                  
                }, error = function(e) {
                  ColStrategy[RunNumber] = strategy_names[as.character(strategy)]
                  ColIVM_field_dep[RunNumber] = IVM_field_dep
                  ColPropIVM[RunNumber] = PropIVM
                  ColNumber_of_cycle[RunNumber] = nc
                  ColIVM_Pregnancy[RunNumber] = IVM_Preg
                  ColGainPrev[RunNumber] = NA
                  print(paste("Erreur run", RunNumber, ":", e$message))
                })
              }
            }
          }
        }
      }
      
      # Creer tableau pour ce p_f
      Tableau_pf = data.frame( IVM_Regimen = ColStrategy,IVM_clearance = ColIVM_field_dep,
                               Prop_IVM = ColPropIVM, n_c = ColNumber_of_cycle,
                               IVM_Pregnancy = ColIVM_Pregnancy, p_f = p_f_val, GainPrev = ColGainPrev)
      
      # Save fichier CSV pour ce p_f
      nom_fichier = paste0("Sensitivity_GainPrev_pf_", gsub("\\.", "_", p_f_val), ".csv")
      write.table(Tableau_pf, file = nom_fichier, sep = ",", row.names = FALSE)
      print(paste("Fichier sauvegarde:", nom_fichier))
    }
    
    # Assembler tous les tableaux
    Tableau_complet = assembler_tableaux_sensitivity(Vectp_f)
    
    return(Tableau_complet)
  }
  
  #assembler les tableaux
  assembler_tableaux_sensitivity <- function(Vectp_f) {
    tous_tableaux = list()
    
    for (p_f_val in Vectp_f) {
      nom_fichier = paste0("Sensitivity_GainPrev_pf_", gsub("\\.", "_", p_f_val), ".csv")
      
      if (file.exists(nom_fichier)) {
        tableau = read.csv(nom_fichier, stringsAsFactors = FALSE)
        tous_tableaux[[as.character(p_f_val)]] = tableau
      }
    }
    
    if (length(tous_tableaux) > 0) {
      Tableau_complet = do.call(rbind, tous_tableaux)
      write.table(Tableau_complet, file = "Sensitivity_GainPrev_complet.csv", sep = ",", row.names = FALSE)
      
      # Exporter aussi en Excel
      wb = createWorkbook()
      
      # Feuille complete
      addWorksheet(wb, "All_data")
      writeData(wb, "All_data", Tableau_complet)
      
      # Feuille par p_f
      for (p_f_val in Vectp_f) {
        sheet_name = paste0("pf_", gsub("\\.", "_", p_f_val))
        addWorksheet(wb, sheet_name)
        df_pf = subset(Tableau_complet, p_f == p_f_val)
        writeData(wb, sheet_name, df_pf)
      }
      
      saveWorkbook(wb, "Sensitivity_GainPrev_Analysis.xlsx", overwrite = TRUE)
      print("Fichier Excel sauvegarde: Sensitivity_GainPrev_Analysis.xlsx")
      
      return(Tableau_complet)
    }
  }
  
  #RUN
  Tableau_sensitivity = run_sensitivity_analysis(Prev_mos10 = 1)
  
  tous_tableaux = list()
  for (p_f_val in Vectp_f) {
    nom_fichier = paste0("Sensitivity_GainPrev_pf_", gsub("\\.", "_", p_f_val), ".csv")
    if (file.exists(nom_fichier)) {
      tous_tableaux[[as.character(p_f_val)]] = read.csv(nom_fichier, stringsAsFactors = FALSE)
    }
  }
}

