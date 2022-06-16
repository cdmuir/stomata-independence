source("r/header.R")

# Source fitted model ----
fit12 = read_rds("objects/fit12.rds")

# Prepare for plotting ----
df_samples12 = fit12 |>
  as_draws_df(pars = c(
    
    # Trait averages
    "b_adaxialstomatallengthum_Intercept",
    "b_abaxialstomatallengthum_Intercept",
    "b_adaxialstomataldensitymm2_Intercept",
    "b_abaxialstomataldensitymm2_Intercept",
    
    # Trait (co)variances
    "b_sigma_adaxialstomatallengthum_Intercept",
    "b_sigma_abaxialstomatallengthum_Intercept",
    "b_sigma_adaxialstomataldensitymm2_Intercept",
    "b_sigma_abaxialstomataldensitymm2_Intercept",

    # Effect pair_age on trait (co)variances
    "b_sigma_adaxialstomatallengthum_pair_age",
    "b_sigma_abaxialstomatallengthum_pair_age",
    "b_sigma_adaxialstomataldensitymm2_pair_age",
    "b_sigma_abaxialstomataldensitymm2_pair_age",
    
    # Trait correlations
    "rescor__adaxialstomatallengthum__adaxialstomataldensitymm2",
    "rescor__adaxialstomatallengthum__abaxialstomatallengthum", 
    "rescor__adaxialstomataldensitymm2__abaxialstomatallengthum",
    "rescor__adaxialstomatallengthum__abaxialstomataldensitymm2",
    "rescor__adaxialstomataldensitymm2__abaxialstomataldensitymm2",
    "rescor__abaxialstomatallengthum__abaxialstomataldensitymm2",
    
    # nu parameter
    "nu"
    
  )) |>
  rename(
    r_adl_add = rescor__adaxialstomatallengthum__adaxialstomataldensitymm2,
    r_adl_abl = rescor__adaxialstomatallengthum__abaxialstomatallengthum, 
    r_add_abl = rescor__adaxialstomataldensitymm2__abaxialstomatallengthum,
    r_adl_abd = rescor__adaxialstomatallengthum__abaxialstomataldensitymm2,
    r_add_abd = rescor__adaxialstomataldensitymm2__abaxialstomataldensitymm2,
    r_abl_abd = rescor__abaxialstomatallengthum__abaxialstomataldensitymm2
  ) |>
  crossing(pair_age = c(0, 10)) |>
  mutate(
    
    # Standard deviations
    sigma_adl = exp(b_sigma_adaxialstomatallengthum_Intercept +
      b_sigma_adaxialstomatallengthum_pair_age * pair_age),
    sigma_abl = exp(b_sigma_abaxialstomatallengthum_Intercept +
      b_sigma_abaxialstomatallengthum_pair_age * pair_age),
    sigma_add = exp(b_sigma_adaxialstomataldensitymm2_Intercept +
      b_sigma_adaxialstomataldensitymm2_pair_age * pair_age),
    sigma_abd = exp(b_sigma_abaxialstomataldensitymm2_Intercept +
      b_sigma_abaxialstomataldensitymm2_pair_age * pair_age),
    
    # Variances
    sigma2_adl = sigma_adl ^ 2,
    sigma2_abl = sigma_abl ^ 2,
    sigma2_add = sigma_add ^ 2,
    sigma2_abd = sigma_abd ^ 2,
    
    # Covariances
    cov_adl_add = r_adl_add * sigma_adl * sigma_add, 
    cov_adl_abl = r_adl_abl * sigma_adl * sigma_abl, 
    cov_add_abl = r_add_abl * sigma_add * sigma_abl, 
    cov_adl_abd = r_adl_abd * sigma_adl * sigma_abd, 
    cov_add_abd = r_add_abd * sigma_add * sigma_abd, 
    cov_abl_abd = r_abl_abd * sigma_abl * sigma_abd
    
  )

write_rds(df_samples12, "objects/df_samples12.rds")
