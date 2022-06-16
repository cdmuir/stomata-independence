source("r/header.R")

# Source fitted model ----
fit3 = read_rds("objects/fit3.rds")

# Prepare for plotting ----

df_draws3 = fit3 |>
  as_draws_df(pars = c(
    
    # Trait averages
    "b_difflog10stomatallengthum_Intercept",
    "b_difflog10stomataldensitymm2_Intercept",
    "b_difflog10stomatallengthum_surfaceadaxial",
    "b_difflog10stomataldensitymm2_surfaceadaxial",

    # Trait (co)variances
    "b_sigma_difflog10stomatallengthum_Intercept",
    "b_sigma_difflog10stomataldensitymm2_Intercept",
    "b_sigma_difflog10stomatallengthum_surfaceadaxial",
    "b_sigma_difflog10stomataldensitymm2_surfaceadaxial",

    # Effect pair_age on trait (co)variances
    "b_sigma_difflog10stomatallengthum_pair_age",
    "b_sigma_difflog10stomataldensitymm2_pair_age",
    "b_sigma_difflog10stomatallengthum_pair_age:surfaceadaxial",
    "b_sigma_difflog10stomataldensitymm2_pair_age:surfaceadaxial",
    
    # Effect focal_fs on trait (co)variances
    "b_sigma_difflog10stomatallengthum_surfaceabaxial:focal_fs",
    "b_sigma_difflog10stomataldensitymm2_surfaceabaxial:focal_fs",
    "b_sigma_difflog10stomatallengthum_surfaceadaxial:focal_fs",
    "b_sigma_difflog10stomataldensitymm2_surfaceadaxial:focal_fs",
    
    # Trait correlation
    "rescor__difflog10stomatallengthum__difflog10stomataldensitymm2",

    # nu parameter
    "nu"
    
  )) |>
  rename(
    r = rescor__difflog10stomatallengthum__difflog10stomataldensitymm2,
  ) 

df_samples3 = df_draws3 |>
  crossing(
    pair_age = 0,
    nesting(
      surface = c("abaxial", "adaxial"),
      surface01 = c(0, 1)
    ),
    focal_fs = seq(0, max(fit3$data$focal_fs), length.out = 25)
  ) |>
  mutate(
    
    # Averages
    mu_density = b_difflog10stomataldensitymm2_Intercept +
      b_difflog10stomataldensitymm2_surfaceadaxial * surface01,
    
    mu_length = b_difflog10stomatallengthum_Intercept +
      b_difflog10stomatallengthum_surfaceadaxial * surface01,
    
    # Standard deviations
    sigma_density = exp(
      b_sigma_difflog10stomataldensitymm2_Intercept +
        b_sigma_difflog10stomataldensitymm2_surfaceadaxial * surface01 +
        pair_age * (
          b_sigma_difflog10stomataldensitymm2_pair_age + 
            `b_sigma_difflog10stomataldensitymm2_pair_age:surfaceadaxial` * surface01
        ) +
        focal_fs * (surface == "abaxial") * 
        `b_sigma_difflog10stomataldensitymm2_surfaceabaxial:focal_fs` + 
        focal_fs * (surface == "adaxial") * 
        `b_sigma_difflog10stomataldensitymm2_surfaceadaxial:focal_fs`
    ),
    
    sigma_length = exp(
      b_sigma_difflog10stomatallengthum_Intercept +
        b_sigma_difflog10stomatallengthum_surfaceadaxial * surface01 +
        pair_age * (
          b_sigma_difflog10stomatallengthum_pair_age + 
            `b_sigma_difflog10stomatallengthum_pair_age:surfaceadaxial` * surface01
        ) +
        focal_fs * (surface == "abaxial") * 
        `b_sigma_difflog10stomatallengthum_surfaceabaxial:focal_fs` + 
        focal_fs * (surface == "adaxial") * 
        `b_sigma_difflog10stomatallengthum_surfaceadaxial:focal_fs`
      
    )
  )

# 95% HDI intervals for difference in coefficient of how focal_fs affects 
# variance in divergence
df_draws3 |>
  mutate(
    diff_length = `b_sigma_difflog10stomatallengthum_surfaceadaxial:focal_fs` - `b_sigma_difflog10stomatallengthum_surfaceabaxial:focal_fs`,
    diff_density = `b_sigma_difflog10stomataldensitymm2_surfaceadaxial:focal_fs` - `b_sigma_difflog10stomataldensitymm2_surfaceabaxial:focal_fs`
  ) |>
  select(.draw,diff_length, diff_density) |>
  pivot_longer(cols = starts_with("diff")) |>
  group_by(name) |>
  point_interval(.width = 0.95, .point = median, .interval = hdi) |>
  write_rds("objects/h3_hpd.rds")
  
write_rds(df_samples3, "objects/df_samples3.rds")
