source("r/header.R")

fit12 = read_rds("objects/fit12.rds")

df = as_draws_df(fit12)

# Average Delta log(trait)
mat_mu = df |>
  select(matches("^b_a[A-Za-z0-9]+_Intercept$")) |>
  pivot_longer(everything()) |>
  mutate(name = case_when(
    name == "b_abaxialstomataldensitymm2_Intercept" ~ "D_ab",
    name == "b_adaxialstomataldensitymm2_Intercept" ~ "D_ad",
    name == "b_abaxialstomatallengthum_Intercept" ~ "L_ab",
    name == "b_adaxialstomatallengthum_Intercept" ~ "L_ad",
  )) |>
  group_by(name) |>
  point_interval(.point = median, .interval = hdi) |>
  select(name, value, .lower, .upper) |>
  mutate(across(where(is.numeric), signif, digits = 2)) |>
  mutate(hdi = glue("$[{.lower},{.upper}]$")) |>
  select(name, value, hdi) |>
  as.matrix() |>
  (function(.x) set_rownames(.x, .x[,"name"]))() 

# Standard deviation of Delta log(trait)
mat_sigma = df |>
  select(matches("^b_sigma_a[A-Za-z0-9]+_Intercept$")) |>
  pivot_longer(everything()) |>
  mutate(name = case_when(
    name == "b_sigma_abaxialstomataldensitymm2_Intercept" ~ "D_ab",
    name == "b_sigma_adaxialstomataldensitymm2_Intercept" ~ "D_ad",
    name == "b_sigma_abaxialstomatallengthum_Intercept" ~ "L_ab",
    name == "b_sigma_adaxialstomatallengthum_Intercept" ~ "L_ad",
  )) |>
  group_by(name) |>
  point_interval(.point = median, .interval = hdi) |>
  select(name, value, .lower, .upper) |>
  mutate(across(where(is.numeric), exp)) |>
  mutate(across(where(is.numeric), signif, digits = 2)) |>
  mutate(hdi = glue("$[{.lower},{.upper}]$")) |>
  select(name, value, hdi) |>
  as.matrix() |>
  (function(.x) set_rownames(.x, .x[,"name"]))() 

# Effect of pair age on standard deviation of Delta log(trait)
mat_age = df |>
  select(matches("^b_sigma_a[A-Za-z0-9]+_pair_age$")) |>
  pivot_longer(everything()) |>
  mutate(name = case_when(
    name == "b_sigma_abaxialstomataldensitymm2_pair_age" ~ "D_ab",
    name == "b_sigma_adaxialstomataldensitymm2_pair_age" ~ "D_ad",
    name == "b_sigma_abaxialstomatallengthum_pair_age" ~ "L_ab",
    name == "b_sigma_adaxialstomatallengthum_pair_age" ~ "L_ad",
  )) |>
  group_by(name) |>
  point_interval(.point = median, .interval = hdi) |>
  select(name, value, .lower, .upper) |>
  mutate(across(where(is.numeric), signif, digits = 2)) |>
  mutate(hdi = glue("$[{.lower},{.upper}]$")) |>
  select(name, value, hdi) |>
  as.matrix() |>
  (function(.x) set_rownames(.x, .x[,"name"]))() 

# Correlation between Delta log(trait)
mat_cor = df |>
  select(matches("^rescor__")) |>
  pivot_longer(everything()) |>
  mutate(name = case_when(
    name == "rescor__adaxialstomataldensitymm2__abaxialstomataldensitymm2" ~ "D_ab-D_ad",
    name == "rescor__abaxialstomatallengthum__abaxialstomataldensitymm2" ~ "D_ab-L_ab",
    name == "rescor__adaxialstomatallengthum__abaxialstomataldensitymm2" ~ "D_ab-L_ad",
    name == "rescor__adaxialstomataldensitymm2__abaxialstomatallengthum" ~ "D_ad-L_ab",
    name == "rescor__adaxialstomatallengthum__adaxialstomataldensitymm2" ~ "D_ad-L_ad",
    name == "rescor__adaxialstomatallengthum__abaxialstomatallengthum" ~ "L_ab-L_ad",
  )) |>
  group_by(name) |>
  point_interval(.point = median, .interval = hdi) |>
  select(name, value, .lower, .upper) |>
  mutate(across(where(is.numeric), signif, digits = 2)) |>
  mutate(hdi = glue("$[{.lower},{.upper}]$")) |>
  select(name, value, hdi) |>
  as.matrix() |>
  (function(.x) set_rownames(.x, .x[,"name"]))() 

# Student t family parameter nu
mat_nu = df |>
  select(nu) |>
  pivot_longer(everything()) |>
  point_interval(value, .point = median, .interval = hdi) |>
  select(value, .lower, .upper) |>
  mutate(across(where(is.numeric), signif, digits = 2)) |>
  mutate(hdi = glue("$[{.lower},{.upper}]$"), name = "nu") |>
  select(name, value, hdi) |>
  as.matrix() |>
  (function(.x) set_rownames(.x, .x[,"name"]))() 

# Combine and write
list(
  mu = mat_mu,
  sigma = mat_sigma,
  age = mat_age,
  cor = mat_cor,
  nu = mat_nu
) |>
 write_rds("objects/h12output.rds")
