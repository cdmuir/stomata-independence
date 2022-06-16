source("r/header.R")

fit3 = read_rds("objects/fit3.rds")

df = as_draws_df(fit3)

# Effect of f_S on standard deviation of Delta log(trait)
mat_fs = df |>
  select(matches("^b_sigma_difflog10stomatal(densitymm2|lengthum)_surfacea(b|d)axial:focal_fs")) |>
  pivot_longer(everything()) |>
  mutate(name = case_when(
    name == "b_sigma_difflog10stomatallengthum_surfaceabaxial:focal_fs" ~ "L_ab",
    name == "b_sigma_difflog10stomatallengthum_surfaceadaxial:focal_fs" ~ "L_ad",
    name == "b_sigma_difflog10stomataldensitymm2_surfaceabaxial:focal_fs" ~ "D_ab",
    name == "b_sigma_difflog10stomataldensitymm2_surfaceadaxial:focal_fs" ~ "D_ad",
  )) |>
  group_by(name) |>
  point_interval(.point = median, .interval = hdi) |>
  select(name, value, .lower, .upper) |>
  mutate(across(where(is.numeric), signif, digits = 2)) |>
  mutate(hdi = glue("$[{.lower},{.upper}]$")) |>
  select(name, value, hdi) |>
  as.matrix() |>
  (function(.x) set_rownames(.x, .x[,"name"]))() 

# Write
list(fs = mat_fs) |>
 write_rds("objects/h3output.rds")
