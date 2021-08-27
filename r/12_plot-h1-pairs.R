source("r/header.R")

# Source saved data and objects ----
pair_div = read_rds("objects/pair_div.rds")
df_samples = read_rds("objects/df_samples.rds")

# Prepare ----
df_medians = df_samples |>
  filter(pair_age == 0) |>
  select(starts_with("b_a"), starts_with("sigma2"), starts_with("cov"), nu) |>
  summarize(across(everything(), median))

# Abaxial panel ----
center_ab = df_medians |>
  select(b_abaxialstomatallengthum_Intercept,
         b_abaxialstomataldensitymm2_Intercept) |>
  as.numeric()

shape_ab = df_medians |>
  select(sigma2_abl, cov_abl_abd, sigma2_abd) |>
  as.numeric() |>
  rep(c(1, 2, 1)) |>
  matrix(ncol = 2)

df_ab_ellipse = my_confidence_ellipse(center_ab, shape_ab, df_medians$nu, 0.95, 100) |>
  as_tibble(.name_repair = "unique") |>
  rename(abaxial_stomatal_length_um = ...1, abaxial_stomatal_density_mm2 = ...2)

gp_ab = ggplot(
  pair_div, 
  aes(abaxial_stomatal_length_um, abaxial_stomatal_density_mm2)
) +
  geom_polygon(data = df_ab_ellipse, fill = "grey", color = "darkgrey") +
  geom_point() +
  scale_x_continuous(limits = c(-0.5, 0.5)) + 
  scale_y_continuous(limits = c(-2, 2)) +
  xlab(expression(paste(Delta, "log(stomatal length [", mu, "m])"))) +
  ylab(expression(paste(Delta, "log(stomatal density [pores m", m^-2, "])"))) +
  ggtitle("abaxial surface") +
  theme_cowplot()

# Adaxial panel ----
center_ad = df_medians |>
  select(b_adaxialstomatallengthum_Intercept,
         b_adaxialstomataldensitymm2_Intercept) |>
  as.numeric()

shape_ad = df_medians |>
  select(sigma2_adl, cov_adl_add, sigma2_add) |>
  as.numeric() |>
  rep(c(1, 2, 1)) |>
  matrix(ncol = 2)

df_ad_ellipse = my_confidence_ellipse(center_ad, shape_ad, df_medians$nu, 0.95, 100) |>
  as_tibble(.name_repair = "unique") |>
  rename(adaxial_stomatal_length_um = ...1, adaxial_stomatal_density_mm2 = ...2)

gp_ad = ggplot(
  pair_div, 
  aes(adaxial_stomatal_length_um, adaxial_stomatal_density_mm2)
) +
  geom_polygon(data = df_ad_ellipse, fill = "grey", color = "darkgrey") +
  geom_point() +
  scale_x_continuous(limits = c(-0.5, 0.5)) + 
  scale_y_continuous(limits = c(-2, 2)) +
  xlab(expression(paste(Delta, "log(stomatal length [", mu, "m])"))) +
  ylab(expression(paste(Delta, "log(stomatal density [pores m", m^-2, "])"))) +
  ggtitle("adaxial surface") +
  theme_cowplot()

h1_pairs = plot_grid(gp_ab, gp_ad)

write_rds(h1_pairs, "objects/h1_pairs.rds")
