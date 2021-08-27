source("r/header.R")

# Source saved data and objects ----
pair_div = read_rds("objects/pair_div.rds")
df_samples = read_rds("objects/df_samples.rds")

# Prepare for plotting ----
df_density = df_samples |>
  filter(pair_age == 0) |>
  select(
    b_abaxialstomataldensitymm2_Intercept,
    b_adaxialstomataldensitymm2_Intercept,
    sigma_add, sigma_abd, r_add_abd
  ) |>
  mutate(
    sma_slope = sign(r_add_abd) * sigma_abd / sigma_add,
    sma_intercept =  b_adaxialstomataldensitymm2_Intercept - 
      b_adaxialstomataldensitymm2_Intercept * sma_slope
  ) 

df_length = df_samples |>
  filter(pair_age == 0) |>
  select(
    b_abaxialstomatallengthum_Intercept,
    b_adaxialstomatallengthum_Intercept,
    sigma_adl, sigma_abl, r_adl_abl
  ) |>
  mutate(
    sma_slope = sign(r_adl_abl) * sigma_abl / sigma_adl,
    sma_intercept =  b_adaxialstomatallengthum_Intercept - 
      b_abaxialstomatallengthum_Intercept * sma_slope
  ) 

write_rds(df_density, "objects/df_density.rds")
write_rds(df_length, "objects/df_length.rds")

# Tibbles for confidence ribbons ----
df_density_ribbon = df_density |>
  crossing(
    x = seq(min(pair_div$abaxial_stomatal_density_mm2), 
            max(pair_div$abaxial_stomatal_density_mm2), 
            length.out = 50)
  ) %>%
  mutate(y = sma_intercept + sma_slope * x) %>%
  group_by(x) %>%
  point_interval(y, .point = median, .interval = hdi)

df_length_ribbon = df_length %>%
  crossing(
    x = seq(min(pair_div$abaxial_stomatal_length_um), 
            max(pair_div$abaxial_stomatal_length_um), 
            length.out = 50)
  ) %>%
  mutate(y = sma_intercept + sma_slope * x) %>%
  group_by(x) %>%
  point_interval(y, .point = median, .interval = hdi)

# Make figures ----
gp_density = ggplot(pair_div, aes(abaxial_stomatal_density_mm2, 
                           adaxial_stomatal_density_mm2)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, size = 1.2, linetype = "dashed",
              color = "grey") +
  geom_ribbon(
    data = df_density_ribbon,
    mapping = aes(x = x, y = y, ymin = .lower, ymax = .upper),
    inherit.aes = FALSE, fill = "grey", alpha = 0.5
    # inherit.aes = FALSE, fill = "tomato", alpha = 0.5
  ) +
  geom_line(
    data = df_density_ribbon,
    mapping = aes(x = x, y = y),
    inherit.aes = FALSE, color = "darkgrey", size = 1.2
    # inherit.aes = FALSE, color = "tomato", size = 1.2
  ) +
  xlab(expression(paste(Delta, "abaxial"))) +
  ylab(expression(paste(Delta, "adaxial"))) +
  ggtitle(expression(paste("log(stomtal density [pores m", m^-2, "])"))) +
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_continuous(limits = c(-2, 2)) +
  coord_equal() +
  theme_cowplot()

gp_length = ggplot(pair_div, aes(abaxial_stomatal_length_um,
                           adaxial_stomatal_length_um)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, size = 1.2, linetype = "dashed",
              color = "grey") +
  geom_ribbon(
    data = df_length_ribbon,
    mapping = aes(x = x, y = y, ymin = .lower, ymax = .upper),
    inherit.aes = FALSE, fill = "grey", alpha = 0.5
    # inherit.aes = FALSE, fill = "steelblue", alpha = 0.5
  ) +
  geom_line(
    data = df_length_ribbon,
    mapping = aes(x = x, y = y),
    inherit.aes = FALSE, color = "darkgrey", size = 1.2
    # inherit.aes = FALSE, color = "steelblue", size = 1.2
  ) +
  xlab(expression(paste(Delta, "abaxial"))) +
  ylab(expression(paste(Delta, "adaxial"))) +
  ggtitle(expression(paste("log(stomtal length [", mu, "m])"))) +
  scale_x_continuous(limits = c(-0.5, 0.5)) +
  scale_y_continuous(limits = c(-0.5, 0.5)) +
  coord_equal() +
  theme_cowplot()

h2_pairs = plot_grid(gp_density, gp_length, align = "hv")

write_rds(h2_pairs, "objects/h2_pairs.rds")
