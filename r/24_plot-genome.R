source("r/header.R")

fit_1_genome = read_rds("objects/fit_1_genome.rds")

pair_div_genome = read_rds("objects/pair_div_genome.rds") |>
  select(dna_amount_2c_pg, matches("^a[b|d]axial_stomatal_[density_mm2|length_um]")) |>
  pivot_longer(-dna_amount_2c_pg) |>
  mutate(
    trait = str_extract(name, "density|length"),
    surface = str_extract(name, "a[b|d]axial")
  ) |>
  mutate(surface = factor(surface, levels = c("adaxial", "abaxial")))

# Plot of coefficients and CIs ----
df_coef = fit_1_genome |>
  as_draws_df() |>
  select(
    b_adaxialstomatallengthum_dna_amount_2c_pg,
    b_abaxialstomatallengthum_dna_amount_2c_pg,
    b_adaxialstomataldensitymm2_dna_amount_2c_pg,
    b_abaxialstomataldensitymm2_dna_amount_2c_pg
  ) |>
  pivot_longer(everything()) |>
  mutate(
    comparison = str_extract(name, "density|length"),
    surface = str_extract(name, "a[b|d]axial")
  )

gp1 = ggplot(df_coef, aes(surface, value)) +
  facet_wrap(comparison ~ .) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(show.legend = FALSE, point_interval = median_hdi) +
  ylab(expression(paste("Effect of ", Delta, "log(2C DNA) on ", Delta, "log(trait)"))) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Plot of data ----
df_samples = fit_1_genome |>
  as_draws_df() |>
  select(
    .draw,
    b_adaxialstomatallengthum_Intercept,
    b_abaxialstomatallengthum_Intercept,
    b_adaxialstomataldensitymm2_Intercept,
    b_abaxialstomataldensitymm2_Intercept,
    b_adaxialstomatallengthum_dna_amount_2c_pg,
    b_abaxialstomatallengthum_dna_amount_2c_pg,
    b_adaxialstomataldensitymm2_dna_amount_2c_pg,
    b_abaxialstomataldensitymm2_dna_amount_2c_pg
  ) |>
  pivot_longer(-.draw) |>
  mutate(
    trait = str_extract(name, "density|length"),
    surface = str_extract(name, "a[b|d]axial"),
    coef = str_extract(name, "Intercept|dna_amount_2c_pg") |>
      str_replace("dna_amount_2c_pg", "slope")
  ) |>
  pivot_wider(id_cols = c(".draw", "trait", "surface"), names_from = "coef") |>
  crossing(dna_amount_2c_pg = seq(
    min(pair_div_genome$dna_amount_2c_pg), 
    max(pair_div_genome$dna_amount_2c_pg), 
    length.out = 1e2
  )) |>
  mutate(value = Intercept + slope * dna_amount_2c_pg) |>
  mutate(surface = factor(surface, levels = c("adaxial", "abaxial")))

gp2 = ggplot(pair_div_genome, aes(dna_amount_2c_pg, value)) +
  facet_wrap(surface ~ trait, scales = "free_y") +
  stat_lineribbon(
    data = df_samples, 
    .width = 0.95, point_interval = "median_hdi", fill = "grey"
  ) +
  geom_point() +
  xlab(expression(paste(Delta, "log(2C DNA)"))) +
  ylab(expression(paste(Delta, "log(trait)"))) +
  theme_cowplot()

plot_grid(gp2, gp1, nrow = 1, labels = "auto", rel_widths = c(0.6, 0.4))

ggsave("figures/genome.pdf", width = 7, height = 5)
ggsave("figures/genome.png", width = 7, height = 5)
