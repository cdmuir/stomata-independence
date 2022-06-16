source("r/header.R")

# Read samples ----
fit3 = read_rds("objects/fit3.rds")
df_samples3 = read_rds("objects/df_samples3.rds")

# Summary for main figure ----
df_main = df_samples3 |>
  select(.draw, surface, focal_fs, starts_with("mu"), starts_with("sigma"), nu) |>
  mutate(nu_density = nu, nu_length = nu) |>
  select(-nu) |>
  pivot_longer(cols = c(starts_with("mu"), starts_with("sigma"), starts_with("nu"))) |>
  separate(name, c("parameter", "trait"), sep = "_", remove = TRUE) |>
  pivot_wider(names_from = "parameter") |>
  mutate(
    ymin = brms::qstudent_t(0.025, nu, mu, sigma),
    ymax = brms::qstudent_t(0.975, nu, mu, sigma)
  ) |>
  group_by(surface, focal_fs, trait) |>
  point_interval(.width = 0.95, .point = median, .interval = qi) |>
  mutate(
    surface = factor(surface, levels = c("adaxial", "abaxial")),
    trait1 = case_when(
      trait == "density" ~ 'paste("stomatal density [pores m", m^-2, "]")',
      trait == "length" ~ 'paste("stomatal length [", mu, "m]")'
    )
  )

# Main figure: f_s versus delta(trait) ----

ggplot(df_main, aes(focal_fs, abs(mu), ymin = 0, ymax = ymax)) +
  facet_grid(
    surface ~ trait1,
    labeller = label_parsed
    ) +
  geom_point(
    data = fit3$data |>
      pivot_longer(starts_with("diff_"), values_to = "mu") |>
      mutate(
        trait = str_extract(name, "density|length"), 
        ymin = numeric(1), ymax =  numeric(1),
        surface = factor(surface, levels = c("adaxial", "abaxial")),
        trait1 = case_when(
          trait == "density" ~ 'paste("stomatal density [pores m", m^-2, "]")',
          trait == "length" ~ 'paste("stomatal length [", mu, "m]")'
        )
      ), alpha = 0.5
  ) +
  geom_ribbon(alpha = 0.5) +
  # geom_line() +
  xlab("Fraction of epidermal area allocated to stomata per surface") +
  ylab(expression(paste("|", Delta, "log(trait)|"))) +
  theme_cowplot()

ggsave("figures/h3.pdf", width = 7, height = 7)
ggsave("figures/h3.png", width = 7, height = 7)

# Summary for supporting figure ----

df_supporting = df_samples3 |>
  select(.draw, surface, focal_fs, starts_with("sigma")) |>
  pivot_longer(cols = starts_with("sigma")) |>
  group_by(surface, focal_fs, name) |>
  point_interval(.width = 0.95, .point = median, .interval = hdi) |>
  mutate(
    surface = factor(surface, levels = c("adaxial", "abaxial")),
    trait1 = case_when(
      name == "sigma_density" ~ 'paste("stomatal density [pores m", m^-2, "]")',
      name == "sigma_length" ~ 'paste("stomatal length [", mu, "m]")'
    )
  )
  
# Supporting figure: f_s versus sigma ----

ggplot(df_supporting, aes(focal_fs, value, ymin = .lower, ymax = .upper, 
                          fill = surface, linetype = surface)) +
  facet_wrap(. ~ trait1, label = label_parsed) +
  geom_ribbon(alpha = 0.5, color = "black") +
  scale_fill_manual(values = c("black", "white")) +
  geom_line() +
  xlab("Fraction of epidermal area allocated to stomata per surface") +
  ylab(expression(sqrt(paste("Var[", Delta, "log(trait)]")))) +
  theme_cowplot()

ggsave("figures/fs-sigma.pdf", width = 7, height = 3.5)
ggsave("figures/fs-sigma.png", width = 7, height = 3.5)
