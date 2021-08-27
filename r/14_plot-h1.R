source("r/header.R")

# Source saved data and objects ----
pair_dov = read_rds("objects/pair_div.rds")
df_samples = read_rds("objects/df_samples.rds")
h1_pairs = read_rds("objects/h1_pairs.rds")

# Prepare for plotting ----
df_pointinterval = df_samples |>
  filter(pair_age == 0) |>
  select(starts_with("sigma2"), cov_ad = cov_adl_add, cov_ab = cov_abl_abd) |>
  mutate(
    diff_sigma2_l = sigma2_adl - sigma2_abl,
    diff_sigma2_d = sigma2_add - sigma2_abd,
    diff_cov = cov_ad - cov_ab,
  ) |>
  pivot_longer(everything()) |>
  mutate(
    surface = str_extract(name, "a[bd]") |>
      str_c("axial"),
    surface = replace_na(surface, "difference\nadaxial-abaxial"),
    surface = factor(surface, levels = c("abaxial", "adaxial", "difference\nadaxial-abaxial")),
    comparison = case_when(
      str_detect(name, "sigma2_[a]*[bd]*l$") ~ "var(length)",
      str_detect(name, "sigma2_[a]*[bd]*d$") ~ "var(density)",
      str_detect(name, "cov") ~ "cov(length,density)"
    ),
    comparison = factor(comparison, levels = c("var(length)", "var(density)",
                                               "cov(length,density)"))
  )

# Hypothesis testing ----
df_pointinterval |>
  filter(str_detect(name, "^diff")) |>
  select(name, value) |>
  group_by(name) |>
  point_interval(.point = median, .interval = hdi) |>
  write_rds("objects/h1_parameters.rds")

# Plot ----
h1_summary = ggplot(df_pointinterval, aes(surface, value, color = surface)) +
  facet_wrap(. ~ comparison, scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(position = position_dodge(), show.legend = FALSE,
                     point_interval = median_hdi) +
  scale_color_manual(values = c("grey", "grey", "black")) +
  xlab("leaf surface") +
  ylab("(co)variance") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

plot_grid(h1_pairs, h1_summary, nrow = 2, labels = "auto", hjust = -4)

ggsave("figures/h1.pdf", width = 7, height = 7)
ggsave("figures/h1.png", width = 7, height = 7)
