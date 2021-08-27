# Plot size-density scaling on each surface
source("r/header.R")

trimmed_data = read_rds("processed-data/trimmed-data.rds") |>
  select(phy_name, starts_with("a")) |>
  pivot_longer(-phy_name) |>
  transmute(
    phy_name = phy_name,
    surface = str_extract(name, "^a[bd]axial"),
    trait = str_extract(name, "stomatal_[a-z]+"),
    value = value
  ) |> 
  pivot_wider(names_from = "trait", values_from = "value")

ggplot(
  trimmed_data, 
  aes(stomatal_length, stomatal_density)
) +
  facet_grid(. ~ surface) +
  geom_hex() +
  scale_x_log10(limits = c(5, 100)) + 
  scale_y_log10(limits = c(0.18, 1500)) +
  scale_fill_gradient(low = "grey20", high = "grey80") +
  xlab(expression(paste("Stomatal length [", mu, "m]"))) +
  ylab(expression(paste("Stomatal density [pores m", m^-2, "]"))) +
  theme_cowplot()

ggsave("figures/h1-raw.pdf", width = 6.5, height = 3.25)
ggsave("figures/h1-raw.png", width = 6.5, height = 3.25)
