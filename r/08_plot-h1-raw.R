# Plot size-density scaling on each surface
source("r/header.R")

trimmed_data = read_rds("processed-data/trimmed-data.rds") |>
  select(phy_name, matches("^a[b|d]axial_stomatal_density|length")) |>
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
  scale_fill_gradient(low = "grey20", high = "grey80", name = "Taxa per\nhex bin") +
  xlab(expression(paste("Stomatal length [", mu, "m]"))) +
  ylab(expression(paste("Stomatal density [pores m", m^-2, "]"))) +
  coord_fixed(ratio =  log10(100 / 5) / log10(1500 / 0.18)) +
  theme_cowplot()

ggsave("figures/h1-raw.pdf", width = 6.5, height = 3.25)
ggsave("figures/h1-raw.png", width = 6.5, height = 3.25)
