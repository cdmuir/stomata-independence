# Plot stomatal density and length data
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
  mutate(
    surface = str_c(surface, " value"),
    trait = case_when(
      trait == "stomatal_density" ~ 'paste("stomatal density [pores m", m^-2, "]")',
      trait == "stomatal_length" ~ 'paste("stomatal length [", mu, "m]")'
    )) |>
  pivot_wider(names_from = "surface", values_from = "value")

my_breaks = function(.x) {
  
  if (max(.x) > 100) {
    return(10 ^ c(0:3))
  } else {
    return(10 ^ c(1:2))
  }
}

my_lims = function(.x) {
  
  if (max(.x) > 100) {
    return(c(0.18, 1500))
  } else {
    return(c(5, 100))
  }
}

ggplot(
  trimmed_data, 
  aes(`abaxial value`, `adaxial value`)
) +
  facet_wrap(vars(trait), scales = "free", labeller = label_parsed) +
  geom_hex() +
  geom_abline(slope = 1, intercept = 0, size = 1.2, linetype = "dashed",
              color = "grey") +
  scale_x_log10(breaks = my_breaks, limits = my_lims, expand=expansion(0,0)) +
  scale_y_log10(breaks = my_breaks, limits = my_lims, expand=expansion(0,0)) +
  scale_fill_gradient(low = "grey20", high = "grey80", name = "Taxa per\nhex bin") +
  theme_cowplot()

ggsave("figures/h2-raw.pdf", width = 6.5, height = 3.25)
ggsave("figures/h2-raw.png", width = 6.5, height = 3.25)
