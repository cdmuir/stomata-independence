source("r/header.R")

# Prepare data ----

## Add absolute trait values for each species
trimmed_dat_sp1 = read_rds("processed-data/trimmed-data.rds") |>
  select(sp1 = phy_name, matches("^a[b|d]axial_fs|gmax")) |>
  rename_with(function(.x) str_c("sp1_", .x), -sp1)

trimmed_dat_sp2 = read_rds("processed-data/trimmed-data.rds") |>
  select(sp2 = phy_name, matches("^a[b|d]axial_fs|gmax")) |>
  rename_with(function(.x) str_c("sp2_", .x), -sp2)

## Pivot pair_div to have column for surface
pair_div = read_rds("objects/pair_div.rds") |>
  pivot_longer(cols = matches("^a[b|d]axial")) |>
  extract(name, c("surface", "trait"), "(a[b|d]axial)_([[:alnum:]_]+)") |>
  mutate(trait = str_c("diff_log10_", trait)) |>
  pivot_wider(names_from = "trait", values_from = "value") |> 
  left_join(trimmed_dat_sp1, by = "sp1") |>
  left_join(trimmed_dat_sp2, by = "sp2") |>
  mutate(
    abaxial_fs = (sp1_abaxial_fs + sp2_abaxial_fs) / 2,
    adaxial_fs = (sp1_adaxial_fs + sp2_adaxial_fs) / 2,
    focal_fs = case_when(
      surface == "abaxial" ~ abaxial_fs,
      surface == "adaxial" ~ adaxial_fs
    )
  )

# Prepare Model ----

## Divergence slows near high surface f_S
bf_3 = bf(mvbind(diff_log10_stomatal_length_um, diff_log10_stomatal_density_mm2) ~
            surface, sigma ~ pair_age * surface + focal_fs:surface)

# Fit model ----
fit3 = brm(
  bf_3 + set_rescor(TRUE),
  data = pair_div,
  family = "student", 
  backend = "cmdstanr", chains = 2L, cores = 2L, 
  seed = 2284290
)

write_rds(fit3, "objects/fit3.rds")
