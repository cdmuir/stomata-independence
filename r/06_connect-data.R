# Connect trait and phylogenetic data

source("r/header.R")

resolved_names = read_rds("processed-data/resolved_names.rds") |>
  rename(source_id = user_supplied_name) |>
  mutate(resolved_name = str_replace_all(resolved_name, " ", "_"))

# this tree is full resolved and ultrametric
phy2 = read_rds("processed-data/phy2.rds") 

# Some tip labels are not in resolved_name. This is because either multiple 
# varieties corresponded to the same tip name or the tip name didn't have the 
# variety. The code below harmonizes names by randomly removing varieties 
# without a resolved name, or reconciling names
tibble(
  phy2_tip = phy2$tip.label,
  in_resolved_names = phy2_tip %in% resolved_names$resolved_name
) |> 
  filter(!in_resolved_names)

resolved_names = resolved_names |>
  # Remove resolved names no longer in use
  filter(!(resolved_name %in% c(
    "Euphorbia_milii_var._hislopii",
    "Euphorbia_celastroides_var._kaenana",
    "Euphorbia_celastroides_var._lorifolia",
    "Euphorbia_celastroides_var._stokesii",
    "Euphorbia_celastroides_var._stokesii",  
    "Euphorbia_celastroides_var._tomentella"
  ))) |>
  mutate(resolved_name = case_when(
    resolved_name == "Euphorbia_milii_var._splendens" ~ "Euphorbia_milii",
    resolved_name == "Euphorbia_celastroides_var._amplectens" ~ "Euphorbia_celastroides",
    resolved_name == "Euphorbia_multiformis_var._microphylla" ~ "Euphorbia_multiformis",
    resolved_name == "Solanum_lycopersicum_var._cerasiforme" ~ "Solanum_lycopersicum",
    TRUE ~ resolved_name
  ))

# Check that no more tips are missing from resolved names
tibble(
  phy2_tip = phy2$tip.label,
  in_resolved_names = phy2_tip %in% resolved_names$resolved_name
) |> 
  filter(!in_resolved_names) |>
  nrow() |>
  equals(0) |>
  assert_true()

# If study reports both pore and guard cell length, use guard cell length for
# greater consistency with other studies in the data set
ex = ropenstomata::stomata_anatomy |>
  filter(str_detect(trait, "length")) |>
  group_by(source, source_id, trait) |>
  mutate(n = n()) |>
  select(source, source_id, trait, n) |>
  pivot_wider(names_from = "trait", values_from = "n") |>
  filter(across(starts_with("a"), ~!is.na(.x))) |>
  select(source, source_id, contains("pore")) |>
  pivot_longer(contains("pore"), names_to = "trait") |>
  select(-value)

trimmed_data = ropenstomata::stomata_anatomy |>
  anti_join(ex, by = c("source_id", "trait", "source")) |>
  left_join(resolved_names, by = "source_id") |>
  filter(!is.na(resolved_name))

# Three species (Lycopersicum esculentum, Abromeitiella brevifolia, and 
# Symphyandra armena) did not make it into Some species did not make it into tree.
# I have data on Lycopersicum esculentum (Solanum lycopersicum) for anothr 
# study, so this does not matter. I am not sure why Abromeitiella brevifolia and 
# Symphyandra armena were excluded, but I am not going to worry since it is only
# two species.

trimmed_data = trimmed_data |>
  filter(resolved_name %in% phy2$tip.label) |>
  select(source, source_id, trait, mu) |>
  mutate(trait = str_replace(trait, "guard_cell", "stomatal")) |>
  filter(!str_detect(trait, "width")) |>
  
  # Assume stomatal_length = 2 * pore_length (Sack and Buckley 2016)
  mutate(
    mu = ifelse(str_detect(trait, "pore_length"), 2 * mu, mu),
    trait = str_replace(trait, "pore_length", "stomatal_length")
  ) |>
  
  pivot_wider(names_from = "trait", values_from = "mu") |>
  filter_at(vars(starts_with("a")), ~ !is.na(.x)) |>
  left_join(resolved_names, by = "source_id") |>
  mutate(phy_name = str_replace_all(resolved_name, " ", "_"))

# Remove cases where different sources measured the same species by choosing
# most recent observation or at random if they are same year. This occurred
# because sometimes two source_ids resolved to same resolved_name, for example
trimmed_data = trimmed_data |>
  mutate(source_year = str_extract(source, "[1-2][0-9]{3}")) |>
  group_by(resolved_name) |>
  summarize(
    # n = n(), 
    # sources = str_c(source, collapse = ","),
    i = last(order(source_year)), 
    source = source[i], 
    source_id = source_id[i],             
    abaxial_stomatal_density_mm2 = abaxial_stomatal_density_mm2[i],
    adaxial_stomatal_density_mm2 = adaxial_stomatal_density_mm2[i],
    abaxial_stomatal_length_um = abaxial_stomatal_length_um[i],
    adaxial_stomatal_length_um = adaxial_stomatal_length_um[i],
    resolved_name = resolved_name[i],
    phy_name = phy_name[i],
    .groups = "drop"
  ) |>
  select(-i)
  
# should be no duplicated resolved_names
trimmed_data |>
  group_by(resolved_name) |>
  summarize(n = n(), .groups = "drop") |>
  filter(n > 1L) |>
  nrow() |>
  equals(0) |>
  assert_true()

# Export trimmed dataset and phylogeny ----

trimmed_data |>
  ungroup() |>
  write_rds("processed-data/trimmed-data.rds")

phy2 |>
  keep.tip(trimmed_data$phy_name) |>
  write_rds("processed-data/trimmed-phylogeny.rds")
