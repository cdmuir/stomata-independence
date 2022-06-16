source("r/header.R")

resolved_names = read_rds("processed-data/resolved_names.rds")

lookup_table = resolved_names |>
  pull(resolved_name) |>
  unique() |>
  lookup_table(by_species = TRUE) |>
  rownames_to_column("species") |>
  # Removing Gymnosperms because there are so few in the data set
  filter(group == "Angiosperms")

# Check for resupinate taxa based on Chitwood et al. (2012)
lookup_table |>
  filter(
    family == "Alstroemeriaceae" |
      genus == "Luzuriaga" |
      genus == "Pharus" |
      genus == "Stylidium" |
      genus == "Nepenthes"
  ) |>
  nrow() |>
  equals(0) |>
  assert_true()

# Takes ~1 minute
phy = phylo.maker(lookup_table)

write_rds(phy, "processed-data/phy.rds")

# Export look up table of plant families to streamline search of Plant DNA C-values
lookup_table |>
  select(family, genus, species) |>
  arrange(family, genus, species) |>
  write_csv("processed-data/lookup_table.csv")
