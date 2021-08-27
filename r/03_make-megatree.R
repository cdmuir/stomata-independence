source("r/header.R")

resolved_names = read_rds("processed-data/resolved_names.rds")

lookup_table = resolved_names |>
  pull(resolved_name) |>
  unique() |>
  lookup_table(by_species = TRUE) |>
  rownames_to_column("species")

# Takes ~1 minute
phy = phylo.maker(lookup_table)

write_rds(phy, "processed-data/phy.rds")
