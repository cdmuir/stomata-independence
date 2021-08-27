source("r/header.R")

distinct_names = stomata_anatomy |>
  select(source_id) |>
  distinct()

# Specify data source(s)
ids = gnr_datasources() |>
  filter(title %in% c("Open Tree of Life Reference Taxonomy",
                      "The International Plant Names Index",
                      "Bishop Museum",
                      "Tropicos - Missouri Botanical Garden",
                      "NCBI",
                      "GRIN Taxonomy for Plants")) |>
  pull(id)

# Divide data into groups so maximum of 100 names are submitted at once
num_groups = ceiling(nrow(distinct_names) / 100)

taxize_output = distinct_names |>
  group_split(group = (row_number() - 1) %/% (n() / num_groups)) |>
  map_dfr(~ {gnr_resolve(sci = .x$source_id, data_source_ids = ids)})

write_rds(taxize_output, "processed-data/taxize_output.rds")
