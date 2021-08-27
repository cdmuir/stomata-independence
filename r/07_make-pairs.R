source("r/header.R")

trimmed_data = read_rds("processed-data/trimmed-data.rds")
trimmed_phy = read_rds("processed-data/trimmed-phylogeny.rds")

set.seed(12340089)
sisters = extract_sisters(multi2di(trimmed_phy), TRUE)
contrasts = extract_contrasts(trimmed_phy) |>
  select(-tree)

# calculate pair_age for contrasts
vcv_phylo = vcv(trimmed_phy)
max_age = max(diag(vcv_phylo))

contrasts = contrasts |>
  rowwise() |>
  mutate(pair_age = max_age - vcv_phylo[sp1, sp2]) |>
  ungroup()

contrasts = bind_rows(sisters, contrasts)

# remove duplicates
contrasts = contrasts |>
  rowwise() |>
  mutate(
    x = str_c(sort(c(sp1, sp2)), collapse = "_"),
    keep = !duplicated(x)
  ) |>
  filter(keep) |>
  select(tree_node, pair_age, sp1, sp2)
 
pair_div = contrasts |>
  pmap_dfr(~ {
    tibble(tree_node = ..1, pair_age = ..2, sp1 = ..3, sp2 = ..4) |>
      bind_cols(trimmed_data |>
                  filter(phy_name == ..3 | phy_name == ..4) |>
                  mutate_if(is.numeric, log10) |>
                  summarize_if(is.numeric, diff))
  })

write_rds(pair_div, "objects/pair_div.rds")
