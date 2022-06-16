# Set of phylogenenetically independent contrasts (PICs)

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
  ungroup() |>
  bind_rows(sisters) |>
  # remove duplicates
  distinct(sp1, sp2, .keep_all = TRUE) |>
  select(tree_node, pair_age, sp1, sp2)
 
pair_div = contrasts |>
  pmap_dfr(~ {
    tibble(tree_node = ..1, pair_age = ..2, sp1 = ..3, sp2 = ..4) |>
      bind_cols(trimmed_data |>
                  filter(phy_name == ..3 | phy_name == ..4) |>
                  mutate(across(where(is.numeric), log10)) |>
                  summarize(across(where(is.numeric), diff)))
  }) |>
  select(-matches("^[a-z]{1}$")) |>
  mutate(pair = as.numeric(as.factor(str_c(sp1, sp2))))

# check that all pairs are unique
assert_true(all(table(pair_div$pair) == 1))

write_rds(pair_div, "objects/pair_div.rds")

# Separate set of PICs to maximize overlap between stomatal anatomy and genome
# size

trimmed_data = read_rds("processed-data/trimmed-data.rds") |>
  filter(!is.na(dna_amount_2c_pg))
trimmed_phy = read_rds("processed-data/trimmed-phylogeny.rds") |>
  keep.tip(trimmed_data$phy_name)

set.seed(760963261)
sisters = extract_sisters(multi2di(trimmed_phy), TRUE)
contrasts = extract_contrasts(trimmed_phy) |>
  select(-tree)

# calculate pair_age for contrasts
vcv_phylo = vcv(trimmed_phy)
max_age = max(diag(vcv_phylo))

contrasts = contrasts |>
  rowwise() |>
  mutate(pair_age = max_age - vcv_phylo[sp1, sp2]) |>
  ungroup() |>
  bind_rows(sisters) |>
  # remove duplicates
  distinct(sp1, sp2, .keep_all = TRUE) |>
  select(tree_node, pair_age, sp1, sp2)

pair_div = contrasts |>
  pmap_dfr(~ {
    tibble(tree_node = ..1, pair_age = ..2, sp1 = ..3, sp2 = ..4) |>
      bind_cols(trimmed_data |>
                  filter(phy_name == ..3 | phy_name == ..4) |>
                  mutate(across(where(is.numeric), log10)) |>
                  summarize(across(where(is.numeric), diff)))
  }) |>
  select(-matches("^[a-z]{1}$")) |>
  mutate(pair = as.numeric(as.factor(str_c(sp1, sp2))))

# check that all pairs are unique
assert_true(all(table(pair_div$pair) == 1))

write_rds(pair_div, "objects/pair_div_genome.rds")
