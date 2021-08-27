source("r/header.R")

set.seed(149288644)
phy1 = phy2 = read_rds("processed-data/phy1.rds")
polytomies1 = polytomies2 = read_rds("objects/polytomies.rds") |> 
  filter(!(genus %in% c("Kobresia"))) |>
  mutate(
    seed = round(runif(1, 0, 1e9)),
    # When multiple clades, pull first genus because that is how files were
    # saved
    clade = str_extract(genus, "^[A-Z][a-z]+")
  )

for (i in seq_len(nrow(polytomies1))) {

  clade1 = polytomies1$clade[i]
  seed1 = polytomies1$seed[i]
  
  node = polytomies2 |>
    filter(clade == clade1) |>
    pull(node)
  
  phy2 = replace_polytomy(phy2, node, clade1, seed1)
  polytomies2 = find_polytomies(phy2, n_tip = 3) |> 
    mutate(clade = str_extract(genus, "^[A-Z][a-z]+"))
  
}

# Randomly remove all but one Kobresia
set.seed(23789534)
phy2 = extract.clade(phy2, polytomies2$node[polytomies2$clade == "Kobresia"]) |>
  use_series(tip.label) |>
  (function(.x) sample(.x, length(.x) - 1))() |>
  (function(.x) drop.tip(phy2, .x))()
  
find_polytomies(phy2, 2) |>
  filter(n > 3) |>
  nrow() |>
  equals(0) |>
  assert_true() # check that major polytomies are resolved

# Resolve 3-way polytomies randomly. If these are internal branches, it won't
# matter. If they are 3 tips, it is identical to randomly dropping one tip
set.seed(220750989)
phy2 = multi2di(phy2)

# Check that tree is approximately ultrametric (within numeric rounding error)
h = diag(vcv(phy2))
d = max(h) - h
assert_true(max(d) < 1e-5)
phy2 = force.ultrametric(phy2, method = "extend")
assert_true(is.ultrametric(phy2))

write_rds(phy2, "processed-data/phy2.rds")
