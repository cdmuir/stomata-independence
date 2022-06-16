source("r/header.R")

phy = read_rds("processed-data/phy.rds")$scenario.3

# Find polytomies ----
polytomies = find_polytomies(phy, 3)

# Note on polytomy at node 1047:
# node: 1047, Zinnia, Verbesina, Pittocaulon, Montanoa, Jaegeria, Ageratina, ...
#
# tl;dr version: I removed 5 taxa that wouldn't have added any contrasts in
# subsequent analysis. Hence, it doesn't affect the current analysis but I 
# should reevaluate if new data are added.
#
# Visual inspection shows the problematic taxa are:
# Jaegeria hirta, Montanoa grandiflora, Pittocaulon praecox, Verbesina virgata, 
# and Zinnia peruviana
#
# I attempted to resolve this polytomy using PyPHLAWD but there were no useful
# clusters that included members of these genera
#
# Next, I attempted to resolve using Open Tree of Life. The OTL tree resolved
# some of this polytomy but in other ways in incongruent with the current tree.
#
# Visual inspection showed that they wouldn't have added any additional 
# contrasts because they were either unresolved or were taxon c in a triplet 
# like: ((a,b),c);
#
# Below is the code download and view OTL v13.4 tree with rotl version 3.0.11
# subtree = extract.clade(phy, 1047)
# taxa = subtree$tip.label |>
#   str_replace_all("_", " ") |>
#   unique() |>
#   tnrs_match_names(context = "Asteraceae")
# 
# taxon_map = structure(taxa$search_string, names = taxa$unique_name)
# 
# x = taxa |>
#   ott_id() |>
#   is_in_tree()
# tr = tol_induced_subtree(ott_id(taxa)[x])
# write_rds(tr, "objects/rotl_aster_tree.rds")

# Similar issue with two taxa in Apocynaceae, therefore I dropped:
# "Ervatamia_divaricata"  "Cascabela_thevetia"
# It shouldn't affect current analyses, but reevaluate if more Apocynaceae are
# added
# subtree = extract.clade(phy, 1290)
# taxa = subtree$tip.label |>
#   str_replace_all("_", " ") |>
#   unique() |>
#   tnrs_match_names(context = "Asterids")
# 
# taxon_map = structure(taxa$search_string, names = taxa$unique_name)
# 
# x = taxa |>
#   ott_id() |>
#   is_in_tree()
# tr = tol_induced_subtree(ott_id(taxa)[x])
# write_rds(tr, "objects/rotl_apocynac_tree.rds")

# Remove Aster polytomy and re-run polytomies ----
phy1 = drop.tip(phy, c("Jaegeria_hirta", "Montanoa_grandiflora",
                      "Pittocaulon_praecox", "Verbesina_virgata",
                      "Zinnia_peruviana", "Ervatamia_divaricata",
                      "Cascabela_thevetia"))

write_rds(phy1, "processed-data/phy1.rds")
rm(phy)

polytomies = find_polytomies(phy1, 3)
write_rds(polytomies, "objects/polytomies.rds")

# For these clades, I used pyPHLAWD to generate multi-locus alignments
# see 'stomata-independence-desktop for scripts

# Prepare RAxML for clades with lots of polytomies ----
file.create("raxml/run-raxml.sh")

# May need to run taxize::use_entrez() to get and set API key
# For some reason, NCBI sometimes fails and command needs to be executed again
# to work.
Lotus = prepare_raxml(filter(polytomies, genus == "Lotus")$node, phy1, TRUE)

Euphorbia = prepare_raxml(
  filter(polytomies, genus == "Euphorbia")$node, phy1, TRUE
)

Astragalus = prepare_raxml(
  filter(polytomies, genus == "Astragalus")$node, phy1, TRUE
)

Mentha = prepare_raxml(filter(polytomies, genus == "Mentha")$node, phy1, TRUE)

Eucalyptus = prepare_raxml(
  filter(polytomies, genus == "Eucalyptus")$node, phy1, TRUE
)

Stipa = prepare_raxml(filter(polytomies, genus == "Stipa")$node, phy1, TRUE)

# no pyPHLAWD clusters
# Kobresia = prepare_raxml(filter(polytomies, genus == "Kobresia")$node, phy1, TRUE) 

Oryza = prepare_raxml(filter(polytomies, genus == "Oryza")$node, phy1, TRUE)

Scabiosa = prepare_raxml(filter(polytomies, genus == "Scabiosa")$node, phy1, TRUE)

Delphinium = prepare_raxml(
  filter(polytomies, genus == "Delphinium")$node, phy1, TRUE
)

Centaurea = prepare_raxml(
  filter(polytomies, genus == "Centaurea")$node, phy1, TRUE
)

Polygonum = prepare_raxml(
  filter(polytomies, genus == "Polygonum, Atraphaxis")$node, phy1, TRUE, 
  clade = "Polygonum"
)

Oxytropis = prepare_raxml(
  filter(polytomies, genus == "Oxytropis")$node, phy1, TRUE
)

Potentilla = prepare_raxml(
  filter(polytomies, genus == "Potentilla")$node, phy1, TRUE
)

Alchemilla  = prepare_raxml(
  filter(polytomies, genus == "Alchemilla")$node, phy1, TRUE
)

Echinochloa = prepare_raxml(
  filter(polytomies, genus == "Echinochloa")$node, phy1, TRUE
)

Dubautia = prepare_raxml(
  filter(polytomies, genus == "Dubautia")$node, phy1, TRUE
)

Solidago = prepare_raxml(
  filter(polytomies, genus == "Solidago")$node, phy1, TRUE
)

Artemisia = prepare_raxml(
  filter(polytomies, genus == "Artemisia")$node, phy1, TRUE
)

Salvia = prepare_raxml(filter(polytomies, genus == "Salvia")$node, phy1, TRUE)

Gentiana = prepare_raxml(
  filter(polytomies, genus == "Gentiana")$node, phy1, TRUE
)

Trigonella = prepare_raxml(
  filter(polytomies, genus == "Trigonella")$node, phy1, TRUE
)

Vicia = prepare_raxml(filter(polytomies, genus == "Vicia")$node, phy1, TRUE)

Onobrychis = prepare_raxml(
  filter(polytomies, genus == "Onobrychis")$node, phy1, TRUE
)

Salix = prepare_raxml(filter(polytomies, genus == "Salix")$node, phy1, TRUE)

Ranunculus = prepare_raxml(
  filter(polytomies, genus == "Ranunculus")$node, phy1, TRUE
)

Thalictrum = prepare_raxml(
  filter(polytomies, genus == "Thalictrum")$node, phy1, TRUE
)

Bromus = prepare_raxml(filter(polytomies, genus == "Bromus")$node, phy1, TRUE)                                                      
# I obtained a lot of clusters for Carex and Solanum during initial June 2020
# analysis. However, I encountered errors in August 2021 when I reran with 
# software and database updates. However, since the number of sequence clusters
# were not limiting for these clades, I am using 2020 data rather than debugging
Carex = prepare_raxml(filter(polytomies, genus == "Carex")$node, phy1, FALSE)

Solanum = prepare_raxml(
  filter(polytomies, genus == "Solanum, Lycopersicon")$node, phy1, FALSE, 
  clade = "Solanum"
)

# Now run 'raxml/run-raxml.sh' to make phylogenies