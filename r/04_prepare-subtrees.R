source("r/header.R")

phy = read_rds("processed-data/phy.rds")$scenario.3

# Find polytomies ----
polytomies = find_polytomies(phy, 3)

# Note on polytomy at node 1181:
# node: 1181, Zinnia, Verbesina, Pittocaulon, Montanoa, Jaegeria, Ageratina, ...
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
# Visual inspection showed that they wouldn't have added any additional contrats
# because they were either unresolved or were taxon c in a triplet like:
# ((a,b),c);
#
# Below is the code download and view OTL v13.4 tree with rotl version 3.0.11
# subtree = extract.clade(phy, 1181)
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
# subtree = extract.clade(phy, 1437)
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

Lotus       = prepare_raxml(1613, phy1, TRUE)
Euphorbia   = prepare_raxml(1739, phy1, TRUE)
Astragalus  = prepare_raxml(1600, phy1, TRUE)
Mentha      = prepare_raxml(1363, phy1, TRUE)
Eucalyptus  = prepare_raxml(1831, phy1, TRUE)
Stipa       = prepare_raxml(2037, phy1, TRUE)
# Kobresia    = prepare_raxml(2060, phy1, TRUE) # no pyPHLAWD clusters
Oryza       = prepare_raxml(2040, phy1, TRUE)
Scabiosa    = prepare_raxml(1343, phy1, TRUE)
Delphinium  = prepare_raxml(1894, phy1, TRUE)
Centaurea   = prepare_raxml(1281, phy1, TRUE)
Polygonum   = prepare_raxml(1536, phy1, TRUE, clade = "Polygonum")
Oxytropis   = prepare_raxml(1603, phy1, TRUE)
Potentilla  = prepare_raxml(1653, phy1, TRUE)
Alchemilla  = prepare_raxml(1659, phy1, TRUE)
Echinochloa = prepare_raxml(1961, phy1, TRUE)
Dubautia    = prepare_raxml(1205, phy1, TRUE)     
Solidago    = prepare_raxml(1235, phy1, TRUE)   
Baccharis   = prepare_raxml(1241, phy1, TRUE)
Artemisia   = prepare_raxml(1248, phy1, TRUE)
Salvia      = prepare_raxml(1373, phy1, TRUE)
Gentiana    = prepare_raxml(1457, phy1, TRUE)
Trigonella  = prepare_raxml(1587, phy1, TRUE)
Vicia       = prepare_raxml(1590, phy1, TRUE)
Onobrychis  = prepare_raxml(1607, phy1, TRUE)
Salix       = prepare_raxml(1719, phy1, TRUE)
Ranunculus  = prepare_raxml(1869, phy1, TRUE)
Thalictrum  = prepare_raxml(1901, phy1, TRUE)
Bromus      = prepare_raxml(2030, phy1, TRUE)                                                      
# I obtained a lot of clusters for Carex and Solanum during initial June 2020
# analysis. However, I encountered errors in August 2021 when I reran with 
# software and database updates. However, since the number of sequence clusters
# were not limiting for these clades, I am using 2020 data rather than debugging
Carex       = prepare_raxml(2054, phy1, FALSE)
Solanum     = prepare_raxml(1463, phy1, FALSE, clade = "Solanum")

# Now run 'raxml/run-raxml.sh' to make phylogenies