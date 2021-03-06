# CRAN packages ----
install.packages("ape")
install.packages("brms")
install.packages("checkmate")
install.packages("cowplot")
install.packages("crayon")
install.packages("diverge")
install.packages("ellipse")
install.packages("dplyr")
install.packages("ggdist")
install.packages("ggplot2")
install.packages("glue")
install.packages("hexbin")
install.packages("magrittr")
install.packages("mvnfast")
install.packages("phytools")
install.packages("purrr")
install.packages("purrrlyr")
install.packages("readr")
install.packages("readxl")
install.packages("RefManageR")
install.packages("seqinr")
install.packages("stringr")
install.packages("taxize")
install.packages("tibble")
install.packages("tidyr")
install.packages("units")

# GitHub packages ----
# Sys.unsetenv("GITHUB_PAT") # GITHUB_PAT was causing an issue
install.packages("devtools")

devtools::install_github("cdmuir/ropenstomata")
devtools::install_github("ropenscilabs/datastorr")
devtools::install_github("wcornwell/taxonlookup")
devtools::install_github("jinyizju/V.PhyloMaker")

# cmdstanr ----
# we recommend running this is a fresh R session or restarting your current session
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
