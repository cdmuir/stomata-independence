source("r/header.R")

trimmed_data = read_rds("processed-data/trimmed-data.rds")
trimmed_phy = read_rds("processed-data/trimmed-phylogeny.rds")
pair_div = read_rds("objects/pair_div.rds")

bf_1 = bf(mvbind(adaxial_stomatal_length_um, adaxial_stomatal_density_mm2,
                    abaxial_stomatal_length_um, abaxial_stomatal_density_mm2) ~
               1, sigma ~ pair_age)

bf_2 = bf(mvbind(adaxial_stomatal_length_um, adaxial_stomatal_density_mm2,
                  abaxial_stomatal_length_um, abaxial_stomatal_density_mm2) ~
             1, sigma ~ 1)

fit_1 = brm(
  bf_1 + set_rescor(TRUE),
  data = pair_div,
  family = "student", 
  backend = "cmdstanr", chains = 2L, cores = 2L, seed = 579515
) |>
  add_criterion("loo")

fit_2 = brm(
  bf_2 + set_rescor(TRUE),
  data = pair_div,
  family = "student", 
  backend = "cmdstanr", chains = 2L, cores = 2L, seed = 152442
) |>
  add_criterion("loo")

write_rds(fit_1, "objects/fit_1.rds")
write_rds(fit_2, "objects/fit_2.rds")

loo_compare = loo_compare(fit_1, fit_2)

write_rds(loo_compare, "objects/loo_compare.rds")

file.copy(glue("objects/{fit}.rds", fit = rownames(loo_compare)[1]),
          "objects/fit.rds", overwrite = TRUE)
