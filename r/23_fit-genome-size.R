source("r/header.R")

pair_div_genome = read_rds("objects/pair_div_genome.rds")

bf_1 = bf(mvbind(adaxial_stomatal_length_um, adaxial_stomatal_density_mm2,
                 abaxial_stomatal_length_um, abaxial_stomatal_density_mm2) ~
            dna_amount_2c_pg, sigma ~ pair_age)

fit_1_genome = brm(
  bf_1 + set_rescor(TRUE),
  data = pair_div_genome,
  family = "student", 
  backend = "cmdstanr", chains = 2L, cores = 2L, seed = 673308
)

write_rds(fit_1_genome, "objects/fit_1_genome.rds")
