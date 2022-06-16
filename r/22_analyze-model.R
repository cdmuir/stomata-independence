source("r/header.R")

# Use TRY to get ballpark range of stomatal index and its variance among species
si = read_excel("raw-data/21413_02062022232051/21413.xlsx") |>
  filter(TraitID == 67) |>
  mutate(stomatal_index = as.numeric(OrigValueStr))

# hist(si$stomatal_index)
# summary(si$stomatal_index)
# sd(si$stomatal_index / 100)

# The commented-out code uses numerical integration. Delete once I am sure that
# I will not use it.
df = crossing(
  mu_i = log(c(0.1, 0.15, 0.2)),
  sigma_i = 0.1, # 10 ^ c(-1.5, -1, -0.5),
  mu_m = log(200), # log( 200 um^3 ), but value does not actually matter
  rel_sigma_m = 10 ^ seq(-1, 1, 0.1),
  A = 5, # Assuming epidermal cell area is ~1000 um^2, but value does not actually matter
  B = 0.25, # actual value has a very small effect on Corr_ds
  n_sim = 1e5
) |>
  rowwise() |>
  mutate(
    sigma_m = sigma_i * rel_sigma_m,
    # sigma2_logitI = Var_logitI(mu_i, sigma_i),
    # Mu is mean of logit_I and log(M)
    # Mu = list(c(E_logitI(mu_i, sigma_i) , mu_m)),
    # Sigma is Covariance matrix of logit_I and log(M)
    # Sigma = list(matrix(c(sigma2_logitI, 0, 0, sigma_m ^ 2), ncol = 2)),
    # Var_d = Var_d(A, B, Mu, Sigma),
    # Var_s = Var_s(A, B, Mu, Sigma),
    # Cov_ds = Cov_ds(A, B, Mu, Sigma),
    # Corr_ds = Cov_ds / (sqrt(Var_d) * sqrt(Var_s)),
    sim = list(tibble(
      A = A, B = B, 
      m = rnorm(n_sim, mu_m, sigma_m),
      # logit_I = rnorm(n_sim, E_logitI(mu_i, sigma_i), sqrt(sigma2_logitI)),
      # I = plogis(logit_I),
      i = rnorm(n_sim, mu_i, sigma_i),
      I = exp(i),
      e = log(A) + m, E = exp(e),
      s = log(B) + e, S = exp(s),
      D = I / (I * S + (1 - I) * E),
      d = log(D)
    ) |>
      filter(I > 0, I < 1)),
    Var_d = var(sim$d),
    Var_s = var(sim$s),
    Cov_ds = cov(sim$d, sim$s),
    Corr_ds = cor(sim$d, sim$s),
    predicted_Corr_ds = - sigma_m ^ 2 / (sigma_m * sqrt(sigma_i ^ 2 + sigma_m ^ 2)) 
  ) |>
  mutate(mu_I = glue('mu[italic(I)]=={I}', I = exp(mu_i)))

ggplot(df, aes(predicted_Corr_ds, Corr_ds)) +
  facet_grid(rows = "mu_I", labeller = label_parsed) +
  geom_abline(intercept = 0, slope = 1, color = "grey", size = 1.1, linetype = "dashed") +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  xlim(-1, 0) +
  ylim(-1, 0) +
  xlab(expression(paste("Approximate Corr[", Delta * italic(d), ",", Delta * italic(s), "]"))) +
  ylab(expression(paste("Simulated Corr[", Delta * italic(d), ",", Delta * italic(s), "]"))) +
  coord_equal() +
  theme_cowplot()

ggsave("figures/check-approximation.pdf", width = 3.5, height = 8)

#### MOST OF THIS STUFF BELOW IS PROBABLY NOT USEFUL. NUMERICAL INTEGRATION DIDN'T WORK AS WELL AS RANDOM SIMULATIONS ----
# Check that functions to numerically calculate variance and covariance are
# correct. These lines can be removed later 
Mu = c(0, 0)
Sigma = diag(2)

# Check that fd integrates to 1
integrate(
  fd,
  lower = -Inf, upper = Inf,
  A = 10, B = 3, Mu = Mu, Sigma = Sigma, log = FALSE
)

# Check that fs integrates to 1
integrate(
  fs,
  lower = -Inf, upper = Inf,
  A = 10, B = 3, Mu = Mu, Sigma = Sigma, log = FALSE
)

# Check that Cov_ds is correct
rmvn(1e6, Mu, Sigma) |>
  as_tibble(.name_repair = "unique") |>
  rename(logit_I = ...1, m = ...2) |>
  mutate(
    A = 10, B = 0.1,
    M = exp(m), I = plogis(logit_I),
    E = A * M, S = B * E,
    D = I / (I * S + (1 - I) * E),
    d = log(D),
    s = log(S)
  ) |>
  select(d, s) |>
  cov()

Var_d(A = 10, B = 0.1, Mu, Sigma) # 1.7746017
Var_s(A = 10, B = 0.1, Mu, Sigma) # 0.9986218
Cov_ds(A = 10, B = 0.1, Mu, Sigma) # -0.9990488

# Points I can make with model:
# Developmental constraint -> integration
# (Relative) variance in I (comapred to m) modulates strength of integration
# If m is primary driver of s-d relationship, slope is ??? relative var is ???
# In order for m to drive major patterns of stomatal evolution:
#   - A and B fixed
#   - Var(logit_I) cannot be much larger than Var(m)
#   
#   
# Questions, if A and B are fixed:
# How does Cov(d,s) behave with different variance in logit_I versus m
# How could we get isometry in size between surfaces, but more variance in density? If I can vary between surfaces

# Reasonable parameter values
hist(plogis(rnorm(1e3, qlogis(0.1), 0.5)), breaks = seq(0, 1, 0.05))
# WHAT'S A REASONABLE VALUE FOR mu_m?
Mu = c(qlogis(0.1), 0)
Sigma = diag(2) * 0.5
Cov_ds(A = 1e3, B = 0.1, Mu, Sigma)
Sigma_ds(A = 1e1, B = 0.5, Mu, Sigma) |>
  cov2cor()
Es(A = 1e3, B = 0.1, Mu, Sigma)
Es(A = 1e2, B = 0.1, Mu, Sigma)

# Scenarios:
Mu = c(qlogis(0.1), 0)

# Var(logit_I) > Var(m)
Sigma = matrix(c(0.5, 0, 0, 0.05), ncol = 2)
Sigma_ds(A = 1e3, B = 0.1, Mu, Sigma) |>
  cov2cor()

# Var(logit_I) = Var(m)
Sigma = matrix(c(0.5, 0, 0, 0.5), ncol = 2)
Sigma_ds(A = 1e3, B = 0.1, Mu, Sigma) |>
  cov2cor()

# Var(logit_I) < Var(m)
Sigma = matrix(c(0.05, 0, 0, 0.5), ncol = 2)
Sigma_ds(A = 1e3, B = 0.1, Mu, Sigma) |>
  cov2cor()

# Equations
# E = A M
# S = B E = A B M
# D = I / (IS + (1 - I) E)
# 
# log-transformed
# e = a + m
# s = a + b + m
# d = i - log(IS + (1 - I) E)

A = 10
B = 0.1
M = 0.1
I = 0.15
E = A * M
S = B * E
D = I / (I * S + (1 - I) * E)
d = log(D)
m = log(M)
s = log(S)
i = exp(d + log(A) + m - log(1 + exp(log(A) + m + d) - exp(d + s)))

A = 100
B = 0.1

mu_I = 0.1
mu_M = 2
mu_logit_I = qlogis(mu_I)
mu_m = log(mu_M)
sigma_logit_I = 0.1
sigma_m = 0.1
n = 1e3
dat = tibble(
  A = A,
  B = B,
  a = log(A),
  b = log(B),
  logit_I = rnorm(n, mu_logit_I, sigma_logit_I),
  I = plogis(I),
  i = log(I),
  m = rnorm(n, mu_m, sigma_m),
  M = exp(m),
  s = a + b + m,
  S = exp(s),
  E = exp(a + m),
  D = I / (I * S + (1 - I) * E),
  d = i - log(I * S + (1 - I) * E),
  d1 = i - a - m - log(I * B + (1 - I))
)

select(dat, d, d1)
qplot(dat$s, dat$d)
cov(select(dat, s, d))
cor(select(dat, s, d))


# inner integration. m is set. integrate {i, -Inf, Inf}
# (value of d) * (probability of d conditional on m)
mu_logit_I = qlogis(mu_I)
Mu = c(mu_logit_I, mu_m)
Sigma = matrix(c(sigma_logit_I ^ 2, 0, 0, sigma_m ^ 2), 2, 2)

# Probability density function of d
I = 0.15
A = 10
B = 0.1
M = 3
E = A * M
S = A * B * M

I / (I * S + (1 - I) * E)
D = I / (M * (I * A * B + (1 - I) * A))

