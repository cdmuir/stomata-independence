source("r/header.R")

# Source saved data and objects ----
pair_div = read_rds("objects/pair_div.rds")
h2_pairs = read_rds("objects/h2_pairs.rds")

df_density = read_rds("objects/df_density.rds") |>
  transmute(r2_density = r_add_abd ^ 2, slope_density = sma_slope)

df_length  = read_rds("objects/df_length.rds") |>
  transmute(r2_length = r_adl_abl ^ 2, slope_length = sma_slope)

# Prepare for plotting ----
df_pointinterval = bind_cols(df_density, df_length) |>
  pivot_longer(everything()) |>
  mutate(
    trait = str_extract(name, "_[a-z]+$") |>
      str_remove("^_"),
    trait = factor(trait, levels = c("density", "length")),
    parameter = str_extract(name, "^[a-z0-9]+_") |>
      str_remove("_$"),
    parameter1 = case_when(
      parameter == "slope" ~ "SMA~slope",
      parameter == "r2" ~ "italic(r)^2",
    ),
    parameter1 = factor(parameter1, levels = c("SMA~slope", "italic(r)^2"))
  )

# Hypothesis testing ----
df_pointinterval |>
  select(name, value) |>
  group_by(name) |>
  point_interval(.point = median, .interval = hdi) |>
  write_rds("objects/h2_parameters.rds")

# Plot ----
h2_summary = ggplot(df_pointinterval, aes(trait, value)) +
  facet_wrap(. ~ parameter1, labeller = label_parsed) +
  geom_hline(yintercept = c(0, 1), linetype = "dashed") +
  stat_pointinterval(position = position_dodge(), show.legend = FALSE,
                     point_interval = median_hdi) +
  xlab("stomatal trait") +
  ylab("parameter estimate") +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_cowplot()

h2_summary1 = plot_grid(
  ggplot()+theme_void(), 
  h2_summary, 
  ggplot()+theme_void(), 
  nrow = 1, rel_widths = c(1, 2, 1)
)

plot_grid(h2_pairs, h2_summary1, nrow = 2, labels = "auto")

ggsave("figures/h2.pdf", width = 7, height = 7)
ggsave("figures/h2.png", width = 7, height = 7)
