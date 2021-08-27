source("r/header.R")

# Values for mean, cov, and range from Liu et al. 2021 
df_constraint = crossing(
  f_s = 1,
  A_s = exp(seq(2.444085, 8.034594, length.out = 101))
) %>%
  mutate(
    D_s = f_s * 1e6 / A_s,
    a_s = log(A_s),
    d_s = log(D_s),
    f_s = set_units(D_s, 1/mm^2) * set_units(A_s, um ^ 2)
  )

df_ellipse = ellipse(
  x = matrix(c(0.5187177, -0.2829251, -0.2829251, 0.5901746), ncol = 2),
  centre = c(5.571961, 5.121915), 
  level = 0.95
) |>
  as_tibble(.name_repair = NULL) |>
  rename(a_s = x, d_s = y)

panel_a = ggplot(df_ellipse, aes(exp(a_s), exp(d_s))) +
  scale_x_continuous(
    trans = "log", 
    limits = exp(5.571961 + c(-4.5, 4.5) * sqrt(0.5187177)),
    breaks = c(10, 100, 1000)
  ) + 
  scale_y_continuous(
    trans = "log", 
    limits = exp(5.121915 + c(-4.5, 5) * sqrt(0.5187177)),
    breaks = c(10, 100, 1000)
  ) + 
  geom_polygon(fill = "grey") +
  # geom_line(data = df_constraint) +
  xlab(expression(paste("Stomatal size [", mu, m ^2, "], log-scale"))) +
  ylab(expression(atop(Stomatal~density, paste(group("[", pores~mm ^-2, "]"), ",",~log-scale)))) +
  ggtitle("Inverse size-density scaling") +
  coord_equal() +
  theme_cowplot() +
  NULL

df_area = tibble(
  sr = seq(0, 1, 1e-2),
  y = 0.5 * dbeta(sr, 1, 10) + 0.5 * dbeta(sr, 10, 10)
)

panel_b = ggplot(df_area, aes(sr, y)) +
  geom_area(fill = "grey") +
  geom_text(data = tibble(sr = c(0, 0.5), y = 6, 
                          label = c("hypostomy", "amphistomy")),
            mapping = aes(label = label), hjust = 0) +
  geom_segment(data = tibble(sr = c(0, 0.5), xend = c(0, 0.5),
                             y = c(5, 2), yend = c(5.75, 5.75)),
            mapping = aes(xend = xend, yend = yend)) +
  xlab("Stomatal ratio [adaxial:abaxial]") +
  ylab("Probability density") +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  ylim(0, 6) +  
  ggtitle("Bimodal stomatal ratio") +
  coord_fixed(ratio = 1/6) +
  theme_cowplot() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ); panel_b

gp = plot_grid(panel_a, panel_b, nrow = 1, align = "hv", axis = "tl",
                rel_widths = 1, labels = "auto", hjust = -2)

ggsave("figures/concepts.pdf", width = 8, height = 4)
