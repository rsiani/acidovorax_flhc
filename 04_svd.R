### SVD
## Roberto Siani
# 15.02.22

# testing media differences

mod <- my_data |>
  filter(Gene %in% c("flhC", "flhD")) |>
  lme(lnRatio ~ flhc * media * Gene,
    weights = varComb(
      varIdent(form = ~ 1 | flhc),
      varFixed(~tech_var)
    ),
    random = ~ 1 | strain,
    data = _
  )


comp <- marginaleffects::comparisons(mod,
  variables = "media",
  by = c("flhc", "Gene")
)


pval_df <- marginaleffects::hypotheses(comp, multcomp = "BH") |>
  tidy() |>
  mutate(
    flhc = c("flhC+", "flhC+", "flhC-", "flhC-"),
    Gene = c("flhC", "flhD", "flhC", "flhD"),
    contrast = c("Lj+Ri - Lj"),
    strain = c("LR140", "LR140", "LR124", "LR124"),
    lnRatio = 5
  )


plot_flhcd <-
  my_data |>
  filter(Gene %in% c("flhC", "flhD")) |>
  ggplot(aes(
    x = strain, y = lnRatio, color = flhc_media,
    fill = after_scale(color)
  )) +
  geom_point(
    aes(shape = strain),
    position = position_jitterdodge(jitter.height = 0)
  ) +
  geom_pointrange(
    stat = "summary",
    color = "#555555",
    position = position_dodge(width = 0.75),
    shape = 95,
    size = 3 / .pt
  ) +
  geom_text(
    data = pval_df,
    color = "#555555",
    nudge_x = 0.5,
    size = 15 / .pt,
    aes(label = format(p.value, scientific = T, digits = 2))
  ) +
  facet_wrap(~Gene) +
  scale_color_manual(
    values = pal_bac4,
    labels = c(
      "flhC+ Lj" = "*flhC*+ *Lj*",
      "flhC- Lj" = "*flhC*- *Lj*",
      "flhC+ Lj+Ri" = "*flhC*+ *Lj*+*Ri*",
      "flhC- Lj+Ri" = "*flhC*- *Lj*+*Ri*"
    )
  ) +
  theme(
    axis.text.x.bottom = marquee::element_marquee(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.justification = "right",
    strip.text = element_text(face = "italic"),
    legend.text = marquee::element_marquee(size = 12)
  ) +
  scale_shape_manual(
    values = pal_shape,
    labels = c(
      "LR140" = "LR140",
      "LR124" = "LR124",
      "Comp_140" = "LR140{.sup *ΔflhC;ΔflhC*}",
      "Delta_140" = "LR140{.sup *ΔflhC*}"
    )
  ) +
  scale_x_discrete(labels = c(
    "LR140" = "LR140",
    "LR124" = "LR124",
    "Comp_140" = "LR140{.sup *ΔflhC;ΔflhC*}",
    "Delta_140" = "LR140{.sup *ΔflhC*}"
  )) +
  scale_y_continuous(name = "Log Ratio") +
  guides(
    colour = guide_legend(override.aes = list(size = 5)),
    shape = guide_legend(override.aes = list(size = 5))
  )

plot_flhcd

my_ggsave("flhcd_diff", 2.5, 5)

# singular value decomposition (wildtype only)

res_svd <-
  my_data |>
  filter(mutant %in% "Wildtype") |>
  pivot_wider(
    id_cols = sample,
    names_from = Gene,
    values_from = lnRatio
  ) |>
  column_to_rownames("sample") |>
  as.matrix() |>
  prcomp()

svd_df <-
  broom::augment(res_svd) |>
  set_names(c("sample", broom::tidy(res_svd, "d") |>
    mutate(name = str_c(
      "Comp. ", PC, " ",
      round(percent, 2)
    )) |>
    pull(name))) |>
  left_join(my_data |>
    filter(Gene %in% "flhC") |>
    select(
      sample, strain, depth, lnRatio,
      flhc, flhc_media, sparsity, media, mutant
    ))

# correlation with flhC, depth and sparsity

ggcorrplot::ggcorrplot(
  corr =
    svd_df |>
      select(where(is.numeric)) |>
      cor() |>
      round(1)
)

plot_svd(svd_df, 1, 2)

set.seed(1)
vegan::adonis2(dist(res_svd$x) ~ flhc * media + depth * sparsity,
  data = svd_df,
  permutations = 9999,
  by = "terms"
) |>
  as.data.frame() |>
  gt::gt(rownames_to_stub = T)

# now also with mutants

res_svd2 <-
  my_data |>
  pivot_wider(
    id_cols = sample,
    names_from = Gene,
    values_from = lnRatio
  ) |>
  column_to_rownames("sample") |>
  select(where(is.numeric)) |>
  as.matrix() |>
  prcomp()

svd_df2 <-
  broom::augment(res_svd2) |>
  set_names(c("sample", broom::tidy(res_svd2, "d") |>
    mutate(name = str_c(
      "Comp. ", PC, " ",
      round(percent, 2)
    )) |>
    pull(name))) |>
  left_join(my_data |>
    filter(Gene %in% "flhC") |>
    select(
      sample, strain, lnRatio, depth, sparsity,
      flhc, flhc_media, media, mutant
    ))

ggcorrplot::ggcorrplot(
  corr =
    svd_df2 |>
      select(where(is.numeric)) |>
      cor() |>
      round(1)
)

# visualization

plot_svd(svd_df, 1, 2) +
  plot_svd(svd_df2, 1, 2) +
  plot_flhcd +
  plot_layout(design = "
              ######CCC
              AAABBBCCC")

my_ggsave("svd", 9, 6)

set.seed(1)
vegan::adonis2(dist(res_svd2$x) ~ flhc / mutant * media + depth * sparsity,
  data = svd_df2,
  permutations = 9999,
  by = "terms"
) |>
  as.data.frame() |>
  gt::gt(rownames_to_stub = T)
