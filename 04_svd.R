### SVD
## Roberto Siani
# 15.02.22


plot_flhcd =
  my_data |>
  filter(Gene %in% c("flhC", "flhD")) |>
  ggplot(aes(x = strain, y = lnRatio, color = flhc_media,
             fill = after_scale(color))) +
  geom_point(
    aes(shape = strain),
    position = position_jitterdodge(jitter.height = 0)) +
  geom_pointrange(stat = "summary",
                  color = "#555555",
                  position = position_dodge(width = 0.75),
                  shape = 95,
                  size = 3/.pt) +
  facet_wrap(~ Gene) +
  scale_color_manual(values = pal_bac4,
                     label = list("flhC+ Lj" = "_flhC_+ _Lj_",
                                  "flhC- Lj" = "_flhC_- _Lj_",
                                  "flhC+ Lj+Ri" = "_flhC_+ _Lj_+_Ri_",
                                  "flhC- Lj+Ri" = "_flhC_- _Lj_+_Ri_")) +
  theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.justification = "right",
        legend.text = ggtext::element_markdown(size = 12)) +
  scale_shape_manual(values = pal_shape,
                     label = list("Comp_140" = "LR140<sup>_ΔflhC;ΔflhC_</sup>",
                          "Delta_140" = "LR140<sup>_ΔflhC_</sup>"))  +
  scale_x_discrete(label = list("Comp_140" = "LR140<sup><i>ΔflhC;ΔflhC</i></sup>",
                                "Delta_140" = "LR140<sup>_ΔflhC_</sup>")) +
  scale_y_continuous(name = "Log Ratio") +
  guides(colour = guide_legend(override.aes = list(size = 5)),
         shape = guide_legend(override.aes = list(size = 5)))

plot_flhcd

my_ggsave("flhcd_diff", 3, 3)

# testing media differences

mod = my_data |>
  filter(Gene %in% c("flhC", "flhD")) |>
  lme(lnRatio ~ flhc * media * Gene,
      weights = varComb(varIdent(form = ~ 1 | flhc),
                        varFixed(~ tech_var)),
      random = ~ 1 | strain,
      data =  _)


marginaleffects::comparisons(mod,
                             variables = "media",
                               by = c("flhc", "Gene"),
                               p_adjust = "BH") |>
  tidy()



# singular value decomposition (wildtype only)

res_svd =
  my_data |>
  filter(mutant %in% "Wildtype") |>
  pivot_wider(id_cols = sample,
              names_from = Gene,
              values_from = lnRatio) |>
  column_to_rownames("sample") |>
  as.matrix() |>
  prcomp()

svd_df =
  broom::augment(res_svd) |>
  set_names(c("sample", broom::tidy(res_svd, "d") |>
                mutate(name = str_c("PC", PC, " ",
                                    round(percent, 2))) |>
                pull(name))) |>
  left_join(my_data |>
              filter(Gene %in% "flhC") |>
              select(sample, strain, depth, lnRatio,
                    flhc, flhc_media, sparsity, media, mutant))

# correlation with flhC, depth and sparsity

ggcorrplot::ggcorrplot(corr =
                         svd_df |>
                         select(where(is.numeric)) |>
                         cor() |>
                         round(1))

plot_svd(svd_df, 1, 2)

set.seed(1)
vegan::adonis2(dist(res_svd$x) ~ flhc * media + depth * sparsity,
               data = svd_df,
               permutations = 9999,
               by = "terms") |>
  as.data.frame() |>
  gt::gt(rownames_to_stub = T)

# now also with mutants

res_svd2 =
  my_data |>
  pivot_wider(id_cols = sample,
              names_from = Gene,
              values_from = lnRatio) |>
  column_to_rownames("sample") |>
  select(where(is.numeric)) |>
  as.matrix() |>
  prcomp()

svd_df2 =
  broom::augment(res_svd2) |>
  set_names(c("sample", broom::tidy(res_svd2, "d") |>
                mutate(name = str_c("PC", PC, " ",
                                    round(percent, 2))) |>
                pull(name))) |>
  left_join(my_data |>
              filter(Gene %in% "flhC") |>
              select(sample, strain, lnRatio, depth, sparsity,
                     flhc, flhc_media, media, mutant))

ggcorrplot::ggcorrplot(corr =
                         svd_df2 |>
                         select(where(is.numeric)) |>
                         cor() |>
                         round(1))

# visualization

plot_svd(svd_df, 1, 2) +
  plot_svd(svd_df2,1, 2) +
  plot_flhcd +
  plot_layout(design = "
              ######CCC
              AAABBBCCC")

my_ggsave("svd", 9, 6)

set.seed(1)
vegan::adonis2(dist(res_svd2$x) ~ flhc/mutant * media + depth * sparsity,
               data = svd_df2,
               permutations = 9999,
               by = "terms") |>
  as.data.frame() |>
  gt::gt(rownames_to_stub = T)

