### Plant Inoculation
### Roberto Siani
### 19.03.25

source("scripts/00_helperFunctions.R")

# bacterial colonization --------------------------------------------------

colonization_bac <-
  read_csv("data/CFU-0316.csv") |>
  mutate(
    strain = factor(strain),
    flhc = case_when(
      str_detect(strain, "LR140|Comp") ~ "flhC+",
      str_detect(strain, "LR124|Delta") ~ "flhC-",
      .default = "control"
    ) |>
      factor(levels = c("control", "flhC+", "flhC-")),
    CFUxg = CFU / `tissue weight`,
    plant_rep = str_c(strain, row_number()),
    AMF = factor(AMF),
    flhc_AMF = str_c(flhc, " ", AMF),
    .by = c(strain, AMF)
  )

mod2 <-
  nlme::lme(CFUxg ~ 0 + flhc + flhc:AMF,
    random = ~ 1 | strain,
    weight = nlme::varIdent(form = ~ 1 | strain * AMF),
    data = colonization_bac |> mutate(flhc = droplevels(flhc))
  )


qqnorm(mod2)

plot(resid(mod2))

augment(mod2) |>
  DescTools::LeveneTest(.resid ~ flhc * AMF, data = _)


df_mod <-
  marginaleffects::avg_predictions(
    mod2,
    variables = c("flhc", "AMF"),
    hypothesis = difference ~ pairwise
  ) |>
  tidy() |>
  separate_wider_delim(hypothesis, delim = " - ", names = c("start", "end")) |>
  mutate(
    start = str_remove_all(start, "\\(|\\)"),
    end = str_remove_all(end, "\\(|\\)"),
    p.adj = p.adjust(p.value, "fdr"),
    pos_y = seq(4500, 5000, by = 100)
  )



p1 <-
  colonization_bac |>
  ggplot(aes(
    x = flhc_AMF, y = CFUxg, color = flhc_AMF, shape = strain,
    fill = after_scale(color)
  )) +
  geom_point(position = position_dodge(width = .5), size = 3, alpha = .8) +
  stat_summary(
    geom = "pointrange",
    fun.min = ~ mean(.x) - sd(.x),
    fun.max = ~ mean(.x) + sd(.x),
    fun = mean,
    shape = 95, size = 5,
    show.legend = F, position = position_dodge(width = .5)
  ) +
  scale_color_manual(values = pal_growth7) +
  scale_shape_manual(
    values = pal_shape,
    label = list(
      "Comp_140" = "LR140<sup>_ΔflhC;flhC_</sup>",
      "Delta_140" = "LR140<sup>_ΔflhC_</sup>"
    )
  ) +
  scale_y_continuous(name = "CFU / Tissue Weight [g]") +
  theme(
    axis.text.x = element_text(face = "italic", angle = 30, hjust = 1),
    legend.text = ggtext::element_markdown(),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  guides(
    color = "none",
    shape = guide_legend(override.aes = list(size = 5))
  ) +
  geomtextpath::geom_textsegment(
    data = df_mod |>
      filter(row_number() %in% c(1, 2, 5, 6)),
    aes(
      x = start,
      xend = end,
      y = pos_y,
      yend = pos_y,
      label = format(p.adj, scientific = T, digits = 2)
    ), inherit.aes = F,
    size = 12 / .pt
  )


p1


# colonization R.irregularis ---------------------------------------------------

colonization_AMF <-
  read_tsv("data/colonization_190325.csv") |>
  mutate(
    flhc = case_when(
      str_detect(strain, "LR140|Comp") ~ "flhC+",
      str_detect(strain, "LR124|Delta") ~ "flhC-",
      .default = "control"
    ) |>
      factor(levels = c("control", "flhC+", "flhC-")),
    plant_rep = str_c(strain, row_number()), .by = strain
  )


mod1 <-
  colonization_AMF |>
  pivot_longer(Vesicle:Hyphae) |>
  nest(.by = name) |>
  mutate(fit = map(
    data,
    ~ lme4::glmer(cbind(value, 100 - value) ~ flhc + (1 | strain),
      family = "binomial",
      data = .x
    )
  ))

pluck(mod1, "fit") |> map(~ qqnorm(resid(.x)))
pluck(mod1, "fit") |> map(~ plot(resid(.x)))

pluck(mod1, "fit") |>
  map(~ DescTools::LeveneTest(.resid ~ flhc, data = augment(.x)))

df_mod <-
  pluck(mod1, "fit") |>
  map(~ marginaleffects::hypotheses(.x, difference ~ pairwise) |>
    tidy()) |>
  set_names(mod1$name) |>
  list_rbind(names_to = "name") |>
  separate_wider_delim(hypothesis, delim = " - ", names = c("start", "end")) |>
  mutate(
    start = case_match(
      start,
      "(flhcflhC+)" ~ "flhC+",
      "(flhcflhC-)" ~ "flhC-"
    ),
    end = case_match(
      end,
      "(flhcflhC+)" ~ "flhC+",
      "((Intercept))" ~ "control"
    ),
    p.adj = p.adjust(p.value, "fdr"),
    pos_y = c(80, 75, 70),
    name = factor(name, levels = c("Hyphae", "Vesicle", "Arbuscule")),
    .by = name
  )


p2 <- colonization_AMF |>
  pivot_longer(Vesicle:Hyphae) |>
  mutate(name = factor(name,
    levels = c("Hyphae", "Vesicle", "Arbuscule")
  )) |>
  ggplot(aes(
    x = flhc, y = value, color = flhc, shape = strain,
    fill = after_scale(color)
  )) +
  geom_point(position = position_dodge(width = .5), size = 3, alpha = .8) +
  stat_summary(
    geom = "pointrange",
    fun.min = ~ mean(.x) - sd(.x),
    fun.max = ~ mean(.x) + sd(.x),
    fun = mean,
    shape = 95, size = 5,
    show.legend = F
  ) +
  scale_color_manual(
    values = c("grey75", "#FC7D0B", "#1170AA"),
    labels = c(
      "flhC+" = "*flhC*+",
      "flhC-" = "*flhC*-"
    )
  ) +
  scale_shape_manual(
    values = pal_shape,
    labels = c(
      "control",
      "LR124",
      "LR140",
      "Comp_140" = "LR140{.sup *ΔflhC;flhC*}",
      "Delta_140" = "LR140{.sup *ΔflhC*}"
    )
  ) +
  scale_y_continuous(name = "Colonization Level (%)") +
  scale_x_discrete(labels = list(
    "flhC+" = "*flhC*+",
    "flhC-" = "*flhC*-"
  )) +
  theme(
    axis.text.x.bottom = marquee::element_marquee(angle = 30, hjust = 1),
    legend.text = marquee::element_marquee(),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    legend.justification.bottom = "right"
  ) +
  guides(
    # color = "none",
    shape = guide_legend(override.aes = list(size = 5))
  ) +
  facet_wrap(vars(name)) +
  geomtextpath::geom_textsegment(
    data = df_mod,
    aes(
      x = start,
      xend = end,
      y = pos_y,
      yend = pos_y,
      label = format(p.adj, scientific = T, digits = 2)
    ), inherit.aes = F,
    size = 12 / .pt
  )

p1 + p2 + plot_layout(design = "ABB") +
  plot_annotation(tag_levels = "a")

my_ggsave("figure7", 9, 5)
