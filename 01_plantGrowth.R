### Plant growth data
## Roberto Siani
# 2024

source("scripts/00_helperFunctions.R")



# swimming ----------------------------------------------------------------

raw_swimming =
  read_csv("data/diameter.csv") |>
  mutate(
    Strain = factor(Strain,
                   levels = c("control",
                              "LR140",
                              "LR124",
                              "Delta_140",
                              "Comp_140")),
    flhc = case_match(
      Strain,
    "LR124" ~ "flhC-",
    "LR140" ~ "flhC+",
    "Delta_140" ~ "flhC-",
    "Comp_140" ~ "flhC+") |>
      factor(levels = c("flhC+", "flhC-")),
    mutant = case_match(
      Strain,
      "LR124" ~ "Wildtype",
      "LR140" ~ "Wildtype",
      "Delta_140" ~ "Mutant",
      "Comp_140" ~ "Mutant") |>
      factor(levels = c("Wildtype", "Mutant"))) |>
  pivot_longer(starts_with("Rep"),
               values_to = "length [mm]")

ggplot(raw_swimming) +
  geom_point(aes(x = Strain,
                 y = `length [mm]`,
                 # ymin = `length [mm]` - std.dev,
                 # ymax = `length [mm]` + std.dev,
                 shape = Strain,
                 fill = after_scale(color),
                 color = flhc),
             size = 3,
             position = ggbeeswarm::position_quasirandom()) +
  scale_color_manual(values = pal_bac2) +
  facet_wrap(~ mutant, scales = "free_x") +
  theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1)) +
  scale_shape_manual(values = c(21, 24, 2, 1))  +
  scale_x_discrete(label = list("Comp_140" = "LR140<sup>_ΔflhC;ΔflhC_</sup>",
                                "Delta_140" = "LR140<sup>_ΔflhC_</sup>"))

my_ggsave("swimming_diameter", 2, 3)

# source data -------------------------------------------------------------

# raw_growth =
#   read_csv("data/growth_data.csv") |>
#   mutate(
#     AMF = factor(AMF, levels = c("Lj", "Lj+Ri")),
#     strain = factor(strain,
#                     levels = c("control",
#                                "LR140",
#                                "Comp_140",
#                                "LR124",
#                                "Delta_140")),
#     flhc = case_match(
#       strain,
#       "control" ~ "control",
#       "LR124" ~ "flhC-",
#       "LR140" ~ "flhC+",
#       "Delta_140" ~ "flhC-",
#       "Comp_140" ~ "flhC+") |>
#       factor(levels = c("control", "flhC+", "flhC-")),
#     flhc_AMF = str_c(flhc, AMF, sep = " ") |>
#       factor(levels = c(
#         "control Lj", "control Lj+Ri",
#         "flhC+ Lj", "flhC+ Lj+Ri",
#         "flhC- Lj", "flhC- Lj+Ri")),
#     mutant = case_match(strain,
#                         "control" ~ "control",
#                         "LR140" ~ "Wildtype",
#                         "LR124" ~ "Wildtype",
#                         .default = "Mutant") |>
#       factor(levels = c("control", "Wildtype", "Mutant")),
#     Plant = str_c(Pot, Plant),
#     Pot2 = str_c(strain, Pot))
#
# skimr::skim(raw_growth)
#
# tidy_growth =
#   raw_growth |>
#   mutate(fresh_g.Root = gross.fresh.Root - alu.Root,
#          fresh_g.Shoot = gross.fresh.Shoot - alu.Shoot,
#          dry_g.Root = gross.dry.Root - alu.Root,
#          dry_g.Shoot = gross.dry.Shoot - alu.Shoot,
#          .keep = "unused") |>
#   mutate(
#     water_g.Root = fresh_g.Root - dry_g.Root,
#     water_g.Shoot = fresh_g.Shoot - dry_g.Shoot) |>
#   pivot_longer(where(is.numeric)) |>
#   separate_wider_delim(name, names = c("metric", "compartment"), delim = ".") |>
#   filter(metric %in% c("dry_g", "water_g", "length_mm")) |>
#   mutate(metric = factor(metric,
#                             levels = c("length_mm",
#                                                "dry_g",
#                                                "water_g"),
#                             labels = c("length [mm]",
#                                        "dry weight [g]",
#                                        "water weight [g]")))
#
# ggplot(tidy_growth,
#        aes(
#          x = strain,
#          y = value,
#          color = flhc_AMF,
#          fill = flhc_AMF
#        )) +
#   # geom_violin(
#   #   data = ~ filter(.x, AMF == "AMF"),
#   #   aes(xmin = after_scale(x)),
#   #   color = "white",
#   #   alpha = .5,
#   #   linewidth = .5,
#   #   draw_quantiles = c(0.25, 0.5, 0.75),
#   #   position = position_nudge(x = 0.1)) +
#   # geom_violin(
#   #   data = ~ filter(.x, AMF == "Lj"),
#   #   aes(xmax = after_scale(x)),
#   #   color = "white",
#   #   alpha = .33,
#   #   linewidth = .5,
#   #   draw_quantiles = c(0.25, 0.5, 0.75),
#   #   position = position_nudge(x = -0.1)) +
#   geom_point(
#     aes(shape = strain),
#     size = 1,
#     alpha = .9,
#     position = ggbeeswarm::position_quasirandom(dodge.width = 0.5)) +
#   geom_point(
#     stat = "summary",
#     fun = mean,
#     shape = 95,
#     size = 3 * .pt,
#     color = "#555555",
#     position = position_dodge(width = 0.5)
#   ) +
#   facet_grid(metric ~ compartment, scale = "free", switch = "y",) +
#   scale_color_manual(values = pal_growth7,
#                      aesthetics = c("color", "fill")) +
#   scale_y_continuous(limits = c(0, NA)) +
#   theme(axis.title.y = element_blank(),
#         strip.placement = "outside",
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_shape_manual(values = c(22, 21, 1, 24, 2))
#
# my_ggsave("growth_data", 4, 6)

# fit a mixed linear model

# res_fit =
#   tidy_growth |>
#   mutate(z.score = (value - mean(value)) / sd(value),
#          .by = c(metric, compartment),
#          compartment = as.factor(compartment),
#          metric = as.factor(metric)) |>
#   lmerTest::lmer(z.score ~ flhc * AMF +
#                    (1 | strain/Pot2/Plant/compartment/metric),
#                  data = _)
#
# broom::augment(res_fit) |>
#   DescTools::LeveneTest(.resid ~ flhc * AMF,
#                         data = _)
#
# qqnorm(resid(res_fit))
#
# anova(res_fit) |> gt::gt(rownames_to_stub = T)
#
# list("Effects on growth" = res_fit) |>
#   modelsummary::modelsummary(
#     fmt =  modelsummary::fmt_statistic(estimate = 2, conf.low = 2, conf.high = 2,
#                                        p.value = modelsummary::fmt_sci(2)),
#     estimate = c("Estimate" = "{estimate} [{conf.low}, {conf.high}]"),
#     statistic = c("p" = "{p.value} {stars}"),
#     coef_omit = "Intercept|SD",
#     coef_rename = c("flhcflhC+" = "flhC+",
#                     "flhcflhC-" = "flhC-",
#                     "flhcflhC+:AMFLj+Ri" = "flhC+ Lj+Ri",
#                     "flhcflhC-:AMFLj+Ri" = "flhC- Lj+Ri",
#                     "AMFLj+Ri" = "Lj+Ri"),
#     gof_omit = 'DF|Deviance|R2|RMSE|BIC',
#     shape = term ~ model + statistic,
#     align = "lll",
#     output = "figures/table_effects_growth.png")
#
# tidied_fit =
#   res_fit |>
#   broom::tidy(effects = "fixed") |>
#   filter(term != "(Intercept)") |>
#   mutate(fdr(p.value),
#          term = case_match(term,
#                            .default = NA,
#                            "flhcflhC+" ~ "flhC+",
#                            "flhcflhC-" ~ "flhC-",
#                            "flhcflhC+:AMFLj+Ri" ~ "flhC+ Lj+Ri",
#                            "flhcflhC-:AMFLj+Ri" ~ "flhC- Lj+Ri",
#                            "AMFLj+Ri" ~ "Lj+Ri") |>
#            factor(levels = c(c("Lj+Ri", "flhC+", "flhC+ Lj+Ri",
#                                "flhC-", "flhC- Lj+Ri"))))
#
# p_estimate =
#   tidied_fit |>
#   ggplot() +
#   geom_pointrange(aes(x = term,
#                       y = estimate,
#                       ymin = estimate - std.error,
#                       ymax = estimate + std.error,
#                       color = term),
#                   size = 1,
#                   linewidth = .5) +
#   scale_color_manual(values = pal_bac7) +
#   geom_text(
#     aes(x = term,
#         y = 1.25,
#         fontface = if_else(FDR < 0.05, "bold", "plain"),
#         label = format(FDR, scientific = T, digits = 2)),
#     color = "#555555",
#     family = "Arial",
#     size = 12/.pt) +
#   theme(
#     axis.title.x = element_blank(),
#     panel.grid.major.y = element_line(linetype = 3,
#                                       color = "#757575")) +
#   scale_y_continuous(breaks = seq(-1, 1, 0.5))
#
# p1 + p2 + p_estimate +
#   plot_layout(guides = "collect",
#               design = "
#               AAB
#               AAB
#               CCC")
#
#
# my_ggsave("growth_data", 5.5, 3.5)


plant_growth =
  read_csv("data/growth_in_plates.csv") |>
  pivot_longer(control:Comp_140,
               names_to = "strain") |>
  mutate(
    strain = factor(strain, levels = c("control", "LR140", "LR124", "Delta_140", "Comp_140")),
    flhc = case_match(
    strain,
    "LR124" ~ "flhC-",
    "LR140" ~ "flhC+",
    "Delta_140" ~ "flhC-",
    "Comp_140" ~ "flhC+",
    .default = "control") |>
      factor(levels = c("control", "flhC+", "flhC-")),
    mutant = case_match(
      strain,
      "LR124" ~ "Wildtype",
      "LR140" ~ "Wildtype",
      "Delta_140" ~ "Mutant",
      "Comp_140" ~ "Mutant",
      .default = "control") |>
      factor(levels = c( "control", "Wildtype", "Mutant")),
    plate_factor = str_c(strain, plate, sep = " "))

plant_growth |>
  ggplot(aes(strain, value, color = flhc)) +
  geom_point(aes(shape = strain,
                 fill = after_scale(color)),
             size = 1, alpha = .9,
             position = ggbeeswarm::position_quasirandom()) +
  geom_point(
    stat = "summary",
    fun = mean,
    shape = 95,
    size = 3 * .pt,
    color = "#555555",
    position = position_dodge(width = 0.5)) +
  facet_grid(metric ~ ., scales = "free", switch = "y", ) +
  scale_color_manual(values = c("grey69", pal_bac2)) +
  scale_shape_manual(values = c(22, 21, 24, 2, 1)) +
  theme(axis.title.y = element_blank(),
        strip.placement = "outside",
        axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1, vjust = 1),
        strip.text.y.left = element_text(hjust = 1)) +
  scale_x_discrete(label = list("Comp_140" = "LR140<sup>_ΔflhC;ΔflhC_</sup>",
                                "Delta_140" = "LR140<sup>_ΔflhC_</sup>"))


my_ggsave("growth_data", 4, 6)

model_growth =
  plant_growth |>
  filter(metric != "CFU") |>
  nest(.by = metric) |>
  mutate(mod = map(data, ~ nlme::lme(value ~ strain,
                                     random = ~ 1 | plate_factor,
                                     weights = varIdent(form = ~ 1 | strain),
                                     data = .x)),
         tidied = map(mod, ~ tidy(.x, effects = "fixed")))

model_growth |>
  pluck("mod", 1) |> anova()

model_growth |>
  unnest(tidied) |>
  filter(term != "(Intercept)") |>
  mutate(fdr(p.value, 0.05), .by = metric)

#
# aug = model_growth |>
#   mutate(augmented = map(mod, augment)) |>
#   unnest(augmented)
#
# qqnorm(aug |> filter(metric %in% "shoot length [mm]") |> pull(.resid))
