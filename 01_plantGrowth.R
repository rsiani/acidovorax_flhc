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


# plant bioassays ------------------------------------------------------------

plant_growth =
  read_csv("data/growth_in_plates.csv") |>
  pivot_longer(control:Comp_140,
               names_to = "strain") |>
  mutate(
    strain = factor(strain,
                     levels = c("control", "LR140", "LR124", "Delta_140", "Comp_140")),
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

# modelling

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
  pluck("mod", 2) |> anova()

model_growth |>
  unnest(tidied) |>
  filter(term != "(Intercept)") |>
  mutate(fdr(p.value, 0.05), .by = metric)
