### Plate motility asssay
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
  scale_color_manual(values = pal_bac2,
                     label = list("flhC+" = "<i>flhC</i>+",
                                  "flhC-" = "<i>flhC</i>-")) +
  facet_wrap(~ mutant, scales = "free_x") +
  theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.direction = "vertical",
        axis.title.x = element_blank(),
        legend.margin = margin(0, 0,0,0),
        legend.text = ggtext::element_markdown(size = 12)) +
  scale_shape_manual(values = c(21, 24, 2, 1),
                     label = list("Comp_140" = "LR140<sup>_ΔflhC;ΔflhC_</sup>",
                                  "Delta_140" = "LR140<sup>_ΔflhC_</sup>"))  +
  scale_x_discrete(label = list("Comp_140" = "LR140<sup>_ΔflhC;ΔflhC_</sup>",
                                "Delta_140" = "LR140<sup>_ΔflhC_</sup>")) +
  guides(colour = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4)))

my_ggsave("swimming_diameter", 7, 8.5)


