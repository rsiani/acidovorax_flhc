### Plate motility asssay
## Roberto Siani
# 2024

source("scripts/00_helperFunctions.R")

# swimming ----------------------------------------------------------------

raw_swimming <-
  read_csv("data/diameter.csv") |>
  mutate(
    Strain = factor(Strain,
      levels = c(
        "control",
        "LR140",
        "LR124",
        "Delta_140",
        "Comp_140"
      )
    ),
    flhc = case_match(
      Strain,
      "LR124" ~ "flhC-",
      "LR140" ~ "flhC+",
      "Delta_140" ~ "flhC-",
      "Comp_140" ~ "flhC+"
    ) |>
      factor(levels = c("flhC+", "flhC-")),
    mutant = case_match(
      Strain,
      "LR124" ~ "Wildtype",
      "LR140" ~ "Wildtype",
      "Delta_140" ~ "Mutant",
      "Comp_140" ~ "Mutant"
    ) |>
      factor(levels = c("Wildtype", "Mutant"))
  ) |>
  pivot_longer(starts_with("Rep"),
    values_to = "diameter [mm]"
  )

ggplot(raw_swimming) +
  geom_point(
    aes(
      x = Strain,
      y = `diameter [mm]`,
      shape = Strain,
      fill = after_scale(color),
      color = flhc
    ),
    size = 3,
    position = ggbeeswarm::position_quasirandom()
  ) +
  scale_color_manual(
    values = pal_bac2,
    label = list(
      "flhC+" = "*flhC*+",
      "flhC-" = "*flhC*-"
    )
  ) +
  facet_wrap(~mutant, scales = "free_x") +
  theme(
    axis.text.x.bottom = marquee::element_marquee(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.direction = "vertical",
    axis.title.x = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.text = marquee::element_marquee(size = 12)
  ) +
  scale_shape_manual(
    values = c(21, 24, 2, 1),
    label = list(
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
  guides(
    colour = guide_legend(override.aes = list(size = 4)),
    shape = guide_legend(override.aes = list(size = 4))
  )

my_ggsave("swimming_diameter", 3, 3.5)
