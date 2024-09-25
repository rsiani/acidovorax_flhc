### Differential Gene/Pathway Expression
## Roberto Siani
# 2024

source("scripts/preProcess.R")

# model --------------------------------------------------------------------

# mixed-effects models

res_lme =
  my_data |>
  nest(.by = c(Gene, target_id, Product,
               Description, InterPro, Pfam,
               starts_with("PGP"))) |>
  mutate(mod =
           map(
             data,
             .progress = T,
             \(.x) lme(
               lnRatio ~ flhc * media,
               weights = varComb(varIdent(form = ~ 1 | flhc),
                                 varFixed(~ tech_var)),
               random = ~ 1 | strain,
               data = .x)))

# for quality check

augmented =
  res_lme |>
  mutate(augmented = map2(mod, data,
                          \(.x, .y) broom.mixed::augment(.x, data = .y))) |>
  unnest(augmented)

qqnorm(augmented$.resid)

ggplot(augmented) +
  geom_boxplot(aes(.resid, str_c(strain, media)))

# extract model summaries

tidied =
  res_lme |>
  mutate(tidied = map(mod,
                      \(.x) tidy(.x, effects = "fixed"))) |>
  unnest(tidied) |>
  filter(!term %in% "(Intercept)") |>
  mutate(adjust_Q(p.value, 0.05),
         term = case_match(term,
                           "flhcflhC-" ~ "Expression",
                           "mediaLj+Ri" ~ "Response",
                           "flhcflhC-:mediaLj+Ri" ~ "Regulation") |>
           factor(levels = c("Expression", "Response", "Regulation")),
         side = case_when(
           term == "Expression" & statistic > 0 ~ "flhC-",
           term == "Response" & statistic > 0 ~ "Lj+Ri",
           term == "Regulation" & statistic > 0 ~ "flhC- Lj+Ri",
           term == "Expression" & statistic < 0 ~ "flhC+",
           term == "Response" & statistic < 0 ~ "Lj",
           term == "Regulation" & statistic < 0 ~ "flhC+ Lj+Ri") |>
           factor(levels = c("flhC+", "flhC-", "Lj", "Lj+Ri",
                             "flhC+ Lj+Ri", "flhC- Lj+Ri")))

tidied |>
  janitor::tabyl(side, significant)

tidied |>
  janitor::tabyl(side)


# upset plot -------------------------------------------------------------

upset_df =
  tidied |>
  filter(significant) |>
  mutate(DEG = 1,
         side = str_c(term, side, sep = ": "),
         PMRs = if_else(is.na(PGP), "", "PMR")) |>
  pivot_wider(id_cols = c(Gene, PMRs),
              names_from = side, values_from = DEG, values_fill = list(DEG = 0)) |>
  as.data.frame()

sets =
  c("Regulation: flhC+ Lj+Ri", "Regulation: flhC- Lj+Ri",
    "Response: Lj", "Response: Lj+Ri",
    "Expression: flhC+", "Expression: flhC-")

meta_for_sets =
  data.frame(
    sets = sets) |>
  mutate(
    type = str_extract(sets, "flhC\\+|flhC\\-"))

pal_sets =
  c(pal_bac4[c(2,4)], "grey80", "grey90", pal_bac4[c(1,3)]) |>
  set_names(sets)

png("upset.png",width = 6, height = 3, res = 300, units = "in")
UpSetR::upset(upset_df,
              nintersects = 15,
              sets = sets,
              text.scale = c(1, 1, 1, 1, 1, 1.5),
              order.by = "freq",
              mb.ratio = c(0.5, 0.5),
              keep.order = T,
              point.size = 2,
              line.size = 1,
              set_size.show= T,
              set_size.scale_max = 650,
              sets.bar.color = pal_sets,
              matrix.color = "grey40",
              main.bar.color = "grey80",
              query.legend = "top",
              queries = list(
                list(query = UpSetR::elements,
                     params = list("PMRs", "PMR"),
                     color = "grey40",
                     active = T,
                     query.name = "PMRs")),
              set.metadata =
                list(
                  data = meta_for_sets,
                  plots = list(
                    list(type = "matrix_rows",
                         column = "sets",
                         alpha = 1,
                         colors = pal_sets))))
dev.off()


# enrichment --------------------------------------------------------------

dataset_enrich =
  tidied |>
  drop_na(PGP) |>
  select(Gene, side, term, starts_with("PGP")) |>
  pivot_longer(cols = c(starts_with("PGP")),
               values_to = "category",
               names_to = "level") |>
  mutate(category = str_to_lower(category) |>
           str_replace_all("\\|", ", ") |>
           str_replace_all("_", " ") |>
           str_replace("prtoeins", "proteins") |>
           str_replace_all("plant vitamin|root vitamin", "vitamin")) |>
  separate_rows(category, sep = "; ") |>
  filter(str_detect(category, "putative", negate = T),
         level %in% c("PGP2", "PGP3", "PGP4")) |>
  distinct() |>
  count(category, side, level, term) |>
  mutate(side_generic = case_when(
    side %in% c("Lj", "flhC+", "flhC+ Lj+Ri") ~ "one",
    .default = "two")) |>
  select(side, side_generic, category, n, level, term) |>
  pivot_wider(names_from = side_generic, values_from = c(n,
                                                         side)) |>
  replace_na(list(`n_one` = 0,
                  `n_two` = 0)) |>
  mutate(diff = abs(`n_one` - `n_two`),
         tot = `n_one` + `n_two`)


# pathway analysis --------------------------------------------------------

pmr_model =
  tidied |>
  mutate(across(starts_with("PGP"), ~ (.x) |>
                  str_replace("prtoeins", "proteins") |>
                  str_replace_all("plant vitamin|root vitamin", "vitamin"))) |>
  separate_rows(PGP5, sep = "; ") |>
  drop_na(PGP5) |>
  filter(n_distinct(Gene) >= 3, .by = PGP5) |>
  mutate(n_of_genes = n(), .by = c(PGP5, term)) |>
  nest(.by = c(PGP5, n_of_genes)) |>
  mutate(mod =
           map(
             data,
             .progress = T,
             \(.x) gls(
               estimate ~ 0 + term,
               weights = ~ std.error ^ 2,
               data = .x)))

pmr_tidied =
  pmr_model |>
  mutate(tidied = map(mod,
                      \(.x) tidy(.x, effects = "fixed"))) |>
  unnest(tidied) |>
  hoist(data, "Gene") |>
  select(-c(data, mod)) |>
  filter(!term %in% "(Intercept)") |>
  mutate(adjust_Q(p.value, 0.05),
         term = str_remove(term, "term") |>
           factor(levels = c("Expression", "Response", "Regulation")),
         side = case_when(
           term == "Expression" & statistic > 0 ~ "flhC-",
           term == "Response" & statistic > 0 ~ "Lj+Ri",
           term == "Regulation" & statistic > 0 ~ "flhC- Lj+Ri",
           term == "Expression" & statistic < 0 ~ "flhC+",
           term == "Response" & statistic < 0 ~ "Lj",
           term == "Regulation" & statistic < 0 ~ "flhC+ Lj+Ri") |>
           factor(levels = c("flhC+", "flhC-", "Lj", "Lj+Ri",
                             "flhC+ Lj+Ri", "flhC- Lj+Ri"))) |>
  left_join(all_pgp |>
              select(PGP1, PGP2, PGP3, PGP4, PGP5) |>
              distinct())


# plotting ----------------------------------------------------------------

(plot_volcano("Expression") /
  plot_enrich("Expression")) |
  plot_pmr("Expression")

my_ggsave("expression", 8, 8)

(plot_volcano("Response") /
    plot_enrich("Response")) |
  plot_pmr("Response")

my_ggsave("Response", 8, 10)

(plot_volcano("Regulation") /
    plot_enrich("Regulation")) |
  plot_pmr("Regulation")

my_ggsave("Regulation", 8, 8)


#
# a = plot_volcano(c("Expression", "Response", "Regulation")) +
#   facet_grid(. ~ term)
# b = plot_pmr(c("Expression", "Response", "Regulation")) +
#   facet_grid(PGP3 ~ term, scales = "free_y", space = "free_y")
#
# a + b +
#   plot_layout(heights = c(2/5, 3/5))
#
# my_ggsave("volcano_pmr", 8, 10)
