### Differential Gene/Pathway Expression
## Roberto Siani
# 2024

# model --------------------------------------------------------------------

# mixed-effects models

res_lme <-
  my_data |>
  nest(.by = c(
    Gene, target_id, Product,
    Description, InterPro, Pfam,
    starts_with("PGP")
  )) |>
  mutate(
    mod =
      map(
        data,
        .progress = T,
        \(.x) lme(
          lnRatio ~ flhc * media,
          weights = varComb(
            varIdent(form = ~ 1 | flhc),
            varFixed(~tech_var)
          ),
          random = ~ 1 | strain,
          data = .x
        )
      )
  )

# for quality check

augmented <-
  res_lme |>
  mutate(augmented = map2(
    mod, data,
    \(.x, .y) broom.mixed::augment(.x, data = .y)
  )) |>
  unnest(augmented)

qqnorm(augmented$.resid)

ggplot(augmented) +
  geom_boxplot(aes(.resid, str_c(strain, media)))

# extract model summaries

tidied <-
  res_lme |>
  mutate(tidied = map(
    mod,
    \(.x) tidy(.x, effects = "fixed")
  )) |>
  unnest(tidied) |>
  filter(!term %in% "(Intercept)") |>
  mutate(adjust_Q(p.value, 0.05),
    term = case_match(
      term,
      "flhcflhC-" ~ "flhC+ vs flhC-",
      "mediaLj+Ri" ~ "Lj+Ri vs Lj",
      "flhcflhC-:mediaLj+Ri" ~ "Lj+Ri: flhC+ vs flhC-"
    ) |>
      factor(levels = c("flhC+ vs flhC-", "Lj+Ri vs Lj", "Lj+Ri: flhC+ vs flhC-")),
    side = case_when(
      term == "flhC+ vs flhC-" & statistic > 0 ~ "flhC-",
      term == "Lj+Ri vs Lj" & statistic > 0 ~ "Lj+Ri",
      term == "Lj+Ri: flhC+ vs flhC-" & statistic > 0 ~ "flhC- Lj+Ri",
      term == "flhC+ vs flhC-" & statistic < 0 ~ "flhC+",
      term == "Lj+Ri vs Lj" & statistic < 0 ~ "Lj",
      term == "Lj+Ri: flhC+ vs flhC-" & statistic < 0 ~ "flhC+ Lj+Ri"
    ) |>
      factor(levels = c(
        "flhC+", "flhC-", "Lj", "Lj+Ri",
        "flhC+ Lj+Ri", "flhC- Lj+Ri"
      ))
  )

tidied |>
  janitor::tabyl(side, significant)

tidied |>
  select(Gene, Product, term, estimate, std.error, df, statistic, p.value, q.value, s.value, lfdr, significant, side) |>
  write_csv("SupplementaryData1.csv")

# enrichment --------------------------------------------------------------

dataset_enrich <-
  tidied |>
  drop_na(PGP) |>
  select(Gene, side, term, starts_with("PGP")) |>
  pivot_longer(
    cols = c(starts_with("PGP")),
    values_to = "category",
    names_to = "level"
  ) |>
  mutate(category = str_to_lower(category) |>
    str_replace_all("\\|", ", ") |>
    str_replace_all("_", " ") |>
    str_replace("prtoeins", "proteins") |>
    str_replace_all("plant vitamin|root vitamin", "vitamin")) |>
  separate_rows(category, sep = "; ") |>
  filter(
    str_detect(category, "putative", negate = T),
    level %in% c("PGP2")
  ) |>
  distinct() |>
  count(category, side, level, term) |>
  mutate(side_generic = case_when(
    side %in% c("Lj", "flhC+", "flhC+ Lj+Ri") ~ "one",
    .default = "two"
  )) |>
  select(side, side_generic, category, n, level, term) |>
  pivot_wider(names_from = side_generic, values_from = c(
    n,
    side
  )) |>
  replace_na(list(
    `n_one` = 0,
    `n_two` = 0
  )) |>
  mutate(
    total = `n_one` + `n_two`,
    rel_one = n_one / total,
    rel_two = n_two / total,
    difference = abs(`rel_one` - `rel_two`)
  )

background_short |>
  filter(level %in% "PGP2") |>
  count(category)


dataset_enrich |>
  write_csv("SupplementaryData2.csv")

# pathway analysis --------------------------------------------------------

pmr_model <-
  tidied |>
  mutate(across(starts_with("PGP"), ~ (.x) |>
    str_replace("prtoeins", "proteins") |>
    str_replace_all("plant vitamin|root vitamin", "vitamin"))) |>
  separate_rows(PGP5, sep = "; ") |>
  drop_na(PGP5) |>
  filter(n_distinct(Gene) >= 3, .by = PGP5) |>
  mutate(n_of_genes = n(), .by = c(PGP5, term)) |>
  nest(.by = c(PGP5, n_of_genes)) |>
  mutate(
    mod =
      map(
        data,
        .progress = T,
        \(.x) gls(
          estimate ~ 0 + term,
          weights = ~ std.error^2,
          data = .x
        )
      )
  )



pmr_tidied <-
  pmr_model |>
  mutate(tidied = map(
    mod,
    \(.x) tidy(.x, effects = "fixed")
  )) |>
  unnest(tidied) |>
  hoist(data, "Gene") |>
  select(-c(data, mod)) |>
  filter(!term %in% "(Intercept)") |>
  mutate(adjust_Q(p.value, 0.05),
    term = str_remove(term, "term") |>
      factor(levels = c("flhC+ vs flhC-", "Lj+Ri vs Lj", "Lj+Ri: flhC+ vs flhC-")),
    side = case_when(
      term == "flhC+ vs flhC-" & statistic > 0 ~ "flhC-",
      term == "Lj+Ri vs Lj" & statistic > 0 ~ "Lj+Ri",
      term == "Lj+Ri: flhC+ vs flhC-" & statistic > 0 ~ "flhC- Lj+Ri",
      term == "flhC+ vs flhC-" & statistic < 0 ~ "flhC+",
      term == "Lj+Ri vs Lj" & statistic < 0 ~ "Lj",
      term == "Lj+Ri: flhC+ vs flhC-" & statistic < 0 ~ "flhC+ Lj+Ri"
    ) |>
      factor(levels = c(
        "flhC+", "flhC-", "Lj", "Lj+Ri",
        "flhC+ Lj+Ri", "flhC- Lj+Ri"
      ))
  ) |>
  left_join(all_pgp |>
    select(PGP1, PGP2, PGP3, PGP4, PGP5) |>
    distinct())

pmr_tidied |>
  write_csv("SupplementaryData3.csv")

# plotting ----------------------------------------------------------------

(free(plot_volcano("flhC+ vs flhC-")) /
  plot_enrich("flhC+ vs flhC-")) |
  plot_pmr("flhC+ vs flhC-")

my_ggsave("expression", 8, 10)

(free(plot_volcano("Lj+Ri vs Lj")) /
  plot_enrich("Lj+Ri vs Lj")) |
  plot_pmr("Lj+Ri vs Lj")

my_ggsave("Response", 8, 10)

(free(plot_volcano("Lj+Ri: flhC+ vs flhC-")) /
  plot_enrich("Lj+Ri: flhC+ vs flhC-")) |
  plot_pmr("Lj+Ri: flhC+ vs flhC-")

my_ggsave("Regulation", 8, 10)
