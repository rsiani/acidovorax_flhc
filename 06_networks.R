# network --------------------------------------------------------------------

source("scripts/preProcess.R")

p_load(propr)

summarise(my_data, n_distinct(sample), .by = c(flhc, strain))

my_data |>
  slice_min(sparsity, n = 14, by = c(Gene, flhc)) |>
  filter(sum(counts >= med) >= 7, .by = c(Gene, flhc)) |>
  summarise(n_distinct(Gene), .by = flhc)

set.seed(37)

nested_propr =
  my_data |>
  slice_min(sparsity, n = 14, by = c(Gene, flhc)) |>
  filter(sum(counts >= med) >= 7, .by = c(Gene, flhc)) |>
  nest(.by  = flhc) |>
  mutate(
    corr = map(data,
                    \(.x)
                    pivot_wider(.x, id_cols = sample,
                                values_from = lnRatio,
                                names_from = Gene) |>
                      column_to_rownames("sample") |>
                      propr(
                        method = "pcor.bshrink",
                        ivar = NA,
                        p = 99,
                        fixseed = T) |>
                      updateCutoffs(custom_cutoffs = seq(-.95, .95, by = .05),
                                    tails = "both",
                                    ncores = 7)))


nested_propr |> pluck("corr", 1, "fdr")
nested_propr |> pluck("corr", 2, "fdr")

flhc_tgr =
  nested_propr |>
  pluck("corr", 2) |>
  getResults() |>
  filter(abs(propr) >= .75) |>
  select(from = Pair,
         to = Partner,
         weights = propr) |>
  as_tbl_graph(directed = F)  |>
  mutate(
    order = graph_order(),
    size = graph_size(),
    group = group_fast_greedy(weights = abs(weights)),
    modularity = graph_modularity(group),
    degree = centrality_degree(),
    n_groups = n_distinct(group),
    efficiency = graph_efficiency())

delta_tgr =
  nested_propr |>
  pluck("corr", 1) |>
  getResults() |>
  filter(abs(propr) >= .75) |>
  select(from = Pair,
         to = Partner,
         weights = propr) |>
  as_tbl_graph(directed = F) |>
  mutate(
    order = graph_order(),
    size = graph_size(),
    group = group_fast_greedy(weights = abs(weights)),
    modularity = graph_modularity(group),
    degree = centrality_degree(),
    n_groups = n_distinct(group),
    efficiency = graph_efficiency())

p1 =
  ggraph(flhc_tgr, "nicely") +
  geom_edge_density(fill = pal_bac4[2]) +
  geom_edge_link(aes(alpha = abs(weights)),
                color = pal_bac2[1]) +
  geom_node_point(
    aes(
      size = log(degree),
      shape = if_else(degree > quantile(degree, .95), I(19), I(1))),
    alpha = 2/3,
    color = pal_bac2[1]) +
  geom_node_text(data = ~filter(.x, name %in% c("flhC", "flhD")),
                 aes(label = name),
                 size = 12/.pt,
                 color = "#555555") +
  theme_graph(plot_margin = margin(0, 0, 0, 0)) +
  theme(legend.position = "none") +
  scale_size(range = c(0.1, 1.5))

p2 = delta_tgr |>
  ggraph(layout = "nicely") +
  geom_edge_density(fill = pal_bac4[4]) +
  geom_edge_link(aes(alpha = abs(weights)),
                 color = pal_bac2[2]) +
  geom_node_point(
    aes(
      size = log(degree),
      shape = if_else(degree > quantile(degree, .95), I(19), I(1))),
    alpha = 2/3,
    color = pal_bac2[2]) +
  geom_node_text(data = ~filter(.x, name %in% c("flhC", "flhD")),
                 aes(label = name),
                 size = 12/.pt,
                 color = "#555555") +
  theme_graph(plot_margin = margin(0, 0, 0, 0)) +
  theme(legend.position = "none") +
  scale_size(range = c(0.1, 1.5))

compare_graphs =
  bind_rows(list("flhC+" = flhc_tgr |> as_tibble(),
                 "flhC-" = delta_tgr |> as_tibble()),
            .id = "flhc") |>
  mutate(flhc = factor(flhc, levels = c("flhC+", "flhC-")),
         group = as.factor(group),
         Gene = name,
         .keep = "unused") |>
  summarise(across(where(is.numeric), mean),
                                        .by = flhc) |>
  mutate(across(where(is.numeric), ~ round(.x, 2)))

p1 + p2

my_ggsave("network", w = 9, h = 4.5)

compare_graphs |>
  gt::gt() |>
  gt::opt_table_font(font = "Arial") |>
  gt::tab_style(style = gt::cell_text(weight = "bold"),
                locations = gt::cells_column_labels()) |>
  gt::gtsave("compare_graphs.png")

total_graph =
  bind_rows(list("plus" = flhc_tgr |> as_tibble(),
                 "minus" = delta_tgr |> as_tibble()),
            .id = "flhc_graph") |>
  pivot_wider(id_cols = c(name),
              values_from = c(group, degree),
              names_from = flhc_graph,
              values_fill = 0) |>
  left_join(my_data |>
              select(Gene, Product, Description, module, starts_with("PGP")) |>
              distinct(),
            join_by(name == Gene))

total_graph |>
  filter(group_plus == 2) |>
  pull(name) -> cluster_flhdc

overrepp =
  clusterProfiler::compareCluster(
      name ~ group,
      data = flhc_tgr |> as_tibble() |>
        filter(n() > 1, .by = group),
      fun = "enricher",
      TERM2GENE = background_short |>
        filter(level %in% "PGP4") |>
        select(category, Gene),
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      pAdjustMethod = "none") |> pluck("compareClusterResult") |>
  select(-c(p.adjust,qvalue)) |>
  mutate(adjust_Q(pvalue, 0.2))

overrepp |>
  filter(group == 2, significant) |> View()
