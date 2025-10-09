### Gene-associations networks
## Roberto Siani
# 2024

# network --------------------------------------------------------------------

summarise(my_data, n_distinct(sample), .by = c(flhc, strain, media))

my_data |>
  slice_min(sparsity, n = 14, by = c(Gene, flhc)) |>
  filter(sum(counts >= med) >= 7, .by = c(Gene, flhc)) |>
  summarise(n_distinct(Gene), .by = flhc)

set.seed(37)

# compute networks on 14 samples, with higher detection stringency
# associations are recovered via proportionality coefficient rho

nested_propr <-
  my_data |>
  slice_min(sparsity, n = 14, by = c(Gene, flhc)) |>
  filter(sum(counts >= med) >= 7, .by = c(Gene, flhc)) |>
  nest(.by = flhc) |>
  mutate(
    corr = map(
      data,
      \(.x)
      pivot_wider(.x,
        id_cols = sample,
        values_from = lnRatio,
        names_from = Gene
      ) |>
        column_to_rownames("sample") |>
        propr(
          metric = "rho",
          ivar = NA,
          p = 999
        ) |>
        updateCutoffs(
          ncores = 7,
          custom_cutoffs = seq(0.05, 0.95, length.out = 19)
        )
    )
  )

# predicted fdr at correlation of .75

nested_propr |> pluck("corr", 1, "fdr")
nested_propr |> pluck("corr", 2, "fdr")

# filter and compute network properties

flhc_tgr <-
  nested_propr |>
  pluck("corr", 2) |>
  getResults() |>
  filter(propr >= .75) |>
  select(
    from = Pair,
    to = Partner,
    weights = propr
  ) |>
  as_tbl_graph(directed = F) |>
  mutate(
    order = graph_order(),
    size = graph_size(),
    group = group_fast_greedy(weights = abs(weights)),
    modularity = graph_modularity(group),
    degree = centrality_degree(),
    n_groups = n_distinct(group),
    efficiency = graph_efficiency()
  )

delta_tgr <-
  nested_propr |>
  pluck("corr", 1) |>
  getResults() |>
  filter(propr >= .75) |>
  select(
    from = Pair,
    to = Partner,
    weights = propr
  ) |>
  as_tbl_graph(directed = F) |>
  mutate(
    order = graph_order(),
    size = graph_size(),
    group = group_fast_greedy(weights = abs(weights)),
    modularity = graph_modularity(group),
    degree = centrality_degree(),
    n_groups = n_distinct(group),
    efficiency = graph_efficiency()
  )

# plot networks

p1 <-
  ggraph(flhc_tgr, "nicely") +
  geom_edge_link(
    aes(alpha = sqrt(weights)),
    color = pal_bac2[1]
  ) +
  geom_node_point(
    alpha = 2 / 3,
    color = pal_bac2[1]
  ) +
  geom_node_text(
    data = ~ filter(.x, name %in% c("flhC")),
    aes(label = name),
    size = 15 / .pt,
    fontface = "italic",
    color = "#555555"
  ) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "lines")
  )

p2 <- delta_tgr |>
  ggraph(layout = "nicely") +
  geom_edge_link(
    aes(alpha = sqrt(weights)),
    color = pal_bac2[2]
  ) +
  geom_node_point(
    alpha = 2 / 3,
    color = pal_bac2[2]
  ) +
  geom_node_text(
    data = ~ filter(.x, name %in% c("flhC")),
    aes(label = name),
    size = 15 / .pt,
    fontface = "italic",
    color = "#555555"
  ) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "lines")
  )

# table of graph properties

compare_graphs <-
  bind_rows(
    list(
      "flhC+" = flhc_tgr |> as_tibble(),
      "flhC-" = delta_tgr |> as_tibble()
    ),
    .id = "flhc"
  ) |>
  mutate(
    flhc = factor(flhc, levels = c("flhC+", "flhC-")),
    group = as.factor(group),
    Gene = name,
    .keep = "unused"
  ) |>
  summarise(across(where(is.numeric), mean),
    .by = flhc
  ) |>
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# plotting

p1 + p2

my_ggsave("network", w = 9, h = 4.5)

compare_graphs |>
  pivot_longer(where(is.numeric),
    names_to = "metric"
  ) |>
  mutate(flhc = str_c("*", flhc, "*")) |>
  pivot_wider(names_from = flhc) |>
  mutate(metric = case_match(
    metric,
    "order" ~ "N. of Genes",
    "size" ~ "N. of Associations",
    "modularity" ~ "Modularity",
    "degree" ~ "Avg. N. of Associations x Gene",
    "n_groups" ~ "N. of Groups",
    "efficiency" ~ "Efficiency"
  )) |>
  gt::gt() |>
  gt::opt_table_font(font = "Arial", size = 20) |>
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_column_labels()
  ) |>
  gt::fmt_number(drop_trailing_zeros = T) |>
  gt::cols_label_with(fn = \(.x) gt::md(.x)) |>
  gt::gtsave("compare_graphs2.png")

# clustering

total_graph <-
  bind_rows(
    list(
      "plus" = flhc_tgr |> as_tibble(),
      "minus" = delta_tgr |> as_tibble()
    ),
    .id = "flhc_graph"
  ) |>
  pivot_wider(
    id_cols = c(name),
    values_from = c(group, degree),
    names_from = flhc_graph,
    values_fill = 0
  ) |>
  left_join(
    my_data |>
      select(Gene, Product, Description, starts_with("PGP")) |>
      distinct(),
    join_by(name == Gene)
  )

# flhdc cluster

total_graph |>
  filter(group_plus == 2) |>
  pull(name) -> cluster_flhdc

# overrepresentation in clusters

overrepp <-
  clusterProfiler::compareCluster(
    name ~ group,
    data = flhc_tgr |> as_tibble() |>
      filter(n() > 1, .by = group),
    fun = "enricher",
    TERM2GENE = background_short |>
      filter(level %in% "PGP3") |>
      select(category, Gene)
  ) |>
  pluck("compareClusterResult")


# we are only really interested in cluster flhDC

overrepp |>
  filter(group == 2)

# how are cluster flhDC genes grouped in flhC-

total_graph |>
  filter(group_plus == 2) |>
  count(group_minus) |>
  print(n = 26)

p_load(ggrepel)

flhc_tgr |>
  filter(group == 2) |>
  ggraph("backbone") +
  geom_edge_bundle_path(
    aes(alpha = sqrt(weights)),
    color = pal_bac2[1]
  ) +
  geom_node_point(
    alpha = 2 / 3,
    color = pal_bac2[1]
  ) +
  geom_node_text(
    data = ~ filter(.x, !is.na(gene_extract(name))),
    aes(label = name),
    size = 20 / .pt,
    fontface = "italic",
    color = "#555555",
    repel = T
  ) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "lines")
  ) +
  scale_size(range = c(1, 2))

my_ggsave("sub_network", 9, 9)


my_data |>
  distinct(Gene, Description, Product, PGP4, ) |>
  filter(Gene %in% cluster_flhdc) |>
  View()
separate_longer_delim(PGP4, delim = "; ") |>
  summarise(tot = n(), genes = str_c(Gene, collapse = ", "), .by = PGP4) |>
  arrange(-tot) |>
  filter(tot > 2, !is.na(PGP4)) |>
  View()

#

nested_propr |>
  pluck("corr", 2) |>
  getResults() |>
  filter(
    propr >= .75,
    Pair %in% c("flhC", "flhD") |
      Partner %in% c("flhC", "flhD")
  ) |>
  select(
    from = Pair,
    to = Partner,
    weights = propr
  ) |>
  arrange(from)


flhc_tgr |>
  to_local_neighborhood(node = name == "flhC", order = 2) |>
  pluck(1) |>
  as_tibble() |>
  left_join(distinct(my_data, Gene, Product), by = c("name" = "Gene")) |>
  View()
