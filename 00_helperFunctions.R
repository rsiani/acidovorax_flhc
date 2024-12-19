### Helper functions and other stuff
## Roberto Siani
# 2024

# # interproscan annotation
#
# ~/00_lib/interproscan-5.64-96.0/interproscan.sh -cpu 7 -i LR140.faa -b pfamLR140 -dra -appl Pfam
# ~/00_lib/interproscan-5.64-96.0/interproscan.sh -cpu 7 -i LR124.faa -b pfamLR124 -dra -appl Pfam
#
# # reciprocal best hit
#
# mmseqs easy-rbh LR140.ffn LR124.ffn -s 7.5 --search-type 3 rbh.tsv tmp --format-output query,target,evalue,pident,alnlen
# mmseqs easy-rbh LR140.ffn LR124.ffn -s 7.5 --search-type 3 rbh.tsv tmp --format-output query,target,evalue,pident,alnlen --translation-table 11



# load libraries

pacman::p_load(pacman, tidyverse, patchwork,
               ggraph, tidygraph, broom,
              broom.mixed, nlme, ggrepel, propr)

# ggplot theme for visualization

theme_set(
  theme_minimal() +
    theme(
      plot.margin = margin(3, 3, 3, 3),
      text = element_text(size = 15,
                          family = "Fira Sans",
                          color = "#555555"),
      panel.grid = element_blank(),
      axis.line = element_line(color = "#555555", linewidth = .5),
      axis.ticks = element_line(color = "#555555", linewidth = .5),
      axis.title = element_text(hjust = 1),
      legend.position = "none",
      legend.title = element_blank(),
      panel.spacing.x = unit(.5, "lines"),
      panel.spacing.y = unit(.5, "lines"),
    )
)

# extract gene ID

gene_extract = \(.x) str_extract(.x,
                                 "\\p{lower}{3}\\d?\\p{upper}\\b")

# preconfigured figure saving

my_ggsave =
  function(.x, w, h) {
    ggsave(filename = str_c("figures/", .x, ".svg", sep = ""),
           width = w,
           height = h,
           device = svglite::svglite,
           units = "in",
    )
  }

# plot single value decomposition

plot_svd =
  function(my.data, PCx, PCy){
    ggplot(
      data =
        my.data,
      aes(
        x = .data[[names(my.data)[PCx + 1]]],
        y = .data[[names(my.data)[PCy + 1]]]
      )) +
      geom_point(
        aes(
          shape = strain,
          color = flhc_media,
          fill = after_scale(color),
        ),
        stroke = 1,
        size = 2,
        show.legend = T
      ) +
      scale_shape_manual(values = pal_shape) +
      scale_color_manual(
        values = pal_bac4,
        aesthetics = c("color", "fill")) +
      ggalt::geom_encircle(
        aes(fill = flhc_media,
            color = flhc_media),
        alpha = .1,
        expand = 0.01,
        show.legend = F) +
      ggside::geom_xsidedensity(aes(fill = flhc_media),
                                alpha = 1/2,
                                color = NA,
                                show.legend = F) +
      ggside::geom_ysidedensity(aes(fill = flhc_media),
                                alpha = 1/2,
                                color = NA,
                                show.legend = F) +
      ggside::scale_xsidey_continuous(labels = NULL) +
      ggside::scale_ysidex_continuous(labels = NULL) +
      theme(text = element_text(size = 20))
  }

# simple fdr adjustment

fdr = \(.x, fdr_level) tibble(
  FDR = p.adjust(.x, "fdr"),
  significant = FDR <= fdr_level,
  s.value = -log2(.x))

# lfdr adjustment for gene level comparisons

adjust_Q = function(.x, fdr_level) {
  q = qvalue::qvalue(.x, fdr.level = fdr_level)
  tibble(
    q.value = q$qvalues,
    s.value = -log2(q.value),
    lfdr = q$lfdr,
    significant = q$significant)
}

# read pfam interpro annotation

read_pfam =
  function(filename) {
    read_tsv(filename,
             col_names = c("target_id", "MD5", "Length", "Analysis",
                           "Signature_accession", "Signature_description",
                           "Start", "End", "score",
                           "Status", "Date", "IP_accession", "IP_description"),
             na = "-") |>
      group_by(target_id) |>
      summarise(
        InterPro = str_c(str_unique(IP_accession, ignore_case = T),
                         collapse = "; ") |>
          str_remove_all("NA(; )?") |>
          na_if(""),
        Pfam = str_c(str_unique(Signature_accession, ignore_case = T),
                     collapse = "; ") |>
          str_remove_all("NA(; )?") |>
          na_if(""),
        Description = str_c(str_unique(Signature_description, ignore_case = T),
                            collapse = "; ") |>
          str_remove_all("NA(; )?") |>
          na_if(""))
    # PREPARADO = any(PREPARADO),
    # PA = any(PA))
  }

# read reciprocal best hits from mmseqs

rbh =
  read_tsv("data/bakta/rbh.tsv",
           col_names = c("LR140", "LR124"),
           col_select = 1:2) |>
  full_join(read_tsv("data/bakta/LR140.tsv",
                     col_select = 6,
                     skip = 5),
            by = c("LR140" = "Locus Tag"),
            na_matches = "never") |>
  full_join(read_tsv("data/bakta/LR124.tsv",
                     col_select = 6,
                     skip = 5),
            by = c("LR124" = "Locus Tag"),
            na_matches = "never") |>
  filter(!(is.na(LR124) & is.na(LR140))) |>
  mutate(IDX = row_number())

# read bakta annotation

read_annotation =
  function(ID) {
    read_tsv(str_c("data/bakta/", {{ID}}, ".tsv"), skip = 5) |>
      select("target_id" = `Locus Tag`, Type, Gene, Product) |>
      filter(!is.na(target_id)) |>
      left_join(read_tsv(str_c("data/bakta/", {{ID}}, "_pgp.csv"),
                         col_names = c("PGP", "target_id"),
                         col_select = 1:2)) |>
      left_join(read_tsv(str_c("data/bakta/", {{ID}}, "_pgp2.txt"),
                         col_names = str_c("PGP", 1:6),
                         col_select = 1:6,
                         skip = 1) |>
                  separate_wider_delim(cols = PGP6, names = c("PGP6", "PGP"), delim = "->")) |>
      mutate(across(starts_with("PGP"), ~ str_to_lower(.x) |>
                      str_replace_all("\\|", ", ") |>
                      str_replace_all("_", " "))) |>
      left_join(read_pfam(str_c("data/bakta/pfam", {{ID}}, ".tsv"))) |>
      summarise(
        across(where(is.character), \(.x) str_c(str_unique(.x, ignore_case = T), collapse = "; ")),
        .by = c(1:4))
  }

# calculate sparsity of a vector

sparsity = \(.x) sum(.x == 0, na.rm = T)/length(.x)

# color and shape palettes

pal_bac2 = c("flhC+" = "#FC7D0B", "flhC-" = "#1170AA")

pal_bac4 =
  c("flhC+ Lj" = "#FC7D0B",
    "flhC+ Lj+Ri" = "#FFBC79",
    "flhC- Lj" = "#1170AA",
    "flhC- Lj+Ri" = "#A3CCE9")

pal_bac7 = c("Lj"= "#8BC34A",
             "Lj+Ri"= "#E1BEE7",
             "flhC+" = "#FC7D0B",
             "flhC+ Lj+Ri" = "#FFBC79",
             "flhC-" = "#1170AA",
             "flhC- Lj+Ri" = "#A3CCE9",
             "ns" = "grey75")

pal_growth7 = c("control Lj"= "#8BC34A",
            "control Lj+Ri"= "#E1BEE7",
            "flhC+ Lj" = "#FC7D0B",
            "flhC+ Lj+Ri" = "#FFBC79",
            "flhC- Lj" = "#1170AA",
            "flhC- Lj+Ri" = "#A3CCE9",
            "ns" = "grey75",
            "Ri" = "#E1BEE7",
            "flhC+" =  "#FC7D0B",
            "flhC+ & Ri" = "#FFBC79",
            "flhC-" = "#1170AA",
            "flhC- & Ri" =  "#A3CCE9")

pal_shape =
  c("LR140" = 21,
    "Comp_140" = 1,
    "LR124" = 24,
    "Delta_140" = 2,
    "media" = 15,
    "control" = 15)

# summarize fastp reports

summary_to_df = function(x) {
  path = x
  tmp = jsonlite::read_json(path = path) |> pluck("summary")
  list_rbind(tmp[3:4] |> map(as.data.frame), names_to = "step") |>
    mutate(sample_lane = {{x}} |>
             str_remove("data/fastp_[0-9]{4}/") |>
             str_remove("_R1_001.fastq.json"))
}

# centered log ratio transformation

centered_log_ratio = function(.x, .pseudo) {
  local_min = min(.x[.x != 0])
  pseudo = if_else(.x == 0, local_min * .pseudo, .x)
  log(pseudo) - mean(log(pseudo))
}


# for cluster profiler

enricher = clusterProfiler::enricher

# for PGP5 pathway analysis

all_pgp =
  bind_rows(read_tsv("data/bakta/LR124_pgp2.txt",
                     col_names = str_c("PGP", 1:6),
                     col_select = 1:6,
                     skip = 1) |>
              separate_wider_delim(cols = PGP6, names = c("PGP6", "PGP"), delim = "->"),
            read_tsv("data/bakta/LR124_pgp2.txt",
                     col_names = str_c("PGP", 1:6),
                     col_select = 1:6,
                     skip = 1) |>
              separate_wider_delim(cols = PGP6, names = c("PGP6", "PGP"), delim = "->")) |>
  mutate(across(starts_with("PGP"), ~ str_to_lower(.x) |>
                  str_replace_all("\\|", ", ") |>
                  str_replace_all("_", " ") |>
                  str_replace("prtoeins", "proteins") |>
                  str_replace_all("plant vitamin|root vitamin", "vitamin")))

# volcano plot

plot_volcano =
  function(.term) {
  tidied |>
  filter(term %in% {{.term}}) |>
  mutate(group = if_else(significant, side, "ns") |>
           as.factor(),
         estimate = estimate,
         s.value = s.value) |>
  ggplot(aes(x = estimate, y = s.value, color = group)) +
  ggpp::geom_quadrant_lines(yintercept = -log2(0.05),
                            linetype = 3,
                            linewidth = 0.5,
                            color = "#555555") +
  geom_point(
    size = 1,
    alpha = .8,
    shape = 15) +
  ggpp::stat_dens2d_labels(
    data = ~ filter(.x,
                    significant,
                    !is.na(gene_extract(Gene))
    ),
    aes(label = Gene),
    geom = "text_repel",
    color = "#555555",
    fontface = "italic",
    keep.fraction = 0.05,
    max.overlaps = 50, force_pull = 0.5) +
  ggpp::stat_group_counts(size = 10/.pt,
                          aes(label = sprintf("%i",
                                              after_stat(count))),
                          parse = T) +
      scale_color_manual(values = pal_bac7,
                         label = list("flhC+" = "_flhC_+",
                                      "flhC-" = "_flhC_-",
                                      "flhC+ Lj+Ri" = "_flhC_+ _Lj_+_Ri_",
                                      "flhC- Lj+Ri" = "_flhC_- _Lj_+_Ri_",
                                      "Lj" = "_Lj_",
                                      "Lj+Ri" = "_Lj_+_Ri_")) +
      theme(legend.position = "left",
            legend.text = ggtext::element_markdown(size = 10)) +
      guides(colour = guide_legend(override.aes = list(size = 5))) +
      scale_x_continuous(name = "Log Fold Change")
  }

# pathway analysis

plot_pmr =
  function(.term) {
  pmr_tidied |>
  filter(significant, term %in% {{.term}}) |>
  ggplot(aes(y = PGP4,
             x = estimate,
             color = side)) +
  geom_point(size = 3,
             shape = 15,
             alpha = .8) +
  scale_color_manual(values = pal_bac7) +
  facet_grid(PGP3 ~ ., scales = "free", space = "free",
             labeller = labeller(PGP3 = label_wrap_gen(30))) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_line(linetype = 3,
                                          color = "grey69",
                                          linewidth = 0.2),
        strip.clip = "off") +
      scale_x_continuous(name = "Log Fold Change")
  }

# representation of terms

plot_enrich  =
  function(.term) {
  dataset_enrich |>
  filter(term %in% {{.term}}, level %in% c("PGP2")) |>
  slice_max(diff, n = 10, by = term) |>
  pivot_longer(contains("_"),
               names_to = c(".value", "side_generic"),
               names_sep = "_") |>
  ggplot() +
  geom_segment(data = ~
                 pivot_wider(.x,
                             id_cols = c(category, level, term, diff),
                             names_from = side_generic, values_from = n),
               aes(x = one,
                   xend = two,
                   color = diff,
                   y = category,
                   alpha = diff,
                   yend = category),
               linewidth = 3, color = "#a1a1a1",
               show.legend = F) +
  geom_point(aes(x = n, y = category, color = side),
             size = 3) +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x = element_blank(),
        panel.grid.major.x = element_line(linewidth = .5,
                                          linetype = 3,
                                          color = "#555555"),
        panel.spacing.y = unit(0, "lines"),
        panel.spacing.x = unit(0, "lines"),
        panel.border = element_rect(linewidth = .5,
                                    linetype = 1,
                                    color = "#555555",
                                    fill = NA),
        strip.text.y = element_text(angle = -90, hjust = 0),
        strip.clip = "off") +
  scale_color_manual(values = pal_bac7) +
  scale_y_discrete(labels = \(.x) str_wrap(.x, width = 30)) +
  scale_alpha(range = c(1/4, 3/4)) +
  scale_x_continuous(name = "N. of differentially expressed genes")
  }

# model testing for the variance parameters

# gene_ten = sample(unique(my_data$Gene), size = 10)
#
# model = lme(lnRatio ~ flhc * media + Gene,
#              weights = ~ tech_var,
#              random = ~ 1 | strain,
#              data = my_data |>
#                filter(Gene %in% gene_ten),
#              method = "ML")
#
# model1 = lme(lnRatio ~ flhc * media + Gene,
#              weights = varComb(varIdent(form = ~ 1 | flhc * media),
#                                varFixed(~ tech_var)),
#              random = ~ 1 | strain,
#              data = my_data |>
#                filter(Gene %in% gene_ten),
#              method = "ML")
#
# model2 = lme(lnRatio ~ flhc * media + Gene,
#              weights = varComb(varIdent(form = ~ 1 | flhc),
#                                varFixed(~ tech_var)),
#              random = ~ 1 | strain,
#              data = my_data |>
#                filter(Gene %in% gene_ten),
#              method = "ML")
#
# model3 = lme(lnRatio ~ flhc * media + Gene,
#              weights = varComb(varIdent(form = ~ 1 | media),
#                                varFixed(~ tech_var)),
#              random = ~ 1 | strain,
#              data = my_data |>
#                filter(Gene %in% gene_ten),
#              method = "ML")
#
#
# model4 = lme(lnRatio ~ flhc * media + Gene,
#              weights = varComb(varIdent(form = ~ 1 | strain),
#                                varFixed(~ tech_var)),
#              random = ~ 1 | strain,
#              data = my_data |>
#                filter(Gene %in% gene_ten),
#              method = "ML")
#
# anova(model, model1, model2,model3, model4)
#
# plot(model, paste(flhc, media) ~ resid(., scaled = T))
# plot(model1, paste(flhc, media) ~ resid(., scaled = T))
# plot(model2, paste(flhc, media) ~ resid(., scaled = T))
# plot(model3, paste(flhc, media) ~ resid(., scaled = T))
# plot(model4, paste(flhc, media) ~ resid(., scaled = T))
#
# qqnorm(resid(model))
# qqnorm(resid(model1))
# qqnorm(resid(model2))
# qqnorm(resid(model3))
# qqnorm(resid(model4))
#
