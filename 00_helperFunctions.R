# my suite of basic packages

pacman::p_load(pacman, tidyverse, patchwork,
               ggraph, tidygraph, broom,
              broom.mixed, nlme, ggrepel)

# ggplot theme for visualization

theme_set(
  theme_minimal() +
    theme(
      plot.margin = margin(3, 3, 3, 3),
      text = element_text(size = 12,
                          family = "Arial",
                          color = "#555555"),
      panel.grid = element_blank(),
      axis.line = element_line(color = "#555555", linewidth = .5),
      axis.ticks = element_line(color = "#555555", linewidth = .5),
      axis.title = element_text(hjust = 1),
      strip.text = element_text(face = "bold"),
      legend.position = "none",
      legend.title = element_blank(),
      panel.spacing.x = unit(.5, "lines"),
      panel.spacing.y = unit(.5, "lines"),
    )
)

# read HMMER tblout, modified from https://github.com/arendsee/rhmmer

read_tblout <- function(file){

  col_types <-
    readr::cols(
      query_name         = readr::col_character(),
      query_accession    = readr::col_character(),
      domain_name          = readr::col_character(),
      domain_accession     = readr::col_character(),
      sequence_evalue     = readr::col_double(),
      sequence_score      = readr::col_double(),
      sequence_bias       = readr::col_double(),
      best_domain_evalue  = readr::col_double(),
      best_domain_score   = readr::col_double(),
      best_domain_bis     = readr::col_double(),
      domain_number_exp   = readr::col_double(),
      domain_number_reg   = readr::col_integer(),
      domain_number_clu   = readr::col_integer(),
      domain_number_ov    = readr::col_integer(),
      domain_number_env   = readr::col_integer(),
      domain_number_dom   = readr::col_integer(),
      domain_number_rep   = readr::col_integer(),
      domain_number_inc   = readr::col_character()
    )

  N <- length(col_types$cols)

  # the line delimiter should always be just "\n", even on Windows
  lines <- readr::read_lines(file, lazy=FALSE, progress=FALSE)

  table <- sub(
    pattern = sprintf("(%s).*", paste0(rep('\\S+', N), collapse=" +")),
    replacement = '\\1',
    x=lines,
    perl = TRUE
  ) %>%
    gsub(pattern="  *", replacement="\t") %>%
    paste0(collapse="\n") %>%
    readr::read_tsv(
      col_names=names(col_types$cols),
      comment='#',
      na='-',
      col_types = col_types,
      lazy=FALSE,
      progress=FALSE
    )

  descriptions <- lines[!grepl("^#", lines, perl=TRUE)] %>%
    sub(
      pattern = sprintf("%s *(.*)", paste0(rep('\\S+', N), collapse=" +")),
      replacement = '\\1',
      perl = TRUE
    )

  table$description <- descriptions[!grepl(" *#", descriptions, perl=TRUE)]

  table
}

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

# plot sing value decomposition results

plot_svd_nice =
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
          shape = str_c(flhc,mutant),
          color = flhc_media,
          fill = after_scale(color),
        ),
        stroke = 1,
        size = 2,
        show.legend = T
      ) +
      scale_shape_manual(values = c(2, 24, 1, 21)) +
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
      # coord_fixed(ratio = str_extract(names(my.data)[PCx + 1], "0.[0-9]*$") |> as.numeric() /
      #               str_extract(names(my.data)[PCx + 1], "0.[0-9]*$") |> as.numeric()) +
      theme(text = element_text(size = 20))
  }

plot_svd_nice2 =
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
          shape = str_c(flhc,mutant),
          color = flhc_media,
          fill = after_scale(color),
        ),
        size = 2,
        stroke = 1,
        show.legend = T
      ) +
      ggalt::geom_encircle(
        aes(fill = flhc_media,
            color = flhc_media),
        alpha = .1,
        expand = 0.01,
        show.legend = F) +
      scale_shape_manual(values = c(24, 21)) +
      scale_color_manual(
        values = pal_bac4,
        aesthetics = c("color", "fill")) +
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
      # coord_fixed(ratio = str_extract(names(my.data)[PCx + 1], "0.[0-9]*$") |> as.numeric() /
      #               str_extract(names(my.data)[PCx + 1], "0.[0-9]*$") |> as.numeric()) +
      theme(text = element_text(size = 20))
  }

# adjust p value

fdr = \(.x, fdr_level) tibble(
  FDR = p.adjust(.x, "fdr"),
  significant = FDR <= fdr_level,
  s.value = -log2(.x))


adjust_Q = function(.x, fdr_level) {
  q = qvalue::qvalue(.x, fdr.level = fdr_level)
  tibble(
    q.value = q$qvalues,
    s.value = -log2(q.value),
    lfdr = q$lfdr,
    significant = q$significant)
}

# read annotation from interpro

read_interpro =
  function(filename) {
    read_tsv(filename,
             col_names = c("LR", "MD5", "Length", "Analysis",
                           "Signature_accession", "Signature_description",
                           "Start", "End", "score",
                           "Status", "Date", "IP_accession", "IP_description"),
             na = "-") |>
      group_by(LR) |>
      summarise(
        InterPro = str_c(str_unique(paste(IP_accession,
                                          IP_description,
                                          sep = ": ")),
                         collapse = "; ") |>
          str_remove_all("NA: NA(; )?") |>
          na_if(""),
        Full = str_c(str_unique(paste(Signature_accession,
                                      Signature_description,
                                      sep = ": ")),
                     collapse = "; ") |>
          str_remove_all("NA: NA(; )?") |>
          na_if(""))
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

# palette for RNAseq experiment

pal_bac4 =
  c("flhC+ Lj" = "#FC7D0B",
    "flhC+ Lj+Ri" = "#FFBC79",
    "flhC- Lj" = "#1170AA",
    "flhC- Lj+Ri" = "#A3CCE9")

pal_bac2 = c("flhC+" = "#FC7D0B", "flhC-" = "#1170AA")

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
            "ns" = "grey75")

# summarize fastp reports

summary_to_df = function(x) {
  path = x
  tmp = jsonlite::read_json(path = path) |> pluck("summary")
  list_rbind(tmp[3:4] |> map(as.data.frame), names_to = "step") |>
    mutate(sample_lane = {{x}} |>
             str_remove("data/fastp_[0-9]{4}/") |>
             str_remove("_R1_001.fastq.json"))
}

centered_log_ratio = function(.x, .pseudo) {
  local_min = min(.x[.x != 0])
  pseudo = if_else(.x == 0, local_min * .pseudo, .x)
  log(pseudo) - mean(log(pseudo))
}


enricher = clusterProfiler::enricher

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
    keep.fraction = 0.05,
    max.overlaps = 50, force_pull = 0.5) +
  ggpp::stat_group_counts(size = 10/.pt,
                          aes(label = sprintf("%i",
                                              after_stat(count))),
                          parse = T) +
  scale_color_manual(values = pal_bac7)}

plot_pmr =
  function(.term) {
  pmr_tidied |>
  filter(significant, term %in% {{.term}}) |>
  ggplot(aes(y = PGP4,
             x = estimate,
             alpha = if_else(significant, I(1), I(1/3)),
             color = if_else(significant, side, "ns"))) +
  geom_point(shape = 15, size = 3) +
  scale_color_manual(values = pal_bac7) +
  facet_grid(PGP3 ~ ., scales = "free", space = "free",
             labeller = labeller(PGP3 = label_wrap_gen(30))) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_line(linetype = 3,
                                          color = "grey69",
                                          linewidth = 0.2),
        strip.clip = "off")
    }

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
        axis.title.x = element_blank(),
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
  facet_grid(level ~ ., scales = "free_y", space = "free") +
  scale_alpha(range = c(1/4, 3/4))
  }


# testing the variance parameters (sequencing depth and technical variance)

# we use a subset of genes for the testing to speed up

# gene_ten = sample(unique(my_data$Gene), size = 2)
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
