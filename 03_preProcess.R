### Pre-processing of read count estimates
## Roberto Siani
# 2024

source("scripts/00_helperFunctions.R")

# LOADING -----------------------------------------------------------------

metadata =
  read_tsv("data/metadata.csv",
           name_repair = "universal") |>
  transmute(
    sample = str_c("R", sample) |> as.factor(),
    strain = factor(strain, levels = c("LR124", "Delta_140", "LR140", "Comp_140")),
    media = case_match(media, "Treatment" ~ "Lj+Ri", "Control" ~ "Lj") |>
      factor(levels = c("Lj", "Lj+Ri")),
    flhc = case_match(
      strain,
      "LR124" ~ "flhC-",
      "LR140" ~ "flhC+",
      "Delta_140" ~ "flhC-",
      "Comp_140" ~ "flhC+") |>
      factor(levels = c("flhC+", "flhC-")),
    flhc_media = paste(flhc, media) |> factor(
      levels = c("flhC+ Lj", "flhC+ Lj+Ri", "flhC- Lj", "flhC- Lj+Ri")),
    mutant = case_match(strain,
                        "LR140" ~ "Wildtype",
                        "LR124" ~ "Wildtype",
                        .default = "Mutant") |>
      factor(levels = c("Wildtype", "Mutant")))

# write_csv(metadata, "metadata.csv")

# palette

metadata |>
  select(flhc, media, strain, mutant) |>
  distinct() |>
  drop_na(flhc) |>
  add_row(flhc = c("control", "control"),
          media = c("Lj", "Lj+Ri"),
          strain = c("media", "media"),
          mutant=  c("media", "media")) |>
  mutate(strain = factor(strain, levels = c("LR140",
                                            "Comp_140",
                                            "LR124",
                                            "Delta_140",
                                            "media")),
         mutant = factor(mutant, levels = c("Wildtype", "Mutant", "media")),
         media = factor(media, levels = c("Lj+Ri", "Lj"))) |>
  ggplot() +
  geom_point(aes(y = strain, x = media, color = str_c(flhc, media, sep = " "),
                 fill = after_scale(color),
                 shape = strain),
             size = 5) +
  scale_color_manual(values = pal_growth7) +
  scale_shape_manual(values = pal_shape) +
  facet_grid(mutant ~ ., scales = "free_y") +
  theme(axis.text.y = ggtext::element_markdown(),
        axis.title = element_blank()) +
  scale_y_discrete(label = list("Comp_140" = "LR140<sup>_ΔflhC;ΔflhC_</sup>",
                                "Delta_140" = "LR140<sup>_ΔflhC_</sup>"))

my_ggsave("legend", 2, 3)

# load annotations (Bakta/Pfam/PLABase)

LR140 =
  read_annotation("LR140") |>
  left_join(rbh |> select(-LR124), join_by(target_id == LR140))

LR124 =
  read_annotation("LR124") |>
  left_join(rbh |> select(-LR140), join_by(target_id == LR124))

# combine annotations (mmseqs rbh). I also complement gene names and description were missing with PFAMs ID and create some shorter and unique labels for the different genes

background =
  bind_rows(LR124, LR140) |>
  mutate(
    Product = na_if(Product,
                    "hypothetical protein") |>
    coalesce(Description) |>
      coalesce(paste0("hp.", ifelse(is.na(IDX), row_number(),
                                    IDX))),
    Gene =
      coalesce(Gene,
               case_when(
                 !is.na(gene_extract(Product)) ~
                   gene_extract(Product),
                 !is.na(gene_extract(Description)) ~
                   gene_extract(Description),
                 str_count(Product, "\\w+") == 1 ~ Product,
                 str_count(Description, "\\w+") == 1 ~ Description,
                 .default = str_remove_all(Product, "[:punct:]") |>
                   str_to_lower() |>
                   abbreviate())),
    descr = str_c(Product, Description) |> str_to_lower()) |>
  arrange(IDX)

write_csv(background, "SupplementaryData4.csv")

# 5 gene difference between LR124 and LR140

uniques =
  rbh |>
  filter(is.na(LR140) | is.na(LR124)) |>
  pivot_longer(1:2) |>
  drop_na(value) |>
  pull(value)

filter(background, target_id %in% uniques) |>
  select(target_id, Gene, Product, InterPro, Pfam, Description) |>
  write_csv("supp_table1.csv")

# to preserve only coding sequences

cds =
  background |>
  filter(Type %in% c("cds"),
         str_detect(descr,
                    "ibosom|trna|rrna|transfer-messenger|elongation|initiation",
                    negate = T),
         ) |>
  pull(Gene) |>
  unique()

# load kallisto transcript counts data and merge everything in one big ol' table

files = fs::dir_ls("data/kall31/", regexp = "abundance.h5", recurse = 2)

# import data, remove non CDS and negative samples

raw_data =
  files |>
  set_names(str_remove_all(files, "data/kall[0-9]*/|_S[0-9]*/abundance.h5")) |>
  map(~ tximport::tximport(.x,
                           type = "kallisto",
                           tx2gene = background |> select(target_id, Gene)) |>
        data.frame() |>
        rownames_to_column("Gene")) |>
  list_rbind(names_to = "sample") |>
  filter(str_detect(sample, "x", negate = T), Gene %in% cds) |>
  complete(Gene, sample) |>
  mutate(across(where(is.numeric), ~ if_else(is.na(.x), 0, ceiling(.x)))) |>
  left_join(metadata)

# check samples quality (n of reads and percentage undetected)

summarise(raw_data,
          depth = sum(counts),
          sparsity = sparsity(counts),
          .by = c(sample, flhc_media, strain)) |> View()

# filter sample with more less than 55 percent genes detected

raw_data |>
  filter(sparsity(counts) < .55, .by = sample) |>
  summarise(n = n_distinct(sample), sparsity(counts), .by = c(flhc, strain, media)) |>
  mutate(sum(n))

# finding a threshold for filtering out poorly covered genes

raw_data |>
  filter(sparsity(counts) < .55, .by = sample) |>
  mutate(
    med = median(if_else(counts == 0, NA, counts), na.rm = T), .by = sample) |>
  summarise(s = sparsity(counts), t = sum(counts >= med), .by = Gene) |>
  ggplot(aes(s, t)) +
  geom_point(shape = "square", alpha = .1, size = 5, color = "#555555") +
  geom_smooth(color = "indianred", se = F) +
  ggpp::geom_quadrant_lines(xintercept = .55, yintercept = 5,
                            color = "#555555", linetype = 3) +
  ggpp::stat_quadrant_counts(xintercept = .55, yintercept = 5,
                             color = "indianred")

# filtering genes detected  (max sparsity at .55, equal or more than median in
# at least 5 samples)

clean_data =
  raw_data |>
  filter(sparsity(counts) < .55, .by = sample) |>
  mutate(med = median(if_else(counts == 0, NA, counts), na.rm = T),
         .by = sample) |>
  filter(
    Gene %in% cds,
    sum(counts >= med) >= 5,
    sparsity(counts) < .55,
    .by = Gene) |>
  mutate(depth = log(sum(counts)),
         sparsity = sparsity(counts), .by = sample)

# sparsity is controlled to ~ 13%

summarise(raw_data, sparsity(counts)) # .429
summarise(clean_data, sparsity(counts)) # .112

# number of transcripts retained

n_distinct(raw_data$Gene) # 2621
n_distinct(clean_data$Gene) # 1605

# check samples quality. We select the best 16 per outer condition (flhc)
# leads to groups of size 3-5 with a slight imbalance towards Lj+Ri

summarise(clean_data,
          sum(counts),
          sparsity(counts),
          local_min = min(if_else(counts == 0, NA, counts), na.rm = T),
          pseudo = local_min * exp(-1),
          .by = c(sample, flhc_media, strain)) |> View()

exp(median(clean_data$depth))

# filtering and centered-log-ratio transformation
# pseudo values is exp(-1) * local minimum detected

transformed_data =
  clean_data |>
  mutate(
    across(c(counts, starts_with("inf")),
                ~ centered_log_ratio(.x, exp(-1)),
                .names = "clr_{.col}"),
         .by = sample)

# variance and depth covariates are calculated

my_data =
  transformed_data |>
  select(!contains("inf")) |>
  left_join(summarise(background,
                      across(where(is.character),
                             ~ str_unique(.x, ignore_case = T) |>
                               str_c(collapse = "; ")),
                      across(where(is.logical),
                             ~ any(.x)),
                      .by = Gene), relationship = "many-to-many") |>
  mutate(lnRatio = Rfast::rowmeans(transformed_data |>
                                    select(starts_with("clr")) |>
                                    as.matrix()),
         tech_var = Rfast::rowVars(transformed_data |>
                                     select(starts_with("clr")) |>
                                     as.matrix()))


my_data |> summarise(depth = mean(depth),
                     tech_var = mean(tech_var),
                      .by = c(sample, flhc_media)) |>
   ggplot() +
  geom_text(aes(depth, tech_var, color = flhc_media, label = sample)) +
  scale_color_manual(values = pal_bac4)

# quality checks

ggplot(my_data) +
  geom_pointrange(stat = "summary", aes(sample, lnRatio))

my_data |>
  summarise(
    vlr = var(lnRatio),
    clr = mean(lnRatio),
    .by = Gene) |>
  ggplot(aes(x = clr, y = vlr, label = Gene)) +
  geom_text() +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 10) +
  ggside::geom_xsidehistogram(bins = 100) +
  ggside::geom_ysidehistogram(bins = 100) +
  scale_color_gradient(low = "slategrey", high = "hotpink2") +
  geom_hline(yintercept = 0, color = "hotpink2", linewidth = .3) +
  geom_vline(xintercept = 0, color = "hotpink2", linewidth = .3)

# clean environment

background_short =
  my_data |>
  drop_na(PGP) |>
  select(Gene, starts_with("PGP")) |>
  pivot_longer(cols = c(starts_with("PGP")),
               values_to = "category",
               names_to = "level") |>
  mutate(category = str_to_lower(category) |>
           str_replace_all("\\|", ", ") |>
           str_replace_all("_", " ") |>
           str_replace("prtoeins", "proteins") |>
           str_replace_all("plant vitamin|root vitamin", "vitamin")) |>
  separate_rows(category, sep = "; ") |>
  filter(str_detect(category, "putative", negate = T)) |>
  distinct()

gc()

