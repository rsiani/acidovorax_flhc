### FlhDC-box
## Roberto Siani
# 2025
#

# -------------------------------------------------------------------------


known_targets <-
  c("flgA", "flgB", "flhB", "fliA", "fliE", "fliF", "fliL")


mod_140 <-
  mutate(LR140, across(c(Gene, Product, Description), ~ str_replace_na(.x)),
    target = str_c(Gene, Product, Description)
  )

map_dfr(known_targets, ~ filter(mod_140, str_detect(toupper(target), toupper(.x)))) |>
  pluck("target_id") |>
  unique() |>
  write_lines("possible_prom.txt")

# region = 150
# bedtools flank -i LR140.bed -g LR140.fai -l $region -r 0 -s > pre_prom.bed
# bedtools slop -i pre_prom.bed -g LR140.fai -l 0 -r 50 -s > promoter.bed
# bedtools getfasta -fi LR140.fna -bed promoter.bed -s > all_promoters.fa
# fasta-get-markov -m 2 all_promoters.fa model_promoters.mm
# cat ../../possible_prom.txt| while read line; do grep $line promoter.bed; done > class_ii.bed
# bedtools getfasta -fi LR140.fna -bed class_ii.bed -s > class_ii.fa
# awk '{print $1","$2","$3","$6","$10}' promoter.bed | sed "s/;.*//g" | sed "s/ID=//g" > fimo_features.csv
# meme class_ii.fa -dna -mod oops -nmotifs 1 -revcomp -bfile model_promoters.mm -V -oc meme_out -cons ANNANNNNNNNNNNNNNNNNNNTNWTNNNNNNWT
# fimo --bgfile model_promoters.mm -no-pgc --thresh 1e-5 --verbosity 5 meme_out/meme.txt all_promoters.fa


# analysis ----------------------------------------------------------------



# validate against ANNANNNNNNNNNNNNNNNNNNTNWTNNNNNNWT

universalmotif::read_meme("data/bakta/meme_out/meme.txt") |>
  universalmotif::view_motifs(normalise.scores = T)

my_ggsave("supp_fig1", 8, 2)

fimo_out <-
  read_tsv("data/bakta/fimo_out/fimo.tsv", comment = "#", name_repair = "universal") |>
  left_join(
    read_csv("data/bakta/fimo_features.csv", col_names = c("A", "B", "C", "D", "target_id")) |>
      mutate(
        sequence_name = str_c(A, ":", B, "-", C, "(", D, ")"),
        .keep = "unused"
      )
  ) |>
  drop_na() |>
  slice_max(score, n = 1, by = sequence_name) |>
  left_join(background |> filter(str_detect(target_id, "LR140"))) |>
  mutate(
    bindingSite =
      case_when(
        q.value <= 0.05 ~ "***",
        q.value <= 0.10 ~ "**",
        q.value <= 0.13 ~ "*",
        .default = NA
      )
  )


fimo_out |>
  ggplot() +
  geom_point(aes(p.value, score, color = bindingSite))

summarise(fimo_out, n(), .by = bindingSite)

fimo_out |>
  select(Gene, bindingSite) |>
  anti_join(my_data)

p2 <-
  my_data |>
  left_join(fimo_out |> select(Gene, bindingSite)) |>
  mutate(bindingSite = if_else(Gene %in% regulators, "Regulator", bindingSite)) |>
  drop_na(bindingSite) |>
  filter(!is.na(gene_extract(Gene))) |>
  ggplot() +
  geom_boxplot(aes(lnRatio, Gene, color = flhc_media),
    outlier.shape = 1,
    linewidth = 1 / 4
  ) +
  scale_color_manual(values = pal_growth7) +
  facet_grid(
    bindingSite ~ .,
    scales = "free",
    space = "free"
  ) +
  labs(x = "Log Ratio") +
  theme(
    axis.text.y.left = element_text(face = "italic"),
    axis.title.y = element_blank()
  )

bound <-
  fimo_out |>
  pluck("Gene")

regulators <-
  c("crp", "crcp", "ompR", "Histone", "qseB", "fur", "tetR", "iclR")


bound[bound %in% cluster_flhdc]

bound[bound %in% regulators]

fimo_out |>
  select(Gene, start, stop, strand, q.value, bindingSite, matched_sequence) |>
  mutate(q.value = format(q.value, scientific = T, digits = 2)) |>
  write_csv("fimo_for_table.csv")
