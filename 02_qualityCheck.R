# quality check after fastp

source("scripts/00_helperFunctions.R")

# load fastp data

fastp_data <-
  fs::dir_ls("data/fastp_2407/", glob = "*.json") |>
  map(function(x) summary_to_df(x)) |>
  list_rbind() |>
  separate_wider_delim(sample_lane, names = c("sample", "lane"), delim = "_L00") |>
  mutate(sample = str_remove(sample, "_S[0-9]*")) |>
  group_by(sample, step) |>
  summarise(across(total_reads:total_bases, sum), across(q20_bases:gc_content, mean))

# write_csv(fastp_data, "fastp_data.csv")

fastp_data |>
  ungroup() |>
  summarise(across(where(is.numeric), median), .by = step)

fastp_data |>
  ggplot() +
  geom_col(aes(y = sample, x = total_reads, fill = step),
    position = "dodge"
  )

# quality check after kallisto and kmer size comparison: 31 is best at the end!

kallisto_check <-
  fs::dir_ls("data/", regexp = "kall.*.json", recurse = 2) |>
  map(~ jsonlite::read_json(.x) |> as.data.frame()) |>
  list_rbind(names_to = "path") |>
  mutate(
    sample = str_remove_all(
      path, "data/kall[0-9]*/|_S[0-9]*/run_info.json"
    ),
    kmer_size = str_extract(path, "kall[0-9]{2}") |>
      str_extract("[0-9]{2}")
  )

# write_csv(kallisto_check |>
#             filter(kmer_size == 31), "kallisto_check.csv")

kallisto_check |>
  filter(str_detect(sample, "data", negate = T)) |>
  select(sample, where(is.numeric), kmer_size) |>
  pivot_wider(names_from = kmer_size, values_from = where(is.numeric))


summarise(kallisto_check |> filter(str_detect(sample, "x", negate = T)),
  across(where(is.numeric), ~ median(.x)),
  .by = kmer_size
)

kallisto_check |>
  ggplot(aes(kmer_size, p_pseudoaligned, group = sample)) +
  geom_point() +
  geom_path() +
  theme(panel.grid.major.y = element_line(
    colour = "#555555",
    linewidth = 0.5,
    linetyp = 3
  )) +
  scale_y_continuous(n.breaks = 6)
