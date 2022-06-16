# Resolving names based on taxize output. Discarded all with score < 0.75 and hybrids.

source("r/header.R")
taxize_output_path = "processed-data/taxize_output.rds"

taxize_output = read_rds(taxize_output_path)

# Tibble of maximum score ----
max_score = taxize_output |>
  group_by(user_supplied_name) |>
  summarize(max_score = max(score, na.rm = TRUE), .groups = "drop")

taxize_output = taxize_output |>
  left_join(max_score, by = "user_supplied_name")

# 1. Select matched_name from those with max_score > 0.75 ----
n = taxize_output |>
  filter(max_score > 0.75, score >= max_score) |>
  mutate(n_words = str_count(matched_name, "\\S+") ) |>
  group_by(user_supplied_name) |>
  summarize(n_words = min(n_words), n_distinct = n_distinct(matched_name),
            .groups = "drop")

taxize_output = taxize_output |>
  filter(max_score > 0.75, score >= max_score) |>
  left_join(n, by = "user_supplied_name")

chosen_names = taxize_output |>
  slice_rows("user_supplied_name") |>
  by_slice(~ {choose_matched_name(.x$matched_name, .x$user_supplied_name)},
           .collate = "rows")

# Manually choose names for n_names > 1
chosen_names |> filter(n_names > 1) 

chosen_names = chosen_names |>
  mutate(
    chosen_name = str_replace(chosen_name, "Caltha palustris var. minor, Caltha palustris subsp. minor", "Caltha palustris subsp. minor"),
    
    chosen_name = str_replace(chosen_name, "Echinochloa crus-galli var. hispidula, Echinochloa crus-galli subsp. hispidula", "Echinochloa crus-galli var. hispidula"),
    
    chosen_name = str_replace(chosen_name, "Lotus corniculatus var. corniculatus, Lotus corniculatus subsp. corniculatus", "Lotus corniculatus var. corniculatus"),
    
    chosen_name = str_replace(chosen_name, "Mentha piperita subsp. citrata, Mentha piperita var. citrata", "Mentha piperita var. citrata")
  )
  
taxize_output1 = taxize_output |>
  left_join(select(chosen_names, user_supplied_name, chosen_name), 
            by = "user_supplied_name")

# 2. Select matched_name from those with max_score == 0.75 ----

taxize_output = read_rds(taxize_output_path) |>
  left_join(max_score, by = "user_supplied_name")

n = taxize_output |>
  filter(max_score <= 0.75, score >= max_score) |>
  mutate(n_words = str_count(matched_name, "\\S+")) |>
  group_by(user_supplied_name) |>
  summarize(n_words = min(n_words), n_distinct = n_distinct(matched_name),
            .groups = "drop")

taxize_output = taxize_output |>
  filter(max_score <= 0.75, score >= max_score) |>
  left_join(n, by = "user_supplied_name")

.x = taxize_output |>
  filter(user_supplied_name == "Nigella damscena")
choose_matched_name(.x$matched_name, .x$user_supplied_name)

chosen_names = taxize_output |>
  slice_rows("user_supplied_name") |>
  by_slice(~ {choose_matched_name(.x$matched_name, .x$user_supplied_name)},
           .collate = "rows")

# Manually choose names for n_names > 1 & n_words == 1 (if needed)
chosen_names |> 
  mutate(n_words = str_count(chosen_name, "\\S+")) |>
  filter(n_names > 1 | n_words == 1)

chosen_names = chosen_names |>
  mutate(
    chosen_name = str_replace(chosen_name, "Nigella damascena, Nigela damascena", "Nigella damascena"),
    chosen_name = str_replace(chosen_name, "Astragalus glycyphyllos, Astragalus glaucophyllus, Astragalus glycyphylloides", "Astragalus glycyphyllus")
  )
    
taxize_output2 = taxize_output |>
  left_join(select(chosen_names, user_supplied_name, chosen_name), 
            by = "user_supplied_name")

# 3. Combine and cleanup ----
# NOTE: Seems like in some cases, extra ',' not removed
resolved_names = bind_rows(
  taxize_output1,
  taxize_output2
) |>
  group_by(user_supplied_name) |>
  summarize(resolved_name = unique(chosen_name), n = n_distinct(chosen_name),
            .groups = "drop")

# Check that all taxa have one resolved name
resolved_names |>
  pull(n) |>
  equals(1) |>
  all() |>
  assert_true()

# remove hybrids
resolved_names = resolved_names |>
  mutate(resolved_name = remove_authority(resolved_name)) |>
  filter(!str_detect(user_supplied_name, " x")) |>
  filter(!str_detect(user_supplied_name, " ×")) |>
  filter(!str_detect(resolved_name, " x")) |>
  filter(!str_detect(resolved_name, " ×")) |>
  select(user_supplied_name, resolved_name)

write_rds(resolved_names, "processed-data/resolved_names.rds")
