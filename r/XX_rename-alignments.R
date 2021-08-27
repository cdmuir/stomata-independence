source("r/header.R")

aln = read_lines("alignments/Astragalus_20400_outaln")

# Odd elements are uid; even elements are alignments.
# Extract just uid and remove '>'
uid = aln[seq(1, length(aln), 2)] |>
  str_remove("^>")

# Check that all have 7 or fewer characters
checkmate::assert_true(all(sapply(uid, nchar) <= 7))

nms = id2name(uid, db = "ncbi")

aln[seq(1, 100, 2)] = purrr::imap_chr(unname(nms), ~{
  name = str_replace_all(.x$name, " ", "_")
  str_replace(aln[[.y * 2 - 1]], .x$id, name)
})

write_lines(aln, "alignments/Astragalus_20400_outaln_names")
