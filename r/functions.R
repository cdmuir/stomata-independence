# Functions for taxon name resolution ----

choose_matched_name <- function(matched_name, user_supplied_name) {
  
  .d <- tibble(matched_name = matched_name) %>%
    mutate(n_distinct = n_distinct(matched_name))
  
  if (any(.d$n_distinct != 1L)) {
    
    .d <- .d %>%
      mutate(n_words = str_count(matched_name, "\\S+")) %>%
      filter(n_words == min(n_words)) %>%
      mutate(n_distinct = n_distinct(matched_name))
    
  }
  
  if (any(.d$n_distinct != 1L)) {
    
    .d <- .d %>%
      mutate(matched_name = remove_authority(matched_name)) %>%
      mutate(n_words = str_count(matched_name, "\\S+")) %>%
      filter(n_words == min(n_words)) %>%
      mutate(n_distinct = n_distinct(matched_name))
    
  }    
  
  if (any(.d$n_distinct != 1L)) {

    .d <- .d %>%
      pull(matched_name) %>%
      unique()
    return(tibble(
      chosen_name = str_c(.d, collapse = ", "),
      n_names = length(.d)
    ))
    
  }
  
  return(tibble(
    chosen_name = first(.d$matched_name),
    n_names = 1L
  ))
  
}

remove_authority <- function(.x) {
  
  .x %>%
    str_split(" ") %>%
    map_chr(~ {
      tibble(
        s = .x,
        capitals = str_detect(s, "[A-Z]")
      ) %>%
        mutate(
          row = row_number(),
          keep = row == 1 | !capitals,
          cutoff = ifelse(any(!keep), min(which(!keep)), Inf)
        ) %>%
        filter(row < cutoff) %>%
        pull(s) %>%
        str_c(collapse = " ")
    })
  
}

# Functions for combining trait and phylogenetic data ----
split_tip <- function(phy, tip.label, n) {
  
  tip.label <- str_replace_all(tip.label, " ", "_")
  
  node <- which(phy$tip.label == tip.label)
  
  for (i in 1:n) {
    phy <- bind.tip(phy, str_c(tip.label, LETTERS[i], sep = "_"), where = node)
    node <- which(phy$tip.label == tip.label)
  }
  
  phy <- drop.tip(phy, tip.label)
  
  phy
  
}

# Find polytomies ----
find_polytomies <- function(phy, n_tip) {
  
  descendent_counts <- rle(sort(phy$edge[, 1]))
  polytomies <- tibble(
    node = descendent_counts$values, 
    n = descendent_counts$lengths
  ) %>%
    filter(n > n_tip) %>%
    arrange(desc(n)) %>%
    rowwise() %>%
    mutate(genus = phy %>%
             extract.clade(node) %>%
             use_series(tip.label) %>%
             str_extract("^[A-Z]{1}[a-z]+") %>%
             unique() %>%
             str_c(collapse = ", "))
  
  polytomies
  
}

get_clade = function(phy) {
  # Take a phylogeny and get unique genera names.
  # Meant for getting genus name from tree that has already been extracted
  # a monophyletic genus
  phy |>
    use_series(tip.label) |>
    str_extract("^[A-Z]{1}[a-z]+") |>
    unique()
}

get_clade_num = function(clade, dir) {
  # Alignments are in format alignments/<clade>_<id>
  # This figures out which alignment to use
  tibble(file = list.files(dir)) |>
    filter(str_detect(file, str_c("^", clade, "_[0-9]+"))) |>
    mutate(clade_num = str_extract(file, str_c("^", clade, "_[0-9]+"))) |>
    pull(clade_num) |>
    unique()
}

get_ncbi_name = function(aln) {
  
  # Alignments from PyPHLAWD are labeled with NCBI uid. 
  # This function uses rotl package to get binomial name from NCBI taxonomy
  
  # Odd elements are uid; even elements are alignments.
  # Extract just uid and remove '>'
  uid = aln[seq(1, length(aln), 2)] |>
    str_remove("^>")
  assert_character(uid, min.chars = 1L, any.missing = FALSE, null.ok = FALSE)
  
  message(bold(yellow(glue("  Getting {n} taxonomic names from NCBI...", 
                           n = length(uid)))))
  if (length(uid) > 1000) {
    
    n = length(uid)
    inc = 1000
    i = seq(1, n, by = inc)
    j = pmin(seq(inc, n + inc - 1, by = inc), n)
    
    nms = map2(i, j, ~ {id2name(uid[.x:.y], db = "ncbi")}) |>
      flatten()
    
  } else {
    nms = id2name(uid, db = "ncbi")
  }
  
  nms
  
}

rename_alignment = function(aln_prefix, aln_names, overwrite) {
  
  # Alignments from PyPHLAWD are labeled with NCBI uid. 
  # This function uses output from get_ncbi_name() to rename taxa in alignment
  # and write a renamed alignment
  
  aln = read_lines(glue("{aln_prefix}_outaln", aln_prefix = aln_prefix))
  nms_file = glue("{aln_prefix}_id2name.rds", aln_prefix = aln_prefix)

  if (!file.exists(nms_file) | overwrite) {
  
    nms = get_ncbi_name(aln)
    write_rds(nms, nms_file)
    message(bold(green(glue("  ...done. Lookup table saved in '{nms_file}'", 
                            nms_file = nms_file))))
  } else {
    message(bold(green(glue("  Using saved NCBI names from '{nms_file}'", 
                            nms_file = nms_file))))
    nms = read_rds(nms_file)
  }

  aln[seq(1, length(aln), 2)] = imap_chr(unname(nms), ~{
    name = str_replace_all(.x$name, " ", "_")
    str_replace(aln[[.y * 2 - 1]], .x$id, name)
  })

  write_lines(aln, aln_names)
  message(bold(green(glue("  Named alignment written to '{aln_names}'", 
                          aln_names = aln_names))))

}

rename_duplicates = function(aln, aln_names) {
  
  # Rename duplicates in named alignment
  if (any(duplicated(names(aln)))) {
    
    aln_table = tibble(
      phy_name = names(aln),
      n_missing = map_int(aln, ~ {length(which(.x == "-"))})
    ) |>
      mutate(i = row_number()) |>
      group_by(phy_name) |>
      summarize(i = i[which.min(n_missing)], .groups = "drop")
    
    write.fasta(
      aln[aln_table$i], 
      names = names(aln)[aln_table$i], 
      str_c(aln_names, "1")
    )
    
  } else {
    write.fasta(aln, names = names(aln), str_c(aln_names, "1"))
  }
  
  message(bold(green(glue("  Modified alignment written under '{aln_names}'", 
                          aln_names = aln_names))))
  
} 

write_raxml_command = function(clade_num, raxml_dir, aln_prefix) {
  
  # Write commands for RAxML in shell script
  glue(
    "cd {raxml_dir}/{clade_num}
  raxmlHPC-PTHREADS -T 8 -M -m GTRGAMMA -p 12345 -q ../../{aln_prefix}_outpart -s ../../{aln_prefix}_outaln_names1 -n T22
    raxmlHPC -f I -m GTRGAMMA -p 12345 -t RAxML_bestTree.T22 -n T23
  cd ../..
\n", raxml_dir = raxml_dir, clade_num = clade_num, aln_prefix = aln_prefix
  ) |>
    cat(file = glue("{raxml_dir}/run-raxml.sh", raxml_dir = raxml_dir), 
        append = TRUE)
  
  message(bold(green(glue("  RAxML commands written in shell script, '{raxml_dir}/run-raxml.sh'",
                          raxml_dir = raxml_dir))))
  
}

write_subtree_info = function(phy, clade, aln, file) {
  
  taxa = phy$tip.label |>
    str_extract(str_c("^", clade, "_[a-z]+")) |>
    unique() |>
    intersect(names(aln))
  
  ret = list(taxa = taxa, subtree = phy)
  
  write_rds(ret, file)
  
  message(bold(green(glue("  Subtree and taxa written to '{file}'",
                          file = file))))
 
  return(ret) 
  
}

prepare_raxml <- function(node, phy, overwrite, clade = NULL) {
  
  alignment_dir = "alignments" # where PyPHLAWD alignments exist
  raxml_dir = "raxml" # where raxml output will go
  
  subtree = extract.clade(phy, node)
  
  # Manual override clade. Some clades have contain multiple genera.
  # To keep naming simple, provide manual override by supplying clade name.
  # Otherwise, extract clade name from subtree tip labels
  if (is.null(clade)) clade = get_clade(subtree)

  clade_num = get_clade_num(clade, alignment_dir)
  aln_prefix = glue("{alignment_dir}/{clade_num}", 
                    alignment_dir = alignment_dir, clade_num = clade_num)
    
  dir = glue("{raxml_dir}/{clade_num}", raxml_dir = raxml_dir, 
             clade_num = clade_num)
  if (!dir.exists(dir) | overwrite) {
    unlink(dir, recursive = TRUE)
    dir.create(dir)
    message(bold(green(glue("Directory {dir} created.", dir = dir))))
  }
  
  aln_names = glue("{aln_prefix}_outaln_names", aln_prefix = aln_prefix)
  
  if (!file.exists(aln_names) | overwrite) {
    rename_alignment(aln_prefix, aln_names, overwrite)
  }
  
  aln = read.fasta(aln_names)

  if (!file.exists(paste0(aln_names, "1")) | overwrite) {
    rename_duplicates(aln, aln_names)
  } else {
    message(bold(green(glue("  Modified alignment read from '{aln_names}'", 
                            aln_names = str_c(aln_names, "1")))))
  }
  aln = read.fasta(str_c(aln_names, "1"))
  
  write_raxml_command(clade_num, raxml_dir, aln_prefix)  
  
  # Save information needed to replace subtree with raxml tree
  subtree_info = glue("{raxml_dir}/{clade_num}.rds", raxml_dir = raxml_dir,
       clade_num = clade_num)
  if (!file.exists(subtree_info) | overwrite) {
    ret = write_subtree_info(subtree, clade, aln, subtree_info)
  } else {
    ret = read_rds(subtree_info)
    message(bold(green(glue("  Subtree and taxa read from '{file}'", 
                            file = subtree_info))))
  }
  
  ret
  
}


replace_polytomy <- function(phy, node, clade, seed) {
  
  set.seed(seed) # ensure reproducibility
  message(bold(green(clade)))
  
  assert_integerish(node, lower = 0, upper = Nnode(phy, FALSE), len = 1L)
  
  clade_num = tibble(dir = dir("raxml", pattern = "^[A-Z]{1}[a-z]+_[0-9]+$")) |>
    filter(str_detect(dir, str_c("^", clade, "_[0-9]+"))) |>
    pull(dir) |>
    unique()
  
  # Information on old subtree
  clade_info = read_rds(glue("raxml/{clade_num}.rds", clade_num = clade_num))
  old_subtree = clade_info$subtree
  tips_to_remove = clade_info$subtree$tip.label
  taxa = clade_info$taxa
  node_age = max(nodeHeights(phy)) - nodeheight(phy, node)
  
  if (length(taxa) <= 1L) {
    
    # If there less than 2 taxa in new subtree that overlap with old subtree,
    # then just choose two taxa randomly for phylo contrast
    taxa = sample(tips_to_remove, 2, replace = FALSE)
    phy = drop.tip(phy, tips_to_remove[tips_to_remove != taxa])
    return(phy)
    
  } else {
    
    # Create new subtree to replace old
    new_subtree = read.tree(glue("raxml/{clade_num}/RAxML_rootedTree.T23",
                                  clade_num = clade_num)) |>
      chronos(quiet = TRUE) |>
      compute.brlen()

    new_subtree1 = keep.tip(new_subtree, clade_info$taxa)
    
    diff_age = max(nodeHeights(new_subtree)) - max(nodeHeights(new_subtree1))
    new_subtree1$root.edge = if (is.null(new_subtree1$root.edge)) { 
      diff_age } else {
        new_subtree1$root.edge + diff_age
      }
    
    # Determine if it's better to used resolved parts of tree rather than grafting
    # RAxML subtree
    n_contrast_new = bind_rows(
      select(extract_contrasts(new_subtree1), -tree), 
      extract_sisters(new_subtree1)
    ) |>
      n_distinct()

    n_contrast_old = if (Nnode(old_subtree) > 1) {
      old_subtree1 = pull_resolved(old_subtree) |>
        compute.brlen()
      bind_rows(
        select(extract_contrasts(old_subtree1), -tree), 
        extract_sisters(old_subtree1)
      ) |>
        n_distinct()      
    } else {0}

    new_subtree = if (n_contrast_new > n_contrast_old) {
      new_subtree1
    } else {
      phy = drop.tip(
        phy, old_subtree$tip.label[!(old_subtree$tip.label %in% 
                                       old_subtree1$tip.label)]
      )
      return(phy)
    }
    
    # rescale branch lengths
    new_subtree$edge.length = new_subtree$edge.length * node_age
    new_subtree$root.edge = new_subtree$root.edge * node_age
    
    # Graft placeholder tip
    placeholder_phy = read.tree(
      text = glue("(placeholder_1:{node_age},placeholder_2:{node_age});", 
                  node_age = node_age)
    )
    
    new_phy = phy |>
      bind.tree(placeholder_phy, where = node) |>
      drop.tip(tips_to_remove)
    
    new_node = getMRCA(new_phy, c("placeholder_1", "placeholder_2"))
    
    new_phy = new_phy |>
      bind.tree(new_subtree, where = new_node) |>
      drop.tip(c("placeholder_1", "placeholder_2"))
    
    # Perform some checks to compare phy and new_phy
    ntip1 = length(old_subtree$tip.label) - length(clade_info$taxa)
    ntip2 = length(phy$tip.label) - length(new_phy$tip.label)
    assert_true(ntip1 == ntip2)
    
    nnode1 = Nnode(old_subtree) - Nnode(new_subtree)
    nnode2 = Nnode(phy) - Nnode(new_phy)
    assert_true(nnode1 == nnode2)
    
    age1 = max(nodeHeights(phy))
    age2 = max(nodeHeights(new_phy))
    assert_true(abs(age1 - age2) < 1e-6)
    
    phy = new_phy
    
  }
  
  phy
  
}

# Functions for PICs ----
my_pic <- function(x, phy) {
  ret <- pic(x, phy, scaled = TRUE, var.contrasts = TRUE, rescaled.tree = FALSE)
  # ret <- pic.ortho(x, phy, var.contrasts = TRUE)
  ret <- ret[, "contrasts"]
  names(ret) <- NULL
  ret
}

# Functions for figures ----
# adapted from https://github.com/tidyverse/ggplot2/blob/master/R/stat-ellipse.R
my_confidence_ellipse <- function(center, shape, nu, level, segments) {
  
  dfn <- 2
  chol_decomp <- chol(shape)
  radius <- sqrt(dfn * stats::qf(level, dfn, nu))
  angles <- (0:segments) * 2 * pi/segments
  unit.circle <- cbind(cos(angles), sin(angles))
  ellipse <- t(center + radius * t(unit.circle %*% chol_decomp))
  ellipse
  
}

# Function for scientific notation ----
scientize <- function(x, threshold = -1L, digits = 2L) {
  
  purrr::map_chr(x, function(.x, threshold, digits) {
    
    oom <- log10(abs(.x))
    if (oom < threshold) {
      x1 <- .x |>
        magrittr::multiply_by(10 ^ -floor(oom)) |>
        round(digits)
      
      x2 <- sprintf(glue::glue("%.{digits}f"), x1) |>
        stringr::str_c(" \\times 10^{", floor(oom), "}")
      
      return(x2)
    } else {
      x1 <- .x |>
        signif(digits + 1L) |>
        as.character()
      return(x1)
    }
    
  }, threshold = threshold, digits = digits)
  
}

# For clade with polytomy, pull resolved portion
pull_resolved = function(phy) {
  n_node = Nnode(phy)
  n_tip = Ntip(phy)
  assert_true(n_tip != (n_node + 1)) # make sure there are polytomies
  assert_true(n_node > 1) # make sure there are resolved portions
  phy_edge = as_tibble(phy$edge, .name_repair = "unique")
  degree = tabulate(phy_edge$`...1`)
  target = which(degree > 2)
  tip2drop = phy_edge |>
    filter(`...1` == target) |>
    filter(degree[`...2`] == 0) |>
    pull(`...2`)
  
  drop.tip(phy, tip2drop, root.edge = 1)
}

# Extracts all 4-tip subtrees 
extract_quartets = function(phy) {
  
  phy1 = write.tree(phy) 
  
  tip_string = "[A-Za-z_.]+:[0-9.]+"
  ab_string = glue("\\({tip},{tip}\\):[0-9.]+", tip = tip_string)

  # "(((a,b),c),d)"
  # "((c,(a,b)),d)"
  # "(d,((a,b),c))"
  # "(d,(c,(a,b)))"
  # "((a,b),(c,d))"
  
  quartet1 = glue("\\(\\({ab},{tip}\\):[0-9.]+,{tip}\\):[0-9.]+", 
                  tip = tip_string, ab = ab_string)
  quartet2 = glue("\\(\\({tip},{ab}\\):[0-9.]+,{tip}\\):[0-9.]+", 
                  tip = tip_string, ab = ab_string)
  quartet3 = glue("\\({tip},\\({ab},{tip}\\):[0-9.]+\\):[0-9.]+", 
                  tip = tip_string, ab = ab_string)
  quartet4 = glue("\\({tip},\\({tip},{ab}\\):[0-9.]+\\):[0-9.]+",
                  tip = tip_string, ab = ab_string)
  quartet5 = glue("\\({ab},{ab}\\)", ab = ab_string)
  
  quartet_trees = c(
    str_extract_all(phy1, quartet1),
    str_extract_all(phy1, quartet2),
    str_extract_all(phy1, quartet3),
    str_extract_all(phy1, quartet4),
    str_extract_all(phy1, quartet5)
  ) |>
    flatten()
  
  if (length(quartet_trees) == 0) {
    return(NULL)
  }
  
  quartet_trees |>
    str_c(";") %>%
    read.tree(text = .)
  
}

# Function to extract contrasts from a quartet
extract_contrasts_quartet = function(quartet) {
  
  assert_class(quartet, "phylo")
  assert_false(inherits(quartet, "multiPhylo"))
  assert_true(Ntip(quartet) == 4L)
  
  sisters = extract_sisters(quartet) |>
    select(sp1, sp2)
  
  if (nrow(sisters) == 2L) {
    # If all sister pairs
    return(sisters)
  } else {
    
    sister_pair = sisters |>
      as.matrix() |>
      t() |>
      as.vector()
    other_pair = drop.tip(quartet, sister_pair)
    others = tibble(
      sp1 = other_pair$tip.label[1],
      sp2 = other_pair$tip.label[2]
    )
    return(bind_rows(sisters, others))
  }
  stop("something went wrong in 'extract_contrasts_quartet'")
  
}

# Function to extract contrasts from single phylogeny
extract_contrasts_inner = function(single_phy) {
  
  assert_class(single_phy, "phylo")
  assert_false(inherits(single_phy, "multiPhylo"))
  
  phy1 = extract_quartets(single_phy)
  
  if (is.null(phy1)) {
    return(tibble(tree_node = numeric(0), sp1 = character(0), sp2 = character(0)))
  }
  
  if (inherits(phy1, "multiPhylo")) {
    ret = map_dfr(phy1, extract_contrasts_quartet)
  } else {
    ret = extract_contrasts_quartet(phy1) 
  }
  
  ret |>
    rowwise() |>
    mutate(tree_node = getMRCA(single_phy, c(sp1, sp2))) |>
    ungroup()
  
}

# Function to extract contrasts from multiple phylogenies
extract_contrasts = function(phy) {
  
  if (inherits(phy, "multiPhylo")) {
    
    ret = imap_dfr(phy, ~{
      
      assert_class(.x, "phylo")
      assert_false(inherits(.x, "multiPhylo"))
      
      extract_contrasts_inner(.x) |>
        mutate(tree = .y)
      
    }) 

  } else {
    ret = extract_contrasts_inner(phy) |>
      mutate(tree = 1)
  }
  
  ret
  
}
