# Calculations of gsmax follow Sack and Buckley 2016 ----
biophysical_constant <- function(D_wv, v) D_wv / v

morphological_constant <- function(c, h, j) {
  (pi * c ^ 2) / (j ^ 0.5 * (4 * h * j + pi * c))
}

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

  # For some reason NCBI doesn't return all names (randomly). Go back and fill in
  # names until all there. This won't work if there >1000 missing, but that seems
  # unlikely
  missing = which(sapply(nms, nrow) == 0)
  while(length(missing) > 0) {
    nms[missing] = id2name(uid[missing], db = "ncbi")
    missing = which(sapply(nms, nrow) == 0)
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
    name = str_replace_all(.x$name, " ", "_") |>
      # remove illegal characters in FASTA that sometimes appear
      str_remove_all("\\[") |>
      str_remove_all("\\]")
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
  raxmlHPC-PTHREADS-SSE3 -T 10 -M -m GTRGAMMA -p 12345 -q ../../{aln_prefix}_outpart -s ../../{aln_prefix}_outaln_names1 -n T22
  raxmlHPC-PTHREADS-SSE3 -T 10 -f I -m GTRGAMMA -p 12345 -t RAxML_bestTree.T22 -n T23
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

# Functions for numerical integration to determine Cov(d,s) ----
# This models d (log stomatal density) and s (log stomatal size) when the
# developmental function is fixed (A and B are constants), but stomatal index 
# (I) and meristem size (M) can evolve. The model assumes that logit(I) and 
# log(M) are bivariate normally distributed with mean vector Mu and covariance
# matrix Sigma.

# Probability density of d
fd = function(
    d,
    A, B, Mu, Sigma,
    log = TRUE
) {
  # integrate over all possible logit_I for each d
  ret = map_dbl(d, function(.x, A, B, Mu, Sigma) {
    log(integrate(
      function(logit_I, D, A, B, Mu, Sigma) {
        # get m from logit_I and d
        I = plogis(logit_I)
        i = log(I)
        a = log(A)
        d = log(D)
        m = i - a - d - log(I * B + (1 - I))
        exp(dmvn(cbind(logit_I, m), Mu, Sigma, log = TRUE))
      },
      lower = -Inf, upper = Inf,
      D = exp(.x), A = A, B = B, Mu = Mu, Sigma = Sigma
    )$value)
  }, A = A, B = B, Mu = Mu, Sigma = Sigma)
  
  if (log) {
    return(ret) 
  } else {
    return(exp(ret))
  }
  
}

# Probability density of s
fs = function(
    s,
    A, B, Mu, Sigma,
    log = TRUE
) {
  # integrate over all possible logit_I for each s
  ret = map_dbl(s, function(.x, A, B, Mu, Sigma) {
    log(integrate(
      function(logit_I, S, A, B, Mu, Sigma) {
        # get m from logit_I and d
        I = plogis(logit_I)
        i = log(I)
        a = log(A)
        b = log(B)
        s = log(S)
        m = s - a - b
        exp(dmvn(cbind(logit_I, m), Mu, Sigma, log = TRUE))
      },
      lower = -Inf, upper = Inf,
      S = exp(.x), A = A, B = B, Mu = Mu, Sigma = Sigma
    )$value)
  }, A = A, B = B, Mu = Mu, Sigma = Sigma)
  
  if (log) {
    return(ret) 
  } else {
    return(exp(ret))
  }
  
}

# Expected value of d
Ed = function(A, B, Mu, Sigma) {
  integrate(
    function(d, A, B, Mu, Sigma) {
      d * exp(fd(d, A = A, B = B, Mu = Mu, Sigma = Sigma, log = TRUE))
    }, 
    lower = -Inf, upper = Inf,
    A = A, B = B, Mu = Mu, Sigma = Sigma
  )$value
}

# Expected value of s
Es = function(A, B, Mu, Sigma) {
  integrate(
    function(s, A, B, Mu, Sigma) {
      s * exp(fs(s, A = A, B = B, Mu = Mu, Sigma = Sigma, log = TRUE))
    }, 
    lower = -Inf, upper = Inf,
    A = A, B = B, Mu = Mu, Sigma = Sigma
  )$value
}

# Expected value of d^2
Ed2 = function(A, B, Mu, Sigma) {
  integrate(
    function(d, A, B, Mu, Sigma) {
      d ^ 2 * exp(fd(d, A = A, B = B, Mu = Mu, Sigma = Sigma, log = TRUE))
    }, 
    lower = -Inf, upper = Inf,
    A = A, B = B, Mu = Mu, Sigma = Sigma
  )$value
}

# Expected value of s^2
Es2 = function(A, B, Mu, Sigma) {
  integrate(
    function(s, A, B, Mu, Sigma) {
      s ^ 2 * exp(fs(s, A = A, B = B, Mu = Mu, Sigma = Sigma, log = TRUE))
    }, 
    lower = -Inf, upper = Inf,
    A = A, B = B, Mu = Mu, Sigma = Sigma
  )$value
}

# Variance of d
Var_d = function(A, B, Mu, Sigma) {
  # Var(X) = E[X^2] - E[X] ^ 2
  Ed2(A, B, Mu, Sigma) - Ed(A, B, Mu, Sigma) ^ 2
}

# Variance of s
Var_s = function(A, B, Mu, Sigma) {
  # Var(X) = E[X^2] - E[X] ^ 2
  Es2(A, B, Mu, Sigma) - Es(A, B, Mu, Sigma) ^ 2
}

# Expected value of d * s
Eds = function(A, B, Mu, Sigma) {
  integrate(
    function(m, A, B, Mu, Sigma) {
      map_dbl(m, function(.x, A, B, Mu, Sigma) {
        
        integrate(function(logit_I, m, A, B, Mu, Sigma) {
          i = plogis(logit_I, log.p = TRUE)
          s = m + log(A) + log(B)
          d = i - log(A) - m - log(exp(i + log(B)) + (1 - exp(i)))
          d * s * exp(dmvn(cbind(logit_I, m), Mu, Sigma, log = TRUE))
        }, 
        lower = -Inf, upper = Inf,
        m = .x, A = A, B = B, Mu = Mu, Sigma = Sigma
        )$value
      }, A = A, B = B, Mu = Mu, Sigma = Sigma)
    }, 
    lower = -Inf, upper = Inf,
    A = A, B = B, Mu = Mu, Sigma = Sigma
  )$value
}

# Covariance of d and s
Cov_ds = function(A, B, Mu, Sigma) {
  # Cov(X,Y) = E[XY] - E[X] E[Y]
  Eds(A, B, Mu, Sigma) - Ed(A, B, Mu, Sigma) * Es(A, B, Mu, Sigma)
}

# Covariance matrix
Sigma_ds = function(A, B, Mu, Sigma) {
  matrix(c(Var_d(A, B, Mu, Sigma), rep(Cov_ds(A, B, Mu, Sigma), 2), 
           Var_s(A, B, Mu, Sigma)), ncol = 2)
}

# I probably don't need these functions for anything, but keeping here in case
# Mu = c(mean of logit_I, mean of m)
# Sigma[1,1] = Var(logit_I)
# Sigma[1,2] = Sigma[2,1] = Cov(logit_I, m)
# Sigma[2,2] = Var(m)

# fi(x|m) := dist of logit_I at a given m
dlogit_I_m = function(logit_I, m, Mu, Sigma, log = TRUE) {
  rho = cov2cor(Sigma)[1, 2]
  dnorm(
    logit_I, 
    Mu[1] + rho * sqrt(Sigma[1, 1] / Sigma[2, 2]) * (m - Mu[2]), 
    (1 - rho ^ 2) * Sigma[1,1],
    log = log
  )
  
}

# fm(x|logit_I) := dist of m at a given logit_I
dm_logit_I = function(m, logit_I, Mu, Sigma, log = TRUE) {
  rho = cov2cor(Sigma)[1, 2]
  dnorm(
    m, 
    Mu[2] + rho * sqrt(Sigma[2, 2] / Sigma[1, 1]) * (logit_I - Mu[1]), 
    (1 - rho ^ 2) * Sigma[2, 2],
    log = log
  )
  
}

# Probability density function of logit_I with i ~ Normal(mu_i, sigma_i)
f_logitI = function(logit_I, mu_i, sigma_i) {
  
  # Normalization scalar
  ns = integrate(function(logit_I, mu_i, sigma_i) {
    i = plogis(logit_I, log.p = TRUE)
    exp(dnorm(i, mu_i, sigma_i, log = TRUE) -
          pnorm(0, mu_i, sigma_i, log = TRUE))
  }, lower = -Inf, upper = Inf, mu_i = mu_i, sigma_i = sigma_i,
  subdivisions = 1000L, abs.tol = 1e-16)$value
  
  i = plogis(logit_I, log.p = TRUE)
  exp(dnorm(i, mu_i, sigma_i, log = TRUE) -
        pnorm(0, mu_i, sigma_i, log = TRUE) - log(ns))
  
}

# Expected value of logit_I with i ~ Normal(mu_i, sigma_i)
# This is slightly off, but I can't figure out why
E_logitI = function(mu_i, sigma_i) {
  integrate(function(logit_I, mu_i, sigma_i) {
    logit_I * f_logitI(logit_I, mu_i, sigma_i)
  }, lower = -Inf, upper = Inf, mu_i = mu_i, sigma_i = sigma_i)$value
}

# Expected value of logit_I ^ 2 with i ~ Normal(mu_i, sigma_i)
# This is slightly off, but I can't figure out why
E_logitI2 = function(mu_i, sigma_i) {
  integrate(function(logit_I, mu_i, sigma_i) {
    logit_I ^ 2 * f_logitI(logit_I, mu_i, sigma_i)
  }, lower = -Inf, upper = Inf, mu_i = mu_i, sigma_i = sigma_i)$value
}

# Variance of logit_I  with i ~ Normal(mu_i, sigma_i)
# This is slightly off, but I can't figure out why
Var_logitI = function(mu_i, sigma_i) {
  E_logitI2(mu_i, sigma_i) - E_logitI(mu_i, sigma_i) ^ 2
}
