test_that("greedy_parsimony() groups only proteins with identical peptide sets", {
  # two unrelated tie groups (A, B) and (C, D) plus a single protein E, all with one peptide
  data <- data.frame(protein_id = c("A", "B", "C", "D", "E"), pep = c("p1", "p1", "p2", "p2", "p3"))

  res <- greedy_parsimony(data, peptide = "pep")

  expect_type(res, "character")
  expect_equal(res[c("p1", "p2", "p3")], c(p1 = "A;B", p2 = "C;D", p3 = "E"))
})

test_that("greedy_parsimony() resolves a tie between overlapping proteins one at a time", {
  # A {p1, p2} and B {p2, p3} tie; A wins p1 and p2, B keeps p3
  data <- data.frame(protein_id = c("A", "A", "B", "B"), pep = c("p1", "p2", "p2", "p3"))

  expect_no_message(res <- greedy_parsimony(data, peptide = "pep"))

  expect_equal(res[c("p1", "p2", "p3")], c(p1 = "A", p2 = "A", p3 = "B"))
})

test_that("greedy_parsimony() follows the level order of factor columns", {
  # A {p1, p2} and B {p2, p3} tie; the protein first in factor(protein_id) wins p2
  b_first <- data.frame(
    protein_id = factor(c("A", "A", "B", "B"), levels = c("B", "A")),
    peptide_seq = c("p1", "p2", "p2", "p3")
  )
  peptides_reversed <- data.frame(
    protein_id = c("A", "A", "B", "B"),
    peptide_seq = factor(c("p1", "p2", "p2", "p3"), levels = c("p3", "p2", "p1"))
  )

  expect_identical(greedy_parsimony(b_first), c(p2 = "B", p3 = "B", p1 = "A"))
  expect_identical(greedy_parsimony(peptides_reversed), c(p2 = "A", p1 = "A", p3 = "B"))
})

test_that("greedy_parsimony() groups proteins on the peptides that remain", {
  data <- data.frame(
    protein_id = c("A", "A", "A", "A", "B", "B", "B", "C", "C", "C"),
    peptide_seq = c("p1", "p2", "p3", "p4", "p1", "p5", "p6", "p2", "p5", "p6")
  )

  res <- greedy_parsimony(data)

  expect_equal(res[c("p1", "p2", "p5", "p6")], c(p1 = "A", p2 = "A", p5 = "B;C", p6 = "B;C"))
})

test_that("greedy_parsimony() prefers the protein explaining most peptides", {
  data <- data.frame(protein_id = c("A", "A", "A", "B", "C"), pep = c("p1", "p2", "p3", "p3", "p4"))

  res <- greedy_parsimony(data, peptide = "pep")

  expect_equal(res[c("p1", "p2", "p3", "p4")], c(p1 = "A", p2 = "A", p3 = "A", p4 = "C"))
})

test_that("greedy_parsimony() ranks proteins by peptide count, however widely the peptides are shared", {
  # A explains 4 peptides, 3 of them shared with X, Y, Z; B explains 3, 1 shared with A
  data <- data.frame(
    protein_id = c("A", "A", "A", "A", "B", "B", "B", "X", "Y", "Z"),
    pep = c("p1", "p2", "p3", "p6", "p3", "p4", "p5", "p1", "p2", "p6")
  )

  res <- greedy_parsimony(data, peptide = "pep")

  expect_equal(
    res[c("p1", "p2", "p3", "p4", "p5", "p6")],
    c(p1 = "A", p2 = "A", p3 = "A", p4 = "B", p5 = "B", p6 = "A")
  )
})

test_that("greedy_parsimony() counts a repeated peptide-protein pair once", {
  # A {p1, p2} and B {p2, p3} tie; counting the three rows of B and p3 (e.g. three offsets) would let B win p2
  data <- data.frame(
    protein_id = c("A", "A", "B", "B", "B", "B"),
    peptide_seq = c("p1", "p2", "p2", "p3", "p3", "p3")
  )

  res <- greedy_parsimony(data)

  expect_identical(res, c(p1 = "A", p2 = "A", p3 = "B"))
  expect_identical(res, greedy_parsimony(unique(data)))
})

test_that("greedy_parsimony() breaks ties between character IDs in C-locale order in any locale", {
  # a {p1, p2} and B {p2, p3} tie; "B" comes first in C-locale order, "a" in en_US collation
  suppressWarnings(withr::local_collate("en_US.UTF-8"))
  skip_if_not(Sys.getlocale("LC_COLLATE") == "en_US.UTF-8", "no en_US.UTF-8 locale")
  data <- data.frame(protein_id = c("a", "a", "B", "B"), peptide_seq = c("p1", "p2", "p2", "p3"))

  expect_identical(greedy_parsimony(data), c(p2 = "B", p3 = "B", p1 = "a"))
})

test_that("greedy_parsimony() gives the same result for a tibble and a data frame", {
  skip_if_not_installed("dplyr")
  data <- data.frame(protein_id = c("A", "B", "B"), peptide_seq = c("p1", "p1", "p2"))

  res <- greedy_parsimony(dplyr::as_tibble(data))

  expect_identical(res, c(p1 = "B", p2 = "B"))
  expect_identical(res, greedy_parsimony(data))
})

test_that("greedy_parsimony() ignores unused factor levels", {
  # most protein and peptide levels are unused in the first rows; the first peptide level is unused without it
  first_rows <- protpepmetashort[1:50, ]
  without_first <- protpepmetashort[protpepmetashort$peptide_seq != levels(protpepmetashort$peptide_seq)[1], ]

  res_first_rows <- greedy_parsimony(first_rows)
  res_without_first <- greedy_parsimony(without_first)

  expect_setequal(names(res_first_rows), as.character(first_rows$peptide_seq))
  expect_setequal(names(res_without_first), as.character(without_first$peptide_seq))
  expect_identical(res_first_rows, greedy_parsimony(droplevels(first_rows)))
  expect_identical(res_without_first, greedy_parsimony(droplevels(without_first)))
})

test_that("greedy_parsimony() stops for input other than a data frame", {
  data <- data.frame(protein_id = c("A", "B"), peptide_seq = c("p1", "p2"))

  expect_error(greedy_parsimony(as.matrix(data)), "`annotation` must be a data frame of .*, not matrix\\.$")
})

test_that("greedy_parsimony() stops, naming them, for missing columns", {
  data <- data.frame(protein = c("A", "B"), peptide_seq = c("p1", "p2"))

  expect_error(
    greedy_parsimony(data),
    "must have the columns `protein_id`, `peptide_seq`; missing: `protein_id`\\.$"
  )
  expect_error(greedy_parsimony(data, peptide = "precursor"), "missing: `protein_id`, `precursor`\\.$")
})

test_that("greedy_parsimony() stops, counting them, for NA IDs", {
  data <- data.frame(protein_id = c("A", NA, NA), peptide_seq = c("p1", "p2", NA))

  expect_error(greedy_parsimony(data), "must not have NA IDs; found 2 in `protein_id` and 1 in `peptide_seq`\\.$")
  expect_error(greedy_parsimony(data[1:2, ]), "found 1 in `protein_id`\\.$")
})

test_that("greedy_parsimony() handles 20000 peptides in under 5 s", {
  skip_on_cran()
  # Genes with 1 + Poisson(1) isoforms that share most of their peptides, a tenth of the peptides also match a second
  # gene: thousands of proteins tie in the late rounds. Comparing tied proteins one at a time took 18 s for 5000
  # peptides, recounting the peptides of all proteins in every round 2.6 s for these 20000.
  withr::local_seed(1)
  n_genes <- 6800
  gene <- sample(n_genes, 20000, replace = TRUE)
  n_isoforms <- 1 + stats::rpois(n_genes, 1)
  isoform_gene <- rep(seq_len(n_genes), n_isoforms)
  isoform <- ifelse(sequence(n_isoforms) == 1, "", paste0("-", sequence(n_isoforms)))
  of_peptide <- split(seq_along(isoform_gene), isoform_gene)[gene]
  data <- data.frame(
    protein_id = paste0("G", isoform_gene, isoform)[unlist(of_peptide)],
    peptide_seq = rep(sprintf("PEP%05d", seq_along(gene)), lengths(of_peptide))
  )
  data <- data[!(grepl("-", data$protein_id) & stats::runif(nrow(data)) < 0.2), ]
  shared <- sample(length(gene), 2000)
  data <- rbind(
    data,
    data.frame(protein_id = paste0("G", sample(n_genes, 2000)), peptide_seq = sprintf("PEP%05d", shared))
  )

  time <- system.time(res <- greedy_parsimony(data))[["elapsed"]]

  expect_lt(time, 5)
  expect_setequal(names(res), data$peptide_seq)
  # every protein of a group matches the peptides assigned to it
  proteins <- strsplit(res, ";", fixed = TRUE)
  pairs <- paste(rep(names(res), lengths(proteins)), unlist(proteins))
  expect_true(all(pairs %in% paste(data$peptide_seq, data$protein_id)))
})

test_that("greedy_parsimony() assigns every precursor of the example data", {
  res <- greedy_parsimony(protpepmetashort, peptide = "peptide_mod_seq")

  expect_length(res, length(unique(protpepmetashort$peptide_mod_seq)))
  expect_setequal(names(res), as.character(protpepmetashort$peptide_mod_seq))
  # 101 protein groups: recorded from the fixed implementation (c416c0c), not computed independently
  expect_length(unique(res), 101)
  # every protein of a group matches the precursors assigned to it
  proteins <- strsplit(res, ";", fixed = TRUE)
  pairs <- paste(rep(names(res), lengths(proteins)), unlist(proteins))
  expect_true(all(pairs %in% paste(protpepmetashort$peptide_mod_seq, protpepmetashort$protein_id)))

  # recorded output, guards against unintended changes; skipped on CRAN
  assignment <- res[sort(names(res))]
  expect_snapshot(cat(paste(names(assignment), assignment, sep = " -> "), sep = "\n"))
})

test_that("greedy_to_table() splits the charge at the last dot", {
  groups <- c("AC[+57.0]K.2" = "P1", "PEPTIDEK.3" = "P2;P3")

  tab <- greedy_to_table(groups)

  expect_named(tab, c("peptide", "precursor_charge", "protein_group"))
  expect_equal(tab$peptide, c("AC[+57.0]K", "PEPTIDEK"))
  expect_identical(tab$precursor_charge, c(2L, 3L))
  expect_identical(tab$protein_group, c("P1", "P2;P3"))
  expect_equal(rownames(tab), names(groups))
})

test_that("greedy_to_table() stops for names without a charge", {
  peptide_level <- greedy_parsimony(protpepmetashort)

  expect_error(
    greedy_to_table(c(PEPTIDEK = "P1", "AC[+57.0]K" = "P2", "PEPTIDER.2" = "P3")),
    "names of `groups` must have the form `peptide.charge`.*found: \"PEPTIDEK\", \"AC\\[\\+57.0\\]K\"\\.$"
  )
  expect_error(greedy_to_table(peptide_level), "names of `groups` must have the form `peptide.charge`")
})

test_that("greedy_to_table() handles modified precursors of the example data", {
  protpepmetashort$precursor <- paste(protpepmetashort$peptide_mod_seq, protpepmetashort$precursor_charge, sep = ".")
  groups <- greedy_parsimony(protpepmetashort, peptide = "precursor")

  tab <- greedy_to_table(groups)

  expect_equal(nrow(tab), length(groups))
  expect_equal(paste(tab$peptide, tab$precursor_charge, sep = "."), names(groups))
  expect_equal(tab$protein_group, unname(groups))
})
