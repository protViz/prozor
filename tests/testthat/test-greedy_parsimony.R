test_that("greedy_parsimony() groups only proteins with identical peptide sets", {
  # two unrelated tie groups (A, B) and (C, D) plus a single protein E, all with one peptide
  pepprot <- prepare_matrix(
    data.frame(protein_id = c("A", "B", "C", "D", "E"), pep = c("p1", "p1", "p2", "p2", "p3")),
    peptide_id = "pep"
  )

  res <- greedy_parsimony(pepprot)

  expect_equal(unlist(res)[c("p1", "p2", "p3")], c(p1 = "A;B", p2 = "C;D", p3 = "E"))
})

test_that("greedy_parsimony() resolves a tie between overlapping proteins one at a time", {
  # A {p1, p2} and B {p2, p3} tie; A wins p1 and p2, B keeps p3
  pepprot <- prepare_matrix(
    data.frame(protein_id = c("A", "A", "B", "B"), pep = c("p1", "p2", "p2", "p3")),
    peptide_id = "pep"
  )

  expect_no_message(res <- greedy_parsimony(pepprot))

  expect_equal(unlist(res)[c("p1", "p2", "p3")], c(p1 = "A", p2 = "A", p3 = "B"))
})

test_that("greedy_parsimony() groups proteins on the peptides that remain", {
  data <- data.frame(
    protein_id = c("A", "A", "A", "A", "B", "B", "B", "C", "C", "C"),
    peptide_seq = c("p1", "p2", "p3", "p4", "p1", "p5", "p6", "p2", "p5", "p6")
  )

  res <- greedy_parsimony(prepare_matrix(data))

  expect_equal(unlist(res)[c("p1", "p2", "p5", "p6")], c(p1 = "A", p2 = "A", p5 = "B;C", p6 = "B;C"))
})

test_that("greedy_parsimony() prefers the protein explaining most peptides", {
  pepprot <- prepare_matrix(
    data.frame(protein_id = c("A", "A", "A", "B", "C"), pep = c("p1", "p2", "p3", "p3", "p4")),
    peptide_id = "pep"
  )

  res <- greedy_parsimony(pepprot)

  expect_equal(unlist(res)[c("p1", "p2", "p3", "p4")], c(p1 = "A", p2 = "A", p3 = "A", p4 = "C"))
})

test_that("greedy_parsimony() ranks proteins by peptide count, not by weight sum", {
  # A explains 4 peptides, B 3; with inverse weights B has the larger weight sum (2.5 vs 2)
  data <- data.frame(
    protein_id = c("A", "A", "A", "A", "B", "B", "B", "X", "Y", "Z"),
    pep = c("p1", "p2", "p3", "p6", "p3", "p4", "p5", "p1", "p2", "p6")
  )
  count <- greedy_parsimony(prepare_matrix(data, peptide_id = "pep"))

  inverse <- greedy_parsimony(prepare_matrix(data, peptide_id = "pep", weighting = "inverse"))

  expect_equal(
    unlist(count)[c("p1", "p2", "p3", "p4", "p5", "p6")],
    c(p1 = "A", p2 = "A", p3 = "A", p4 = "B", p5 = "B", p6 = "A")
  )
  expect_identical(inverse, count)
})

test_that("greedy_parsimony() leaves out peptides without a protein", {
  # A and B tie on p1 and are removed together; p2 matches no protein
  pepprot <- Matrix::sparseMatrix(
    i = c(1, 1),
    j = c(1, 2),
    x = 1,
    dims = c(2, 2),
    dimnames = list(c("p1", "p2"), c("A", "B"))
  )

  expect_no_warning(res <- greedy_parsimony(pepprot))

  expect_equal(res, list(p1 = "A;B"))
})

test_that("greedy_parsimony() stops for a matrix without peptide or protein names", {
  pepprot <- Matrix::sparseMatrix(i = c(1, 2), j = c(1, 2), x = 1)

  expect_error(
    greedy_parsimony(pepprot),
    "`pepprot` must have peptide row names and protein column names; it has no row and column names"
  )
  rownames(pepprot) <- c("p1", "p2")
  expect_error(greedy_parsimony(pepprot), "it has no column names")
})

test_that("greedy_parsimony() handles 5000 peptides in a few seconds", {
  skip_on_cran()
  # Genes with 1 + Poisson(1) isoforms that share most of their peptides, a tenth of the peptides also match a second
  # gene: thousands of proteins tie in the late rounds. Comparing tied proteins one at a time took 18 s here.
  withr::local_seed(1)
  n_genes <- 1700
  gene <- sample(n_genes, 5000, replace = TRUE)
  n_isoforms <- 1 + stats::rpois(n_genes, 1)
  isoform_gene <- rep(seq_len(n_genes), n_isoforms)
  isoform <- ifelse(sequence(n_isoforms) == 1, "", paste0("-", sequence(n_isoforms)))
  of_peptide <- lapply(gene, function(g) which(isoform_gene == g))
  data <- data.frame(
    protein_id = paste0("G", isoform_gene, isoform)[unlist(of_peptide)],
    peptide_seq = rep(sprintf("PEP%04d", seq_along(gene)), lengths(of_peptide))
  )
  data <- data[!(grepl("-", data$protein_id) & stats::runif(nrow(data)) < 0.2), ]
  shared <- sample(length(gene), 500)
  data <- rbind(
    data,
    data.frame(protein_id = paste0("G", sample(n_genes, 500)), peptide_seq = sprintf("PEP%04d", shared))
  )
  pepprot <- prepare_matrix(data)

  time <- system.time(res <- greedy_parsimony(pepprot))[["elapsed"]]

  expect_lt(time, 5)
  expect_setequal(names(res), rownames(pepprot))
  matches <- vapply(names(res), function(p) all(pepprot[p, strsplit(res[[p]], ";")[[1]]] > 0), logical(1))
  expect_true(all(matches))
})

test_that("greedy_parsimony() assigns every precursor of the example data", {
  data(protpepmetashort, package = "prozor", envir = environment())
  pepprot <- prepare_matrix(protpepmetashort, peptide_id = "peptide_mod_seq")

  res <- greedy_parsimony(pepprot)

  expect_length(res, nrow(pepprot))
  expect_setequal(names(res), rownames(pepprot))
  # 101 protein groups: recorded from the fixed implementation (c416c0c), not computed independently
  expect_length(unique(unlist(res)), 101)
  # every protein of a group matches the precursors assigned to it
  matches <- vapply(names(res), function(p) all(pepprot[p, strsplit(res[[p]], ";")[[1]]] > 0), logical(1))
  expect_true(all(matches))

  # recorded output, guards against unintended changes; skipped on CRAN
  assignment <- unlist(res)[sort(names(res))]
  expect_snapshot(cat(paste(names(assignment), assignment, sep = " -> "), sep = "\n"))
})

test_that("greedy_to_table() splits the charge at the last dot", {
  res <- list("AC[+57.0]K.2" = "P1", "PEPTIDEK.3" = "P2;P3")

  mat <- greedy_to_table(res)

  expect_equal(mat$peptide, c("AC[+57.0]K", "PEPTIDEK"))
  expect_identical(mat$precursor_charge, c(2L, 3L))
  expect_equal(mat$protein, c("P1", "P2;P3"))
  expect_equal(rownames(mat), names(res))
})

test_that("greedy_to_table() stops for names without a charge", {
  data(protpepmetashort, package = "prozor", envir = environment())
  peptide_level <- greedy_parsimony(prepare_matrix(protpepmetashort))

  expect_error(
    greedy_to_table(list(PEPTIDEK = "P1", "AC[+57.0]K" = "P2", "PEPTIDER.2" = "P3")),
    "names of `res` must have the form `peptide.charge`.*found: \"PEPTIDEK\", \"AC\\[\\+57.0\\]K\"\\.$"
  )
  expect_error(greedy_to_table(peptide_level), "names of `res` must have the form `peptide.charge`")
})

test_that("greedy_to_table() handles modified precursors of the example data", {
  data(protpepmetashort, package = "prozor", envir = environment())
  protpepmetashort$precursor <- paste(protpepmetashort$peptide_mod_seq, protpepmetashort$precursor_charge, sep = ".")
  pepprot <- prepare_matrix(protpepmetashort, peptide_id = "precursor")
  res <- greedy_parsimony(pepprot)

  mat <- greedy_to_table(res)

  expect_equal(nrow(mat), length(res))
  expect_equal(paste(mat$peptide, mat$precursor_charge, sep = "."), names(res))
  expect_equal(mat$protein, unlist(res, use.names = FALSE))
})
