test_that("greedy_parsimony() groups only proteins with identical peptide sets", {
  # two unrelated tie groups (A, B) and (C, D) plus a single protein E, all with one peptide
  pepprot <- prepare_matrix(
    data.frame(proteinID = c("A", "B", "C", "D", "E"), pep = c("p1", "p1", "p2", "p2", "p3")),
    peptide_id = "pep"
  )

  res <- greedy_parsimony(pepprot)

  expect_equal(unlist(res)[c("p1", "p2", "p3")], c(p1 = "A;B", p2 = "C;D", p3 = "E"))
})

test_that("greedy_parsimony() resolves a tie between overlapping proteins one at a time", {
  # A {p1, p2} and B {p2, p3} tie; A wins p1 and p2, B keeps p3
  pepprot <- prepare_matrix(
    data.frame(proteinID = c("A", "A", "B", "B"), pep = c("p1", "p2", "p2", "p3")),
    peptide_id = "pep"
  )

  expect_no_message(res <- greedy_parsimony(pepprot))

  expect_equal(unlist(res)[c("p1", "p2", "p3")], c(p1 = "A", p2 = "A", p3 = "B"))
})

test_that("greedy_parsimony() prefers the protein explaining most peptides", {
  pepprot <- prepare_matrix(
    data.frame(proteinID = c("A", "A", "A", "B", "C"), pep = c("p1", "p2", "p3", "p3", "p4")),
    peptide_id = "pep"
  )

  res <- greedy_parsimony(pepprot)

  expect_equal(unlist(res)[c("p1", "p2", "p3", "p4")], c(p1 = "A", p2 = "A", p3 = "A", p4 = "C"))
})

test_that("greedy_parsimony() assigns every precursor of the example data", {
  data(protpepmetashort, package = "prozor", envir = environment())
  pepprot <- prepare_matrix(protpepmetashort, peptide_id = "peptideModSeq")

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

  expect_equal(mat$Peptide, c("AC[+57.0]K", "PEPTIDEK"))
  expect_equal(mat$precursorCharge, c("2", "3"))
  expect_equal(mat$Protein, c("P1", "P2;P3"))
  expect_equal(rownames(mat), names(res))
})

test_that("greedy_to_table() handles modified precursors of the example data", {
  data(protpepmetashort, package = "prozor", envir = environment())
  protpepmetashort$precursor <- paste(protpepmetashort$peptideModSeq, protpepmetashort$precursorCharge, sep = ".")
  pepprot <- prepare_matrix(protpepmetashort, peptide_id = "precursor")
  res <- greedy_parsimony(pepprot)

  mat <- greedy_to_table(res)

  expect_equal(nrow(mat), length(res))
  expect_equal(paste(mat$Peptide, mat$precursorCharge, sep = "."), names(res))
  expect_equal(mat$Protein, unlist(res, use.names = FALSE))
})
