test_that("greedy_parsimony() groups only proteins with identical peptide sets", {
  # two unrelated tie groups (A, B) and (C, D) plus a single protein E, all with one peptide
  pepprot <- prepareMatrix(
    data.frame(proteinID = c("A", "B", "C", "D", "E"), pep = c("p1", "p1", "p2", "p2", "p3")),
    peptideID = "pep"
  )

  # .greedy2() reports its run time as a message
  res <- suppressMessages(greedy_parsimony(pepprot))

  expect_equal(unlist(res)[c("p1", "p2", "p3")], c(p1 = "A;B", p2 = "C;D", p3 = "E"))
})

test_that("greedy_parsimony() prefers the protein explaining most peptides", {
  pepprot <- prepareMatrix(
    data.frame(proteinID = c("A", "A", "A", "B", "C"), pep = c("p1", "p2", "p3", "p3", "p4")),
    peptideID = "pep"
  )

  res <- suppressMessages(greedy_parsimony(pepprot))

  expect_equal(unlist(res)[c("p1", "p2", "p3", "p4")], c(p1 = "A", p2 = "A", p3 = "A", p4 = "C"))
})

test_that("greedy_parsimony() assigns every precursor of the example data", {
  data(protpepmetashort, package = "prozor", envir = environment())
  pepprot <- prepareMatrix(protpepmetashort, peptideID = "peptideModSeq")

  res <- suppressMessages(greedy_parsimony(pepprot))

  expect_length(res, nrow(pepprot))
  expect_setequal(names(res), rownames(pepprot))
  expect_length(unique(unlist(res)), 101)
  # every protein of a group matches the precursors assigned to it
  matches <- vapply(names(res), function(p) all(pepprot[p, strsplit(res[[p]], ";")[[1]]] > 0), logical(1))
  expect_true(all(matches))

  assignment <- unlist(res)[sort(names(res))]
  expect_snapshot(cat(paste(names(assignment), assignment, sep = " -> "), sep = "\n"))
})

test_that("greedyRes2Matrix() splits peptide and charge", {
  data(protpepmetashort, package = "prozor", envir = environment())
  protpepmetashort$precursor <- paste(protpepmetashort$peptideSeq,
                                      protpepmetashort$precursorCharge, sep = ".")
  pepprot <- prepareMatrix(protpepmetashort, peptideID = "precursor")
  res <- suppressMessages(greedy_parsimony(pepprot))

  mat <- greedyRes2Matrix(res)

  expect_s3_class(mat, "data.frame")
  expect_named(mat, c("Peptide", "precursorCharge", "Protein"))
  expect_equal(nrow(mat), length(res))
  expect_equal(paste(mat$Peptide, mat$precursorCharge, sep = "."), names(res))
})
