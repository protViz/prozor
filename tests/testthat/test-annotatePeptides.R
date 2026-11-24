test_that("annotateAHO() finds every occurrence with 1-based offsets", {
  expect_no_message(res <- annotateAHO(c("PEPTIDEK", "GGGR"), toy_fasta()))
  res <- res[order(res$peptideSeq, res$proteinID, res$Offset), ]

  expect_equal(res$peptideSeq, c("GGGR", rep("PEPTIDEK", 6)))
  expect_equal(res$proteinID, c("P1", "P1", "P2", "P3", "P4", "P5", "P5"))
  expect_equal(res$Offset, c(14, 6, 4, 1, 2, 1, 10))
  expect_equal(res$proteinSequence[res$proteinID == "P2"], "MSSPEPTIDEKLLL")
})

test_that("annotatePeptides() flags each occurrence by the residue before it", {
  expect_no_warning(expect_no_message(res <- annotatePeptides("PEPTIDEK", toy_fasta())))
  res <- res[order(res$proteinID, res$Offset), ]

  expect_equal(res$proteinID, c("P1", "P2", "P3", "P4", "P5", "P5"))
  expect_equal(res$Offset, c(6, 4, 1, 2, 1, 10))
  # after K, after S, N-terminus, after initial M, N-terminus, after S
  expect_equal(res$matched, c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE))
  expect_equal(unique(res$lengthPeptide), 8)
})

test_that("annotatePeptides() uses prefix and suffix", {
  res <- annotatePeptides("PEPTIDEK", toy_fasta(), prefix = "", suffix = "A")
  res <- res[order(res$proteinID, res$Offset), ]

  # only P3 and P4 continue with A after the peptide
  expect_equal(res$matched, c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE))
})

test_that("annotatePeptides() keeps one row per peptide row and protein occurrence", {
  pepinfo <- data.frame(peptideSeq = c("PEPTIDEK", "PEPTIDEK", "GGGR"), charge = c(2, 3, 3))

  res <- annotatePeptides(pepinfo, toy_fasta())

  expect_equal(nrow(res), 2 * 6 + 1)
  expect_equal(res$charge[res$peptideSeq == "GGGR"], 3)
  expect_equal(as.vector(table(res$charge[res$peptideSeq == "PEPTIDEK"])), c(6, 6))
})

test_that("annotatePeptides() reads the peptide column given by `peptide`", {
  pepinfo <- data.frame(seq = c("PEPTIDEK", "GGGR"))

  res <- annotatePeptides(pepinfo, toy_fasta(), peptide = "seq")

  expect_setequal(res$seq, c("PEPTIDEK", "GGGR"))
  expect_equal(nrow(res), 7)
})

test_that("annotatePeptides() warns and returns NULL without matches", {
  expect_warning(res <- annotatePeptides("WWWWW", toy_fasta()), "no matches found")
  expect_null(res)
})
