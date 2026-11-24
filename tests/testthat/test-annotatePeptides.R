toy_fasta <- function() {
  list(
    P1 = seqinr::as.SeqFastaAA("MAAAKPEPTIDEKGGGR", name = "P1", Annot = ">P1"),
    P2 = seqinr::as.SeqFastaAA("MSSPEPTIDEKLLL", name = "P2", Annot = ">P2"),
    P3 = seqinr::as.SeqFastaAA("PEPTIDEKAAA", name = "P3", Annot = ">P3")
  )
}

test_that("annotateAHO() finds every occurrence with 1-based offsets", {
  # dplyr::inner_join() reports the join column as a message
  res <- suppressMessages(annotateAHO(c("PEPTIDEK", "GGGR"), toy_fasta()))
  res <- res[order(res$peptideSeq, res$proteinID), ]

  expect_equal(res$proteinID, c("P1", "P1", "P2", "P3"))
  expect_equal(res$peptideSeq, c("GGGR", "PEPTIDEK", "PEPTIDEK", "PEPTIDEK"))
  expect_equal(res$Offset, c(14, 6, 4, 1))
  expect_equal(res$proteinSequence[res$proteinID == "P2"], "MSSPEPTIDEKLLL")
})

test_that("annotatePeptides() flags tryptic matches and adds peptide length", {
  res <- suppressMessages(annotatePeptides(c("PEPTIDEK", "GGGR"), toy_fasta()))
  res <- res[order(res$peptideSeq, res$proteinID), ]

  expect_equal(res$proteinID, c("P1", "P1", "P2", "P3"))
  # P2 has S before the peptide; P3 starts with it; P1 has K before both
  expect_equal(res$matched, c(TRUE, TRUE, FALSE, TRUE))
  expect_equal(res$lengthPeptide, c(4, 8, 8, 8))
})

test_that("annotatePeptides() keeps extra peptide columns", {
  pepinfo <- data.frame(peptideSeq = c("PEPTIDEK", "GGGR"), charge = c(2, 3))

  res <- suppressMessages(annotatePeptides(pepinfo, toy_fasta()))

  expect_equal(res$charge[res$peptideSeq == "GGGR"], 3)
  expect_equal(unique(res$charge[res$peptideSeq == "PEPTIDEK"]), 2)
})

test_that("annotatePeptides() warns and returns NULL without matches", {
  expect_warning(res <- annotatePeptides("WWWWW", toy_fasta()), "no matches found")
  expect_null(res)
})
