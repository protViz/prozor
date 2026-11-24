test_that("annotate_aho() finds every occurrence with 1-based offsets", {
  expect_no_message(res <- annotate_aho(c("PEPTIDEK", "GGGR"), toy_fasta()))
  res <- res[order(res$peptide_seq, res$protein_id, res$offset), ]

  expect_equal(res$peptide_seq, c("GGGR", rep("PEPTIDEK", 6)))
  expect_equal(res$protein_id, c("P1", "P1", "P2", "P3", "P4", "P5", "P5"))
  expect_equal(res$offset, c(14, 6, 4, 1, 2, 1, 10))
  expect_equal(res$protein_sequence[res$protein_id == "P2"], "MSSPEPTIDEKLLL")
})

test_that("annotate_aho() rejects peptides with characters other than upper-case letters", {
  # AhoCorasickTrie 0.1.3 crashes R on such peptides
  expect_error(annotate_aho(c("PEPTIDEK", "peptidek"), toy_fasta()), "upper-case letters only.*\"peptidek\"")
  expect_error(annotate_aho("AC[+57.0]K", toy_fasta()), "\"AC\\[\\+57\\.0\\]K\"")
})

test_that("annotate_aho() rejects duplicated protein names", {
  fasta <- c(toy_fasta(), toy_fasta()["P1"])

  expect_error(annotate_aho("PEPTIDEK", fasta), "must be unique; duplicated: P1")
})

test_that("annotate_peptides() flags each occurrence by the residue before it", {
  expect_no_warning(expect_no_message(res <- annotate_peptides("PEPTIDEK", toy_fasta())))
  res <- res[order(res$protein_id, res$offset), ]

  expect_equal(res$protein_id, c("P1", "P2", "P3", "P4", "P5", "P5"))
  expect_equal(res$offset, c(6, 4, 1, 2, 1, 10))
  # after K, after S, N-terminus, after initial M, N-terminus, after S
  expect_equal(res$matched, c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE))
  expect_equal(unique(res$peptide_length), 8)
})

test_that("annotate_peptides() uses prefix and suffix", {
  res <- annotate_peptides("PEPTIDEK", toy_fasta(), prefix = "", suffix = "A")
  res <- res[order(res$protein_id, res$offset), ]

  # only P3 and P4 continue with A after the peptide
  expect_equal(res$matched, c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE))
})

test_that("annotate_peptides() keeps one row per peptide row and protein occurrence", {
  pepinfo <- data.frame(peptide_seq = c("PEPTIDEK", "PEPTIDEK", "GGGR"), charge = c(2, 3, 3))

  res <- annotate_peptides(pepinfo, toy_fasta())

  expect_equal(nrow(res), 2 * 6 + 1)
  expect_equal(res$charge[res$peptide_seq == "GGGR"], 3)
  expect_equal(as.vector(table(res$charge[res$peptide_seq == "PEPTIDEK"])), c(6, 6))
})

test_that("annotate_peptides() gives the same result for a tibble and a data frame", {
  pepinfo <- data.frame(peptide_seq = c("PEPTIDEK", "GGGR"), charge = c(2, 3))

  from_df <- annotate_peptides(pepinfo, toy_fasta())
  from_tibble <- annotate_peptides(dplyr::as_tibble(pepinfo), toy_fasta())

  expect_equal(as.data.frame(from_tibble), from_df)
})

test_that("annotate_peptides() reads the peptide column given by `peptide`", {
  pepinfo <- data.frame(seq = c("PEPTIDEK", "GGGR"))

  res <- annotate_peptides(pepinfo, toy_fasta(), peptide = "seq")

  expect_setequal(res$seq, c("PEPTIDEK", "GGGR"))
  expect_equal(nrow(res), 7)
})

test_that("annotate_peptides() warns and returns NULL without matches", {
  expect_warning(res <- annotate_peptides("WWWWW", toy_fasta()), "None of the peptides")
  expect_null(res)
})
