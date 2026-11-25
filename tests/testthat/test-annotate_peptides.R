test_that("annotate_aho() finds every occurrence with 1-based offsets", {
  expect_no_message(res <- annotate_aho(c("PEPTIDEK", "GGGR"), toy_fasta()))
  res <- res[order(res$peptide_seq, res$protein_id, res$offset), ]

  expect_equal(res$peptide_seq, c("GGGR", rep("PEPTIDEK", 6)))
  expect_equal(res$protein_id, c("P1", "P1", "P2", "P3", "P4", "P5", "P5"))
  expect_equal(res$offset, c(14, 6, 4, 1, 2, 1, 10))
  expect_equal(res$protein_sequence[res$protein_id == "P2"], "MSSPEPTIDEKLLL")
})

test_that("annotate_aho() rejects peptides with characters other than upper-case letters", {
  # AhoCorasickTrie 0.1.3 crashes R on such peptides, 0.1.4 gives a less specific error
  expect_error(
    annotate_aho(c("PEPTIDEK", "peptidek"), toy_fasta()),
    "`peptides` must consist of upper-case letters only.*\"peptidek\""
  )
  expect_error(annotate_aho("AC[+57.0]K", toy_fasta()), "\"AC\\[\\+57\\.0\\]K\"")
})

test_that("annotate_aho() rejects empty and NA peptides", {
  # AhoCorasickTrie gives wrong hits for them: "" at most positions, NA as the string "NA"
  expect_error(annotate_aho(c("PEPTIDEK", ""), toy_fasta()), "upper-case letters only.*; found: \"\"\\.")
  expect_error(annotate_aho(c("PEPTIDEK", NA_character_), toy_fasta()), "upper-case letters only.*; found: NA\\.")
})

test_that("annotate_aho() gives offsets in the sequences of `fasta`, also after blanks", {
  res <- annotate_aho("PEPTIDEK", list(P1 = "  MKPEPTIDEK"))

  expect_equal(res$offset, 5)
  expect_equal(res$protein_sequence, "XXMKPEPTIDEK")
})

test_that("annotate_aho() stops for a `fasta` without names", {
  expect_error(annotate_aho("PEPTIDEK", unname(toy_fasta())), "`fasta` must be a named list .*; it has no names\\.")
})

test_that("annotate_aho() rejects peptides with whitespace and shows it escaped", {
  expect_error(annotate_aho(c("PEPTIDEK", "GGGR\r"), toy_fasta()), "without whitespace.*\"GGGR\\\\r\"")
  expect_error(annotate_aho(" PEPTIDEK", toy_fasta()), "\" PEPTIDEK\"")
})

test_that("annotate_aho() searches proteins with lower-case residues or characters other than letters", {
  # a stop codon (*) of a translated database, soft-masked (lower-case) residues
  fasta <- list(P1 = "M*KPEPTIDEKaaa*", P2 = "MRpeptidekR")

  res <- annotate_aho("PEPTIDEK", fasta)

  expect_equal(res$protein_id, c("P1", "P2"))
  expect_equal(res$offset, c(4, 3))
  expect_equal(res$protein_sequence, c("MXKPEPTIDEKAAAX", "MRPEPTIDEKR"))
  expect_equal(annotate_peptides("PEPTIDEK", fasta)$matched, c(TRUE, TRUE))
})

test_that("annotate_aho() rejects duplicated protein names", {
  fasta <- c(toy_fasta(), toy_fasta()["P1"])

  expect_error(annotate_aho("PEPTIDEK", fasta), "must be unique; duplicated: \"P1\"\\.")
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
  peptides <- data.frame(peptide_seq = c("PEPTIDEK", "PEPTIDEK", "GGGR"), charge = c(2, 3, 3))

  res <- annotate_peptides(peptides, toy_fasta())

  expect_equal(nrow(res), 2 * 6 + 1)
  expect_equal(res$charge[res$peptide_seq == "GGGR"], 3)
  expect_equal(as.vector(table(res$charge[res$peptide_seq == "PEPTIDEK"])), c(6, 6))
})

test_that("annotate_peptides() trims whitespace around the peptides", {
  # e.g. "\r" at the end of the last column of a file with Windows line endings
  peptides <- data.frame(peptide_seq = c("PEPTIDEK\r", " GGGR"), charge = c(2, 3))

  res <- annotate_peptides(peptides, toy_fasta())

  expect_equal(nrow(res), 6 + 1)
  expect_setequal(res$peptide_seq, c("PEPTIDEK", "GGGR"))
  expect_equal(res$peptide_length[res$peptide_seq == "GGGR"], 4)
  expect_equal(res$charge[res$peptide_seq == "GGGR"], 3)
})

test_that("annotate_peptides() rejects `peptides` columns named like result columns", {
  # e.g. search engine output with a protein column
  psm <- data.frame(peptide_seq = "PEPTIDEK", protein_id = "P1", offset = 99, matched = TRUE, peptide_length = 8)

  expect_error(
    annotate_peptides(psm, toy_fasta()),
    "`peptides` must not have columns named like result .*; rename `protein_id`, `offset`, `matched`, `peptide_length`\\."
  )
  expect_error(
    annotate_peptides("PEPTIDEK", toy_fasta(), peptide = "offset"),
    "`peptides` must not have columns .*; rename `offset`\\."
  )
})

test_that("annotate_peptides() gives the same result for a tibble and a data frame", {
  skip_if_not_installed("dplyr")
  peptides <- data.frame(peptide_seq = c("PEPTIDEK", "GGGR"), charge = c(2, 3))

  from_df <- annotate_peptides(peptides, toy_fasta())
  from_tibble <- annotate_peptides(dplyr::as_tibble(peptides), toy_fasta())

  expect_equal(as.data.frame(from_tibble), from_df)
})

test_that("annotate_peptides() stops when `peptides` lacks the column `peptide`", {
  peptides <- data.frame(seq = c("PEPTIDEK", "GGGR"))

  expect_error(annotate_peptides(peptides, toy_fasta()), "`peptides` must have the columns `peptide_seq`; missing:")
})

test_that("annotate_peptides() reads the peptide column given by `peptide`", {
  peptides <- data.frame(seq = c("PEPTIDEK", "GGGR"))

  res <- annotate_peptides(peptides, toy_fasta(), peptide = "seq")

  expect_setequal(res$seq, c("PEPTIDEK", "GGGR"))
  expect_equal(nrow(res), 7)
})

test_that("annotate_aho() returns 0 rows with the result columns without matches", {
  expect_no_warning(res <- annotate_aho("WWWWW", toy_fasta(), peptide = "sequence"))

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0)
  expect_named(res, c("protein_id", "sequence", "offset", "protein_sequence"))
})

test_that("annotate_peptides() warns and returns 0 rows with the result columns without matches", {
  peptides <- data.frame(peptide_seq = c("WWWWW", "YYYYY"), charge = c(2, 3))

  expect_warning(res <- annotate_peptides(peptides, toy_fasta()), "None of the peptides")

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0)
  expect_named(res, c("peptide_seq", "protein_id", "offset", "protein_sequence", "matched", "charge", "peptide_length"))
  # downstream functions take the empty result
  expect_equal(nrow(residue_coverage(res)), 0)
  expect_equal(nrow(sequence_coverage(res)), 0)
})
