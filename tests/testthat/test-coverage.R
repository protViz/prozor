test_that("residue_coverage() counts the peptides covering each residue", {
  res <- residue_coverage(toy_annotation())

  expect_named(res, c("protein_id", "position", "residue", "n_peptides"))
  expect_equal(res$position, 1:17)
  expect_equal(paste(res$residue, collapse = ""), "MAAAKLLLRCGGGCWWW")
  expect_equal(res$n_peptides, c(0L, 1L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L))
})

test_that("sequence_coverage() summarizes the coverage per protein", {
  annotation <- rbind(
    toy_annotation(),
    data.frame(protein_id = "Q", peptide_seq = "GGK", offset = c(1, 8), protein_sequence = "GGKAAAAGGK")
  )

  res <- sequence_coverage(annotation)

  expect_equal(res$protein_id, c("P", "Q"))
  expect_equal(res$protein_length, c(17L, 10L))
  expect_equal(res$n_peptides, c(3L, 1L))
  expect_equal(res$n_covered, c(11L, 6L))
  expect_equal(res$coverage, c(11 / 17, 6 / 10))
})

test_that("residue_coverage() counts a peptide occurring twice in a protein at both places", {
  annotation <- data.frame(protein_id = "Q", peptide_seq = "GGK", offset = c(1, 8), protein_sequence = "GGKAAAAGGK")

  res <- residue_coverage(annotation)

  expect_equal(res$n_peptides, c(1L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 1L, 1L))
})

test_that("residue_coverage() and sequence_coverage() count overlapping occurrences of one peptide once", {
  fasta <- list(R = seq_fasta_aa("MGGGGGK", "R"))
  hits <- annotate_aho("GGG", fasta)

  expect_equal(hits$offset, c(2, 3, 4))
  expect_equal(residue_coverage(hits)$n_peptides, c(0L, 1L, 1L, 1L, 1L, 1L, 0L))
  expect_equal(sequence_coverage(hits)$n_peptides, 1L)
})

test_that("residue_coverage() and sequence_coverage() count repeated rows of one occurrence once", {
  fasta <- list(P = seq_fasta_aa("MAAAKLLLRCGGGCWWW", "P"))
  # AAAK with two charge states
  peptides <- data.frame(peptide_seq = c("AAAK", "AAAK", "GGG"), precursor_charge = c(2L, 3L, 2L))

  annotation <- annotate_peptides(peptides, fasta)

  expect_equal(nrow(annotation), 3)
  expect_equal(
    residue_coverage(annotation)$n_peptides,
    c(0L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L)
  )
  expect_equal(sequence_coverage(annotation)$n_peptides, 2L)
})

test_that("residue_coverage() works on annotate_aho() results with another peptide column name", {
  fasta <- list(P = seq_fasta_aa("MAAAKLLLRCGGGCWWW", "P"))
  hits <- annotate_aho(c("AAAK", "AAKLLLR", "GGG"), fasta, peptide = "sequence")

  expect_equal(
    residue_coverage(hits, peptide = "sequence")$n_peptides,
    residue_coverage(toy_annotation())$n_peptides
  )
})

test_that("residue_coverage() and sequence_coverage() return 0 rows with the result columns for 0 rows", {
  empty <- toy_annotation()[0, ]

  residues <- residue_coverage(empty)
  coverage <- sequence_coverage(empty)

  expect_equal(nrow(residues), 0)
  expect_named(residues, c("protein_id", "position", "residue", "n_peptides"))
  expect_identical(lapply(residues, class), lapply(residue_coverage(toy_annotation()), class))
  expect_equal(nrow(coverage), 0)
  expect_named(coverage, c("protein_id", "protein_length", "n_peptides", "n_covered", "coverage"))
  expect_identical(lapply(coverage, class), lapply(sequence_coverage(toy_annotation()), class))
})

test_that("residue_coverage(), sequence_coverage() and plot_coverage() stop for missing columns and unknown proteins", {
  expect_error(residue_coverage(toy_annotation()[, -3]), "missing: `offset`\\.")
  expect_error(sequence_coverage(toy_annotation(), peptide = "sequence"), "missing: `sequence`\\.")
  expect_error(plot_coverage(residue_coverage(toy_annotation()), "X"), "`protein_id` must be one protein")
  expect_error(plot_coverage(toy_annotation(), "P"), "missing columns: `position`, `residue`, `n_peptides`\\.")
})

test_that("residue_coverage() errors do not show the internal check as the call", {
  err <- expect_error(residue_coverage(toy_annotation()[, -3]), "`annotation` must have the columns")

  expect_null(conditionCall(err))
})

test_that("plot_coverage() stops unless width is a positive whole number", {
  residues <- residue_coverage(toy_annotation())

  expect_error(plot_coverage(residues, "P", width = 0), "`width` must be a single positive whole number, not 0\\.")
  expect_error(plot_coverage(residues, "P", width = 2.5), "`width` must be a single positive whole number, not 2.5\\.")
})

test_that("plot_coverage() draws one tile per residue in rows of `width`", {
  residues <- residue_coverage(toy_annotation())

  p <- plot_coverage(residues, "P", width = 5)
  tiles <- ggplot2::layer_data(p, 1)

  expect_s3_class(p, "ggplot")
  expect_equal(nrow(tiles), 17)
  expect_length(unique(tiles$y), 4)
  expect_equal(max(tiles$x), 5)
  # uncovered residues are grey, the doubly covered ones get the second blue
  expect_equal(unique(tiles$fill[residues$n_peptides == 0]), "grey90")
  expect_equal(unique(tiles$fill[residues$n_peptides == 2]), "#6BAED6")
  expect_match(p$labels$title, "P: 65 % covered")
})
