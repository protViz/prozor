toy_annotation <- function() {
  # P: MAAAKLLLRCGGGCWWW; AAAK covers 2-5, AAKLLLR 3-9, GGG 11-13
  data.frame(
    protein_id = "P",
    peptide_seq = c("AAAK", "AAKLLLR", "GGG"),
    offset = c(2, 3, 11),
    protein_sequence = "MAAAKLLLRCGGGCWWW"
  )
}

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

test_that("a peptide occurring twice in a protein covers both places", {
  annotation <- data.frame(protein_id = "Q", peptide_seq = "GGK", offset = c(1, 8), protein_sequence = "GGKAAAAGGK")

  res <- residue_coverage(annotation)

  expect_equal(res$n_peptides, c(1L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 1L, 1L))
})

test_that("overlapping occurrences of one peptide count once", {
  fasta <- list(R = seq_fasta_aa("MGGGGGK", "R"))
  hits <- annotate_aho("GGG", fasta)

  expect_equal(hits$offset, c(2, 3, 4))
  expect_equal(residue_coverage(hits)$n_peptides, c(0L, 1L, 1L, 1L, 1L, 1L, 0L))
  expect_equal(sequence_coverage(hits)$n_peptides, 1L)
})

test_that("repeated rows of one occurrence, e.g. charge states, count once", {
  fasta <- list(P = seq_fasta_aa("MAAAKLLLRCGGGCWWW", "P"))
  pepinfo <- data.frame(peptide_seq = c("AAAK", "AAAK", "GGG"), precursor_charge = c(2L, 3L, 2L))

  annotation <- annotate_peptides(pepinfo, fasta)

  expect_equal(nrow(annotation), 3)
  expect_equal(
    residue_coverage(annotation)$n_peptides,
    c(0L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L)
  )
  expect_equal(sequence_coverage(annotation)$n_peptides, 2L)
})

test_that("annotate_aho() results and other peptide column names work", {
  fasta <- list(P = seq_fasta_aa("MAAAKLLLRCGGGCWWW", "P"))
  hits <- annotate_aho(c("AAAK", "AAKLLLR", "GGG"), fasta, peptide = "sequence")

  expect_equal(
    residue_coverage(hits, peptide = "sequence")$n_peptides,
    residue_coverage(toy_annotation())$n_peptides
  )
})

test_that("missing columns and unknown proteins are errors", {
  expect_error(residue_coverage(toy_annotation()[, -3]), "missing: offset")
  expect_error(sequence_coverage(toy_annotation(), peptide = "sequence"), "missing: sequence")
  expect_error(plot_coverage(residue_coverage(toy_annotation()), "X"), "`protein_id` must be one protein")
  expect_error(plot_coverage(toy_annotation(), "P"), "missing columns: position, residue, n_peptides")
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
