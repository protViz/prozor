test_that("create_fgcz_fasta_db() writes a dated decoy database to output_dir", {
  outdir <- withr::local_tempdir()
  example <- system.file("p1000_db1_example", package = "prozor")

  # progress messages for each file read and written
  res <- suppressMessages(create_fgcz_fasta_db(example, output_dir = outdir, make_summary = FALSE))
  back <- read_peptide_fasta(res$filepath)

  expect_equal(dirname(res$filepath), outdir)
  expect_match(basename(res$filepath), "^p1000_db1_example_d_[0-9]{8}\\.fasta$")
  expect_equal(res$dbname, "p1000_db1_example_d")
  expect_length(res$resDB, 3155)
  expect_length(back, 3155)
  expect_false(anyDuplicated(names(back)) > 0)
})

test_that("create_fgcz_fasta_db() without rev_lab writes targets and contaminants only", {
  outdir <- withr::local_tempdir()
  example <- system.file("p1000_db1_example", package = "prozor")

  res <- suppressMessages(create_fgcz_fasta_db(example, rev_lab = NULL, output_dir = outdir, make_summary = FALSE))

  expect_match(basename(res$filepath), "^p1000_db1_example_[0-9]{8}\\.fasta$")
  expect_length(res$resDB, 1578)
  expect_false(any(grepl("^REV_", names(res$resDB))))
})

test_that("create_fgcz_fasta_db() appends the FASTA summary when make_summary = TRUE", {
  outdir <- withr::local_tempdir()
  example <- system.file("p1000_db1_example", package = "prozor")

  res <- suppressMessages(create_fgcz_fasta_db(example, output_dir = outdir, make_summary = TRUE))

  expect_true(any(grepl("AA frequencies", res$summary)))
})

test_that("create_fgcz_fasta_db() requires output_dir, an existing directory and an annotation file", {
  outdir <- withr::local_tempdir()
  example <- system.file("p1000_db1_example", package = "prozor")
  no_annotation <- withr::local_tempdir()
  file.copy(file.path(example, "Annotation_canSeq.fasta.gz"), no_annotation)

  expect_error(suppressMessages(create_fgcz_fasta_db(example)), "output_dir")
  expect_error(create_fgcz_fasta_db(file.path(outdir, "missing"), output_dir = outdir), "must be an existing directory")
  expect_error(create_fgcz_fasta_db(no_annotation, output_dir = outdir), "Can't find an annotation file")
})

test_that("make_fasta_summary() counts sequences and amino acids", {
  fasta <- list(
    seqinr::as.SeqFastaAA("AAC", name = "a", Annot = ">a"),
    seqinr::as.SeqFastaAA("CD", name = "b", Annot = ">b")
  )

  res <- make_fasta_summary(fasta, as_string = FALSE)
  old <- make_fasta_summary(fasta, old = TRUE, as_string = FALSE)

  expect_equal(res$nrSequences, 2)
  expect_equal(res$aafreq, c(A = 2, C = 2, D = 1))
  expect_equal(c(old$aafreq), c(A = 2L, C = 2L, D = 1L))
  expect_equal(unname(res$lengthSummary[c("Min.", "Max.")]), c(2, 3))
  expect_type(make_fasta_summary(fasta), "character")
})
