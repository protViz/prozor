test_that("create_fgcz_fasta_db() writes a dated decoy database to output_dir", {
  outdir <- withr::local_tempdir()
  example <- system.file("p1000_db1_example", package = "prozor")

  # progress messages for each file read and written
  res <- suppressMessages(create_fgcz_fasta_db(example, output_dir = outdir, make_summary = FALSE))
  back <- read_peptide_fasta(res$file_path)

  expect_equal(normalizePath(dirname(res$file_path)), normalizePath(outdir))
  expect_match(basename(res$file_path), "^p1000_db1_example_d_[0-9]{8}\\.fasta$")
  expect_equal(res$db_name, "p1000_db1_example_d")
  expect_length(res$db, 2941)
  expect_length(back, 2941)
  expect_false(anyDuplicated(names(back)) > 0)
})

test_that("create_fgcz_fasta_db() without rev_lab writes targets and contaminants only", {
  outdir <- withr::local_tempdir()
  example <- system.file("p1000_db1_example", package = "prozor")

  res <- suppressMessages(create_fgcz_fasta_db(example, rev_lab = NULL, output_dir = outdir, make_summary = FALSE))

  expect_match(basename(res$file_path), "^p1000_db1_example_[0-9]{8}\\.fasta$")
  expect_length(res$db, 1471)
  expect_false(any(grepl("^REV_", names(res$db))))
})

test_that("create_fgcz_fasta_db() appends the FASTA summary when make_summary = TRUE", {
  outdir <- withr::local_tempdir()
  example <- system.file("p1000_db1_example", package = "prozor")

  res <- suppressMessages(create_fgcz_fasta_db(example, output_dir = outdir, make_summary = TRUE))

  expect_true(any(grepl("AA frequencies", res$summary)))
})

test_that("create_fgcz_fasta_db() records the argument values in the summary", {
  database_dir <- toy_database_dir()
  outdir <- withr::local_tempdir()
  contaminants <- toy_contaminants()

  res <- suppressMessages(
    create_fgcz_fasta_db(
      database_dir,
      contaminants = contaminants,
      rev_lab = NULL,
      output_dir = outdir,
      make_summary = FALSE
    )
  )

  expect_match(res$summary, paste0("database_dir : ", database_dir, "\n"), fixed = TRUE)
  expect_match(res$summary, "contaminants : 1 sequences\n", fixed = TRUE)
  expect_match(res$summary, "rev_lab : NULL\n", fixed = TRUE)
  expect_match(res$summary, paste0("output_dir : ", outdir, "\n"), fixed = TRUE)
  expect_match(res$summary, "nr sequences : 4\n", fixed = TRUE)
  expect_no_match(res$summary, "contaminants = contaminants", fixed = TRUE)
})

test_that("create_fgcz_fasta_db() returns the entries as read_peptide_fasta() reads them from the file", {
  database_dir <- toy_database_dir()
  outdir <- withr::local_tempdir()

  res <- suppressMessages(
    create_fgcz_fasta_db(database_dir, contaminants = toy_contaminants(), output_dir = outdir, make_summary = FALSE)
  )

  expect_identical(unname(res$db), unname(read_peptide_fasta(res$file_path)))
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

test_that("create_fgcz_fasta_db() reads *.fa, *.faa and gzipped FASTA files", {
  example <- system.file("p1000_db1_example", package = "prozor")
  database_dir <- toy_database_dir("extra.faa")
  file.copy(file.path(example, "Annotation_canSeq.fasta.gz"), file.path(database_dir, "uniprot_mouse.fa.gz"))
  outdir <- withr::local_tempdir()

  res <- suppressMessages(
    create_fgcz_fasta_db(database_dir, contaminants = NULL, rev_lab = NULL, output_dir = outdir, make_summary = FALSE)
  )

  targets <- read_peptide_fasta(file.path(database_dir, "uniprot_mouse.fa.gz"))
  expect_equal(names(res$db)[-1], c("x", "y", names(targets)))
})

test_that("create_fgcz_fasta_db() reads only FASTA files and the exact annotation.txt", {
  database_dir <- toy_database_dir()
  writeLines(c("How this fasta was made:", "downloaded from UniProt"), file.path(database_dir, "fasta_README.txt"))
  writeLines("zz|old|annotation", file.path(database_dir, "annotation.txt.bak"))
  outdir <- withr::local_tempdir()

  res <- suppressMessages(
    create_fgcz_fasta_db(database_dir, contaminants = NULL, rev_lab = NULL, output_dir = outdir, make_summary = FALSE)
  )

  expect_equal(names(res$db)[-1], c("x", "y"))
  expect_equal(attr(res$db[[1]], "Annot"), paste0(">", readLines(file.path(database_dir, "annotation.txt"))))
})

test_that("create_fgcz_fasta_db() stops when database_dir has no FASTA file", {
  database_dir <- toy_database_dir()
  file.remove(file.path(database_dir, "targets.fasta"))
  writeLines(c("How this fasta was made:", "downloaded from UniProt"), file.path(database_dir, "fasta_README.txt"))
  outdir <- withr::local_tempdir()

  expect_error(
    create_fgcz_fasta_db(database_dir, output_dir = outdir),
    "`database_dir` must contain FASTA files .* none in \"[^\"]*toydb\""
  )
  expect_length(list.files(outdir), 0)
})

test_that("create_fgcz_fasta_db() refuses to write into database_dir", {
  database_dir <- toy_database_dir()

  expect_error(
    create_fgcz_fasta_db(database_dir, output_dir = database_dir),
    "`output_dir` must not be `database_dir`"
  )
  expect_error(
    create_fgcz_fasta_db(database_dir, output_dir = file.path(database_dir, ".")),
    "`output_dir` must not be `database_dir`"
  )
  expect_setequal(list.files(database_dir), c("annotation.txt", "targets.fasta"))
})

test_that("make_fasta_summary() counts sequences and amino acids", {
  fasta <- list(
    seq_fasta_aa("AAC", "a"),
    seq_fasta_aa("CD", "b")
  )

  res <- make_fasta_summary(fasta, as_string = FALSE)

  expect_equal(res$nr_sequences, 2)
  expect_equal(res$aa_freq, c(A = 2, C = 2, D = 1))
  expect_equal(unname(res$length_summary[c("Min.", "Max.")]), c(2, 3))
  expect_type(make_fasta_summary(fasta), "character")
})

test_that("make_fasta_summary() counts every character and empty sequences", {
  fasta <- list(seq_fasta_aa("MKm*", "a"), seq_fasta_aa("", "b"), seq_fasta_aa("KK", "c"))

  res <- make_fasta_summary(fasta, as_string = FALSE)

  expect_equal(res$aa_freq[c("M", "K", "m", "*")], c(M = 1, K = 3, m = 1, "*" = 1))
  expect_length(res$aa_freq, 4)
  expect_equal(unname(res$length_summary[c("Min.", "Max.")]), c(0, 4))
})
