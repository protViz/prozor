test_that("fgcz_create_fasta.R adds the defaults: REV_ decoys and the special contaminants, with --summary", {
  skip_on_cran()
  skip_if_not_installed("docopt")
  skip_if_not(is_checking(), "needs the installed package (R CMD check)")
  database_dir <- toy_database_dir()
  outdir <- withr::local_tempdir()

  out <- run_fgcz_create_fasta(c(database_dir, paste0("--output_dir=", outdir), "--summary"))

  expect_null(attr(out, "status"))
  db <- read_fasta(list.files(outdir, "^toydb_d_[0-9]{8}\\.fasta$", full.names = TRUE))
  expect_true(all(c("x", "y", "REV_x", "REV_y") %in% names(db)))
  expect_true(all(names(load_contaminants("fgcz2024")) %in% names(db)))
  # the parameter log is written next to the FASTA file, with the same name
  log <- readLines(list.files(outdir, "^toydb_d_[0-9]{8}\\.txt$", full.names = TRUE))
  expect_true(any(grepl("revLab: REV_", log, fixed = TRUE)))
  expect_true(any(grepl("contamin: special", log, fixed = TRUE)))
  expect_true(any(grepl("AA frequencies:", log, fixed = TRUE)))
})

test_that("fgcz_create_fasta.R takes a set of load_contaminants() and logs nodecoy without revLab", {
  skip_on_cran()
  skip_if_not_installed("docopt")
  skip_if_not(is_checking(), "needs the installed package (R CMD check)")
  database_dir <- toy_database_dir()
  outdir <- withr::local_tempdir()

  out <- run_fgcz_create_fasta(c(database_dir, "nodecoy", "--contamin=universal2024", paste0("--output_dir=", outdir)))

  expect_null(attr(out, "status"))
  db <- read_fasta(list.files(outdir, "^toydb_[0-9]{8}\\.fasta$", full.names = TRUE))
  expect_length(db, 1 + 2 + length(load_contaminants("universal2024")))
  expect_false(any(grepl("^REV_", names(db))))
  log <- readLines(list.files(outdir, "^toydb_[0-9]{8}\\.txt$", full.names = TRUE))
  expect_true(any(grepl("revLab: NULL", log, fixed = TRUE)))
  expect_false(any(grepl("revLab: REV_", log, fixed = TRUE)))
})

test_that("fgcz_create_fasta.R takes the contaminants from a FASTA file", {
  skip_on_cran()
  skip_if_not_installed("docopt")
  skip_if_not(is_checking(), "needs the installed package (R CMD check)")
  database_dir <- toy_database_dir()
  outdir <- withr::local_tempdir()
  contaminant_file <- withr::local_tempfile(fileext = ".fasta")
  write_fasta(toy_contaminants(), contaminant_file)

  out <- run_fgcz_create_fasta(
    c(database_dir, "nodecoy", paste0("--contamin=", contaminant_file), paste0("--output_dir=", outdir))
  )

  expect_null(attr(out, "status"))
  db <- read_fasta(list.files(outdir, "^toydb_[0-9]{8}\\.fasta$", full.names = TRUE))
  expect_named(db, c("zz|p1000_toydb|20261124", "x", "y", "c1"))
})

test_that("fgcz_create_fasta.R names the contaminant sets for an unknown --contamin", {
  skip_on_cran()
  skip_if_not_installed("docopt")
  skip_if_not(is_checking(), "needs the installed package (R CMD check)")
  database_dir <- toy_database_dir()
  outdir <- withr::local_tempdir()

  out <- suppressWarnings(
    run_fgcz_create_fasta(c(database_dir, "--contamin=fgcz2023", paste0("--output_dir=", outdir)))
  )

  expect_equal(attr(out, "status"), 1L)
  expect_match(paste(out, collapse = "\n"), "must be one of .*\"universal2024\", \"special2024\", not \"fgcz2023\"")
  expect_length(list.files(outdir, "\\.fasta$"), 0)
})
