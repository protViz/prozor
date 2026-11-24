# inst/script/fgcz_create_fasta.R loads the installed prozor, which is the working tree only under R CMD check
run_fgcz_create_fasta <- function(args) {
  script <- system.file("script", "fgcz_create_fasta.R", package = "prozor")
  system2(file.path(R.home("bin"), "Rscript"), shQuote(c(script, args)), stdout = TRUE, stderr = TRUE)
}

test_that("fgcz_create_fasta.R takes every bundled contaminant set and logs nodecoy without revLab", {
  skip_on_cran()
  skip_if_not_installed("docopt")
  skip_if_not(is_checking(), "needs the installed package (R CMD check)")
  database_dir <- toy_database_dir()
  outdir <- withr::local_tempdir()
  # the script writes its parameter log <db_name>.txt into the working directory
  withr::local_dir(outdir)

  out <- run_fgcz_create_fasta(c(database_dir, "nodecoy", "--contamin=universal2024", paste0("--output_dir=", outdir)))

  expect_null(attr(out, "status"))
  db <- read_peptide_fasta(list.files(outdir, "^toydb_[0-9]{8}\\.fasta$", full.names = TRUE))
  expect_length(db, 1 + 2 + length(load_contaminants("universal2024")))
  expect_false(any(grepl("^REV_", names(db))))
  log <- readLines(file.path(outdir, "toydb.txt"))
  expect_true(any(grepl("revLab: NULL", log, fixed = TRUE)))
  expect_false(any(grepl("revLab: REV_", log, fixed = TRUE)))
})

test_that("fgcz_create_fasta.R names the contaminant sets for an unknown --contamin", {
  skip_on_cran()
  skip_if_not_installed("docopt")
  skip_if_not(is_checking(), "needs the installed package (R CMD check)")
  database_dir <- toy_database_dir()
  outdir <- withr::local_tempdir()
  withr::local_dir(outdir)

  out <- suppressWarnings(
    run_fgcz_create_fasta(c(database_dir, "--contamin=fgcz2023", paste0("--output_dir=", outdir)))
  )

  expect_equal(attr(out, "status"), 1L)
  expect_match(paste(out, collapse = "\n"), "must be one of .*\"universal2024\", \"special2024\", not \"fgcz2023\"")
  expect_length(list.files(outdir, "\\.fasta$"), 0)
})
