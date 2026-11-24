test_that("createDecoyDB() adds contaminants, decoys and one source entry", {
  file <- system.file("extdata/fgcz_contaminants2021_20210929.fasta.gz", package = "prozor")
  n_db <- length(readPeptideFasta(file))
  cont <- loadContaminantsFGCZ2022()

  # createDecoyDB() reports each file it reads as a message
  full <- suppressMessages(createDecoyDB(file, useContaminants = cont))
  no_decoy <- suppressMessages(createDecoyDB(file, useContaminants = cont, revLab = NULL))
  db_only <- suppressMessages(createDecoyDB(file, useContaminants = NULL, revLab = NULL))

  expect_length(full, 1 + 2 * (n_db + length(cont)))
  expect_equal(sum(grepl("^REV_", names(full))), n_db + length(cont))
  expect_length(no_decoy, 1 + n_db + length(cont))
  expect_length(db_only, 1 + n_db)
})
