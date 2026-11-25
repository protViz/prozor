test_that("create_decoy_db() puts the source entry first, then targets, contaminants and decoys", {
  file <- toy_db_file()

  expect_message(
    res <- create_decoy_db(file, contaminants = toy_contaminants(), header = "zz|src|x"),
    "Reading database"
  )

  expect_equal(as.character(res[[1]]), "CRAPCRAPCRAP")
  expect_equal(attributes(res[[1]])[c("name", "Annot")], list(name = "zz|src|x", Annot = ">zz|src|x"))
  expect_equal(names(res)[-1], c("x", "y", "c1", "REV_x", "REV_y", "REV_c1"))
  expect_equal(as.character(res$REV_x), "DCBA")
  expect_equal(as.character(res$REV_c1), "KIH")
})

test_that("create_decoy_db() labels decoys with rev_lab", {
  file <- toy_db_file()

  res <- suppressMessages(create_decoy_db(file, contaminants = NULL, rev_lab = "DECOY_"))

  expect_equal(names(res)[-1], c("x", "y", "DECOY_x", "DECOY_y"))
})

test_that("create_decoy_db() skips decoys with rev_lab = NULL", {
  file <- toy_db_file()

  res <- suppressMessages(create_decoy_db(file, contaminants = toy_contaminants(), rev_lab = NULL))

  expect_equal(names(res)[-1], c("x", "y", "c1"))
})

test_that("create_decoy_db() skips contaminants with contaminants = NULL", {
  file <- toy_db_file()

  res <- suppressMessages(create_decoy_db(file, contaminants = NULL, rev_lab = NULL))

  expect_equal(names(res)[-1], c("x", "y"))
})

test_that("create_decoy_db() reads every file given", {
  file <- toy_db_file()
  other <- withr::local_tempfile(fileext = ".fasta")
  write_fasta(list(seq_fasta_aa("LMN", "z")), other)

  res <- suppressMessages(create_decoy_db(c(file, other), contaminants = NULL, rev_lab = NULL))

  expect_equal(names(res)[-1], c("x", "y", "z"))
})

test_that("create_decoy_db() keeps only the first entry of each name", {
  file <- toy_db_file()
  cont <- list(x = seq_fasta_aa("HIK", "x"))

  res <- suppressMessages(create_decoy_db(c(file, file), contaminants = cont))

  expect_equal(names(res)[-1], c("x", "y", "REV_x", "REV_y"))
  expect_equal(as.character(res$x), "ABCD")
  expect_equal(as.character(res$REV_x), "DCBA")
})
