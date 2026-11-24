test_that("createDecoyDB() puts the source entry first, then targets, contaminants and decoys", {
  file <- toy_db_file()

  expect_message(res <- createDecoyDB(file, useContaminants = toy_contaminants(), annot = "zz|src|x"), "reading db")

  expect_equal(as.character(res[[1]]), "CRAPCRAPCRAP")
  expect_equal(seqinr::getAnnot(res[[1]]), "zz|src|x")
  expect_equal(names(res)[-1], c("x", "y", "c1", "REV_x", "REV_y", "REV_c1"))
  expect_equal(as.character(res$REV_x), "DCBA")
  expect_equal(as.character(res$REV_c1), "KIH")
})

test_that("createDecoyDB() labels decoys with revLab", {
  file <- toy_db_file()

  res <- suppressMessages(createDecoyDB(file, useContaminants = NULL, revLab = "DECOY_"))

  expect_equal(names(res)[-1], c("x", "y", "DECOY_x", "DECOY_y"))
})

test_that("createDecoyDB() skips decoys with revLab = NULL", {
  file <- toy_db_file()

  res <- suppressMessages(createDecoyDB(file, useContaminants = toy_contaminants(), revLab = NULL))

  expect_equal(names(res)[-1], c("x", "y", "c1"))
})

test_that("createDecoyDB() skips contaminants with useContaminants = NULL", {
  file <- toy_db_file()

  res <- suppressMessages(createDecoyDB(file, useContaminants = NULL, revLab = NULL))

  expect_equal(names(res)[-1], c("x", "y"))
})

test_that("createDecoyDB() reads every file given", {
  file <- toy_db_file()

  res <- suppressMessages(createDecoyDB(c(file, file), useContaminants = NULL, revLab = NULL))

  expect_equal(names(res)[-1], c("x", "y", "x", "y"))
})
