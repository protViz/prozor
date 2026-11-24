test_that("removeSignalPeptide() cuts the signal peptide and renames the chain", {
  file <- system.file("extdata/fgcz_contaminants2021_20210929.fasta.gz", package = "prozor")
  db <- readPeptideFasta(file)[2:3]
  entry <- vapply(strsplit(names(db), split = "|", fixed = TRUE), function(x) x[2], character(1))
  signal <- data.frame(Entry = entry, Signal.peptide = c("SIGNAL 1 20", ""))

  res <- removeSignalPeptide(db, signal)

  expect_equal(nchar(res[[1]]), nchar(db[[1]]) - 20)
  expect_equal(as.character(res[[1]]), substring(as.character(db[[1]]), 21))
  expect_equal(attr(res[[1]], "name"), paste0(attr(db[[1]], "name"), "s21"))
  expect_equal(res[[2]], db[[2]])
})

test_that("removeSignalPeptide() works on a one-entry database", {
  file <- system.file("extdata/fgcz_contaminants2021_20210929.fasta.gz", package = "prozor")
  db <- readPeptideFasta(file)[2]
  entry <- strsplit(names(db), split = "|", fixed = TRUE)[[1]][2]
  signal <- data.frame(Entry = entry, Signal.peptide = "SIGNAL 1 20")

  res <- removeSignalPeptide(db, signal)

  expect_equal(nchar(res[[1]]), nchar(db[[1]]) - 20)
})

test_that("makeID() and makeIDUnip() build chain IDs", {
  expect_equal(makeID("ABCDE", "sp|P30443|1A01_HUMAN", 24), "sp|P30443|1A01_HUMANs24")
  expect_equal(makeIDUnip("ABCDE", "sp|P30443|1A01_HUMAN", 24), "sp|P30443|24-5")
})
