test_that("read_fasta() returns one SeqFastaAA entry per header", {
  file <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c("text before the first header", ">sp|P1|A desc", "MKV", ";comment", "LLL", ">P2\ttab", "PP"), file)

  res <- read_fasta(file)

  expect_equal(res[["sp|P1|A"]], seq_fasta_aa("MKVLLL", "sp|P1|A", ">sp|P1|A desc"))
  # a tab ends the name too
  expect_equal(res[["P2"]], seq_fasta_aa("PP", "P2", ">P2\ttab"))
  expect_named(res, c("sp|P1|A", "P2"))
})

test_that("read_fasta() reads a file without a final newline without a warning", {
  file <- withr::local_tempfile(fileext = ".fasta")
  cat(">a", "MKV", sep = "\n", file = file)

  expect_no_warning(res <- read_fasta(file))

  expect_equal(res, list(a = seq_fasta_aa("MKV", "a")))
})

test_that("read_fasta() gives an entry without sequence lines the empty string", {
  file <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">a", ">b", "PP", ">c"), file)

  res <- read_fasta(file)

  expect_equal(unlist(res, use.names = FALSE), c("", "PP", ""))
  expect_named(res, c("a", "b", "c"))
})

test_that("read_fasta() stops for a file without a header", {
  file <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c("MKV", "LLL"), file)

  expect_error(read_fasta(file), "`file` must be a FASTA file with a line starting with `>`")
})
