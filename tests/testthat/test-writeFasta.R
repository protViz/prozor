test_that("writeFasta() output reads back unchanged", {
  fasta <- loadContaminantsFGCZ2022()[2:6]
  file <- withr::local_tempfile(fileext = ".fasta")

  writeFasta(file, fasta)
  back <- readPeptideFasta(file)

  expect_equal(names(back), names(fasta))
  expect_equal(as.character(unlist(back)), as.character(unlist(fasta)))
  expect_equal(vapply(back, seqinr::getAnnot, character(1)),
               vapply(fasta, seqinr::getAnnot, character(1)))
})

test_that("writeFasta() writes the Annot header and 60 residues per line", {
  fasta <- list(seqinr::as.SeqFastaAA(strrep("A", 70), name = "x", Annot = ">x some protein"))
  file <- withr::local_tempfile(fileext = ".fasta")

  writeFasta(file, fasta)

  expect_equal(readLines(file), c(">x some protein", strrep("A", 60), strrep("A", 10)))
})

test_that("writeFasta() concatenates several FASTA lists", {
  a <- list(seqinr::as.SeqFastaAA("AAA", name = "a", Annot = ">a"))
  b <- list(seqinr::as.SeqFastaAA("CCC", name = "b", Annot = ">b"))
  file <- withr::local_tempfile(fileext = ".fasta")

  writeFasta(file, a, b)

  expect_equal(readLines(file), c(">a", "AAA", ">b", "CCC"))
})
