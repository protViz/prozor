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
