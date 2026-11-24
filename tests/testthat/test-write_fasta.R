test_that("write_fasta() output reads back unchanged", {
  fasta <- load_contaminants("fgcz2022")[2:6]
  file <- withr::local_tempfile(fileext = ".fasta")

  write_fasta(fasta, file)
  back <- read_peptide_fasta(file)

  expect_equal(names(back), names(fasta))
  expect_equal(as.character(unlist(back)), as.character(unlist(fasta)))
  expect_equal(vapply(back, seqinr::getAnnot, character(1)), vapply(fasta, seqinr::getAnnot, character(1)))
})

test_that("write_fasta() writes the Annot header and 60 residues per line", {
  fasta <- list(seqinr::as.SeqFastaAA(strrep("A", 70), name = "x", Annot = ">x some protein"))
  file <- withr::local_tempfile(fileext = ".fasta")

  write_fasta(fasta, file)

  expect_equal(readLines(file), c(">x some protein", strrep("A", 60), strrep("A", 10)))
})

test_that("write_fasta() keeps a '>' inside the header", {
  fasta <- list(seqinr::as.SeqFastaAA("AAA", name = "x", Annot = ">x ADP>ATP converter"))
  file <- withr::local_tempfile(fileext = ".fasta")

  write_fasta(fasta, file)

  expect_equal(readLines(file)[1], ">x ADP>ATP converter")
})
