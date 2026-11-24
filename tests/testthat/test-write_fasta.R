test_that("write_fasta() output reads back unchanged", {
  fasta <- load_contaminants("fgcz2022")[2:6]
  file <- withr::local_tempfile(fileext = ".fasta")

  write_fasta(fasta, file)
  back <- read_peptide_fasta(file)

  expect_equal(names(back), names(fasta))
  expect_equal(as.character(unlist(back)), as.character(unlist(fasta)))
  expect_equal(vapply(back, attr, character(1), "Annot"), vapply(fasta, attr, character(1), "Annot"))
})

test_that("write_fasta() writes the Annot header and 60 residues per line", {
  fasta <- list(seq_fasta_aa(strrep("A", 70), "x", ">x some protein"))
  file <- withr::local_tempfile(fileext = ".fasta")

  write_fasta(fasta, file)

  expect_equal(readLines(file), c(">x some protein", strrep("A", 60), strrep("A", 10)))
})

test_that("write_fasta() keeps a '>' inside the header", {
  fasta <- list(seq_fasta_aa("AAA", "x", ">x ADP>ATP converter"))
  file <- withr::local_tempfile(fileext = ".fasta")

  write_fasta(fasta, file)

  expect_equal(readLines(file)[1], ">x ADP>ATP converter")
})

test_that("write_fasta() writes sequences of 0, 60 and 121 residues", {
  fasta <- list(seq_fasta_aa("", "e"), seq_fasta_aa(strrep("A", 60), "s60"), seq_fasta_aa(strrep("C", 121), "s121"))
  file <- withr::local_tempfile(fileext = ".fasta")

  write_fasta(fasta, file)

  expect_equal(readLines(file), c(">e", ">s60", strrep("A", 60), ">s121", strrep("C", 60), strrep("C", 60), "C"))
  expect_equal(unlist(read_peptide_fasta(file), use.names = FALSE), unlist(fasta, use.names = FALSE))
})

test_that("write_fasta() output reads back unchanged across write chunks", {
  # write_fasta() writes 20,000 entries at a time
  n <- 20001
  fasta <- lapply(seq_len(n), function(i) seq_fasta_aa(strrep("K", i %% 130), paste0("p", i)))
  file <- withr::local_tempfile(fileext = ".fasta")

  write_fasta(fasta, file)

  expect_equal(unname(read_peptide_fasta(file)), fasta)
})
