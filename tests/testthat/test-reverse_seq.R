test_that("reverse_seq() reverses sequences and labels names and headers", {
  fasta <- list(seqinr::as.SeqFastaAA("DYKDDDDK", name = "Flag|FLAG|p2079", Annot = ">Flag|FLAG|p2079 flag tag"))

  res <- reverse_seq(fasta)

  expect_named(res, "REV_Flag|FLAG|p2079")
  expect_equal(as.character(res[[1]]), "KDDDDKYD")
  expect_equal(seqinr::getAnnot(res[[1]]), ">REV_Flag|FLAG|p2079 flag tag")
})

test_that("reverse_seq() uses the given label", {
  fasta <- list(seqinr::as.SeqFastaAA("ABC", name = "x", Annot = ">x"))

  res <- reverse_seq(fasta, rev_lab = "DECOY_")

  expect_named(res, "DECOY_x")
  expect_equal(as.character(res[[1]]), "CBA")
})
