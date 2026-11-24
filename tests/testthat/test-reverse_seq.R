test_that("reverse_seq() reverses sequences and labels names and headers", {
  fasta <- list(seq_fasta_aa("DYKDDDDK", "Flag|FLAG|p2079", ">Flag|FLAG|p2079 flag tag"))

  res <- reverse_seq(fasta)

  expect_named(res, "REV_Flag|FLAG|p2079")
  expect_equal(as.character(res[[1]]), "KDDDDKYD")
  expect_equal(attr(res[[1]], "Annot"), ">REV_Flag|FLAG|p2079 flag tag")
})

test_that("reverse_seq() uses the given label", {
  fasta <- list(seq_fasta_aa("ABC", "x"))

  res <- reverse_seq(fasta, rev_lab = "DECOY_")

  expect_named(res, "DECOY_x")
  expect_equal(as.character(res[[1]]), "CBA")
})
