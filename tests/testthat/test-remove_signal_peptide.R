test_that("remove_signal_peptide() cuts the signal peptide and renames the chain", {
  db <- contaminant_entries(2:3)
  signal <- data.frame(Entry = accessions(db), Signal.peptide = c("SIGNAL 1 20", ""))

  res <- remove_signal_peptide(db, signal)

  expect_equal(nchar(res[[1]]), nchar(db[[1]]) - 20)
  expect_equal(as.character(res[[1]]), substring(as.character(db[[1]]), 21))
  expect_equal(attr(res[[1]], "name"), paste0(attr(db[[1]], "name"), "s21"))
  expect_equal(res[[2]], db[[2]])
})

test_that("remove_signal_peptide() reads the current UniProt range format", {
  db <- contaminant_entries(2)
  signal <- data.frame(Entry = accessions(db), Signal.peptide = "SIGNAL 1..20; /evidence=\"ECO:0000255\"")

  expect_no_warning(res <- remove_signal_peptide(db, signal))

  expect_equal(nchar(res[[1]]), nchar(db[[1]]) - 20)
})

test_that("remove_signal_peptide() works on a one-entry database", {
  db <- contaminant_entries(2)
  signal <- data.frame(Entry = accessions(db), Signal.peptide = "SIGNAL 1 20")

  res <- remove_signal_peptide(db, signal)

  expect_equal(nchar(res[[1]]), nchar(db[[1]]) - 20)
})

test_that("remove_signal_peptide() names the first entry that does not match", {
  db <- contaminant_entries(2:3)
  signal <- data.frame(Entry = c(accessions(db)[1], "P99999"), Signal.peptide = c("", ""))

  expect_error(remove_signal_peptide(db, signal), "first mismatch in row 2: \"P99999\"")
  expect_error(remove_signal_peptide(db, signal[1, ]), "one row per entry of `db` \\(2\\), not 1")
})

test_that("make_chain_id() and make_chain_id_uniprot() build chain IDs", {
  expect_equal(make_chain_id("ABCDE", "sp|P30443|1A01_HUMAN", 24), "sp|P30443|1A01_HUMANs24")
  expect_equal(make_chain_id_uniprot("ABCDE", "sp|P30443|1A01_HUMAN", 24), "sp|P30443|24-5")
})
