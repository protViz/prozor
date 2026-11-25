test_that("remove_signal_peptide() cuts the signal peptide and renames the chain", {
  fasta <- contaminant_entries(2:3)
  signal <- data.frame(Entry = accessions(fasta), Signal.peptide = c("SIGNAL 1 20", ""))

  res <- remove_signal_peptide(fasta, signal)

  expect_equal(nchar(res[[1]]), nchar(fasta[[1]]) - 20)
  expect_equal(as.character(res[[1]]), substring(as.character(fasta[[1]]), 21))
  expect_equal(attr(res[[1]], "name"), paste0(attr(fasta[[1]], "name"), "s21"))
  expect_equal(res[[2]], fasta[[2]])
})

test_that("remove_signal_peptide() renames the list entry and keeps the header description", {
  fasta <- contaminant_entries(2:3)
  signal <- data.frame(Entry = accessions(fasta), Signal.peptide = c("SIGNAL 1 20", ""))
  file <- withr::local_tempfile(fileext = ".fasta")

  res <- remove_signal_peptide(fasta, signal)
  write_fasta(res, file)
  back <- read_fasta(file)

  expect_equal(names(res), c("zz|Y-FGCZCont00001|s21", "zz|Y-FGCZCont00002|"))
  expect_equal(
    attr(res[[1]], "Annot"),
    ">zz|Y-FGCZCont00001|s21  zz|FGCZCont0000_P61626|LYSC_HUMAN blastpHomologue_5.0e-107"
  )
  expect_equal(names(back), names(res))
  expect_equal(attr(back[[1]], "Annot"), attr(res[[1]], "Annot"))
})

test_that("remove_signal_peptide() names the chain with `id_fun`", {
  fasta <- contaminant_entries(2:3)
  signal <- data.frame(Entry = accessions(fasta), Signal.peptide = c("SIGNAL 1 20", ""))

  res <- remove_signal_peptide(fasta, signal, id_fun = make_chain_id_uniprot)

  # the protein is 148 residues long; the chain starts at 21
  expect_equal(names(res), c("zz|Y-FGCZCont00001|21-148", "zz|Y-FGCZCont00002|"))
  expect_equal(attr(res[[1]], "name"), "zz|Y-FGCZCont00001|21-148")
})

test_that("remove_signal_peptide() reads the current UniProt range format", {
  fasta <- contaminant_entries(2)
  signal <- data.frame(Entry = accessions(fasta), Signal.peptide = "SIGNAL 1..20; /evidence=\"ECO:0000255\"")

  expect_no_warning(res <- remove_signal_peptide(fasta, signal))

  expect_equal(nchar(res[[1]]), nchar(fasta[[1]]) - 20)
})

test_that("remove_signal_peptide() works on a one-entry database", {
  fasta <- contaminant_entries(2)
  signal <- data.frame(Entry = accessions(fasta), Signal.peptide = "SIGNAL 1 20")

  res <- remove_signal_peptide(fasta, signal)

  expect_equal(nchar(res[[1]]), nchar(fasta[[1]]) - 20)
})

test_that("remove_signal_peptide() names the first entry that does not match", {
  fasta <- contaminant_entries(2:3)
  signal <- data.frame(Entry = c(accessions(fasta)[1], "P99999"), Signal.peptide = c("", ""))

  expect_error(remove_signal_peptide(fasta, signal), "first mismatch in row 2: \"P99999\"")
  expect_error(remove_signal_peptide(fasta, signal[1, ]), "one row per entry of `fasta` \\(2\\), not 1")
})

test_that("make_chain_id() and make_chain_id_uniprot() build chain IDs", {
  expect_equal(make_chain_id("ABCDE", "sp|P30443|1A01_HUMAN", start = 24), "sp|P30443|1A01_HUMANs24")
  # the number after "-" is the end of the chain, the length of the protein
  expect_equal(make_chain_id_uniprot(strrep("A", 148), "sp|P30443|1A01_HUMAN", start = 21), "sp|P30443|21-148")
})

test_that("remove_signal_peptide() stops for names without an accession and for NA entries", {
  fasta <- contaminant_entries(2:3)
  signal <- data.frame(Entry = accessions(fasta), Signal.peptide = c("", ""))
  no_accession <- fasta
  names(no_accession)[2] <- "P61626"
  na_entry <- signal
  na_entry$Entry[2] <- NA

  expect_error(remove_signal_peptide(no_accession, signal), "accession after the first `|`; \"P61626\" has none")
  expect_error(remove_signal_peptide(fasta, na_entry), "first mismatch in row 2: \"NA\"")
})

test_that("remove_signal_peptide() stops when `signal` lacks the column `Entry` or `Signal.peptide`", {
  fasta <- contaminant_entries(2:3)
  # read.delim(check.names = FALSE) and readr keep UniProt's column name "Signal peptide"
  unchecked_names <- data.frame(Entry = accessions(fasta), `Signal peptide` = c("SIGNAL 1 20", ""), check.names = FALSE)
  no_entry <- data.frame(Signal.peptide = c("SIGNAL 1 20", ""))

  expect_error(
    remove_signal_peptide(fasta, unchecked_names),
    "`signal` must have the columns `Entry`, `Signal.peptide`; missing: `Signal.peptide`\\."
  )
  expect_error(remove_signal_peptide(fasta, no_entry), "missing: `Entry`\\.")
})
