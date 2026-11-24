# Toy proteins for peptide annotation; PEPTIDEK occurs
# - after K (P1, offset 6), after S (P2, offset 4), at the N-terminus (P3, offset 1),
# - after the initial M (P4, offset 2), and twice in P5 (offset 1, and offset 10 after S)
toy_fasta <- function() {
  list(
    P1 = seqinr::as.SeqFastaAA("MAAAKPEPTIDEKGGGR", name = "P1", Annot = ">P1"),
    P2 = seqinr::as.SeqFastaAA("MSSPEPTIDEKLLL", name = "P2", Annot = ">P2"),
    P3 = seqinr::as.SeqFastaAA("PEPTIDEKAAA", name = "P3", Annot = ">P3"),
    P4 = seqinr::as.SeqFastaAA("MPEPTIDEKAAA", name = "P4", Annot = ">P4"),
    P5 = seqinr::as.SeqFastaAA("PEPTIDEKSPEPTIDEK", name = "P5", Annot = ">P5")
  )
}

bundled_fasta_files <- function() {
  c(
    list.files(system.file("extdata", package = "prozor"), pattern = "\\.fasta\\.gz$", full.names = TRUE),
    list.files(system.file("p1000_db1_example", package = "prozor"), pattern = "\\.fasta\\.gz$", full.names = TRUE)
  )
}

# Toy target database (x, y) written to a temporary FASTA file for create_decoy_db()
toy_db_file <- function(env = parent.frame()) {
  file <- withr::local_tempfile(fileext = ".fasta", .local_envir = env)
  write_fasta(
    list(
      seqinr::as.SeqFastaAA("ABCD", name = "x", Annot = ">x"),
      seqinr::as.SeqFastaAA("EFG", name = "y", Annot = ">y")
    ),
    file
  )
  file
}

toy_contaminants <- function() {
  list(c1 = seqinr::as.SeqFastaAA("HIK", name = "c1", Annot = ">c1"))
}

# Bundled contaminant entries and their UniProt accessions for remove_signal_peptide()
contaminant_entries <- function(rows) {
  file <- system.file("extdata/fgcz_contaminants2021_20210929.fasta.gz", package = "prozor")
  read_peptide_fasta(file)[rows]
}

accessions <- function(db) {
  vapply(strsplit(names(db), split = "|", fixed = TRUE), function(x) x[2], character(1))
}
