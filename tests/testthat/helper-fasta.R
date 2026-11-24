# A FASTA entry as read_peptide_fasta() returns it
seq_fasta_aa <- function(sequence, name, annot = paste0(">", name)) {
  structure(sequence, name = name, Annot = annot, class = "SeqFastaAA")
}

# Toy proteins for peptide annotation; PEPTIDEK occurs
# - after K (P1, offset 6), after S (P2, offset 4), at the N-terminus (P3, offset 1),
# - after the initial M (P4, offset 2), and twice in P5 (offset 1, and offset 10 after S)
toy_fasta <- function() {
  list(
    P1 = seq_fasta_aa("MAAAKPEPTIDEKGGGR", "P1"),
    P2 = seq_fasta_aa("MSSPEPTIDEKLLL", "P2"),
    P3 = seq_fasta_aa("PEPTIDEKAAA", "P3"),
    P4 = seq_fasta_aa("MPEPTIDEKAAA", "P4"),
    P5 = seq_fasta_aa("PEPTIDEKSPEPTIDEK", "P5")
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
      seq_fasta_aa("ABCD", "x"),
      seq_fasta_aa("EFG", "y")
    ),
    file
  )
  file
}

# Toy database directory "toydb" for create_fgcz_fasta_db(): annotation.txt and the toy targets (x, y) in `fasta_file`
toy_database_dir <- function(fasta_file = "targets.fasta", env = parent.frame()) {
  database_dir <- file.path(withr::local_tempdir(.local_envir = env), "toydb")
  dir.create(database_dir)
  writeLines("zz|p1000_toydb|20261124 toy database", file.path(database_dir, "annotation.txt"))
  file.copy(toy_db_file(env), file.path(database_dir, fasta_file))
  database_dir
}

toy_contaminants <- function() {
  list(c1 = seq_fasta_aa("HIK", "c1"))
}

# Bundled contaminant entries and their UniProt accessions for remove_signal_peptide()
contaminant_entries <- function(rows) {
  file <- system.file("extdata/fgcz_contaminants2021_20210929.fasta.gz", package = "prozor")
  read_peptide_fasta(file)[rows]
}

accessions <- function(db) {
  vapply(strsplit(names(db), split = "|", fixed = TRUE), function(x) x[2], character(1))
}
