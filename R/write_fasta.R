#' Write protein sequences to a FASTA file
#'
#' Writes each sequence with its `Annot` header and 60 residues per line.
#'
#' @param fasta A list of `SeqFastaAA` objects, e.g. from [read_peptide_fasta()] or [create_decoy_db()].
#' @param file Path of the FASTA file to write.
#' @return No return value, called for side effects (writes a FASTA file).
#' @family FASTA databases
#' @export
#' @examples
#' file <- system.file("extdata/fgcz_contaminants2021_20210929.fasta.gz", package = "prozor")
#' fasta <- read_peptide_fasta(file)
#' decoy_db <- c(fasta, reverse_seq(fasta))
#' out <- tempfile(fileext = ".fasta")
#' write_fasta(decoy_db, out)
#' length(read_peptide_fasta(out))
#' unlink(out)
write_fasta <- function(fasta, file) {
  con <- file(file, "w")
  on.exit(close(con))
  # one writeLines() per chunk of 20,000 entries; only the lines of one chunk are held in memory
  for (chunk in split(seq_along(fasta), ceiling(seq_along(fasta) / 20000))) {
    entries <- fasta[chunk]
    headers <- paste0(">", sub("^>", "", vapply(entries, attr, character(1), "Annot", USE.NAMES = FALSE)))
    sequences <- unlist(entries, use.names = FALSE)
    n_lines <- ceiling(nchar(sequences) / 60)
    starts <- sequence(n_lines, from = 1L, by = 60L)
    lines <- character(length(entries) + sum(n_lines))
    header_pos <- cumsum(c(1, n_lines[-length(n_lines)] + 1))
    lines[header_pos] <- headers
    lines[-header_pos] <- substring(rep.int(sequences, n_lines), starts, starts + 59L)
    writeLines(lines, con)
  }
}
