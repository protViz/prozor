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
  headers <- vapply(fasta, function(x) attributes(x)$Annot, character(1))
  seqinr::write.fasta(fasta, sub("^>", "", headers), file.out = file, nbchar = 60, as.string = TRUE)
}
