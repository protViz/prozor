#' Read a FASTA file of protein sequences
#'
#' Reads amino-acid sequences with [seqinr::read.fasta()], each sequence as one string.
#'
#' @param file Path to a FASTA file, optionally gzip-compressed.
#' @return A named list of `SeqFastaAA` objects, one per sequence.
#' @family FASTA databases
#' @export
#' @examples
#' file <- system.file("extdata/fgcz_contaminants2021_20210929.fasta.gz", package = "prozor")
#' fasta <- read_peptide_fasta(file)
#' length(fasta)
read_peptide_fasta <- function(file) {
  seqinr::read.fasta(file = file, as.string = TRUE, seqtype = "AA")
}
