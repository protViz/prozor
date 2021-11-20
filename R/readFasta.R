#' wrapper setting the correct parameters
#' seqinr::read.fasta for reading peptide sequences
#'
#' peptides which do not have protein assignment drop out
#' @param file - fasta file
#' @export
#' @return list with sequences
#' @examples
#' library(seqinr)
#'
#' file = system.file("extdata/fgcz_contaminants_20150123.fasta.gz",package = "prozor")
#' fasta = readPeptideFasta(file)
#'
readPeptideFasta <- function(file){
    seqinr::read.fasta(file = file, as.string = TRUE, seqtype = "AA")
}
