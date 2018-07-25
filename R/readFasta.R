#' wrapper setting the correct parameters
#'
#' peptides which do not have protein assignment drop out
#' @param file - fasta file
#' @export
#' @examples
#' library(seqinr)
#' library(prozor)
#' file = system.file("extdata/fgcz_contaminants_20150123.fasta.gz",package = "prozor")
#' fasta = readPeptideFasta(file)
#'
readPeptideFasta <- function(file){
  read.fasta(file = file, as.string = TRUE, seqtype="AA")
}
