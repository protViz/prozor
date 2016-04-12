#' wrapper setting the correct parameters
#'
#' peptides which do not have protein assignment drop out
#' @param file - fasta file
#' @export
#' @examples
#' library(seqinr)
#' library(prozor)
#' file = file.path(path.package("prozor"),"extdata/fgcz_contaminants_20150123.fasta")
#' fasta = readPeptideFasta(file)
#'
readPeptideFasta <- function(file){
  read.fasta(file = file, as.string = TRUE, seqtype="AA")
}
