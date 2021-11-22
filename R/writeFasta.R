.writeFasta <- function(fasta , file=NULL  ) {
  namesres <- vapply(fasta, function(x) {attributes(x)$Annot}, character(1))
  namesres <- gsub(">", "", namesres)
  write.fasta(fasta, namesres , file.out = file, nbchar = 60, as.string = TRUE)
}
#' write fasta lists into file
#'
#' peptides which do not have protein assignment drop out
#' @param ... fasta list or single file
#' @param file where to write
#' @export
#' @return writes a file.
#' @examples
#' #example how to create a protein db with decoy sequences
#' library(seqinr)
#' #library(prozor)
#' file = system.file("extdata/fgcz_contaminants2021_20210929.fasta.gz",package = "prozor")
#' fasta = readPeptideFasta(file = file)
#' revfasta <- reverseSeq(fasta)
#' decoyDB <- c(fasta,revfasta)
#' stopifnot(length(decoyDB) == 2 * length(fasta))
#' \donttest{
#' writeFasta(decoyDB, file="test.fasta")
#' }
writeFasta <- function(file, ...  ) {
   fasta <- c(...)
  .writeFasta(fasta,file = file)
}
