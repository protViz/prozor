.writeFasta <- function(fasta , file=NULL  ) {
  namesres <- sapply(fasta, function(x) attributes(x)$Annot)
  namesres <- gsub(">","",namesres)
  write.fasta(fasta, namesres , file.out = file, nbchar = 60, as.string =TRUE)
}
#' write fasta lists into file
#'
#' peptides which do not have protein assignment drop out
#' @param ... fasta list or single file
#' @param file where to write
#' @export
#' @examples
#' #example how to create a protein db with decoy sequences
#' library(seqinr)
#' library(prozor)
#' file = file.path(path.package("prozor"),"extdata/fgcz_contaminants_20150123.fasta")
#' fasta = readPeptideFasta(file = file)
#' revfasta <- reverseSeq(fasta)
#' decoyDB <- c(fasta,revfasta)
#' stopifnot(length(decoyDB) == 2 * length(fasta))
#' writeFasta(decoyDB, file="test.fasta")
#'
writeFasta <- function( file, ...  ) {
   fasta <- c(...)
  .writeFasta(fasta,file=file)
}
