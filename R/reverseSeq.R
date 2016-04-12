.strReverse <- function(x) {
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse = "")
}
.reverseSingleSeq <-function(fasta, revLab="REV_"){
  name <- attributes(fasta)$name
  Annot <- attributes(fasta)$Annot
  revseq <- .strReverse(fasta)
  return(as.SeqFastaAA(revseq, Annot=paste(revLab, Annot, sep=""), name = paste(revLab, name, sep="")))

}
#' create rev sequences to fasta list
#'
#' peptides which do not have protein assignment drop out
#' @param fasta - an r list with SeqFastaAA
#' @param revLab - how to label reverse sequences, default = REV_
#' @export
#' @examples
#' library(seqinr)
#' library(prozor)
#'
#' file = file.path(path.package("prozor"),"extdata/fgcz_contaminants_20150123.fasta")
#' file
#' fasta = readPeptideFasta(file = file)
#' x <- reverseSeq(fasta)
#'
#'
#' revseq <- reverseSeq(fasta ,revLab = "REV_")
#' stopifnot(length(revseq) == length(fasta))
#' stopifnot(grep("^REV_","REV_zz|ZZ_FGCZCont0000|")==1)
#'
#' tmp <- list(as.SeqFastaAA(("DYKDDDDK"),name="Flag|FLAG|p2079",Annot=""))
#'
#' reverseSeq(tmp)
#'
reverseSeq<- function(fasta, revLab="REV_"){
  res <- lapply(fasta, .reverseSingleSeq ,revLab=revLab )
  revnames <- sapply(res ,  function(x){attributes(x)$name})
  names(res) <- revnames
  return(res)
}

