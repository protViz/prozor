.strReverse <- function(x) {
  bb <- lapply(strsplit(x, NULL), rev)
  vapply(bb,  paste, character(1), collapse = "")
}

.reverseSingleSeq <- function(fasta, revLab = "REV_"){
  name <- getName(fasta)
  Annot <- getAnnot(fasta)
  revseq <- .strReverse(fasta)
  return(as.SeqFastaAA(revseq, Annot = paste(">", revLab, gsub("^>","",Annot), sep = ""),
                       name = paste(revLab, name, sep = "")))

}
#' create rev sequences to fasta list
#'
#' peptides which do not have protein assignment drop out
#' @param fasta - an r list with SeqFastaAA
#' @param revLab - how to label reverse sequences, default = REV_
#' @export
#' @return string with reversed sequence
#' @examples
#' library(seqinr)
#' #library(prozor)
#'
#' #file = file.path(path.package("prozor"),"extdata/fgcz_contaminants_20150123.fasta.gz")
#' file = system.file("extdata/fgcz_contaminants_20150123.fasta.gz",package = "prozor")
#' fasta = readPeptideFasta(file = file)
#' getAnnot(fasta[[1]])
#' x <- reverseSeq(fasta)
#'
#' revseq <- reverseSeq(fasta ,revLab = "REV_")
#' stopifnot(length(revseq) == length(fasta))
#' stopifnot(grep("^REV_","REV_zz|ZZ_FGCZCont0000|")==1)
#'
#' tmp <- list(as.SeqFastaAA(("DYKDDDDK"),name="Flag|FLAG|p2079",Annot=""))
#'
#' reverseSeq(tmp)
#'
reverseSeq <- function(fasta, revLab = "REV_"){
  res <- lapply(fasta, .reverseSingleSeq ,revLab = revLab )
  revnames <- vapply(res ,  function(x){attributes(x)$name}, character(1))
  names(res) <- revnames
  return(res)
}

