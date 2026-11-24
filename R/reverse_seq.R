.str_reverse <- function(x) {
  vapply(strsplit(x, NULL), function(chars) paste(rev(chars), collapse = ""), character(1))
}

.reverse_single_seq <- function(fasta, rev_lab) {
  seqinr::as.SeqFastaAA(
    .str_reverse(fasta),
    Annot = paste0(">", rev_lab, sub("^>", "", seqinr::getAnnot(fasta))),
    name = paste0(rev_lab, seqinr::getName(fasta))
  )
}

#' Reverse protein sequences to create decoys
#'
#' @param fasta A list of `SeqFastaAA` objects.
#' @param rev_lab Prefix added to the names and headers of the reversed sequences.
#' @return A list of reversed `SeqFastaAA` objects, named with the prefixed protein names.
#' @family FASTA databases
#' @export
#' @examples
#' file <- system.file("extdata/fgcz_contaminants2021_20210929.fasta.gz", package = "prozor")
#' fasta <- read_peptide_fasta(file)
#' decoys <- reverse_seq(fasta)
#' names(decoys)[1:3]
#'
#' tag <- list(seqinr::as.SeqFastaAA("DYKDDDDK", name = "Flag|FLAG|p2079", Annot = ">Flag|FLAG|p2079"))
#' reverse_seq(tag, rev_lab = "DECOY_")
reverse_seq <- function(fasta, rev_lab = "REV_") {
  res <- lapply(fasta, .reverse_single_seq, rev_lab = rev_lab)
  names(res) <- vapply(res, function(x) attributes(x)$name, character(1))
  res
}
