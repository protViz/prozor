.reverse_single_seq <- function(fasta, rev_lab) {
  .as_seq_fasta_aa(
    intToUtf8(rev(utf8ToInt(fasta))),
    name = paste0(rev_lab, attr(fasta, "name")),
    annot = paste0(">", rev_lab, sub("^>", "", attr(fasta, "Annot")))
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
#' # an entry as read_peptide_fasta() returns it
#' tag <- structure(
#'   "DYKDDDDK",
#'   name = "Flag|FLAG|p2079",
#'   Annot = ">Flag|FLAG|p2079",
#'   class = "SeqFastaAA"
#' )
#' reverse_seq(list(tag), rev_lab = "DECOY_")
reverse_seq <- function(fasta, rev_lab = "REV_") {
  res <- lapply(fasta, .reverse_single_seq, rev_lab = rev_lab)
  names(res) <- vapply(res, function(x) attributes(x)$name, character(1))
  res
}
