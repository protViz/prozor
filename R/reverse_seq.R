#' Reverse protein sequences to create decoys
#'
#' @param fasta A list of `SeqFastaAA` objects.
#' @param rev_lab Prefix added to the names and headers of the reversed sequences.
#' @return A list of reversed `SeqFastaAA` objects, named with the prefixed protein names.
#' @family FASTA databases
#' @export
#' @examples
#' file <- system.file("extdata/fgcz_contaminants2021_20210929.fasta.gz", package = "prozor")
#' fasta <- read_fasta(file)
#' decoys <- reverse_seq(fasta)
#' names(decoys)[1:3]
#'
#' # an entry as read_fasta() returns it
#' tag <- structure(
#'   "DYKDDDDK",
#'   name = "Flag|FLAG|p2079",
#'   Annot = ">Flag|FLAG|p2079",
#'   class = "SeqFastaAA"
#' )
#' reverse_seq(list(tag), rev_lab = "DECOY_")
reverse_seq <- function(fasta, rev_lab = "REV_") {
  sequences <- vapply(
    unlist(fasta, use.names = FALSE),
    function(sequence) intToUtf8(rev(utf8ToInt(sequence))),
    character(1),
    USE.NAMES = FALSE
  )
  annots <- vapply(fasta, attr, character(1), "Annot", USE.NAMES = FALSE)
  .as_seq_fasta_aa(sequences, paste0(">", rev_lab, sub("^>", "", annots), recycle0 = TRUE))
}
