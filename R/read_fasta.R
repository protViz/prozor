#' Read a FASTA file of protein sequences
#'
#' Reads amino-acid sequences, each sequence as one string. Lines starting with `;` and lines before the first
#' header are skipped; an entry without sequence lines gets the empty string.
#'
#' @param file Path to a FASTA file, optionally gzip-compressed.
#' @return A named list with one `SeqFastaAA` object per sequence: the sequence as a string, with the attributes
#'   `name` (the first word of the header, up to the first whitespace, without `>`; also the list name) and `Annot`
#'   (the header line).
#' @family FASTA databases
#' @export
#' @examples
#' file <- system.file("extdata/fgcz_contaminants2021_20210929.fasta.gz", package = "prozor")
#' fasta <- read_fasta(file)
#' length(fasta)
#' attributes(fasta[[1]])
read_fasta <- function(file) {
  fasta_lines <- readLines(file, warn = FALSE)
  fasta_lines <- fasta_lines[!startsWith(fasta_lines, ";")]
  is_header <- startsWith(fasta_lines, ">")
  if (!any(is_header)) {
    stop("`file` must be a FASTA file with a line starting with `>`, but \"", file, "\" has none.")
  }
  entry <- cumsum(is_header)
  in_sequence <- !is_header & entry > 0
  groups <- factor(entry[in_sequence], levels = seq_len(sum(is_header)))
  sequences <- vapply(split(fasta_lines[in_sequence], groups), paste, character(1), collapse = "", USE.NAMES = FALSE)
  .as_seq_fasta_aa(sequences, fasta_lines[is_header])
}

# FASTA entries as read_fasta() returns them, with the class and attributes of seqinr's as.SeqFastaAA(): the
# `sequences` with their header lines `annots`, which start with ">". The `name` attribute and the list name are the
# first word of the header, up to the first whitespace, without ">". write_fasta() writes `Annot` with one leading
# ">", so it also accepts hand-built entries without it.
.as_seq_fasta_aa <- function(sequences, annots) {
  ids <- substring(sub("[[:space:]].*", "", annots), 2)
  res <- mapply(
    function(sequence, name, annot) structure(sequence, name = name, Annot = annot, class = "SeqFastaAA"),
    sequences,
    ids,
    annots,
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )
  names(res) <- ids
  res
}
