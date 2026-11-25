#' Write protein sequences to a FASTA file
#'
#' Writes each sequence with its `Annot` header and 60 residues per line. A header without the leading `>` gets one.
#'
#' @param fasta A list of `SeqFastaAA` objects, e.g. from [read_fasta()] or [create_decoy_db()]; each needs the
#'   attribute `Annot`, a single string.
#' @param file Path of the FASTA file to write.
#' @return `file`, invisibly.
#' @family FASTA databases
#' @export
#' @examples
#' file <- system.file("extdata/fgcz_contaminants2021_20210929.fasta.gz", package = "prozor")
#' fasta <- read_fasta(file)
#' decoy_db <- c(fasta, reverse_seq(fasta))
#' out <- tempfile(fileext = ".fasta")
#' write_fasta(decoy_db, out)
#' length(read_fasta(out))
#' unlink(out)
write_fasta <- function(fasta, file) {
  annots <- lapply(fasta, attr, which = "Annot")
  no_annot <- which(!vapply(annots, function(annot) is.character(annot) && length(annot) == 1, logical(1)))
  if (length(no_annot) > 0) {
    stop(
      "Each entry of `fasta` must have the attribute `Annot`, a single string, as from read_fasta(); entry ",
      no_annot[1],
      " does not."
    )
  }
  con <- file(file, "w")
  on.exit(close(con), add = TRUE)
  # one writeLines() per chunk of 20,000 entries; only the lines of one chunk are held in memory
  for (chunk in split(seq_along(fasta), ceiling(seq_along(fasta) / 20000))) {
    headers <- paste0(">", sub("^>", "", unlist(annots[chunk], use.names = FALSE)))
    sequences <- unlist(fasta[chunk], use.names = FALSE)
    n_lines <- ceiling(nchar(sequences) / 60)
    starts <- sequence(n_lines, from = 1L, by = 60L)
    fasta_lines <- character(length(chunk) + sum(n_lines))
    header_pos <- cumsum(c(1, n_lines[-length(n_lines)] + 1))
    fasta_lines[header_pos] <- headers
    fasta_lines[-header_pos] <- substring(rep.int(sequences, n_lines), starts, starts + 59L)
    writeLines(fasta_lines, con)
  }
  invisible(file)
}
