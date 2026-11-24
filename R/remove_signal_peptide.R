#' Make a UniProt-style ID for a protein chain
#'
#' @param sequence The chain sequence, as one string.
#' @param id The UniProt ID of the protein, e.g. `"sp|P30443|1A01_HUMAN"`.
#' @param sp Start position of the chain.
#' @return The first two fields of `id`, followed by `sp`, `-` and the length of `sequence`, e.g.
#'   `"sp|P30443|24-5"`.
#' @family FASTA databases
#' @export
#' @examples
#' make_chain_id_uniprot("ABCDE", "sp|P30443|1A01_HUMAN", 24)
make_chain_id_uniprot <- function(sequence, id, sp) {
  fields <- strsplit(id, split = "|", fixed = TRUE)[[1]][1:2]
  paste0(paste0(fields, "|", collapse = ""), sp, "-", nchar(sequence))
}

#' Make an ID for a protein chain
#'
#' @param sequence The chain sequence, as one string; not used.
#' @param id The ID of the protein, e.g. `"sp|P30443|1A01_HUMAN"`.
#' @param sp Start position of the chain.
#' @return `id` followed by `s` and `sp`, e.g. `"sp|P30443|1A01_HUMANs24"`.
#' @family FASTA databases
#' @export
#' @examples
#' make_chain_id("ABCDE", "sp|P30443|1A01_HUMAN", 24)
make_chain_id <- function(sequence, id, sp) {
  paste0(id, "s", sp)
}

.make_sig <- function(sequence, chain_start, id_fun) {
  attrs <- attributes(sequence)
  new_id <- id_fun(sequence, attrs$name, chain_start)
  attrs$Annot <- new_id
  attrs$name <- new_id
  chain <- substring(sequence, first = chain_start)
  attributes(chain) <- attrs
  chain
}

#' Remove signal peptides from protein sequences
#'
#' @param db A list of UniProt protein sequences, as returned by [read_peptide_fasta()].
#' @param signal A data frame with one row per entry of `db`, in the same order, and the columns `Entry` (the
#'   UniProt accession) and `Signal.peptide` (UniProt's signal peptide feature, e.g. `"SIGNAL 1 20"` or
#'   `"SIGNAL 1..20; /evidence=..."`; empty if there is none).
#' @param id_fun Function that makes the ID of the remaining chain, e.g. [make_chain_id()] or
#'   [make_chain_id_uniprot()].
#' @return `db` with the signal peptides removed and the IDs of the shortened sequences changed by `id_fun`.
#' @family FASTA databases
#' @export
#' @examples
#' file <- system.file("extdata/fgcz_contaminants2021_20210929.fasta.gz", package = "prozor")
#' db <- read_peptide_fasta(file)[2:3]
#' entry <- vapply(strsplit(names(db), split = "|", fixed = TRUE), function(x) x[2], character(1))
#' signal <- data.frame(Entry = entry, Signal.peptide = c("SIGNAL 1 20", ""))
#' res <- remove_signal_peptide(db, signal)
#' nchar(db[[1]])
#' nchar(res[[1]])
remove_signal_peptide <- function(db, signal, id_fun = make_chain_id) {
  entries <- vapply(strsplit(names(db), split = "|", fixed = TRUE), function(x) x[2], character(1))
  if (length(entries) != nrow(signal)) {
    stop("`signal` must have one row per entry of `db` (", length(entries), "), not ", nrow(signal), ".")
  }
  mismatch <- which(entries != signal$Entry)
  if (length(mismatch) > 0) {
    stop(
      "`signal$Entry` must list the accessions of `db` in the same order; first mismatch in row ",
      mismatch[1],
      ": \"",
      signal$Entry[mismatch[1]],
      "\", expected \"",
      entries[mismatch[1]],
      "\"."
    )
  }
  parts <- regmatches(signal$Signal.peptide, regexec("^SIGNAL [0-9]+(\\.\\.| )([0-9]+)", signal$Signal.peptide))
  signal_end <- vapply(parts, function(p) if (length(p) == 3) as.numeric(p[3]) else NA_real_, numeric(1))

  res <- db
  for (i in seq_along(signal_end)) {
    if (!is.na(signal_end[i])) {
      res[[i]] <- .make_sig(db[[i]], signal_end[i] + 1, id_fun = id_fun)
    }
  }
  res
}
