#' Make a UniProt-style ID for a protein chain
#'
#' @param sequence The chain sequence, as one string.
#' @param id The UniProt ID of the protein, e.g. `"sp|P30443|1A01_HUMAN"`.
#' @param start Start position of the chain.
#' @return The first two fields of `id`, followed by `start`, `-` and the length of `sequence`, e.g.
#'   `"sp|P30443|24-5"`.
#' @family FASTA databases
#' @export
#' @examples
#' make_chain_id_uniprot("ABCDE", "sp|P30443|1A01_HUMAN", 24)
make_chain_id_uniprot <- function(sequence, id, start) {
  fields <- strsplit(id, split = "|", fixed = TRUE)[[1]][1:2]
  paste0(paste0(fields, "|", collapse = ""), start, "-", nchar(sequence))
}

#' Make an ID for a protein chain
#'
#' @param sequence The chain sequence, as one string; not used.
#' @param id The ID of the protein, e.g. `"sp|P30443|1A01_HUMAN"`.
#' @inheritParams make_chain_id_uniprot
#' @return `id` followed by `s` and `start`, e.g. `"sp|P30443|1A01_HUMANs24"`.
#' @family FASTA databases
#' @export
#' @examples
#' make_chain_id("ABCDE", "sp|P30443|1A01_HUMAN", 24)
make_chain_id <- function(sequence, id, start) {
  paste0(id, "s", start)
}

.make_sig <- function(sequence, chain_start, id_fun) {
  new_id <- id_fun(sequence, attr(sequence, "name"), chain_start)
  # the header is ">" and the ID, then the description, which is kept
  .as_seq_fasta_aa(
    substring(sequence, first = chain_start),
    name = new_id,
    annot = paste0(">", new_id, sub("^>?\\S*", "", attr(sequence, "Annot")))
  )
}

#' Remove signal peptides from protein sequences
#'
#' @param fasta A list of UniProt protein sequences, as returned by [read_fasta()], named with UniProt IDs
#'   such as `"sp|P30443|1A01_HUMAN"` (the accession follows the first `|`).
#' @param signal A data frame with one row per entry of `fasta`, in the same order, and the columns `Entry` (the
#'   UniProt accession) and `Signal.peptide` (UniProt's signal peptide feature, e.g. `"SIGNAL 1 20"` or
#'   `"SIGNAL 1..20; /evidence=..."`; empty if there is none).
#' @param id_fun Function that makes the ID of the remaining chain, e.g. [make_chain_id()] or
#'   [make_chain_id_uniprot()].
#' @return `fasta` with the signal peptides removed. A shortened sequence gets the ID made by `id_fun` as its list
#'   name, its `name` and the first word of its header; the rest of the header is kept.
#' @family FASTA databases
#' @export
#' @examples
#' file <- system.file("extdata/fgcz_contaminants2021_20210929.fasta.gz", package = "prozor")
#' fasta <- read_fasta(file)[2:3]
#' entry <- vapply(strsplit(names(fasta), split = "|", fixed = TRUE), function(x) x[2], character(1))
#' signal <- data.frame(Entry = entry, Signal.peptide = c("SIGNAL 1 20", ""))
#' res <- remove_signal_peptide(fasta, signal)
#' nchar(fasta[[1]])
#' nchar(res[[1]])
remove_signal_peptide <- function(fasta, signal, id_fun = make_chain_id) {
  entries <- vapply(strsplit(names(fasta), split = "|", fixed = TRUE), function(x) x[2], character(1))
  no_accession <- which(is.na(entries))
  if (length(no_accession) > 0) {
    stop(
      "The names of `fasta` must be UniProt IDs such as \"sp|P30443|1A01_HUMAN\", with the accession after the first ",
      "`|`; \"",
      names(fasta)[no_accession[1]],
      "\" has none."
    )
  }
  if (length(entries) != nrow(signal)) {
    stop("`signal` must have one row per entry of `fasta` (", length(entries), "), not ", nrow(signal), ".")
  }
  mismatch <- which(is.na(signal$Entry) | entries != signal$Entry)
  if (length(mismatch) > 0) {
    stop(
      "`signal$Entry` must list the accessions of `fasta` in the same order; first mismatch in row ",
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

  res <- fasta
  for (i in seq_along(signal_end)) {
    if (!is.na(signal_end[i])) {
      res[[i]] <- .make_sig(fasta[[i]], signal_end[i] + 1, id_fun = id_fun)
      names(res)[i] <- attr(res[[i]], "name")
    }
  }
  res
}
