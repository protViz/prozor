#' Annotate peptides with the proteins that contain them
#'
#' Finds every occurrence of each peptide in `fasta` with [annotate_aho()] and flags whether it is preceded by a
#' tryptic cleavage site. Peptides that occur in no protein are dropped.
#'
#' Each occurrence of a peptide in a protein is one row. `matched` is `TRUE` when the residues before that
#' occurrence match `prefix` and the residues after it match `suffix`. With the defaults, `matched` means: preceded
#' by K or R or the protein N-terminus (or by the Met at position 1). The end of the peptide is not checked, neither
#' its last residue nor a following proline; pass a `suffix` to check the residues after it.
#'
#' Whitespace around the peptide sequences, such as the `"\r"` of a file with Windows line endings, is removed
#' before the search, and the result holds the trimmed sequences.
#'
#' @param peptides A character vector of peptide sequences, or a data frame (or tibble) with the column named by
#'   `peptide` and any further columns, such as the charge state, which are kept. It must not have the columns
#'   `protein_id`, `offset`, `protein_sequence`, `matched` or `peptide_length`, which the result adds.
#' @param fasta A list of protein sequences with unique names, as returned by [read_fasta()].
#' @param peptide Name of the column with the peptide sequences.
#' @param prefix Regular expression the residues before a match must end with.
#' @param suffix Regular expression the residues after a match must start with.
#' @return A data frame with one row per row of `peptides` and protein occurrence: the columns of `peptides`,
#'   `protein_id`, `offset` (1-based start of the peptide), `protein_sequence` (as searched, see [annotate_aho()]),
#'   `matched` and `peptide_length`. With 0 rows and a warning if no peptide occurs in `fasta`.
#' @family protein inference
#' @export
#' @examples
#' file <- system.file("extdata/IDResults.txt.gz", package = "prozor")
#' peptides <- unique(utils::read.delim(file)$peptideSeq)
#' fasta_file <- system.file("p1000_db1_example/Annotation_canSeq.fasta.gz", package = "prozor")
#' fasta <- read_fasta(fasta_file)
#' annot <- annotate_peptides(peptides[1:20], fasta)
#' head(annot[, c("protein_id", "peptide_seq", "offset", "peptide_length", "matched")])
annotate_peptides <- function(peptides, fasta, peptide = "peptide_seq", prefix = "(([RK])|(^)|(^M))", suffix = "") {
  if (is.null(dim(peptides))) {
    peptides <- data.frame(peptides)
    colnames(peptides) <- peptide
  }
  .check_columns(peptides, peptide, "peptides")
  clashes <- intersect(colnames(peptides), c("protein_id", "offset", "protein_sequence", "matched", "peptide_length"))
  if (length(clashes) > 0) {
    stop(
      "`peptides` must not have columns named like result columns; rename ",
      .quoted(clashes, "`"),
      "."
    )
  }
  peptides[[peptide]] <- trimws(as.character(peptides[[peptide]]))
  peptides$peptide_length <- nchar(peptides[[peptide]])

  restab <- annotate_aho(unique(peptides[[peptide]]), fasta, peptide = peptide)
  if (nrow(restab) == 0) {
    warning("None of the peptides in `peptides` occurs in `fasta`.")
  }
  before <- substr(restab$protein_sequence, 1, restab$offset - 1)
  after <- substring(restab$protein_sequence, restab$offset + nchar(restab[[peptide]]))
  restab$matched <- grepl(paste0("(", prefix, ")$"), before) & grepl(paste0("^(", suffix, ")"), after)
  merge(restab, peptides, by.x = peptide, by.y = peptide)
}


#' Find peptides in protein sequences with the Aho-Corasick algorithm
#'
#' Searches all peptides at once with [AhoCorasickTrie::AhoCorasickSearch()]. Peptides that occur in no protein
#' are dropped.
#'
#' The protein sequences are searched in upper case, with each character other than a letter, such as the stop
#' codon `*` of a translated database or a blank, replaced by `X`. The offsets are those of the sequences in `fasta`.
#'
#' @param peptides Peptide sequences: upper-case letters only, without whitespace or modifications.
#'   [annotate_peptides()] removes surrounding whitespace before it calls this function.
#' @inheritParams annotate_peptides
#' @param peptide Name of the result column that holds the peptide sequences.
#' @return A data frame with one row per occurrence and the columns `protein_id`, `peptide` (named by the
#'   argument), `offset` (1-based start of the peptide) and `protein_sequence` (the sequence as searched); with 0
#'   rows if no peptide occurs.
#' @family protein inference
#' @export
#' @examples
#' file <- system.file("extdata/IDResults.txt.gz", package = "prozor")
#' peptides <- unique(utils::read.delim(file)$peptideSeq)
#' fasta_file <- system.file("p1000_db1_example/Annotation_canSeq.fasta.gz", package = "prozor")
#' fasta <- read_fasta(fasta_file)
#' hits <- annotate_aho(peptides[1:20], fasta)
#' dim(hits)
annotate_aho <- function(peptides, fasta, peptide = "peptide_seq") {
  peptides <- unique(peptides)
  invalid <- peptides[!grepl("^[A-Z]+$", peptides)]
  if (length(invalid) > 0) {
    # rejects modifications, whitespace and lower case (AhoCorasickTrie 0.1.3 crashes R on them, 0.1.4 gives a less
    # specific error), and "" and NA (wrong hits in both versions)
    stop(
      "`peptides` must consist of upper-case letters only, without whitespace or modifications; found: ",
      .quoted(utils::head(invalid, 5), "\""),
      "."
    )
  }
  protein_ids <- names(fasta)
  if (is.null(protein_ids)) {
    stop("`fasta` must be a named list of protein sequences, as returned by read_fasta(); it has no names.")
  }
  duplicated_ids <- unique(protein_ids[duplicated(protein_ids)])
  if (length(duplicated_ids) > 0) {
    stop(
      "The names of `fasta` must be unique; duplicated: ",
      .quoted(utils::head(duplicated_ids, 5), "\""),
      "."
    )
  }
  # AhoCorasickTrie stops on characters outside A-Z, such as the stop codon `*`; one X per character keeps the offsets
  sequences <- gsub("[^A-Z]", "X", toupper(unlist(fasta, use.names = FALSE)), perl = TRUE)
  names(sequences) <- protein_ids

  res <- AhoCorasickTrie::AhoCorasickSearch(peptides, sequences, alphabet = "aminoacid")
  matches <- unlist(res, recursive = FALSE, use.names = FALSE)
  protein_id <- rep(names(res), lengths(res))
  hits <- data.frame(
    protein_id = protein_id,
    peptide = vapply(matches, `[[`, character(1), "Keyword"),
    offset = vapply(matches, `[[`, numeric(1), "Offset"),
    protein_sequence = unname(sequences[protein_id])
  )
  names(hits)[2] <- peptide
  hits
}
