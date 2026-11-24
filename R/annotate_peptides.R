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
#' @param pepinfo A character vector of peptide sequences, or a data frame (or tibble) with the column `peptide`
#'   and any further columns, such as the charge state, which are kept. It must not have the columns `protein_id`,
#'   `offset`, `protein_sequence`, `matched` or `peptide_length`, which the result adds.
#' @param fasta A list of protein sequences with unique names, as returned by [read_peptide_fasta()].
#' @param peptide Name of the column with the peptide sequences.
#' @param prefix Regular expression the residues before a match must end with.
#' @param suffix Regular expression the residues after a match must start with.
#' @return A data frame with one row per row of `pepinfo` and protein occurrence: the columns of `pepinfo`,
#'   `protein_id`, `offset` (1-based start of the peptide), `protein_sequence` (as searched, see [annotate_aho()]),
#'   `matched` and `peptide_length`. `NULL`, with a warning, if no peptide occurs in `fasta`.
#' @family protein inference
#' @export
#' @examples
#' file <- system.file("extdata/IDResults.txt.gz", package = "prozor")
#' peptides <- unique(utils::read.delim(file)$peptideSeq)
#' fasta_file <- system.file("p1000_db1_example/Annotation_canSeq.fasta.gz", package = "prozor")
#' fasta <- read_peptide_fasta(fasta_file)
#' annot <- annotate_peptides(peptides[1:20], fasta)
#' head(annot[, c("protein_id", "peptide_seq", "offset", "peptide_length", "matched")])
annotate_peptides <- function(pepinfo, fasta, peptide = "peptide_seq", prefix = "(([RK])|(^)|(^M))", suffix = "") {
  if (is.null(dim(pepinfo))) {
    pepinfo <- data.frame(pepinfo)
    colnames(pepinfo) <- peptide
  }
  clashes <- intersect(colnames(pepinfo), c("protein_id", "offset", "protein_sequence", "matched", "peptide_length"))
  if (length(clashes) > 0) {
    stop("`pepinfo` has columns named like result columns; rename them: ", paste(clashes, collapse = ", "), ".")
  }
  pepinfo[[peptide]] <- trimws(as.character(pepinfo[[peptide]]))
  pepinfo$peptide_length <- nchar(pepinfo[[peptide]])

  restab <- annotate_aho(unique(pepinfo[[peptide]]), fasta, peptide = peptide)
  if (is.null(restab)) {
    warning("None of the peptides in `pepinfo` occurs in `fasta`.")
    return(NULL)
  }
  before <- substr(restab$protein_sequence, 1, restab$offset - 1)
  after <- substring(restab$protein_sequence, restab$offset + nchar(restab[[peptide]]))
  restab$matched <- grepl(paste0("(", prefix, ")$"), before) & grepl(paste0("^(", suffix, ")"), after)
  merge(restab, pepinfo, by.x = peptide, by.y = peptide)
}


#' Find peptides in protein sequences with the Aho-Corasick algorithm
#'
#' Searches all peptides at once with [AhoCorasickTrie::AhoCorasickSearch()]. Peptides that occur in no protein
#' are dropped.
#'
#' The protein sequences are searched in upper case, with each character other than a letter, such as the stop
#' codon `*` of a translated database, replaced by `X`. The offsets are those of the sequences in `fasta`.
#'
#' @param pepseq Peptide sequences: upper-case letters only, without whitespace or modifications.
#'   [annotate_peptides()] removes surrounding whitespace before it calls this function.
#' @param fasta A list of protein sequences with unique names, as returned by [read_peptide_fasta()].
#' @param peptide Name of the result column that holds the peptide sequences.
#' @return A data frame with one row per occurrence and the columns `protein_id`, `peptide` (named by the
#'   argument), `offset` (1-based start of the peptide) and `protein_sequence` (the sequence as searched);
#'   `NULL` if no peptide occurs.
#' @family protein inference
#' @export
#' @examples
#' file <- system.file("extdata/IDResults.txt.gz", package = "prozor")
#' peptides <- unique(utils::read.delim(file)$peptideSeq)
#' fasta_file <- system.file("p1000_db1_example/Annotation_canSeq.fasta.gz", package = "prozor")
#' fasta <- read_peptide_fasta(fasta_file)
#' hits <- annotate_aho(peptides[1:20], fasta)
#' dim(hits)
annotate_aho <- function(pepseq, fasta, peptide = "peptide_seq") {
  pepseq <- unique(pepseq)
  invalid <- pepseq[!grepl("^[A-Z]+$", pepseq)]
  if (length(invalid) > 0) {
    # AhoCorasickTrie 0.1.3 crashes R on keywords outside A-Z (https://github.com/chambm/AhoCorasickTrie/pull/5)
    stop(
      "Peptide sequences must consist of upper-case letters only, without whitespace or modifications; found: ",
      paste(encodeString(utils::head(invalid, 5), quote = "\""), collapse = ", "),
      "."
    )
  }
  protein_ids <- names(fasta)
  duplicated_ids <- unique(protein_ids[duplicated(protein_ids)])
  if (length(duplicated_ids) > 0) {
    stop(
      "Protein names in `fasta` must be unique; duplicated: ",
      paste(utils::head(duplicated_ids, 5), collapse = ", "),
      "."
    )
  }
  # AhoCorasickTrie stops on characters outside A-Z, such as the stop codon `*`; one X per character keeps the offsets
  sequences <- gsub("[^A-Z]", "X", toupper(trimws(unlist(fasta, use.names = FALSE))), perl = TRUE)
  names(sequences) <- protein_ids

  res <- AhoCorasickTrie::AhoCorasickSearch(pepseq, sequences, alphabet = "aminoacid")
  if (length(res) == 0) {
    return(NULL)
  }
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
