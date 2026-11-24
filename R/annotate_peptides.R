#' Annotate peptides with the proteins that contain them
#'
#' Finds every occurrence of each peptide in `fasta` with [annotate_aho()] and flags whether it is a tryptic
#' match. Peptides that occur in no protein are dropped.
#'
#' Each occurrence of a peptide in a protein is one row. `matched` is `TRUE` when the residues before that
#' occurrence match `prefix` and the residues after it match `suffix`; by default, when the occurrence
#' follows K or R, starts the protein, or follows its initial methionine.
#'
#' @param pepinfo A character vector of peptide sequences, or a data frame (or tibble) with the column `peptide`
#'   and any further columns, such as the charge state, which are kept.
#' @param fasta A list of protein sequences with unique names, as returned by [read_peptide_fasta()].
#' @param peptide Name of the column with the peptide sequences.
#' @param prefix Regular expression the residues before a match must end with.
#' @param suffix Regular expression the residues after a match must start with.
#' @return A data frame with one row per row of `pepinfo` and protein occurrence: the columns of `pepinfo`,
#'   `proteinID`, `Offset` (1-based start of the peptide), `proteinSequence`, `matched` and `lengthPeptide`.
#'   `NULL`, with a warning, if no peptide occurs in `fasta`.
#' @family protein inference
#' @export
#' @examples
#' file <- system.file("extdata/IDResults.txt.gz", package = "prozor")
#' peptides <- unique(utils::read.delim(file)$peptideSeq)
#' fasta_file <- system.file("p1000_db1_example/Annotation_canSeq.fasta.gz", package = "prozor")
#' fasta <- read_peptide_fasta(fasta_file)
#' annot <- annotate_peptides(peptides[1:20], fasta)
#' head(annot[, c("proteinID", "peptideSeq", "Offset", "lengthPeptide", "matched")])
annotate_peptides <- function(pepinfo, fasta, peptide = "peptideSeq", prefix = "(([RK])|(^)|(^M))", suffix = "") {
  if (is.null(dim(pepinfo))) {
    pepinfo <- data.frame(pepinfo)
    colnames(pepinfo) <- peptide
  }
  pepinfo[[peptide]] <- as.character(pepinfo[[peptide]])
  pepinfo$lengthPeptide <- nchar(pepinfo[[peptide]])

  restab <- annotate_aho(unique(pepinfo[[peptide]]), fasta, peptide = peptide)
  if (is.null(restab)) {
    warning("None of the peptides in `pepinfo` occurs in `fasta`.")
    return(NULL)
  }
  before <- substr(restab$proteinSequence, 1, restab$Offset - 1)
  after <- substring(restab$proteinSequence, restab$Offset + nchar(restab[[peptide]]))
  restab$matched <- grepl(paste0("(", prefix, ")$"), before) & grepl(paste0("^(", suffix, ")"), after)
  merge(restab, pepinfo, by.x = peptide, by.y = peptide)
}


#' Find peptides in protein sequences with the Aho-Corasick algorithm
#'
#' Searches all peptides at once with [AhoCorasickTrie::AhoCorasickSearch()]. Peptides that occur in no protein
#' are dropped.
#'
#' @param pepseq Peptide sequences, upper-case letters only.
#' @param fasta A list of protein sequences with unique names, as returned by [read_peptide_fasta()].
#' @param peptide Name of the result column that holds the peptide sequences.
#' @return A data frame with one row per occurrence and the columns `proteinID`, `peptide` (named by the
#'   argument), `Offset` (1-based start of the peptide) and `proteinSequence`; `NULL` if no peptide occurs.
#' @family protein inference
#' @export
#' @examples
#' file <- system.file("extdata/IDResults.txt.gz", package = "prozor")
#' peptides <- unique(utils::read.delim(file)$peptideSeq)
#' fasta_file <- system.file("p1000_db1_example/Annotation_canSeq.fasta.gz", package = "prozor")
#' fasta <- read_peptide_fasta(fasta_file)
#' hits <- annotate_aho(peptides[1:20], fasta)
#' dim(hits)
annotate_aho <- function(pepseq, fasta, peptide = "peptideSeq") {
  pepseq <- unique(trimws(pepseq))
  invalid <- pepseq[!grepl("^[A-Z]+$", pepseq)]
  if (length(invalid) > 0) {
    # AhoCorasickTrie 0.1.3 crashes R on keywords outside A-Z (https://github.com/chambm/AhoCorasickTrie/pull/5)
    stop(
      "Peptide sequences must consist of upper-case letters only; found: ",
      paste0("\"", utils::head(invalid, 5), "\"", collapse = ", "),
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
  sequences <- trimws(unlist(fasta, use.names = FALSE))
  names(sequences) <- protein_ids

  res <- AhoCorasickTrie::AhoCorasickSearch(pepseq, sequences, alphabet = "aminoacid")
  if (length(res) == 0) {
    return(NULL)
  }
  hits <- as.data.frame(dplyr::bind_rows(lapply(res, dplyr::bind_rows), .id = "proteinID"))
  colnames(hits)[colnames(hits) == "Keyword"] <- peptide
  proteins <- data.frame(proteinID = protein_ids, proteinSequence = unname(sequences))
  dplyr::inner_join(hits, proteins, by = "proteinID")
}
