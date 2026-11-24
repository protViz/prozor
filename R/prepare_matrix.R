#' Build a sparse peptide-by-protein matrix
#'
#' @param data A data frame (or tibble) of peptide-protein pairs, e.g. from [annotate_peptides()].
#' @param protein_id Name of the protein ID column (character or factor; unused factor levels are dropped).
#' @param peptide_id Name of the peptide or precursor ID column (character or factor; unused factor levels are
#'   dropped).
#' @param weighting `NULL` (default) gives each peptide-protein pair the weight 1. `"inverse"` divides the weights
#'   of each peptide by the number of proteins it matches. A numeric vector gives one weight per row of `data`.
#' @param use.last.ij Passed to [Matrix::sparseMatrix()]; if `TRUE`, a repeated peptide-protein pair keeps the
#'   weight of its last row, if `FALSE` the sum of the weights of its rows. Applies to every `weighting`.
#' @return A sparse matrix with peptides in rows and proteins in columns.
#' @family protein inference
#' @export
#' @examples
#' data(protpepmetashort)
#' count <- prepare_matrix(protpepmetashort)
#' dim(count)
#' inverse <- prepare_matrix(protpepmetashort, weighting = "inverse")
#' # a peptide shared by n proteins gets the weight 1/n in each of them
#' table(round(as.vector(inverse[inverse > 0]), 2))
prepare_matrix <- function(
  data,
  protein_id = "protein_id",
  peptide_id = "peptide_seq",
  weighting = NULL,
  use.last.ij = TRUE
) {
  weighted <- is.numeric(weighting) && length(weighting) == nrow(data)
  if (!(is.null(weighting) || weighted || identical(weighting, "inverse"))) {
    stop(
      "`weighting` must be NULL, \"inverse\" or a numeric vector with one weight per row of `data`, not a ",
      class(weighting)[1],
      " of length ",
      length(weighting),
      "."
    )
  }
  # factor() drops unused levels of factor columns and keeps the order of the used ones
  proteins <- factor(data[[protein_id]])
  peptides <- factor(data[[peptide_id]])
  pep_prot <- Matrix::sparseMatrix(
    as.integer(peptides),
    as.integer(proteins),
    x = if (weighted) weighting else 1,
    use.last.ij = use.last.ij,
    dimnames = list(levels(peptides), levels(proteins))
  )
  if (identical(weighting, "inverse")) {
    # divide the stored entries by their row sums; sweep() would build a dense peptide x protein array
    pep_prot@x <- pep_prot@x / unname(rowSums(pep_prot))[pep_prot@i + 1L]
  }
  pep_prot
}
