#' Build a sparse peptide-by-protein matrix
#'
#' @param data A data frame (or tibble) of peptide-protein pairs, e.g. from [annotate_peptides()].
#' @param protein_id Name of the protein ID column.
#' @param peptide_id Name of the peptide or precursor ID column.
#' @param weighting `NULL` (default) gives each peptide-protein pair the weight 1. `"inverse"` divides the weights
#'   of each peptide by the number of proteins it matches. A numeric vector gives one weight per row of `data`.
#' @param sep Not used.
#' @param use.last.ij Passed to [Matrix::sparseMatrix()]; if `TRUE`, a repeated peptide-protein pair keeps a single
#'   entry instead of the sum of its entries.
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
  sep = "|",
  use.last.ij = TRUE
) {
  proteins <- as.factor(data[[protein_id]])
  peptides <- as.factor(data[[peptide_id]])
  i <- as.integer(peptides)
  j <- as.integer(proteins)

  if (is.null(weighting)) {
    pep_prot <- Matrix::sparseMatrix(i, j, x = 1, use.last.ij = use.last.ij)
  } else if (is.numeric(weighting) && length(weighting) == nrow(data)) {
    pep_prot <- Matrix::sparseMatrix(i, j, x = weighting)
  } else if (identical(weighting, "inverse")) {
    pep_prot <- Matrix::sparseMatrix(i, j, x = 1, use.last.ij = use.last.ij)
    pep_prot <- sweep(pep_prot, 1, rowSums(pep_prot), "/")
  } else {
    stop(
      "`weighting` must be NULL, \"inverse\" or a numeric vector with one weight per row of `data`, not a ",
      class(weighting)[1],
      " of length ",
      length(weighting),
      "."
    )
  }
  colnames(pep_prot) <- levels(proteins)
  rownames(pep_prot) <- levels(peptides)
  pep_prot
}
