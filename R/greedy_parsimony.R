.greedy2 <- function(pepprot) {
  res <- vector(ncol(pepprot), mode = "list")
  for (i in seq_len(ncol(pepprot))) {
    if (nrow(pepprot) == 0) {
      return(res[seq_len(i - 1)])
    }
    peps_per_prot <- colSums(pepprot)
    if (max(peps_per_prot) == 0) {
      return(res[seq_len(i - 1)])
    }
    idx <- which(peps_per_prot == max(peps_per_prot))
    # Tie: proteins explaining the same peptides as the first tied protein are indistinguishable and reported
    # together; the other tied proteins compete in the next rounds.
    if (length(idx) > 1) {
      first <- pepprot[, idx[1]] > 0
      same <- vapply(idx, function(j) identical(pepprot[, j] > 0, first), logical(1))
      idx <- idx[same]
    }
    explained <- pepprot[, idx[1]] > 0
    res[[i]] <- list(prot = paste(colnames(pepprot)[idx], collapse = ";"), peps = rownames(pepprot)[explained])
    pepprot <- pepprot[!explained, -idx, drop = FALSE]
  }
  res
}


#' Find a minimal protein set with a greedy algorithm
#'
#' Repeatedly assigns all unexplained peptides of the protein that explains the most of them. Proteins that
#' explain exactly the same peptides are reported together, separated by `;`.
#'
#' @param pepprot A peptide-by-protein matrix, as returned by [prepare_matrix()].
#' @return A named list with one element per peptide (row of `pepprot`): the protein, or `;`-separated protein
#'   group, it is assigned to.
#' @family protein inference
#' @export
#' @examples
#' data(protpepmetashort)
#' pepprot <- prepare_matrix(protpepmetashort, peptide_id = "peptideModSeq")
#' groups <- greedy_parsimony(pepprot)
#' # every precursor is assigned to exactly one protein group
#' stopifnot(setequal(names(groups), rownames(pepprot)))
#' head(unlist(groups))
greedy_parsimony <- function(pepprot) {
  groups <- .greedy2(pepprot)
  peps <- lapply(groups, function(g) g$peps)
  prots <- vapply(groups, function(g) g$prot, character(1))
  stats::setNames(as.list(rep(prots, lengths(peps))), unlist(peps))
}


#' Convert a greedy_parsimony() result to a table
#'
#' @param res A result of [greedy_parsimony()] whose names have the form `peptide.charge`. The charge is taken
#'   after the last dot, so modified peptides such as `"AC[+57.0]K.2"` work.
#' @return A data frame with the columns `Peptide`, `precursorCharge` and `Protein`, and the precursor IDs as row
#'   names.
#' @family protein inference
#' @export
#' @examples
#' data(protpepmetashort)
#' protpepmetashort$precursor <- paste(
#'   protpepmetashort$peptideModSeq,
#'   protpepmetashort$precursorCharge,
#'   sep = "."
#' )
#' pepprot <- prepare_matrix(protpepmetashort, peptide_id = "precursor")
#' head(greedy_to_table(greedy_parsimony(pepprot)))
greedy_to_table <- function(res) {
  precursor <- names(res)
  data.frame(
    Peptide = sub("\\.[^.]*$", "", precursor),
    precursorCharge = sub("^.*\\.", "", precursor),
    Protein = unlist(res, use.names = FALSE),
    row.names = precursor
  )
}
