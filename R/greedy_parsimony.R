.greedy2 <- function(pepprot) {
  # a protein explains the peptides with a positive entry; the weights themselves do not count
  pepprot <- pepprot > 0
  res <- vector(ncol(pepprot), mode = "list")
  n_groups <- 0
  # stops when no peptide-protein pair is left: no peptide, no protein, or only peptides without a protein
  while (any(pepprot)) {
    peps_per_prot <- colSums(pepprot)
    idx <- which(peps_per_prot == max(peps_per_prot))
    explained <- pepprot[, idx[1]]
    # Tie: proteins explaining the same peptides as the first tied protein are indistinguishable and reported
    # together; the other tied proteins compete in the next rounds. Tied proteins explain equally many peptides,
    # so those explaining all peptides of the first one explain exactly them.
    if (length(idx) > 1) {
      idx <- idx[colSums(pepprot[explained, idx, drop = FALSE]) == sum(explained)]
    }
    n_groups <- n_groups + 1
    res[[n_groups]] <- list(prot = paste(colnames(pepprot)[idx], collapse = ";"), peps = rownames(pepprot)[explained])
    pepprot <- pepprot[!explained, -idx, drop = FALSE]
  }
  res[seq_len(n_groups)]
}


#' Find a minimal protein set with a greedy algorithm
#'
#' Repeatedly assigns all unexplained peptides of the protein that explains the most of them. Proteins that
#' explain exactly the same unexplained peptides are reported together, separated by `;`. The groups are formed
#' on the peptides that remain in each round: proteins B (peptides p1, p5, p6) and C (p2, p5, p6) form the group
#' `"B;C"` for p5 and p6 after protein A has taken p1 and p2.
#'
#' A protein explains a peptide when their entry in `pepprot` is positive. Proteins are ranked by the number of
#' peptides they explain, not by the sum of their weights, so a matrix from `prepare_matrix(weighting = "inverse")`
#' gives the same result as one with weights 1.
#'
#' @param pepprot A peptide-by-protein matrix with peptide row names and protein column names, as returned by
#'   [prepare_matrix()].
#' @return A named list with one element per peptide that matches a protein (row of `pepprot` with a positive
#'   entry): the protein, or `;`-separated protein group, it is assigned to. Peptides without a protein are left
#'   out.
#' @family protein inference
#' @export
#' @examples
#' data(protpepmetashort)
#' pepprot <- prepare_matrix(protpepmetashort, peptide_id = "peptide_mod_seq")
#' groups <- greedy_parsimony(pepprot)
#' # every precursor is assigned to exactly one protein group
#' stopifnot(setequal(names(groups), rownames(pepprot)))
#' head(unlist(groups))
greedy_parsimony <- function(pepprot) {
  no_names <- c("row", "column")[c(is.null(rownames(pepprot)), is.null(colnames(pepprot)))]
  if (length(no_names) > 0) {
    stop(
      "`pepprot` must have peptide row names and protein column names; it has no ",
      paste(no_names, collapse = " and "),
      " names."
    )
  }
  groups <- .greedy2(pepprot)
  peps <- lapply(groups, function(g) g$peps)
  prots <- vapply(groups, function(g) g$prot, character(1))
  stats::setNames(as.list(rep(prots, lengths(peps))), unlist(peps))
}


#' Convert a greedy_parsimony() result to a table
#'
#' @param res A result of [greedy_parsimony()] whose names have the form `peptide.charge`, i.e. end in a dot and
#'   an integer charge; an error otherwise. The charge is taken after the last dot, so modified peptides such as
#'   `"AC[+57.0]K.2"` work.
#' @return A data frame with the columns `peptide`, `precursor_charge` (integer) and `protein`, and the precursor
#'   IDs as row names.
#' @family protein inference
#' @export
#' @examples
#' data(protpepmetashort)
#' protpepmetashort$precursor <- paste(
#'   protpepmetashort$peptide_mod_seq,
#'   protpepmetashort$precursor_charge,
#'   sep = "."
#' )
#' pepprot <- prepare_matrix(protpepmetashort, peptide_id = "precursor")
#' head(greedy_to_table(greedy_parsimony(pepprot)))
greedy_to_table <- function(res) {
  precursor <- names(res)
  no_charge <- precursor[!grepl("\\.[0-9]+$", precursor)]
  if (length(no_charge) > 0) {
    stop(
      "The names of `res` must have the form `peptide.charge`, e.g. \"PEPTIDEK.2\"; found: ",
      paste0("\"", utils::head(no_charge, 5), "\"", collapse = ", "),
      "."
    )
  }
  data.frame(
    peptide = sub("\\.[0-9]+$", "", precursor),
    precursor_charge = as.integer(sub("^.*\\.", "", precursor)),
    protein = unlist(res, use.names = FALSE),
    row.names = precursor
  )
}
