.greedy_counts <- function(peptides, proteins) {
  n_pep <- nlevels(peptides)
  n_prot <- nlevels(proteins)
  # one number per distinct peptide-protein pair, sorted by protein and then by peptide
  pair <- sort(unique((as.integer(proteins) - 1) * n_pep + as.integer(peptides)))
  pep <- as.integer((pair - 1) %% n_pep) + 1L
  prot <- as.integer((pair - 1) %/% n_pep) + 1L
  # peptides of each protein and proteins of each peptide, both in ascending order
  peps_of <- split(pep, factor(prot, levels = seq_len(n_prot)))
  prots_of <- split(prot, factor(pep, levels = seq_len(n_pep)))
  # unexplained peptides per protein; -1 once the protein is assigned
  counts <- lengths(peps_of, use.names = FALSE)
  alive <- rep(TRUE, n_pep)
  mark <- logical(n_pep)
  res <- vector(n_prot, mode = "list")
  n_groups <- 0
  repeat {
    first <- which.max(counts)
    # stop when no protein explains an unexplained peptide
    if (length(first) == 0 || counts[first] <= 0) {
      break
    }
    explained <- peps_of[[first]][alive[peps_of[[first]]]]
    # the proteins of each explained peptide
    affected <- unlist(prots_of[explained], use.names = FALSE)
    sharing <- unique(affected)
    # Tie: proteins explaining the same peptides as the first tied protein are indistinguishable and reported
    # together; the other tied proteins compete in the next rounds. Such proteins share the explained peptides and
    # explain equally many, so those whose unexplained peptides are all explained by the first one explain exactly
    # them.
    mark[explained] <- TRUE
    tied <- sharing[counts[sharing] == counts[first]]
    same <- vapply(tied, function(j) all(mark[peps_of[[j]][alive[peps_of[[j]]]]]), logical(1))
    group <- sort(tied[same])
    n_groups <- n_groups + 1
    res[[n_groups]] <- list(prot = paste(levels(proteins)[group], collapse = ";"), peps = levels(peptides)[explained])
    # the explained peptides no longer count for the proteins that contain them
    counts[sharing] <- counts[sharing] - tabulate(match(affected, sharing), length(sharing))
    counts[group] <- -1L
    alive[explained] <- FALSE
    mark[explained] <- FALSE
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
#' Proteins are ranked by the number of distinct peptides they explain; a peptide-protein pair that occurs in
#' several rows counts once. Of the proteins that tie, the one first in `factor(annotation$protein_id)` is taken
#' first: alphabetical order for a character column, the level order for a factor column.
#'
#' @param annotation A data frame (or tibble) with one row per peptide-protein match, e.g. from
#'   [annotate_peptides()]: the protein IDs in the column `protein_id` and the peptide IDs in the column named by
#'   `peptide` (character or factor columns, without `NA`; unused factor levels are ignored).
#' @param peptide Name of the column with the peptide IDs. For precursors, add a column such as
#'   `paste(peptide_seq, precursor_charge, sep = ".")` and pass its name.
#' @return A named character vector with one element per distinct peptide ID, named by the peptide ID: the protein,
#'   or `;`-separated protein group, it is assigned to. The elements are in the order in which the groups are chosen,
#'   and the peptides of a group in the order of `factor(annotation[[peptide]])`.
#' @family protein inference
#' @export
#' @examples
#' data(protpepmetashort)
#' groups <- greedy_parsimony(protpepmetashort, peptide = "peptide_mod_seq")
#' # every modified peptide is assigned to exactly one protein group
#' stopifnot(setequal(names(groups), protpepmetashort$peptide_mod_seq))
#' head(groups)
greedy_parsimony <- function(annotation, peptide = "peptide_seq") {
  if (!is.data.frame(annotation)) {
    stop(
      "`annotation` must be a data frame of peptide-protein matches, e.g. from annotate_peptides(), not ",
      class(annotation)[1],
      "."
    )
  }
  columns <- c("protein_id", peptide)
  missing_columns <- setdiff(columns, names(annotation))
  if (length(missing_columns) > 0) {
    stop(
      "`annotation` must have the columns `protein_id` and `",
      peptide,
      "` (argument `peptide`); missing: ",
      paste0("`", missing_columns, "`", collapse = ", "),
      "."
    )
  }
  n_na <- vapply(columns, function(column) sum(is.na(annotation[[column]])), integer(1))
  if (any(n_na > 0)) {
    stop(
      "`annotation` must not have NA IDs; found ",
      paste0(n_na[n_na > 0], " in `", columns[n_na > 0], "`", collapse = " and "),
      "."
    )
  }
  # factor() drops unused levels of factor columns and keeps the order of the used ones
  groups <- .greedy_counts(factor(annotation[[peptide]]), factor(annotation[["protein_id"]]))
  peps <- lapply(groups, function(g) g$peps)
  prots <- vapply(groups, function(g) g$prot, character(1))
  stats::setNames(rep(prots, lengths(peps)), unlist(peps))
}


#' Convert a greedy_parsimony() result to a table
#'
#' @param groups A result of [greedy_parsimony()] whose names have the form `peptide.charge`, i.e. end in a dot and
#'   an integer charge; an error otherwise. The charge is taken after the last dot, so modified peptides such as
#'   `"AC[+57.0]K.2"` work.
#' @return A data frame with the columns `peptide`, `precursor_charge` (integer) and `protein_group`, and the
#'   precursor IDs as row names.
#' @family protein inference
#' @export
#' @examples
#' data(protpepmetashort)
#' protpepmetashort$precursor <- paste(
#'   protpepmetashort$peptide_mod_seq,
#'   protpepmetashort$precursor_charge,
#'   sep = "."
#' )
#' head(greedy_to_table(greedy_parsimony(protpepmetashort, peptide = "precursor")))
greedy_to_table <- function(groups) {
  precursor <- names(groups)
  no_charge <- precursor[!grepl("\\.[0-9]+$", precursor)]
  if (length(no_charge) > 0) {
    stop(
      "The names of `groups` must have the form `peptide.charge`, e.g. \"PEPTIDEK.2\"; found: ",
      paste0("\"", utils::head(no_charge, 5), "\"", collapse = ", "),
      "."
    )
  }
  data.frame(
    peptide = sub("\\.[0-9]+$", "", precursor),
    precursor_charge = as.integer(sub("^.*\\.", "", precursor)),
    protein_group = unname(groups),
    row.names = precursor
  )
}
