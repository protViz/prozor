#' Peptide-protein assignments
#'
#' Peptides with the columns `peptide_seq`, `peptide_mod_seq`, `precursor_charge` and `m` (precursor m/z).
#'
#' @name pepprot
#' @docType data
#' @keywords data
NULL

#' Peptide-protein assignments, small version
#'
#' Peptide-protein pairs with the columns `peptide_seq`, `protein_length`, `protein_id`, `peptide_mod_seq`,
#' `precursor_charge`, `m` (precursor m/z) and `peptide_length`; small, so the examples stay fast.
#'
#' @name protpepmetashort
#' @docType data
#' @keywords data
NULL

#' MS1 precursor masses
#'
#' About 150,000 MS1 precursor masses (m/z).
#'
#' @name masses
#' @docType data
#' @keywords data
NULL

#' Search engine scores with protein IDs
#'
#' A data frame with the columns `protein_id` (decoys start with `REV_`), `score` (smaller is better) and
#' `score2` (larger is better).
#'
#' @name fdr_sample
#' @docType data
#' @keywords data
NULL

#' COMET search results
#'
#' Peptide-spectrum matches as produced by the COMET search engine.
#'
#' @name bcd
#' @docType data
#' @keywords data
NULL
