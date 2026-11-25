#' @section FASTA databases:
#' [read_fasta()] and [write_fasta()] read and write protein sequences. [create_fgcz_fasta_db()] builds a database
#' from a directory of FASTA files with [create_decoy_db()], which adds contaminants from [load_contaminants()] and
#' reversed decoys from [reverse_seq()]; [make_fasta_summary()] summarizes it. [remove_signal_peptide()] cuts signal
#' peptides off UniProt sequences. See `vignette("CreatingDecoyDB", package = "prozor")`.
#'
#' @section Peptide annotation and protein inference:
#' [annotate_peptides()] finds the proteins that contain each peptide, with [annotate_aho()].
#' [greedy_parsimony()] finds the minimal set of proteins that explains the peptides, and [greedy_to_table()] turns
#' its result into a table. See `vignette("PeptideAnnotationwithProzor", package = "prozor")`.
#'
#' @section Sequence coverage:
#' [residue_coverage()] and [sequence_coverage()] compute the peptide coverage of each residue and each protein from
#' the result of [annotate_peptides()]; [plot_coverage()] plots it for one protein. See
#' `vignette("ProteinCoverage", package = "prozor")`.
#'
#' @section Target-decoy FDR:
#' [TargetDecoyFDR] estimates the false discovery rate and the q-values of target-decoy search results, finds the
#' score cut-off for a given FDR, and plots the score distributions and q-values. See
#' `vignette("TargetDecoyFDR_Example", package = "prozor")` and `vignette("TargetDecoyComet", package = "prozor")`.
#'
#' @section DIA windows:
#' [Cdsw] computes precursor isolation windows with about equal numbers of precursors from the MS1 masses;
#' [readjust_windows()] moves their boundaries to regions with few precursors, and [plot_windows()] shows them. See
#' `vignette("cdsw", package = "prozor")`.
#'
#' @keywords internal
#' @importFrom ggplot2 .data
#' @importFrom R6 R6Class
"_PACKAGE"
