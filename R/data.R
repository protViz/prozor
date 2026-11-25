#' Precursors of identified peptides
#'
#' Sequence, charge and m/z of identified precursors, without proteins; [protpepmetashort] holds the proteins of
#' some of them.
#'
#' @format A data frame with 4,938 rows and 4 columns:
#' * `peptide_seq`: peptide sequence (factor).
#' * `peptide_mod_seq`: modified peptide sequence (factor).
#' * `precursor_charge`: precursor charge (integer).
#' * `precursor_mz`: precursor m/z.
#' @family protein inference
"pepprot"

#' Peptide-protein assignments
#'
#' Peptide-protein pairs for 241 mouse proteins and 9 FGCZ contaminants; small, so the examples stay fast.
#'
#' @format A data frame with 423 rows and 7 columns:
#' * `peptide_seq`: peptide sequence.
#' * `protein_length`: length of the protein sequence.
#' * `protein_id`: UniProt protein ID, e.g. `sp|Q62009|POSTN_MOUSE`.
#' * `peptide_mod_seq`: modified peptide sequence.
#' * `precursor_charge`: precursor charge (integer).
#' * `precursor_mz`: precursor m/z.
#' * `peptide_length`: length of the peptide sequence.
#' @source Precursors of [pepprot], annotated with [annotate_peptides()] against mouse UniProt sequences and FGCZ
#'   contaminants.
#' @family protein inference
#' @name protpepmetashort
#' @docType data
#' @keywords datasets
NULL

#' MS1 precursor masses
#'
#' Precursor masses, the input of [Cdsw] for computing DIA windows.
#'
#' @format A numeric vector of 146,009 precursor m/z values.
#' @family DIA windows
"masses"

#' Search engine scores with protein IDs
#'
#' Target and decoy hits with two scores, the input of [TargetDecoyFDR].
#'
#' @format A data frame with 40,000 rows and 3 columns:
#' * `protein_id`: protein ID (factor); decoys start with `REV_`.
#' * `score`: score where smaller is better.
#' * `neg_log10_score`: `-log10(score)`, a score where larger is better.
#' @family FDR
"fdr_sample"

#' COMET search results
#'
#' Peptide-spectrum matches as produced by the COMET search engine, with its column names.
#'
#' @format A data frame with 4,557 rows and 18 columns:
#' * `scan`: scan number.
#' * `num`: rank of the match for the scan.
#' * `charge`: precursor charge.
#' * `exp_neutral_mass`, `calc_neutral_mass`: experimental and calculated neutral mass.
#' * `e.value`: E-value (smaller is better).
#' * `xcorr`: cross-correlation score.
#' * `delta_cn`: difference to the next best xcorr, normalized.
#' * `sp_score`: preliminary score (Sp).
#' * `ions_matched`, `ions_total`: matched and theoretical fragment ions.
#' * `plain_peptide`: peptide sequence.
#' * `modified_peptide`: peptide with flanking residues and modification masses.
#' * `prev_aa`, `next_aa`: residues before and after the peptide.
#' * `protein`: comma-separated protein IDs; decoys start with `REV_`.
#' * `protein_count`: number of proteins.
#' * `modifications`: modifications as position, residue and mass.
#' @source COMET search (low-resolution HCD settings) of the FGCZ quality-control run
#'   `FUSION_1.20180713_005_autoQC4L` against human UniProt sequences with reversed decoys.
#' @family FDR
"bcd"
