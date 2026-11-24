.contaminant_files <- c(
  fgcz2019 = "fgcz_contaminants2019_20190708.fasta.gz",
  fgcz2021 = "fgcz_contaminants2021_20210929.fasta.gz",
  fgcz2022 = "fgcz_contaminants2022_20220405.fasta.gz",
  universal2024 = "fgcz_universal_contaminants_github_20241112.fasta.gz",
  special2024 = "fgcz_special_prot_20241112.fasta.gz"
)

#' Load a bundled set of contaminant sequences
#'
#' The sets:
#' * `"fgcz2019"`, `"fgcz2021"`, `"fgcz2022"`: FGCZ contaminant lists of 2019, 2021 and 2022.
#' * `"universal2024"`: universal contaminants from
#'   <https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics>.
#' * `"special2024"`: FGCZ special proteins of 2024.
#'
#' @param set Name of the contaminant set: `"fgcz2019"`, `"fgcz2021"`, `"fgcz2022"`, `"universal2024"` or
#'   `"special2024"`.
#' @param no_human If `TRUE`, drop the sequences whose header contains `HUMAN`.
#' @return A named list of `SeqFastaAA` objects.
#' @family FASTA databases
#' @export
#' @examples
#' cont <- load_contaminants("fgcz2022")
#' length(cont)
#' length(load_contaminants("fgcz2022", no_human = TRUE))
load_contaminants <- function(set, no_human = FALSE) {
  if (length(set) != 1 || !set %in% names(.contaminant_files)) {
    stop(
      "`set` must be one of ",
      paste0("\"", names(.contaminant_files), "\"", collapse = ", "),
      ", not ",
      paste0("\"", set, "\"", collapse = ", "),
      "."
    )
  }
  contaminants <- read_peptide_fasta(system.file("extdata", .contaminant_files[[set]], package = "prozor"))
  if (no_human) {
    annot <- vapply(contaminants, seqinr::getAnnot, character(1))
    contaminants <- contaminants[!grepl("HUMAN", annot)]
  }
  contaminants
}
