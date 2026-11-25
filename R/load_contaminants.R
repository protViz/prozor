# The bundled FASTA files of each set; the first set is the default of load_contaminants()
.contaminant_files <- list(
  fgcz2024 = c("fgcz_special_prot_20241112.fasta.gz", "fgcz_universal_contaminants_github_20241112.fasta.gz"),
  fgcz2019 = "fgcz_contaminants2019_20190708.fasta.gz",
  fgcz2021 = "fgcz_contaminants2021_20210929.fasta.gz",
  fgcz2022 = "fgcz_contaminants2022_20220405.fasta.gz",
  universal2024 = "fgcz_universal_contaminants_github_20241112.fasta.gz",
  special2024 = "fgcz_special_prot_20241112.fasta.gz"
)

#' Load a bundled set of contaminant sequences
#'
#' The sets:
#' * `"fgcz2024"` (default): the FGCZ default, the sets `"special2024"` and `"universal2024"` together; used by
#'   [create_decoy_db()], [create_fgcz_fasta_db()] and `inst/script/fgcz_create_fasta.R`.
#' * `"fgcz2019"`, `"fgcz2021"`, `"fgcz2022"`: FGCZ contaminant lists of 2019, 2021 and 2022.
#' * `"universal2024"`: universal contaminants from
#'   <https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics>,
#'   described by Frankenfield et al. (2022). See `system.file("COPYRIGHTS", package = "prozor")`.
#' * `"special2024"`: FGCZ special proteins of 2024.
#'
#' @param set Name of the contaminant set: `"fgcz2024"` (default), `"fgcz2019"`, `"fgcz2021"`, `"fgcz2022"`,
#'   `"universal2024"` or `"special2024"`.
#' @param no_human If `TRUE`, drop the sequences whose header contains `HUMAN` or `Homo sapiens`.
#' @return A named list of `SeqFastaAA` objects.
#' @references Frankenfield AM, Ni J, Ahmed M, Hao L (2022). Protein Contaminants Matter: Building Universal
#'   Protein Contaminant Libraries for DDA and DIA Proteomics. *Journal of Proteome Research*, 21(9), 2104-2113.
#'   \doi{10.1021/acs.jproteome.2c00145}
#' @family FASTA databases
#' @export
#' @examples
#' length(load_contaminants())
#' cont <- load_contaminants("fgcz2022")
#' length(cont)
#' length(load_contaminants("fgcz2022", no_human = TRUE))
load_contaminants <- function(
  set = c("fgcz2024", "fgcz2019", "fgcz2021", "fgcz2022", "universal2024", "special2024"),
  no_human = FALSE
) {
  set <- .match_method(set, names(.contaminant_files), arg = "set")
  files <- system.file("extdata", .contaminant_files[[set]], package = "prozor")
  contaminants <- do.call(c, lapply(files, read_fasta))
  if (no_human) {
    annot <- vapply(contaminants, attr, character(1), "Annot")
    contaminants <- contaminants[!grepl("HUMAN|Homo sapiens", annot)]
  }
  contaminants
}
