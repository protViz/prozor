#' Create a database with contaminants and decoys
#'
#' Reads one or more FASTA files, adds contaminants and reversed decoy sequences, and puts first a dummy
#' entry `CRAPCRAPCRAP` whose header records the source of the database. Of entries with the same name, only
#' the first is kept. See `vignette("CreatingDecoyDB", package = "prozor")`.
#'
#' @param dbs Paths to one or more FASTA files.
#' @param contaminants A list of contaminant sequences, e.g. from [load_contaminants()], or `NULL` for none.
#' @param rev_lab Prefix of the reversed decoy sequences, or `NULL` for no decoys.
#' @param annot Header of the dummy entry that records the source of the database, with or without the leading `>`.
#' @return A list of `SeqFastaAA` objects with unique names: the dummy entry, the targets, the contaminants and,
#'   unless `rev_lab` is `NULL`, the reversed decoys of the targets and contaminants.
#' @family FASTA databases
#' @export
#' @examples
#' file <- system.file("extdata/fgcz_contaminants2022_20220405.fasta.gz", package = "prozor")
#' targets <- read_peptide_fasta(file)
#' cont <- load_contaminants("fgcz2022")
#'
#' # the file, read twice, holds the same entries as the default contaminants: each is kept once
#' db <- create_decoy_db(c(file, file), contaminants = cont)
#' stopifnot(length(db) == 1 + 2 * length(targets), identical(names(targets), names(cont)))
#'
#' db <- create_decoy_db(file, rev_lab = NULL, contaminants = NULL)
#' stopifnot(length(db) == 1 + length(targets))
create_decoy_db <- function(
  dbs,
  contaminants = c(load_contaminants("special2024"), load_contaminants("universal2024")),
  rev_lab = "REV_",
  annot = "zz|sourceOf|database"
) {
  header <- paste0(">", sub("^>", "", annot))
  dummy <- .as_seq_fasta_aa("CRAPCRAPCRAP", name = substring(sub(" .*", "", header), 2), annot = header)
  db <- do.call(
    c,
    lapply(dbs, function(file) {
      message("Reading database ", file, ".")
      read_peptide_fasta(file)
    })
  )
  db <- c(db, contaminants)
  if (!is.null(rev_lab)) {
    db <- c(db, reverse_seq(db, rev_lab = rev_lab))
  }
  db <- c(list(dummy), db)
  db[!duplicated(names(db))]
}
