#' Summarize a FASTA database
#'
#' Counts the sequences and amino acids and summarizes the sequence lengths.
#'
#' @param db A list of sequences, as returned by [read_peptide_fasta()].
#' @param old If `TRUE`, count the amino acids by concatenating all sequences (the older, slower method).
#' @param as_string If `TRUE`, return text to pass to [cat()]; if `FALSE`, a list.
#' @return With `as_string = TRUE`, a character vector to pass to [cat()]. Otherwise a list with the elements
#'   `nr_sequences`, `length_summary` and `aa_freq` (counts per amino acid).
#' @family FASTA databases
#' @export
#' @examples
#' file <- system.file("extdata/fgcz_contaminants2022_20220405.fasta.gz", package = "prozor")
#' db <- read_peptide_fasta(file)
#' cat(make_fasta_summary(db))
#' make_fasta_summary(db, as_string = FALSE)
make_fasta_summary <- function(db, old = FALSE, as_string = TRUE) {
  if (old) {
    aa_freq <- table(strsplit(paste(db, collapse = ""), split = "")[[1]])
  } else {
    counts <- unlist(unname(lapply(db, function(x) table(strsplit(x, split = "")[[1]]))))
    aa_freq <- vapply(split(as.numeric(counts), names(counts)), sum, numeric(1))
  }
  length_summary <- summary(vapply(db, seqinr::getLength, numeric(1)))

  if (as_string) {
    c(
      "nr sequences:\n",
      length(db),
      "\n length summary:\n",
      paste0(utils::capture.output(length_summary), "\n"),
      "AA frequencies:\n",
      paste0(utils::capture.output(as.matrix(aa_freq)), "\n")
    )
  } else {
    list(nr_sequences = length(db), length_summary = length_summary, aa_freq = aa_freq)
  }
}


#' Create an FGCZ FASTA database from a directory of FASTA files
#'
#' Reads all FASTA files in `database_dir`, adds contaminants and decoys with [create_decoy_db()], drops
#' duplicated entries and writes the result to `output_dir` as `<name>[_d]_<date>.fasta`, where `_d` marks a
#' database with decoys. The first line of `annotation.txt` in `database_dir` becomes the header of the
#' dummy source entry.
#'
#' @param database_dir Directory with the FASTA files and `annotation.txt`.
#' @param contaminants Contaminant sequences to add, e.g. from [load_contaminants()], or `NULL` for none.
#' @param rev_lab Prefix of the reversed decoy sequences, or `NULL` for no decoys.
#' @param output_dir Directory the FASTA file is written to.
#' @param make_summary If `TRUE`, append the output of [make_fasta_summary()] to the summary.
#' @return A list with the elements `db` (the sequences written), `file_path`, `summary` (text), `call` and
#'   `db_name`.
#' @family FASTA databases
#' @export
#' @examples
#' res <- create_fgcz_fasta_db(
#'   system.file("p1000_db1_example", package = "prozor"),
#'   output_dir = tempdir(),
#'   make_summary = FALSE
#' )
#' basename(res$file_path)
#' unlink(res$file_path)
create_fgcz_fasta_db <- function(
  database_dir,
  contaminants = load_contaminants("fgcz2022"),
  rev_lab = "REV_",
  output_dir,
  make_summary = TRUE
) {
  mcall <- match.call()

  if (!dir.exists(database_dir)) {
    stop("`database_dir` must be an existing directory, not \"", database_dir, "\".")
  }
  db_name <- basename(database_dir)
  files <- file.path(database_dir, grep("fasta", dir(database_dir), value = TRUE))
  annot_file <- grep("annotation.txt", dir(database_dir), value = TRUE)
  if (length(annot_file) == 0) {
    stop("Can't find an annotation file (annotation.txt) in `database_dir` \"", database_dir, "\".")
  }
  annotation <- readLines(file.path(database_dir, annot_file))[1]

  db <- create_decoy_db(files, contaminants = contaminants, annot = annotation, rev_lab = rev_lab)
  db <- db[!duplicated(names(db))]

  if (!is.null(rev_lab)) {
    db_name <- paste0(db_name, "_d")
  }
  file_path <- file.path(output_dir, paste0(db_name, "_", format(Sys.time(), "%Y%m%d"), ".fasta"))
  message("Writing database to ", file_path, ".")
  write_fasta(db, file = file_path)
  summary_str <- paste0(
    "Database created with prozor: ",
    paste(utils::capture.output(mcall), collapse = "\n"),
    "\n",
    "\n      FASTA name : ",
    db_name,
    "\n FASTA file name : ",
    basename(file_path),
    "\n  written to dir : ",
    output_dir,
    "\n       annotation:\n",
    annotation,
    "\n",
    "\n      nr sequences: ",
    length(db),
    "\n"
  )

  if (make_summary) {
    summary_str <- c(summary_str, make_fasta_summary(db, as_string = TRUE))
  }
  list(db = db, file_path = file_path, summary = summary_str, call = mcall, db_name = db_name)
}
