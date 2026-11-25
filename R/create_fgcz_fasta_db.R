#' Summarize a FASTA database
#'
#' Counts the sequences and amino acids and summarizes the sequence lengths.
#'
#' @param fasta A list of sequences, as returned by [read_fasta()].
#' @return A list with the elements `n_sequences`, `length_summary` and `aa_freq` (counts per amino acid).
#' @family FASTA databases
#' @export
#' @examples
#' file <- system.file("extdata/fgcz_contaminants2022_20220405.fasta.gz", package = "prozor")
#' fasta <- read_fasta(file)
#' make_fasta_summary(fasta)
make_fasta_summary <- function(fasta) {
  sequences <- unlist(fasta, use.names = FALSE)
  # residue counts per byte, pasting 20,000 sequences at a time
  counts <- numeric(255)
  for (chunk in split(seq_along(sequences), ceiling(seq_along(sequences) / 20000))) {
    counts <- counts + tabulate(as.integer(charToRaw(paste(sequences[chunk], collapse = ""))), 255)
  }
  present <- which(counts > 0)
  aa_freq <- stats::setNames(counts[present], vapply(as.raw(present), rawToChar, character(1)))
  aa_freq <- aa_freq[order(names(aa_freq))]
  length_summary <- summary(as.numeric(nchar(sequences)))
  list(n_sequences = length(fasta), length_summary = length_summary, aa_freq = aa_freq)
}

# The sequence lengths and amino-acid counts of a make_fasta_summary() result as text, one value per line
.format_fasta_summary <- function(fasta_summary) {
  length_summary <- fasta_summary$length_summary
  aa_freq <- fasta_summary$aa_freq
  paste0(
    "length summary:\n",
    paste(sprintf("%15s : %s\n", names(length_summary), signif(unclass(length_summary), 6)), collapse = ""),
    "AA frequencies:\n",
    paste(sprintf("%15s : %s\n", names(aa_freq), aa_freq), collapse = "")
  )
}

# The FASTA files in `database_dir` and the first line of its `annotation.txt`, as list(files, annotation); stops
# unless `database_dir` and `output_dir` are existing, different directories and `database_dir` holds both
.database_files <- function(database_dir, output_dir) {
  if (!dir.exists(database_dir)) {
    stop("`database_dir` must be an existing directory, not \"", database_dir, "\".", call. = FALSE)
  }
  if (!dir.exists(output_dir)) {
    stop("`output_dir` must be an existing directory, not \"", output_dir, "\".", call. = FALSE)
  }
  if (normalizePath(output_dir) == normalizePath(database_dir)) {
    stop(
      "`output_dir` must not be `database_dir` (\"",
      database_dir,
      "\"): the next run would read the database written there as targets.",
      call. = FALSE
    )
  }
  files <- list.files(database_dir, pattern = "\\.(fasta|fa|faa)(\\.gz)?$", full.names = TRUE)
  if (length(files) == 0) {
    stop(
      "`database_dir` must contain FASTA files (*.fasta, *.fa or *.faa, optionally gzipped), found none in \"",
      database_dir,
      "\".",
      call. = FALSE
    )
  }
  annot_file <- file.path(database_dir, "annotation.txt")
  if (!file.exists(annot_file)) {
    stop("`database_dir` must contain the file `annotation.txt`; \"", database_dir, "\" has none.", call. = FALSE)
  }
  list(files = files, annotation = readLines(annot_file, warn = FALSE)[1])
}

#' Create an FGCZ FASTA database from a directory of FASTA files
#'
#' Reads all FASTA files in `database_dir` (`*.fasta`, `*.fa` or `*.faa`, optionally gzipped), adds
#' contaminants and decoys with [create_decoy_db()], which drops duplicated entries, and writes the result to
#' `output_dir` as `<name>[_d]_<date>.fasta`, where `_d` marks a database with decoys. The first line of
#' `annotation.txt` in `database_dir` becomes the header of the dummy source entry.
#'
#' The whole database, decoys included, is held in memory. Peak memory grows linearly with the database, by about
#' 4 KB per written sequence or 6 to 7 times the size of the written file: 250,000 target proteins (500,000
#' sequences with decoys, a 320 MB file) take 2.2 GB.
#'
#' @param database_dir Directory with the FASTA files and `annotation.txt`.
#' @param output_dir Existing directory the FASTA file is written to; must not be `database_dir`, where the next run
#'   would read the written database as targets.
#' @param contaminants Contaminant sequences to add, e.g. from [load_contaminants()], or `NULL` for none.
#' @inheritParams create_decoy_db
#' @param make_summary If `TRUE`, append the sequence lengths and amino-acid counts of [make_fasta_summary()] to
#'   the summary.
#' @return A list with the elements `fasta` (the sequences written), `file_path`, `summary` (text with the argument
#'   values, the file name and the number of sequences), `call` and `db_name`.
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
  output_dir,
  contaminants = load_contaminants(),
  rev_lab = "REV_",
  make_summary = TRUE
) {
  mcall <- match.call()
  input <- .database_files(database_dir, output_dir)
  annotation <- input$annotation
  fasta <- create_decoy_db(input$files, contaminants = contaminants, rev_lab = rev_lab, header = annotation)

  db_name <- basename(database_dir)
  if (!is.null(rev_lab)) {
    db_name <- paste0(db_name, "_d")
  }
  file_path <- file.path(output_dir, paste0(db_name, "_", format(Sys.time(), "%Y%m%d"), ".fasta"))
  message("Writing database to ", file_path, ".")
  write_fasta(fasta, file = file_path)
  values <- c(
    database_dir = database_dir,
    contaminants = if (is.null(contaminants)) "NULL" else paste(length(contaminants), "sequences"),
    rev_lab = if (is.null(rev_lab)) "NULL" else rev_lab,
    output_dir = output_dir,
    `FASTA name` = db_name,
    `FASTA file name` = basename(file_path),
    annotation = annotation,
    n_sequences = length(fasta)
  )
  summary_str <- paste0(
    "Database created with prozor ",
    utils::packageVersion("prozor"),
    "\n",
    paste(sprintf("%15s : %s\n", names(values), values), collapse = "")
  )

  if (make_summary) {
    summary_str <- paste0(summary_str, .format_fasta_summary(make_fasta_summary(fasta)))
  }
  list(fasta = fasta, file_path = file_path, summary = summary_str, call = mcall, db_name = db_name)
}
