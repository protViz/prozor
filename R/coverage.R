# Stop unless `annotation` has the columns the coverage functions read.
.check_annotation <- function(annotation, peptide) {
  .check_columns(annotation, c("protein_id", "protein_sequence", "offset", peptide), "annotation")
}

# The residues covered by each distinct peptide of each protein. Returns the proteins (`protein_id`, in the order of
# factor(protein_id)), their sequences, their numbers of distinct peptides and, in `residue`, the covered residues: one
# entry per peptide and residue it covers, the residue numbered along the concatenated protein sequences, and in
# `protein` the protein of each entry. Repeated rows (e.g. charge states) and overlapping occurrences of a peptide
# count once; positions outside the protein are dropped.
.covered_residues <- function(annotation, peptide) {
  protein_id <- as.character(annotation$protein_id)
  proteins <- sort(unique(protein_id))
  protein_sequence <- as.character(annotation$protein_sequence)[match(proteins, protein_id)]
  protein_length <- nchar(protein_sequence)
  peptides <- as.character(annotation[[peptide]])
  protein <- match(protein_id, proteins)
  peptide_index <- match(peptides, peptides)
  from <- as.integer(annotation$offset)
  to <- from + nchar(peptides) - 1L
  ordered <- order(protein, peptide_index, from, method = "radix")
  protein <- protein[ordered]
  peptide_index <- peptide_index[ordered]
  from <- from[ordered]
  to <- to[ordered]
  # In this order, the occurrences of a peptide in a protein follow each other by offset; as they have the same
  # length, the previous one ends last, so each one adds only the residues after its end.
  same_pair <- protein == c(0L, utils::head(protein, -1)) & peptide_index == c(0L, utils::head(peptide_index, -1))
  from[same_pair] <- pmax(from[same_pair], to[which(same_pair) - 1L] + 1L)
  from <- pmax(from, 1L)
  to <- pmin(to, protein_length[protein])
  n_residues <- pmax(to - from + 1L, 0L)
  list(
    protein_id = proteins,
    protein_sequence = protein_sequence,
    n_peptides = tabulate(protein[!same_pair], length(proteins)),
    protein = rep.int(protein, n_residues),
    residue = sequence(n_residues, from = c(0L, cumsum(protein_length))[protein] + from)
  )
}


#' Peptide coverage of each residue
#'
#' Counts, for every residue of the proteins in `annotation`, how many distinct peptide sequences cover it.
#' Repeated rows of the same peptide, such as several charge states from [annotate_peptides()], count once. A
#' peptide that occurs twice in a protein covers both places, and where its occurrences overlap it counts once.
#'
#' @param annotation A data frame with the columns `protein_id`, `protein_sequence`, `offset` and the column named
#'   by `peptide`, as returned by [annotate_peptides()] or [annotate_aho()].
#' @inheritParams annotate_peptides
#' @return A data frame with one row per residue of each protein in `annotation`, ordered by `protein_id` and
#'   `position`: `protein_id`, `position`, `residue` and `n_peptides`, the number of peptides covering it; 0 rows
#'   if `annotation` has none.
#' @family sequence coverage
#' @export
#' @examples
#' file <- system.file("extdata/IDResults.txt.gz", package = "prozor")
#' peptides <- unique(utils::read.delim(file)$peptideSeq)
#' fasta_file <- system.file("p1000_db1_example/Annotation_canSeq.fasta.gz", package = "prozor")
#' annot <- annotate_peptides(peptides, read_fasta(fasta_file))
#' residues <- residue_coverage(annot)
#' head(residues)
residue_coverage <- function(annotation, peptide = "peptide_seq") {
  .check_annotation(annotation, peptide)
  covered <- .covered_residues(annotation, peptide)
  protein_length <- nchar(covered$protein_sequence)
  data.frame(
    protein_id = rep.int(covered$protein_id, protein_length),
    position = sequence(protein_length),
    residue = strsplit(paste(covered$protein_sequence, collapse = ""), "", fixed = TRUE)[[1]],
    n_peptides = tabulate(covered$residue, sum(protein_length))
  )
}


#' Sequence coverage of each protein
#'
#' Summarizes [residue_coverage()] per protein: how many distinct peptides match it and which fraction of its
#' residues they cover.
#'
#' @inheritParams residue_coverage
#' @return A data frame with one row per protein, ordered by `protein_id`: `protein_id`, `protein_length`,
#'   `n_peptides` (distinct peptides), `n_covered` (residues covered by at least one peptide) and `coverage`
#'   (`n_covered / protein_length`); 0 rows if `annotation` has none.
#' @family sequence coverage
#' @export
#' @examples
#' file <- system.file("extdata/IDResults.txt.gz", package = "prozor")
#' peptides <- unique(utils::read.delim(file)$peptideSeq)
#' fasta_file <- system.file("p1000_db1_example/Annotation_canSeq.fasta.gz", package = "prozor")
#' annot <- annotate_peptides(peptides, read_fasta(fasta_file))
#' coverage <- sequence_coverage(annot)
#' head(coverage[order(coverage$coverage, decreasing = TRUE), ])
sequence_coverage <- function(annotation, peptide = "peptide_seq") {
  .check_annotation(annotation, peptide)
  covered <- .covered_residues(annotation, peptide)
  protein_length <- nchar(covered$protein_sequence)
  n_covered <- tabulate(covered$protein[!duplicated(covered$residue)], length(protein_length))
  data.frame(
    protein_id = covered$protein_id,
    protein_length = protein_length,
    n_peptides = covered$n_peptides,
    n_covered = n_covered,
    coverage = n_covered / protein_length
  )
}


#' Plot the peptide coverage of a protein sequence
#'
#' Writes the protein sequence in rows of `width` residues. Each residue is a tile colored by the number of
#' peptides that cover it; uncovered residues are gray.
#'
#' @param residues A result of [residue_coverage()].
#' @param protein_id The protein to plot: one protein ID, a value of `residues$protein_id`.
#' @param width Residues per row: a positive whole number.
#' @return A ggplot object; the y axis labels the position of the first residue of each row.
#' @family sequence coverage
#' @export
#' @examples
#' file <- system.file("extdata/IDResults.txt.gz", package = "prozor")
#' peptides <- unique(utils::read.delim(file)$peptideSeq)
#' fasta_file <- system.file("p1000_db1_example/Annotation_canSeq.fasta.gz", package = "prozor")
#' annot <- annotate_peptides(peptides, read_fasta(fasta_file))
#' coverage <- sequence_coverage(annot)
#' best <- coverage$protein_id[which.max(coverage$coverage)]
#' plot_coverage(residue_coverage(annot), best)
plot_coverage <- function(residues, protein_id, width = 50) {
  .check_columns(residues, c("protein_id", "position", "residue", "n_peptides"), "residues")
  if (length(protein_id) != 1 || !protein_id %in% residues$protein_id) {
    stop("`protein_id` must be one protein of `residues`, not ", .quoted(protein_id, "\""), ".")
  }
  .check_number(width, "width", positive = TRUE, whole = TRUE)
  x <- residues[residues$protein_id == protein_id, ]
  x <- x[order(x$position), ]
  x$column <- (x$position - 1L) %% width + 1L
  x$row <- (x$position - 1L) %/% width
  bins <- c("0", "1", "2", "3", "4+")
  x$depth <- factor(bins[pmin(x$n_peptides, 4L) + 1L], levels = bins)
  x$text_colour <- ifelse(x$n_peptides >= 3, "white", "black")
  rows <- seq_len(max(x$row) + 1L) - 1L

  ggplot2::ggplot(x, ggplot2::aes(x = .data$column, y = .data$row)) +
    # show.legend = TRUE draws the keys of bins no residue falls into
    ggplot2::geom_tile(ggplot2::aes(fill = .data$depth), colour = "white", linewidth = 0.3, show.legend = TRUE) +
    ggplot2::geom_text(ggplot2::aes(label = .data$residue, colour = .data$text_colour), size = 2.8) +
    ggplot2::scale_fill_manual(
      values = c("0" = "grey90", "1" = "#C6DBEF", "2" = "#6BAED6", "3" = "#2171B5", "4+" = "#08306B"),
      drop = FALSE,
      name = "peptides"
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_y_reverse(breaks = rows, labels = rows * width + 1L, expand = c(0, 0)) +
    ggplot2::scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      x = NULL,
      y = "position",
      title = sprintf("%s: %.0f %% covered", protein_id, 100 * mean(x$n_peptides > 0))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}
