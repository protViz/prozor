# Stop unless `annotation` has the columns the coverage functions read.
.check_annotation <- function(annotation, peptide) {
  needed <- c("protein_id", "protein_sequence", "offset", peptide)
  missing_cols <- setdiff(needed, colnames(annotation))
  if (length(missing_cols) > 0) {
    stop(
      "`annotation` must have the columns ",
      paste0("`", needed, "`", collapse = ", "),
      " (as returned by annotate_peptides()); missing: ",
      paste0("`", missing_cols, "`", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
}

# One row per peptide occurrence: repeated rows of the same occurrence (e.g. charge states) count once.
.peptide_occurrences <- function(annotation, peptide) {
  occurrences <- data.frame(
    protein_id = as.character(annotation$protein_id),
    peptide = as.character(annotation[[peptide]]),
    offset = as.integer(annotation$offset),
    protein_sequence = as.character(annotation$protein_sequence)
  )
  occurrences[!duplicated(occurrences[, c("protein_id", "peptide", "offset")]), ]
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
  occurrences <- .peptide_occurrences(annotation, peptide)
  per_protein <- lapply(split(occurrences, occurrences$protein_id), function(x) {
    protein_seq <- x$protein_sequence[1]
    len <- nchar(protein_seq)
    ends <- x$offset + nchar(x$peptide) - 1L
    # residues covered by each distinct peptide, so overlapping occurrences of one peptide count once
    covered <- lapply(split(seq_len(nrow(x)), x$peptide), function(i) {
      unique(unlist(Map(seq.int, x$offset[i], ends[i]), use.names = FALSE))
    })
    depth <- tabulate(unlist(covered, use.names = FALSE), len)
    data.frame(
      protein_id = x$protein_id[1],
      position = seq_len(len),
      residue = strsplit(protein_seq, "", fixed = TRUE)[[1]],
      n_peptides = as.integer(depth)
    )
  })
  # rbind() drops the 0-row template unless there are no proteins, which gives a 0-row result
  template <- data.frame(protein_id = character(), position = integer(), residue = character(), n_peptides = integer())
  res <- do.call(rbind, c(list(template), unname(per_protein)))
  rownames(res) <- NULL
  res
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
  residues <- residue_coverage(annotation, peptide = peptide)
  occurrences <- .peptide_occurrences(annotation, peptide)
  depth <- split(residues$n_peptides, residues$protein_id)
  n_peptides <- tapply(occurrences$peptide, occurrences$protein_id, function(p) length(unique(p)))
  n_covered <- vapply(depth, function(d) sum(d > 0), integer(1), USE.NAMES = FALSE)
  protein_length <- lengths(depth, use.names = FALSE)
  data.frame(
    protein_id = names(depth),
    protein_length = protein_length,
    n_peptides = as.integer(n_peptides[names(depth)]),
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
  missing_cols <- setdiff(c("protein_id", "position", "residue", "n_peptides"), colnames(residues))
  if (length(missing_cols) > 0) {
    stop(
      "`residues` must be a result of residue_coverage(); missing columns: ",
      paste0("`", missing_cols, "`", collapse = ", "),
      "."
    )
  }
  if (length(protein_id) != 1 || !protein_id %in% residues$protein_id) {
    stop("`protein_id` must be one protein of `residues`, not \"", paste(protein_id, collapse = "\", \""), "\".")
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
