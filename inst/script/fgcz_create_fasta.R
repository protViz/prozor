#!/usr/bin/env Rscript

suppressMessages(library(docopt))
suppressMessages(library(prozor))

doc <- "Create fasta file with prozor

Usage:
  fgcz_create_fasta.R <fasta_dir>  [--contamin=<contamin>] [--revLab=<revLab>] [--output_dir=<output_dir>] [--summary]
  fgcz_create_fasta.R <fasta_dir> nodecoy  [--contamin=<contamin>] [--output_dir=<output_dir>] [--summary]


Options:
  -o --output_dir=<output_dir> output directory, default next to fasta_dir
  -c --contamin=<contamin> add contaminants: special (special2024 and universal2024), `none`, a set of
                           load_contaminants() such as fgcz2022, or a FASTA file [default: special]
  -r --revLab=<revLab> create reverse sequences with prefix [default: REV_].


Arguments:
  fasta_dir input file
"

opt <- docopt(doc)

fasta_dir <- opt$fasta_dir
nodecoy <- opt$nodecoy
make_summary <- opt$summary
contamin <- opt[["--contamin"]]
rev_lab <- opt[["--revLab"]]
output_dir <- if (is.null(opt[["--output_dir"]])) dirname(fasta_dir) else opt[["--output_dir"]]
if (nodecoy) {
  rev_lab <- NULL
}

params <- c(
  fasta_dir = fasta_dir,
  nodecoy = nodecoy,
  summary = make_summary,
  contamin = contamin,
  revLab = if (is.null(rev_lab)) "NULL" else rev_lab,
  output_dir = output_dir
)
params_text <- paste0(
  "\nParameters used:\n",
  paste0(sprintf("\t%10s: %s\n", names(params), params), collapse = ""),
  "\n"
)
cat(params_text)

contaminants <- switch(
  contamin,
  special = c(load_contaminants("special2024"), load_contaminants("universal2024")),
  none = NULL,
  if (file.exists(contamin)) read_fasta(contamin) else load_contaminants(contamin)
)

res_db <- create_fgcz_fasta_db(
  fasta_dir,
  contaminants = contaminants,
  rev_lab = rev_lab,
  output_dir = output_dir,
  make_summary = make_summary
)

cat(res_db$summary)

writeLines(
  c(paste0("fgcz_create_fasta.R run with params\n", params_text), paste(res_db$summary, collapse = "")),
  con = paste0(res_db$db_name, ".txt")
)
