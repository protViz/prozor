#!/usr/bin/env Rscript

suppressMessages(library(docopt))
suppressMessages(library(prozor))

doc <- "Create fasta file with prozor

Usage:
  fgcz_create_fasta.R <fasta_dir>  [--contamin=<contamin>] [--revLab=<revLab>] [--output_dir=<output_dir>] [--summary]
  fgcz_create_fasta.R <fasta_dir> nodecoy  [--contamin=<contamin>] [--output_dir=<output_dir>] [--summary]


Options:
  -o --output_dir=<output_dir> output directory, default next to fasta_dir
  -c --contamin=<contamin> add contaminants [default: special] or `none`, or path to fasta file with contaminants
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

params <- c(
  "\nParameters used:\n\t",
  " fasta_dir: ",
  fasta_dir,
  "\n\t",
  "   nodecoy: ",
  nodecoy,
  "\n\t",
  "   summary: ",
  make_summary,
  "\n\t",
  "  contamin: ",
  contamin,
  "\n\t",
  "    revLab: ",
  rev_lab,
  "\n\t",
  "output_dir: ",
  output_dir,
  "\n\n\n\t"
)
cat(paste(params, collapse = ""))

contaminants <- switch(
  contamin,
  special = c(load_contaminants("special2024"), load_contaminants("universal2024")),
  fgcz2022 = load_contaminants("fgcz2022"),
  fgcz2021 = load_contaminants("fgcz2021"),
  fgcz2019 = load_contaminants("fgcz2019"),
  none = NULL,
  if (file.exists(contamin)) {
    read_peptide_fasta(contamin)
  } else {
    stop("--contamin must be special, fgcz2022, fgcz2021, fgcz2019, none or a FASTA file, not ", contamin)
  }
)

if (nodecoy) {
  rev_lab <- NULL
}

res_db <- create_fgcz_fasta_db(
  fasta_dir,
  contaminants = contaminants,
  rev_lab = rev_lab,
  output_dir = output_dir,
  make_summary = make_summary
)

cat(res_db$summary)

res <- paste(c("fgcz_create_fasta.R run with params\n", params), collapse = "")
writeLines(c(res, paste(res_db$summary, collapse = "")), con = paste0(res_db$dbname, ".txt"))
