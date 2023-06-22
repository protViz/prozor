#!/usr/bin/Rscript

suppressMessages(library(docopt))
suppressMessages(library(prozor))

doc <- "Create fasta file with prozor

Usage:
  fgcz_create_fasta.R <fasta_dir>  [--contamin=<contamin>] [--revLab=<revLab>] [--output_dir=<output_dir>] [--summary]
  fgcz_create_fasta.R <fasta_dir> nodecoy  [--contamin=<contamin>] [--output_dir=<output_dir>]


Options:
  -o --output_dir=<output_dir> output directory, default next to fasta_dir
  -c --contamin=<contamin> add contaminants [default: fgcz2022] or `none`, or path to fasta file with contaminants
  -r --revLab=<revLab> create reverse sequences with prefix [default: REV_].


Arguments:
  fasta_dir input file
"

if (FALSE) {
  args <- c("C:\\Users\\wewol\\Dropbox\\DataAnalysis\\p65\\fgcz_9606_SARS_CoV_2_reviewed_cnl", "-o" ,"c:/users/wewol")
  #args <- c("C:\\Users\\wewol\\Dropbox\\DataAnalysis\\p65\\fgcz_10116_RattusNor_reviewed_cnl")
  args <- c("Z:/p65/Proteomics/fasta_db/o24206_db1_Synechococcus_sp_PCC7336", "nodecoy", "-o", "c:/users/wewol")
  args <- c("/Users/witoldwolski/p2799_d_db10_Vibrio_O1")
  args <- c("/Users/witoldwolski/p2799_d_db10_Vibrio_O1", "--summary")
  opt <- docopt(doc,args = args)
}else{
  opt <- docopt(doc)
}


params <- c("\nParameters used:\n\t",
            " fasta_dir: ", fasta_dir <- opt$fasta_dir, "\n\t",
            "   nodecoy: ",   nodecoy <- opt$nodecoy, "\n\t",
            "   summary: ",   summary <- opt$summary, "\n\t",
            "  contamin: ", contamin <- opt[["--contamin"]], "\n\t",
            "    revLab: ", revLab <- opt[["--revLab"]], "\n\t",
            "output_dir: ", output_dir <- opt[["--output_dir"]], "\n\n\n\t")

cat(paste(params, collapse = ""))

if (contamin == "fgcz2022") {
  contamin <- loadContaminantsFGCZ2022()
} else if (contamin == "fgcz2021") {
  contamin <- loadContaminantsFasta2021()
} else if (contamin == "fgcz2019") {
  contamin <- loadContaminantsFasta2019()
}else if (contamin == "none") {
  contamin <- NULL
}else{
  if (file.exists(contamin)) {
    contamin <- prozor::readPeptideFasta(contamin)
  }
}

if (nodecoy) {
  revLab <- NULL
}

resDB <-  create_fgcz_fasta_db(fasta_dir,
                               useContaminants = contamin,
                               revLab = revLab,
                               outputdir = output_dir,
                               summary = opt$summary)

cat(resDB$summary)

res <- paste(c("fgcz_create_fasta.R run with params\n",params),collapse = "")
writeLines( c(res, paste(resDB$summary, collapse = "")), con = paste0(resDB$dbname,".txt"))

