#' load list of contaminant sequences
#'
#' @export
#' @examples
#' library(prozor)
#' cont <- loadContaminantsFasta()
#' cont[[1]]
#' #example how to create a protein db with decoy sequences
loadContaminantsFasta <- function(){
  file = file.path(path.package("prozor"),"extdata/fgcz_contaminants_20150123.fasta")
  contaminants <- readPeptideFasta(file)
}
#' load human
#' @export
#' @examples
#' library(prozor)
#' fastaMusMusc <- loadHomoSapiensFasta()
#'
loadHomoSapiensFasta <- function(){
  file = file.path(path.package("prozor"),"extdata/uniprot_homo_sapiens_9609_23032016.fasta.gz")
  fasta <- readPeptideFasta(file)
}
#'
#' load human signal peptides
#' @export
#' @examples
#' library(prozor)
#' signal <- loadHomoSapiensSignalPeptides()
#'
loadHomoSapiensSignalPeptides <- function(){
  file = file.path(path.package("prozor"),"extdata/uniprot_signal_homo_sapiens_9606_23032016.tab.gz")
  fasta <- read.csv(file,sep="\t", stringsAsFactors = FALSE)
}

#'
#' load mus musculus signal peptides
#' @export
#' @examples
#' library(prozor)
#' signal <- loadMusMusculusSignalPeptides()
#' head(signal)
loadMusMusculusSignalPeptides <- function(){
  file = file.path(path.package("prozor"),"extdata/uniprot_signal_mus_musculus_29032016.tab.gz")
  fasta <- read.csv(file,sep="\t", stringsAsFactors = FALSE)
}
#'
#' load mus musculus fasta
#' @export
#' @examples
#' library(prozor)
#' fastaMusMusc <- loadMusMusculusFasta()
#'
loadMusMusculusFasta <- function(){
  file = file.path(path.package("prozor"),"extdata/uniprot_mus_musculus_29032016.fasta.gz")
  fasta <- readPeptideFasta(file)
}


