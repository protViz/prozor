#' load list of contaminant sequences
#'
#' @export
#' @examples
#' library(prozor)
#' cont <- loadContaminantsFasta()
#' cont[[1]]
#' #example how to create a protein db with decoy sequences
loadContaminantsFasta <- function(){
  file = file.path(path.package("prozor"),"extdata/fgcz_ContaminantsWithAnnotation.fasta")
  contaminants <- readPeptideFasta(file)
}
#' load list of contaminant without human sequences
#'
#' @export
#' @examples
#' library(prozor)
#' cont <- loadContaminantsNoHumanFasta()
#' cont[[1]]
#' #example how to create a protein db with decoy sequences
loadContaminantsNoHumanFasta <- function(){
    file = file.path(path.package("prozor"),"extdata/fgcz_ContaminantsWithAnnotationNoHuman.fasta")
    contaminants <- readPeptideFasta(file)
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
  fasta <- utils::read.csv(file,sep="\t", stringsAsFactors = FALSE)
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
  fasta <- utils::read.csv(file,sep="\t", stringsAsFactors = FALSE)
}

