#' load list of contaminant sequences FGCZ 2019
#'
#' @export
#' @param noHuman should human contaminants be excluded? default FALSE
#' @return list with contaminant sequences
#' @examples
#' #library(prozor)
#' cont <- loadContaminantsFasta2019()
#' length(cont)
#' contNH <- loadContaminantsFasta2019()
#' length(contNH)
#' #example how to create a protein db with decoy sequences
loadContaminantsFasta2019 <- function(noHuman = FALSE){
  file = system.file("extdata/fgcz_contaminants2019_20190708.fasta.gz",package = "prozor")
  contaminants <- readPeptideFasta(file)
  if (noHuman) {
    annot <- vapply(contaminants, seqinr::getAnnot, character(1))
    contaminants <- contaminants[!grepl("HUMAN",annot)]
  }
  invisible(contaminants)
}


#' load list of contaminant sequences FGCZ 2021
#'
#' @export
#' @param noHuman should human contaminants be excluded? default FALSE
#' @return list with contaminant sequences
#' @examples
#' #library(prozor)
#' cont <- loadContaminantsFasta2021()
#' length(cont)
#' contNH <- loadContaminantsFasta2021()
#' length(contNH)
#' #example how to create a protein db with decoy sequences
loadContaminantsFasta2021 <- function(noHuman = FALSE){
  file = system.file("extdata/fgcz_contaminants2021_20210929.fasta.gz",package = "prozor")
  contaminants <- readPeptideFasta(file)
  if (noHuman) {
    annot <- vapply(contaminants, seqinr::getAnnot,character(1))
    contaminants <- contaminants[!grepl("HUMAN",annot)]
  }
  invisible(contaminants)
}
