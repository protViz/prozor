#' load list of contaminant sequences
#'
#' @export
#' @examples
#' #library(prozor)
#' cont <- loadContaminantsFasta()
#' cont[[1]]
#' #example how to create a protein db with decoy sequences
loadContaminantsFasta <- function(){
  file = system.file("extdata/fgcz_ContaminantsWithAnnotation.fasta.gz",package = "prozor")
  #file = file.path(path.package("prozor"),"extdata/fgcz_ContaminantsWithAnnotation.fasta.gz")
  contaminants <- readPeptideFasta(file)
}
#' load list of contaminant without human sequences
#'
#' @export
#' @examples
#' #library(prozor)
#' cont <- loadContaminantsNoHumanFasta()
#' cont[[1]]
#' #example how to create a protein db with decoy sequences
loadContaminantsNoHumanFasta <- function(){
    file = system.file("extdata/fgcz_ContaminantsWithAnnotationNoHuman.fasta.gz",package = "prozor")
    #file = file.path(path.package("prozor"),"extdata/fgcz_ContaminantsWithAnnotationNoHuman.fasta.gz")
    contaminants <- readPeptideFasta(file)
    invisible(contaminants)
}

#' load list of contaminant sequences FGCZ 2019
#'
#' @export
#' @param noHuman should human contaminants be excluded? default FALSE
#' @examples
#' #library(prozor)
#' cont <- loadContaminantsFasta2019()
#' length(cont)
#' contNH <- loadContaminantsFasta2019()
#' length(contNH)
#' #example how to create a protein db with decoy sequences
loadContaminantsFasta2019 <- function(noHuman = FALSE){
  file = system.file("extdata/fgcz_contaminants2019_20190708.fasta.gz",package = "prozor")
  #file = file.path(path.package("prozor"),"extdata/fgcz_ContaminantsWithAnnotation.fasta.gz")
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
