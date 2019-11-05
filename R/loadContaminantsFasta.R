#' load list of contaminant sequences
#'
#' @export
#' @examples
#' library(prozor)
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
#' library(prozor)
#' cont <- loadContaminantsNoHumanFasta()
#' cont[[1]]
#' #example how to create a protein db with decoy sequences
loadContaminantsNoHumanFasta <- function(){
    file = system.file("extdata/fgcz_ContaminantsWithAnnotationNoHuman.fasta.gz",package = "prozor")
    #file = file.path(path.package("prozor"),"extdata/fgcz_ContaminantsWithAnnotationNoHuman.fasta.gz")
    contaminants <- readPeptideFasta(file)
}

#' load list of contaminant sequences FGCZ 2019
#'
#' @export
#' @examples
#' library(prozor)
#' cont <- loadContaminantsFasta2019()
#' cont[[1]]
#' #example how to create a protein db with decoy sequences
loadContaminantsFasta2019 <- function(){
  file = system.file("extdata/fgcz_contaminants2019_20190708.fasta.gz",package = "prozor")
  #file = file.path(path.package("prozor"),"extdata/fgcz_ContaminantsWithAnnotation.fasta.gz")
  contaminants <- readPeptideFasta(file)
}
