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
#' contNH <- loadContaminantsFasta2021(noHuman = TRUE)
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


#' load list of contaminant sequences FGCZ 2022
#'
#' @export
#' @param noHuman should human contaminants be excluded? default FALSE
#' @return list with contaminant sequences
#' @examples
#' #library(prozor)
#' cont <- loadContaminantsFGCZ2022()
#' length(cont)
#' contNH <- loadContaminantsFGCZ2022(noHuman = TRUE)
#' length(contNH)
#' #example how to create a protein db with decoy sequences
loadContaminantsFGCZ2022 <- function(noHuman = FALSE){
    file = system.file("extdata/fgcz_contaminants2022_20220405.fasta.gz",package = "prozor")
    contaminants <- readPeptideFasta(file)
    if (noHuman) {
        annot <- vapply(contaminants, seqinr::getAnnot, character(1))
        contaminants <- contaminants[!grepl("HUMAN",annot)]
    }
    invisible(contaminants)
}

#' load universal contaminants
#'
#' These sequences are downloaded from here https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics
#'
#' @export
#' @param noHuman should human contaminants be excluded? default FALSE
#' @return list with contaminant sequences
#' @examples
#' #library(prozor)
#' cont <- load_universal_contaminants_2024()
#' length(cont)
#' contNH <- load_universal_contaminants_2024(noHuman = TRUE)
#' length(contNH)
#' #example how to create a protein db with decoy sequences
load_universal_contaminants_2024 <- function(noHuman = FALSE){
  file = system.file("extdata/fgcz_universal_contaminants_github_20241112.fasta.gz",package = "prozor")
  contaminants <- readPeptideFasta(file)
  if (noHuman) {
    annot <- vapply(contaminants, seqinr::getAnnot, character(1))
    contaminants <- contaminants[!grepl("HUMAN",annot)]
  }
  invisible(contaminants)
}


#' load special prot
#'
#' @export
#' @param noHuman should human contaminants be excluded? default FALSE
#' @return list with contaminant sequences
#' @examples
#' #library(prozor)
#' cont <- load_special_prot_2024()
#' length(cont)
#' contNH <- load_special_prot_2024(noHuman = TRUE)
#' length(contNH)
#' #example how to create a protein db with decoy sequences
load_special_prot_2024 <- function(noHuman = FALSE){
  file = system.file("extdata/fgcz_special_prot_20241112.fasta.gz", package = "prozor")
  contaminants <- readPeptideFasta(file)
  if (noHuman) {
    annot <- vapply(contaminants, seqinr::getAnnot, character(1))
    contaminants <- contaminants[!grepl("HUMAN",annot)]
  }
  invisible(contaminants)
}
