#' Create db with decoys and contaminants
#'
#' For more details and references see package vignette
#' \code{vignette("CreateDecoyDB", package = "prozor")}
#'
#' @param dbs a path to a fasta file or an array of files
#' @param useContaminants list with contaminant sequences
#' @param revLab label for reversed peptides (if NULL do not generate decoys)
#' @param annot source of database
#' @export
#' @examples
#' #file = file.path(path.package("prozor"),"extdata/shortfasta.fasta.gz")
#' file = system.file("extdata/fgcz_contaminants_20150123.fasta.gz",package = "prozor")
#' cont <- loadContaminantsFasta()
#' rabbit <-readPeptideFasta(file)
#' tmp <- 2*(2*length(rabbit)+length(cont)) + 1
#'
#' res <- createDecoyDB(c(file,file))
#' length(res)
#' tmp
#' stopifnot(length(res) == tmp)
#'
#' res <- createDecoyDB(c(file,file), revLab=NULL)
#' stopifnot(length(res) == (2*length(rabbit)+length(cont) + 1))
#' res <- createDecoyDB(c(file,file), revLab=NULL, useContaminants = NULL)
#' stopifnot(length(res) == (2*length(rabbit) + 1) )
#'
createDecoyDB <- function(dbs ,
                          useContaminants = loadContaminantsFasta(),
                          revLab= "REV_",
                          annot="zz|sourceOf|database")
{
    dummy <-as.SeqFastaAA("CRAPCRAPCRAP", Annot=annot, name= annot)
    dbsfasta <- NULL
    for(db in dbs){
        message( "reading db :" , db )
        dbsfasta <- c(dbsfasta,readPeptideFasta(db) )
    }
    if(!is.null(useContaminants)){
        dbsfasta = c(dbsfasta,useContaminants)
    }
    if(!is.null(revLab)){
        dbsfasta <- c(dbsfasta,reverseSeq(dbsfasta ,revLab = revLab))
    }
    return(c(list(dummy),dbsfasta))
}
