#' create db with decoys and contaminants
#'
#' @param dbs a path to a fasta file or an array of files
#' @param useContaminants add fgcz contaminants
#' @param revLab label for reversed peptides (if NULL do not generate decoys)
#' @param annot source of database
#' @export
#' @examples
#' file = file.path(path.package("prozor"),"extdata/uniprot_taxonomy_Oryctolagus_cuniculus.fasta.gz")
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
#' res <- createDecoyDB(c(file,file), revLab=NULL, useContaminants = FALSE)
#' stopifnot(length(res) == (2*length(rabbit) + 1) )
#'
createDecoyDB <- function(dbs ,
                          useContaminants = TRUE,
                          revLab= "REV_",
                          annot="zz|sourceOf|database")
{
    dummy <-as.SeqFastaAA("CRAPCRAPCRAP", Annot=annot, name= annot)

    dbsfasta <- NULL
    if(useContaminants){
        dbsfasta = loadContaminantsFasta()
    }
    for(db in dbs){
        message( "reading db :" , db )
        dbsfasta <- c(dbsfasta,readPeptideFasta(db) )
    }
    if(!is.null(revLab)){
        dbsfasta <- c(dbsfasta,reverseSeq(dbsfasta ,revLab = revLab))
    }
    return(c(list(dummy),dbsfasta))
}
