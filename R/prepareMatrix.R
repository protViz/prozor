.makePeptideID <- function(peprot){
    peprot[,c("peptideModSequence","z")]

}
#' given table of peptide protein assigments generate matrix
#' @param data  column 1 proteins, column 2 peptides
#' @value sparse matrix
#' @export
prepareMatrix <- function(data){

    fprots = as.factor(data[,"proteinID"])
    prots = as.integer(fprots)
    fpeps = as.factor(data[,"peptideSequence"])
    peps = as.integer(fpeps)

    pepProt =sparseMatrix(peps , prots,x = 1 )

    colnames(pepProt) <- levels(fprots)
    rownames(pepProt) <- levels(fpeps)

    dim(pepProt)
    return(pepProt)
}
