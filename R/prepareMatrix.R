.makePeptideID <- function(peprot){
    tmp = peprot[,c("peptideModSequence","z")]
    tmp = apply(tmp,2,as.character)
    apply(tmp, 1, paste, collapse = '.')
}
#' given table of peptide protein assigments generate matrix
#' @param data  column 1 proteins, column 2 peptides
#' @value sparse matrix
#' @export
prepareMatrix <- function(data, weight = c("count", "AA", "coverage")){
    fprots = as.factor(data[,"proteinID"])
    prots = as.integer(fprots)
    fpeps = as.factor(.makePeptideID(data))
    peps = as.integer(fpeps)

    pepProt =sparseMatrix(peps , prots,x = 1 )

    colnames(pepProt) <- levels(fprots)
    rownames(pepProt) <- levels(fpeps)

    dim(pepProt)
    return(pepProt)
}

