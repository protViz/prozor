#' given table of peptide protein assigments generate matrix
#' @param data  column 1 proteins, column 2 peptides
#' @value sparse matrix
#' @export

prepareMatrix <- function(data){
    fprots = as.factor(data[,"prots"])
    prots = as.integer(fprots)
    fpeps = as.factor(data[,"peps"])
    peps = as.integer(fpeps)


    pepProt =sparseMatrix(peps, prots,x = 1 )
    dim(pepProt)

    sum(duplicated(levels(fpeps)))
    sum(duplicated(levels(fprots)))

    colnames(pepProt) <- levels(fprots)
    rownames(pepProt) <- levels(fpeps)

    dim(pepProt)
    return(pepProt)
}
