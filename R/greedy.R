.removeZeroRows <- function(x){
    tmp <- rowSums(x)
    y = x[tmp > 0,, drop = FALSE ]
    return(y)
}

.greedy2 <- function(pepprot){
    ncolX = ncol(pepprot)
    res <- vector(ncolX , mode = "list")
    oldtime <- Sys.time()
    for (i in seq_len(ncolX))
    {
        if (nrow(pepprot) == 0) {
            return(res[seq_len(i - 1)])
        }
        pepsPerProt <- colSums(pepprot)
        if (max(pepsPerProt) == 0) {
            return(res[seq_len(i - 1)])
        }
        mymax <- function(x){ which(max(x) == x) }
        idx <- mymax(pepsPerProt)
        # Tie: proteins explaining the same peptides as the first tied protein are
        # indistinguishable and reported together; the other tied proteins compete
        # in the next rounds.
        if (length(idx) > 1) {
            first <- pepprot[, idx[1]] > 0
            same <- vapply(idx, function(j) identical(pepprot[, j] > 0, first), logical(1))
            idx <- idx[same]
        }
        dele <- pepprot[,idx[1]] > 0
        tmpRes = list(prot = paste(colnames(pepprot)[idx] , collapse = ";"), peps = rownames(pepprot)[dele])
        pepprot <- pepprot[!dele,-idx,drop = FALSE]
        res[[i]] <- tmpRes
    }
    newtime <- Sys.time()
    message("time : ", newtime - oldtime)
    return(res)
}



#' given matrix (columns protein rows peptides), compute minimal protein set using greedy algorithm
#' @param pepprot matrix as returned by prepareMatrix
#' @return list of peptide protein assignment
#' @export
#' @examples
#' data(protpepmetashort)
#' colnames(protpepmetashort)
#' xx = prepareMatrix(protpepmetashort, peptideID = "peptideModSeq")
#' dim(xx)
#' stopifnot(nrow(xx) == length(unique(protpepmetashort$peptideModSeq)))
#'
#' es = greedy_parsimony(xx)
#' # every precursor is explained by exactly one protein (group)
#' stopifnot(setequal(names(es), rownames(xx)))
#'
greedy_parsimony <- function(pepprot) {
    protPepAssingments <- .greedy2(pepprot)
    matrixlist <- lapply(protPepAssingments,function(x){ t(cbind(x$peps, rep(x$prot,length(x$peps)))) })
    res = matrix(unlist(matrixlist), ncol = 2, byrow = TRUE)
    ltmp = as.list(res[,2])
    names(ltmp) <- res[,1]
    return(ltmp)
}
#' converts result of greedy_parsimony function to a matrix with 3 columns - peptide - charge and protein
#' @return data.frame of peptide protein assignments with columns Peptide, precursorCharge and Protein
#' @param res result of function prozor::greedy_parsimony; its names must have the form peptide.charge
#' @export
#' @examples
#' data(protpepmetashort)
#' protpepmetashort$precursor <- paste(protpepmetashort$peptideSeq,
#'                                     protpepmetashort$precursorCharge, sep = ".")
#' xx <- prepareMatrix(protpepmetashort, peptideID = "precursor")
#' es <- greedy_parsimony(xx)
#' head(greedyRes2Matrix(es))
#'
greedyRes2Matrix <- function(res){
    res <- (cbind(names(res),unlist(res)))
    cnamessplit <- strsplit(as.character(res[,1]),split = "\\.")
    protnam <- do.call("rbind",cnamessplit)
    res <- cbind(protnam,res[,2])
    colnames(res) <- c("Peptide", "precursorCharge", "Protein")
    res <- as.data.frame(res)
    return(res)
}

