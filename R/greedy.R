.removeZeroRows <- function(x){
    tmp <- rowSums(x)
    y = x[tmp > 0,, drop = FALSE ]
    return(y)
}

#'
if(FALSE){
    data(protpepmetashort)
    protpepmetashort
    xx = prepareMatrix(protpepmetashort, peptideID = "peptideModSeq")
    debug(prozor:::.greedy2)
    prozor:::.greedy2(xx)
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
        # if there is tie you need to resolve it.
        # if peptides are disjoint, return arbitrary since it will win in the next round .
        # if peptides are not disjoint:
        # - either they share all peptides then merge/group and remove
        # - or there are overlapping N and O disjoint peptides then sample winner and warn.
        if (length(idx) > 1) {
            # check if disjoint
            mm <- pepprot[,idx]
            find_overlap <- function(mm){
                ii_x <- mdist(t(as.matrix(mm)), .func = prozor:::inter_x)
                ii_y <- mdist(t(as.matrix(mm)), .func = prozor:::inter_y)
                if (!all(na.omit(ii_x == ii_y))) { warning("ii_x and ii_y are not identical.") }
                cS <- mean(colSums(mm))
                ib <- which( ii_x == cS, arr.ind = TRUE)
                if (nrow(ib) > 1) {
                    return(unique(as.integer(ib)))
                }
                return(1)
            }
            idx <- idx[find_overlap(mm)]
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
#' #library(prozor)
#'
#' data(protpepmetashort)
#' colnames(protpepmetashort)
#' dim(unique(protpepmetashort[,4]))
#' xx = prepareMatrix(protpepmetashort, peptideID = "peptideModSeq")
#' dim(xx)
#' stopifnot(dim(xx)[1] == dim(unique(protpepmetashort[,4]))[1])
#'
#' es = greedy_parsimony(xx)
#' debug(prozor:::.greedy2)
#' stopifnot(length(unique(names(es))) == dim(unique(protpepmetashort[,4]))[1])
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
#' @return matrix of peptide protein assignments
#' @param res result of function prozor::greedy_parsimony
#' @export
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

