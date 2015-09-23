.removeZeroRows <-function(x){
    tmp <- rowSums(x)
    y = x[tmp > 0,, drop=FALSE ]
    return(y)
}

.greedy2 <- function( pepprot ){
    ncolX = ncol(pepprot)
    res<-vector(ncolX , mode="list")
    oldtime <- Sys.time()
    for(i in 1:ncolX)
    {
        if(nrow(pepprot) == 0){
            return(res[1:(i-1)])
        }
        pepsPerProt <- colSums(pepprot)
        if(max(pepsPerProt) == 0){
            return(res[1:(i-1)])
        }
        idx <- which.max(pepsPerProt)
        if(length(idx) > 1){
            idx<-idx[1]
        }
        dele <- pepprot[,idx] > 0
        tmpRes = list(prot = colnames(pepprot)[idx], peps = rownames(pepprot)[dele])
        #message(paste(i, " protein : ", colnames(pepprot)[idx], " matched by peptides: ", sum(dele)))
        pepprot <- pepprot[!dele,-idx,drop=FALSE]
        res[[i]] <- tmpRes
    }
    newtime <- Sys.time()
    message(paste("time ",newtime - oldtime))
    return(res)
}


.greedy <- function( pepprot ){
    ncolX = ncol(pepprot)
    res<-vector(ncolX , mode="list")
    idxx <-NULL
    for(i in 1:ncolX)
    {
        if(i %% 10 == 0){
            pepprot <- .removeZeroRows(pepprot)
            pepprot <- pepprot[,-idxx,drop=FALSE]
            idxx <-NULL
            message(paste("new dim" , dim(pepprot)))
        }
        if(nrow(pepprot) == 0){
            return(res[1:(i-1)])
        }
        oldtime <- Sys.time()
        pepsPerProt <- colSums(pepprot)
        if(max(pepsPerProt) == 0){
            return(res[1:(i-1)])
        }
        idx <- which.max(pepsPerProt)
        if(length(idx) > 1){
            idx<-idx[1]
        }
        idxx <- c(idxx, idx )
        dele <- pepprot[,idx]
        tmpRes = list(prot = colnames(pepprot)[idx], peps = rownames(pepprot)[dele > 0])
        res[[i]] <- tmpRes
        message(paste(i, " protein : ", colnames(pepprot)[idx], " matched by peptides: ", sum(dele)))
        if(sum(dele) > 0){
            set = cbind(rep(dele > 0, ncol(pepprot)))
            pepprot[set] <- 0
        }
        newtime <- Sys.time()
        message(paste("time ",newtime - oldtime))
    }
    return(res)
}
#' given matrix (columns protein rows peptides), compute minimal protein set using greedy algorithm
#' @param pepprot matrix as returned by prepareMatrix
#' @return list of peptide protein assignment
#' @export
#' @examples
#' library(prozor)
#'
#' data(protpepmetashort)
#' head(protpepmetashort)
#' dim(unique(protpepmetashort[,4:5]))
#' xx = prepareMatrix(protpepmetashort, weight= "one")
#' dim(xx)
#' stopifnot(dim(xx)[1] == dim(unique(protpepmetashort[,4:5]))[1])
#' es = greedy(as.matrix(xx))
#' stopifnot(length(unique(names(es))) == dim(unique(protpepmetashort[,4:5]))[1])
#'
greedy <- function( pepprot ){
    protPepAssingments <- .greedy2(pepprot)
    matrixlist <- lapply(protPepAssingments,function(x){ t(cbind(x$peps, rep(x$prot,length(x$peps)))) })
    res = matrix(unlist(matrixlist), ncol=2, byrow = TRUE)
    ltmp = as.list(res[,2])
    names(ltmp) <- res[,1]
    return(ltmp)
}


