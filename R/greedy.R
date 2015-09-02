.removeZeroRows <-function(x){
    tmp <- rowSums(x)
    cat("sum" , sum(tmp == 0) )
    y = x[tmp > 0,, drop=FALSE ]
    return(y)
}
#' given matrix
#' @param pepprot matrix as returned by prepareMatrix
#' @export
#' @examples
#' library(prozor)
#'
#' data(protpepmeta)
#' xx = prepareMatrix(protpepmeta, weight= "count")
#' dim(xx)
#' es = occam(as.matrix(xx))
occam <- function(pepprot ){
    ncolX = ncol(pepprot)
    res<-vector(ncolX , mode="list")
    idxx <-NULL
    for(i in 1:ncolX)
    {
        if(i %% 10 == 0){
            pepprot <- .removeZeroRows(pepprot)
            #drumm <<-pepprot
            #idxxx <<- idxx
            message(paste("length(idxx)" , length(idxx),sep= " "))
            pepprot <- pepprot[,-idxx,drop=FALSE]

            idxx <-NULL
            message(paste("new dim" , dim(pepprot)))
        }
        if(nrow(pepprot) == 0)
            break()
        oldtime <- Sys.time()
        pepsPerProt <- colSums(pepprot)
        idx <- which.max(pepsPerProt)
        if(length(idx) > 1){
            idx<-idx[1]
        }
        idxx <- c( idxx, idx )
        dele <- pepprot[,idx]
        tmpRes = list(prot = colnames(pepprot)[idx], peps = rownames(pepprot)[dele>0])
        res[[i]] <- tmpRes
        message(paste(i, " ", idx, " ", sum(dele)))
        if(sum(dele) > 0){
            set = cbind(rep(dele > 0, ncol(pepprot)))
            pepprot[set] <- 0
        }
        newtime <- Sys.time()
        message(paste("time ",newtime - oldtime))

    }
    return(list(res = res))
}
