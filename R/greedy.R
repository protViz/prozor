#' Given matrix
#'
#'
occam <- function(pepprot, ncolX = ncol(pepprot)){
    res<-vector(ncolX , mode="list")
    idxx <-NULL
    for(i in 1:ncolX)
    {

        if(i %% 10 == 0){
            pepprot <- removeZeroRows(pepprot)
            drumm <<-pepprot
            idxxx <<- idxx
            cat("length(idxx)" , length(idxx), "\n")
            pepprot <- pepprot[,-idxx,drop=FALSE]

            idxx <-NULL
            cat("new dim" , dim(pepprot), "\n")
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
        cat(i, " ", idx, " ", sum(dele), "\n")
        if(sum(dele) > 0){
            set = cbind(rep(dele > 0, ncol(pepprot)))
            #for(j in 1:ncol(pepprot)){
            pepprot[set] <- 0
            #}
        }
        newtime <- Sys.time()
        cat("time ",newtime - oldtime, "\n")

    }
    return(list(res = res))
}
