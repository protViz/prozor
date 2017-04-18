.annotateProteinIDGrep <- function(x , fasta, digestPattern="(([RK])|(^))"){
    sequence = x
    idx <- grep (sequence,  fasta, fixed = TRUE)
    if(length(idx) > 1){
        pattern = paste(digestPattern, sequence, sep='')
        selected <- fasta[idx]
        idx2 <- grep(pattern, selected, fixed=FALSE)
        idx<-idx[idx2]
    }
    return(idx)
}


.getMatchingProteinIDX <- function(data,
                                   fasta,
                                   digestPattern = "(([RK])|(^))",
                                   mcCores=NULL
){
    timeStart <- Sys.time();
    if(is.null(mcCores)){
        mcCores <- min(6,parallel::detectCores(logical=FALSE))
    }
    if( length(data) > 100 & mcCores > 1){
        message(paste("going to use : " , mcCores ," cores."))
        registerDoParallel(mcCores)
        res <- foreach(i = data ) %dopar% .annotateProteinIDGrep(i, fasta, digestPattern)
        stopImplicitCluster()
    }else{
        res <- lapply(data, .annotateProteinIDGrep, fasta, digestPattern)
    }
    names(res) = data
    timeEnd <- Sys.time();
    message(paste("time taken: ", difftime(timeEnd, timeStart, units='mins'),  "minutes"))
    return(res)
}
#' annotate vector of petpide sequences against fasta file (Deprecated)
#'
#' @param pepseq peptide sequences
#' @param fasta fasta file
#' @param digestPattern digest pattern as regex
#' @param mcCores nr of cores to use
#' @examples
#'
#' library(prozor)
#' file = file.path(path.package("prozor"),"extdata/shortfasta.fasta" )
#' fasta = readPeptideFasta(file = file)
#'
#' res = annotateVec(pepprot[1:20,"peptideSeq"],fasta)
#' head(res)
#' @export
annotateVec <- function(pepseq, fasta,digestPattern = "(([RK])|(^)|(^M))",mcCores=NULL ){
    res = .getMatchingProteinIDX(pepseq, fasta,digestPattern,mcCores)
    lengthFasta  = sapply(fasta,nchar)
    namesFasta = names(fasta)
    protLength = vector(length(res),mode="list")
    for(i in 1:length(res)){
        protLength[[i]] =rbind("lengthProtein"=lengthFasta[res[[i]]],
                               "proteinID"=namesFasta[res[[i]]],
                               "peptideSeq"=names(res)[i])
    }

    checkdim <- sapply(protLength, function(x){dim(x)[1]})
    which2remove <- which(checkdim == 1)
    if( length(which2remove) > 0 ){
        protLength <- protLength[-which2remove]
    }
    restab = matrix(unlist(protLength),ncol=3,byrow=TRUE)
    colnames(restab) = c("lengthProtein","proteinID","peptideSeq")
    return(restab)
}
