#R
.annotateProteinIDGrep <- function(x , fasta, digestPattern="(([RK])|(^)|(^M))"){
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

#' annotate peptides with protein ids
#' @param data - list of peptides - sequence, optional modified sequence, charge state.
#' @param fasta - object as created by read.fasta in pacakge seqinr
#' @param digestPatter - default "(([RK])|(^)|(^M))"
#' @export
#' @examples
#' library(doParallel)
#' library(foreach)
#' data(pepprot)
#' dim(pepprot)
#' pepprot<-apply(pepprot,2,as.character)
#' pepprot <- pepprot[!is.na(pepprot[,1]),]
#' pepprot <- unique(pepprot[,2:5])
#' dim(pepprot)
#' library(seqinr)
#' file = file.path(path.package("protrazor"),"extdata/fgcz_10090_20140715.fasta" )
#' fasta = read.fasta(file = file, as.string = TRUE, seqtype="AA")
#' pepseq  = unique(as.character(pepprot[,"peptideSequence"]))
#' res = annotatePeptides(pepseq, fasta)
#' length(res)
#' length(pepseq)
#' restab = NULL
#' for(i in res){
#'  restab = rbind(restab,i)
#' }
#' dim(pepprot)
#' res2 = merge(restab,pepprot,by.x="peptideSequence",by.y="peptideSequence")
#' dim(res2)
#'
#'
annotatePeptides <- function(data,
                             fasta,
                             digestPattern = "(([RK])|(^)|(^M))"
){
    timeStart <- Sys.time();
    if( length(data) > 200 & parallel::detectCores(logical=FALSE) > 1){
        mcCores <- min(6,parallel::detectCores(logical=FALSE))
        message(paste("going to use : " , mcCores ," cores."))
        registerDoParallel(mcCores)
        data <- foreach(i=data ) %dopar% protrazor:::.annotateProteinIDGrep(i, fasta, digestPattern)
        stopImplicitCluster()
    }else{
        data <- lapply(data, .annotateProteinIDGrep, fasta, digestPattern)
    }
    timeEnd <- Sys.time();
    message(paste("time taken: ", difftime(timeEnd, timeStart, units='mins'),  "minutes"))
    return(data)
}


