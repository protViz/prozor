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

.getMatchingProteinIDX <- function(data,
                                   fasta,
                                   digestPattern = "(([RK])|(^)|(^M))"
){
    timeStart <- Sys.time();
    if( length(data) > 200 & parallel::detectCores(logical=FALSE) > 1){
        mcCores <- min(6,parallel::detectCores(logical=FALSE))
        message(paste("going to use : " , mcCores ," cores."))
        registerDoParallel(mcCores)
        res <- foreach(i=data ) %dopar% protrazor:::.annotateProteinIDGrep(i, fasta, digestPattern)
        stopImplicitCluster()
    }else{
        res <- lapply(data, .annotateProteinIDGrep, fasta, digestPattern)
    }
    names(res) = data
    timeEnd <- Sys.time();
    message(paste("time taken: ", difftime(timeEnd, timeStart, units='mins'),  "minutes"))
    return(res)
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
#' head(pepprot)
#' file = file.path(path.package("protrazor"),"extdata/fgcz_10090_20140715.fasta" )
#' fasta = read.fasta(file = file, as.string = TRUE, seqtype="AA")
#' res = annotatePeptides(pepprot, fasta)
#' head(res2)
#'
annotatePeptides <- function(pepinfo,
                                fasta,
                                digestPattern = "(([RK])|(^)|(^M))"
){
    pepinfo = apply(pepinfo,2,as.character)
    lengthPeptide = sapply(pepinfo[,"peptideSequence"],nchar)
    pepinfo = cbind(pepinfo,"lengthPeptide"=lengthPeptide)

    pepseq  = unique(as.character(pepinfo[,"peptideSequence"]))
    res = protrazor:::.getMatchingProteinIDX(pepseq, fasta,digestPattern)
    lengthFasta  = sapply(fasta,nchar)
    namesFasta = names(fasta)
    protLength = vector(length(res),mode="list")
    for(i in 1:length(res)){
        protLength[[i]] =rbind("lengthProtein"=lengthFasta[res[[i]]],"proteinID"=namesFasta[res[[i]]],"peptideSequence"=names(res)[i])
    }
    restab = matrix(unlist(protLength),ncol=3,byrow=TRUE)
    colnames(restab) = c("lengthProtein","proteinSequence","peptideSequence")
    restab <<- restab
    pepinfo <<- pepinfo
    res = merge(restab,pepinfo,by.x="peptideSequence",by.y="peptideSequence")
    return(res)
}


