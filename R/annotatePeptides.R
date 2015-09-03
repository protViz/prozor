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
                                   digestPattern = "(([RK])|(^)|(^M))",mcCores=NULL
){
    timeStart <- Sys.time();
    if( length(data) > 100 & parallel::detectCores(logical=FALSE) > 1){
        if(is.null(mcCores)){
            mcCores <- min(6,parallel::detectCores(logical=FALSE))
        }
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


#' annotate peptides with protein ids
#' @param pepinfo - list of peptides - sequence, optional modified sequence, charge state.
#' @param fasta - object as created by read.fasta in pacakge seqinr
#' @param digestPattern - default "(([RK])|(^)|(^M))"
#' @param mcCores number of cores to use
#' @export
#' @examples
#' library(prozor)
#' library(doParallel)
#' library(foreach)
#' library(seqinr)
#' data(pepdata)
#' head(pepdata)
#'
#' file = file.path(path.package("prozor"),"extdata/shortfasta.fasta" )
#' fasta = read.fasta(file = file, as.string = TRUE, seqtype="AA")
#' # we use a subset of the data to speedup the computation
#' #res = annotatePeptides(pepdata, fasta)
#' res = annotatePeptides(pepdata[1:20,], fasta,mcCores=1)
#' head(res)
#'

annotatePeptides <- function(pepinfo,
                                fasta,
                                digestPattern = "(([RK])|(^)|(^M))",mcCores=NULL
){
    pepinfo = apply(pepinfo,2,as.character)
    lengthPeptide = sapply(pepinfo[,"peptideSequence"],nchar)
    pepinfo = cbind(pepinfo,"lengthPeptide"=lengthPeptide)

    pepseq  = unique(as.character(pepinfo[,"peptideSequence"]))
    res = .getMatchingProteinIDX(pepseq, fasta,digestPattern,mcCores)
    lengthFasta  = sapply(fasta,nchar)
    namesFasta = names(fasta)
    protLength = vector(length(res),mode="list")
    for(i in 1:length(res)){
        protLength[[i]] =rbind("lengthProtein"=lengthFasta[res[[i]]],"proteinID"=namesFasta[res[[i]]],"peptideSequence"=names(res)[i])
    }
    restab = matrix(unlist(protLength),ncol=3,byrow=TRUE)
    colnames(restab) = c("lengthProtein","proteinID","peptideSequence")
    restab <<- restab
    pepinfo <<- pepinfo
    res = merge(restab,pepinfo,by.x="peptideSequence",by.y="peptideSequence")
    return(res)
}


