#R
.annotateProteinIDGrep <- function(x , fasta, digestPattern ){
    sequence = x$peptideSequence
    idx <- grep (sequence,  fasta, fixed = TRUE)

    pattern = paste(digestPattern, sequence, sep='')
    selected <- fasta[idx]
    idx2 <- grep(pattern, selected, fixed=FALSE)
    idx<-idx[idx2]
    x$proteinInformation = names(fasta)[idx]
    return(x)
}

#' annotate peptides with protein ids
#' @param data - list of peptides - sequence, optional modified sequence, charge state.
#' @param fasta - object as created by read.fasta in pacakge seqinr
#' @param digestPatter - default "(([RK])|(^)|(^M))"
#' @export
#' @examples
#' data(peptides)
#' fasta = read.fasta(file = file, as.string = TRUE, seqtype="AA")
#' annotatePeptides(data, fasta)
annotatePeptides <- function(data,
                             fasta,
                             digestPattern = "(([RK])|(^)|(^M))"
){
    timeStart <- Sys.time();

    if( length(data) > 200 & parallel::detectCores(logical=FALSE) > 1){
        data <- parallel::mclapply(data, .annotateProteinIDGrep, fasta, digestPattern,
                                   mc.cores = min(4,parallel::detectCores(logical=FALSE) ))

    }else{
        data <- lapply(data, .annotateProteinIDGrep, fasta, digestPattern)
    }
    class(data) <- "psmSet"
    timeEnd <- Sys.time();
    message(paste("time taken: ", difftime(timeEnd, timeStart, units='mins'),  "minutes"))
    return(data)
}


