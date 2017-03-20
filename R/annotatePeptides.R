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


#' Annotate peptides with protein ids
#'
#' peptides which do not have protein assignment drop out
#' @param pepinfo - list of peptides - sequence, optional modified sequence, charge state.
#' @param fasta - object as created by readPeptideFasta
#' @param digestPattern - default "(([RK])|(^))"
#' @import stringr
#' @export
#' @examples
#' library(prozor)
#' data(pepdata)
#'
#' file = file.path(path.package("prozor"),"extdata/shortfasta.fasta" )
#'
#' fasta = readPeptideFasta(file = file)
#' res = annotatePeptides(pepdata[1:20,], fasta)
#' head(res)
#' res = annotatePeptides(pepdata[1:20,"peptideSequence"],fasta)
#' length(res)
annotatePeptides <- function(pepinfo,
                             fasta,
                             digestPattern = c("R","K","")
){

    if(is.null(dim(pepinfo))){
        pepinfo = matrix(pepinfo,ncol=1)
        colnames(pepinfo) = "peptideSequence"
    }
    pepinfo = pepinfo[,"proteinID" != colnames(pepinfo),drop=FALSE]

    pepinfo = apply(pepinfo,2,as.character)
    lengthPeptide = sapply(pepinfo[,"peptideSequence"],nchar)
    pepinfo = cbind(pepinfo,"lengthPeptide"=lengthPeptide)
    pepseq  = unique(as.character(pepinfo[,"peptideSequence"]))
    restab <- annotateAHO(pepseq, fasta)
    restab <- filterSequences(restab, digestPattern = digestPattern)
    res = merge(restab,pepinfo,by.x="peptideSequence",by.y="peptideSequence")
    res[,"peptideSequence"] <- as.character( res[,"peptideSequence"])
    res[,"proteinID"]<- as.character(res[,"proteinID"])

    return(res)
}

#' annotate vector of petpide sequences against fasta file
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
#' res = annotateVec(pepdata[1:20,"peptideSequence"],fasta)
#' head(res)
#' @export
annotateVec <- function(pepseq, fasta,digestPattern = "(([RK])|(^))",mcCores=NULL ){
    res = .getMatchingProteinIDX(pepseq, fasta,digestPattern,mcCores)
    lengthFasta  = sapply(fasta,nchar)
    namesFasta = names(fasta)
    protLength = vector(length(res),mode="list")
    for(i in 1:length(res)){
        protLength[[i]] =rbind("lengthProtein"=lengthFasta[res[[i]]],"proteinID"=namesFasta[res[[i]]],"peptideSequence"=names(res)[i])
    }

    checkdim <- sapply(protLength, function(x){dim(x)[1]})
    which2remove <- which(checkdim == 1)
    if( length(which2remove) > 0 ){
        protLength <- protLength[-which2remove]
    }
    restab = matrix(unlist(protLength),ncol=3,byrow=TRUE)
    colnames(restab) = c("lengthProtein","proteinID","peptideSequence")
    return(restab)
}

#'
#' annotate peptides using AhoCorasickTrie
#'
#' peptides which do not have protein assignment drop out
#' @param pepseq - list of peptides - sequence, optional modified sequence, charge state.
#' @param fasta - object as created by readPeptideFasta
#' @import AhoCorasickTrie
#' @import stringr
#' @examples
#'
#' library(prozor)
#' file = file.path(path.package("prozor"),"extdata/shortfasta.fasta" )
#' fasta = readPeptideFasta(file = file)
#' #res = annotateVec(pepdata[1:20,"peptideSequence"],fasta)
#' system.time(res2 <- annotateAHO(pepdata[1:20,"peptideSequence"],fasta))
#' colnames(res2)
#' @export
annotateAHO <- function(pepseq,fasta){
    #100_000 peptides
    #40_000 Proteine

    pepseq <-stringr::str_trim(unique(pepseq))

    proteinIDS <- names(unlist(fasta))
    fasta <- stringr::str_trim(unlist(fasta))
    names(fasta) <- proteinIDS

    system.time(res <- AhoCorasickSearch(unique(pepseq) , unlist(fasta), alphabet = "aminoacid"))

    simplifyAhoCorasickResult <- function(x, name){t <- as.data.frame(do.call("rbind",(x))); t$proteinID <- name; return(t)}
    tmp <- mapply(simplifyAhoCorasickResult, res, names(res), SIMPLIFY=FALSE)

    xx <- plyr::rbind.fill(tmp)
    colnames(xx)[colnames(xx)=="Keyword"]<-"peptideSequence"
    xx$peptideSequence <- as.character(xx$peptideSequence)
    xx$Offset <- as.numeric(xx$Offset)
    dbframe <- data.frame(proteinID = names(fasta), proteinSequence = as.character(unlist(fasta)),stringsAsFactors = FALSE)
    matches <- merge(xx, dbframe )
    return(matches)
}
#' Filter for specific residues
#'
#' Will check if AA at Offset is a valid cleavage site
#'
#' @param matches must have 2 columns proteinSequnce and Offset
#' @param digestPattern - list of N terminal amino acids including empty string (protein start) default tryptic = c("","K","R")
#' @export
#'
filterSequences <- function(matches,digestPattern = c("","K","R") ){
    matches$predcessor <- substr(matches$proteinSequence,matches$Offset-1 , matches$Offset-1 )
    finmat <- matches[matches$predcessor %in% digestPattern,]
}

