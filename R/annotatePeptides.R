#' Annotate peptides with protein ids
#'
#' peptides which do not have protein assignment drop out
#'
#' @param pepinfo - list of peptides - sequence, optional modified sequence, charge state.
#' @param fasta - object as created by readPeptideFasta
#' @param prefix - default "(([RK])|(^)|(^M))"
#' @param suffix - default ""
#' @import stringr
#' @export
#' @examples
#' library(prozor)
#' data(pepprot)
#' file = system.file("extdata/shortfasta.fasta.gz",package = "prozor")
#'
#' fasta = readPeptideFasta(file = file)
#' res = annotatePeptides(pepprot[1:20,], fasta)
#' head(res)
#' res = annotatePeptides(pepprot[1:20,"peptideSeq"],fasta)
#' length(res)
annotatePeptides <- function(pepinfo,
                             fasta,
                             prefix = "(([RK])|(^)|(^M))",
                             suffix = ""

){

    if(is.null(dim(pepinfo))){
        pepinfo = matrix(pepinfo,ncol=1)
        colnames(pepinfo) = "peptideSeq"
    }
    pepinfo = pepinfo[,"proteinID" != colnames(pepinfo),drop=FALSE]

    pepinfo = apply(pepinfo,2,as.character)
    lengthPeptide = sapply(pepinfo[,"peptideSeq"],nchar)
    pepinfo = cbind(pepinfo,"lengthPeptide"=lengthPeptide)
    pepseq  = unique(as.character(pepinfo[,"peptideSeq"]))
    restab <- annotateAHO(pepseq, fasta)
    restab <- filterSequences(restab, prefix = prefix, suffix = suffix)
    res = merge(restab,pepinfo,by.x="peptideSeq",by.y="peptideSeq")
    res[,"peptideSeq"] <- as.character( res[,"peptideSeq"])
    res[,"proteinID"]<- as.character(res[,"proteinID"])

    return(res)
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
#' library(AhoCorasickTrie)
#' file = system.file("extdata/shortfasta.fasta.gz",package = "prozor")
#' fasta = readPeptideFasta(file = file)
#' pepprot <- get(data("pepprot", package = "prozor"))
#' system.time( res2 <- annotateAHO( pepprot[1:20,"peptideSeq"], fasta))
#' colnames(res2)
#'
#' @export
annotateAHO <- function(pepseq,fasta){
    #100_000 peptides
    #40_000 Proteine

    pepseq <-stringr::str_trim(unique(pepseq))
    proteinIDS <- names(unlist(fasta))
    fasta <- stringr::str_trim(unlist(fasta))
    names(fasta) <- proteinIDS

    res <- AhoCorasickSearch(unique(pepseq) , unlist(fasta), alphabet = "aminoacid")
    if(length(res) == 0)
    {
        return(NULL)
    }
    xx <- purrr::map2_df(names(res),res,
                          .f=function(name,x){data.frame(proteinID=name,
                                                         map_df(x,.f=function(x){x}), stringsAsFactors = FALSE)
                          })
    colnames(xx)[colnames(xx)=="Keyword"]<-"peptideSeq"
    dbframe <- data.frame(proteinID = names(fasta), proteinSequence = as.character(unlist(fasta)),stringsAsFactors = FALSE)
    matches <- dplyr::inner_join(xx, dbframe )

    return(matches)
}

.matchPepsequence <- function(matches, prefix= "(([RK])|(^)|(^M))", suffix =""){

    seqpattern <-paste(prefix, matches$peptideSeq[1], suffix, sep="")
    idx2 <- grep(seqpattern, matches$proteinSequence, fixed=FALSE)
    if(length(idx2) > 0){
        matchesres <- matches[idx2,]
        matchesres$pattern <- seqpattern
        return(matchesres)
    }else{
        return(NULL)
    }

}

#'
#' Filter for specific residues
#'
#' Will check if AA at Offset is a valid cleavage site
#'
#' @param matches must have 2 columns proteinSequnce and Offset
#' @param prefix - regular expression describing the prefix of the peptide sequence e.g. (([RK])|(^)|(^M))
#' @param suffix - regular expression describing the suffix of the peptide sequence
#' @export
#'
filterSequences <- function(matches,prefix = "(([RK])|(^)|(^M))", suffix="" ){
    x <- plyr::ddply(matches, ~peptideSeq, .matchPepsequence, prefix = prefix, suffix = suffix)
}


