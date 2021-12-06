


.matchPepsequence <- function(peptideSeq,
                              proteinSeq ,
                              proteinID,
                              prefix = "(([RK])|(^)|(^M))",
                              suffix = "")
{
    seqpattern <- paste(prefix, peptideSeq, suffix, sep = "")
    idx2 <- grepl(seqpattern, proteinSeq, fixed = FALSE)
    return(idx2)
}



#' Annotate peptides with protein ids
#'
#' peptides which do not have protein assignment drop out
#'
#' @param pepinfo - list of peptides - sequence, optional modified sequence, charge state.
#' @param fasta - object as created by readPeptideFasta
#' @param peptide - name of column containing peptide sequences default "peptideSeq"
#' @param prefix - default "(([RK])|(^)|(^M))"
#' @param suffix - default ""
#' @import stringr
#' @return data.frame with columns "peptideSeq", "proteinID","Offset","proteinSequence","matched", "lengthPeptide","proteinlength"
#' @export
#' @examples
#'
#' library(dplyr)
#'
#' file = system.file("extdata/IDResults.txt.gz" , package = "prozor")
#' specMeta <- readr::read_tsv(file)
#' upeptide <- unique(specMeta$peptideSeq)
#' resCan <-
#'    prozor::readPeptideFasta(
#'        system.file("p1000_db1_example/Annotation_canSeq.fasta.gz" , package = "prozor"))
#'
#' annotAll = prozor::annotatePeptides(upeptide[seq_len(20)], resCan)
#' dim(annotAll)
#'
#' res <-  mutate(annotAll, proteinlength = nchar(proteinSequence))
#' res <-  select(res, proteinID, peptideSeq, proteinlength, Offset, lengthPeptide)
#' head(res)
annotatePeptides <- function(pepinfo,
                             fasta,
                             peptide = "peptideSeq",
                             prefix = "(([RK])|(^)|(^M))",
                             suffix = "") {
    if (is.null(dim(pepinfo))) {
        pepinfo = data.frame(pepinfo, stringsAsFactors = FALSE)
        colnames(pepinfo) = peptide
    }
    pepinfo <- dplyr::mutate_at(pepinfo, peptide,
                                dplyr::funs(as.character))
    pepinfo <- dplyr::mutate_at(pepinfo, peptide,
                                .funs = dplyr::funs("lengthPeptide" := nchar))

    pepseq  = unique(as.character(pepinfo[, peptide]))
    restab <- annotateAHO(pepseq, fasta)
    if (!is.null(restab)) {

        restab <-
            dplyr::group_by_at(restab, peptide)
        restab <-
            dplyr::mutate(restab,
                matched = .matchPepsequence(
                    dplyr::first(!!dplyr::sym(peptide)),
                    !!dplyr::sym("proteinSequence") ,
                    !!dplyr::sym("proteinID"),
                    prefix = prefix,
                    suffix = suffix
                )
            )
        res = merge(restab, pepinfo, by.x = peptide, by.y = peptide)
        return(res)

    } else {
        warning("no matches found")
        return(NULL)
    }

}


#' annotate peptides using AhoCorasickTrie
#'
#' peptides which do not have protein assignment drop out
#' @param pepseq - list of peptides - sequence, optional modified sequence, charge state.
#' @param fasta - object as created by \code{readPeptideFasta}
#' @import AhoCorasickTrie
#' @import stringr
#' @importFrom purrr map2_df map_df
#' @return A data.frame with proteinID, peptideSeq, Offset and proteinSequence
#' @examples
#'
#' library(dplyr)
#'
#' file = system.file("extdata/IDResults.txt.gz" , package = "prozor")
#' specMeta <- readr::read_tsv(file)
#' upeptide <- unique(specMeta$peptideSeq)
#' resCan <-
#'    prozor::readPeptideFasta(
#'        system.file("p1000_db1_example/Annotation_canSeq.fasta.gz" , package = "prozor"))
#' resCanU <- resCan[!duplicated(unlist(resCan))]
#' annotAll = annotateAHO(upeptide[seq_len(20)], resCanU)
#' dim(annotAll)
#'
#' @export
annotateAHO <- function(pepseq, fasta) {
    #100_000 peptides
    #40_000 Proteine

    pepseq <- stringr::str_trim(unique(pepseq))
    proteinIDS <- names(unlist(fasta))
    fasta <- stringr::str_trim(unlist(fasta))
    names(fasta) <- proteinIDS

    res <-
        AhoCorasickSearch(unique(pepseq) , unlist(fasta), alphabet = "aminoacid")
    if (length(res) == 0)
    {
        return(NULL)
    }
    xx <- purrr::map2_df(
        names(res),
        res,
        .f = function(name, x) {
            data.frame(
                proteinID = name,
                map_df(
                    x,
                    .f = function(x) {
                        x
                    }
                ),
                stringsAsFactors = FALSE
            )
        }
    )
    colnames(xx)[colnames(xx) == "Keyword"] <- "peptideSeq"
    dbframe <- data.frame(
        proteinID = names(fasta),
        proteinSequence = as.character(unlist(fasta)),
        stringsAsFactors = FALSE
    )
    matches <- dplyr::inner_join(xx, dbframe)
    return(matches)
}
