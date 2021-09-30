

.matchPepsequence <- function(
    peptideSeq,
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
#' @export
#' @examples
#'
#' library(prozor)
#' library(dplyr)
#' file = system.file("extdata/shortfasta.fasta.gz",package = "prozor")
#'
#' fasta = readPeptideFasta(file = file)
#' res = annotatePeptides(pepprot[1:20,], fasta)
#' res = annotatePeptides(pepprot[1:20,"peptideSeq"],fasta)
#' str(res)
#' res %>% mutate(proteinlength = nchar(proteinSequence)) -> res
#' res %>% select(proteinID, peptideSeq, proteinlength, Offset, lengthPeptide)
#'
annotatePeptides <- function(pepinfo,
                             fasta,
                             peptide = "peptideSeq",
                             prefix = "(([RK])|(^)|(^M))",
                             suffix = "") {
    if (is.null(dim(pepinfo))) {
        pepinfo = data.frame(pepinfo, stringsAsFactors = FALSE)
        colnames(pepinfo) = peptide
    }
    pepinfo <-
        pepinfo %>% dplyr::mutate_at(peptide, dplyr::funs(as.character)) %>%
        dplyr::mutate_at(peptide, .funs = dplyr::funs("lengthPeptide" := nchar))
    pepseq  = unique(as.character(pepinfo[, peptide]))
    restab <- annotateAHO(pepseq, fasta)

    restab <- restab %>%
        dplyr::group_by_at(peptide) %>%
        dplyr::mutate(
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
}

#'
#' annotate peptides using AhoCorasickTrie
#'
#' peptides which do not have protein assignment drop out
#' @param pepseq - list of peptides - sequence, optional modified sequence, charge state.
#' @param fasta - object as created by readPeptideFasta
#' @import AhoCorasickTrie
#' @import stringr
#' @importFrom purrr map2_df map_df
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
