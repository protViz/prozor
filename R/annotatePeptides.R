
#' Annotate peptides with protein ids
#'
#' peptides which do not have protein assignment drop out
#'
#' Each occurrence of a peptide in a protein is one row. \code{matched} is TRUE when the residues before
#' that occurrence match \code{prefix} and the residues after it match \code{suffix}; by default, when the
#' occurrence follows K or R, starts the protein, or follows its initial methionine.
#'
#' @param pepinfo - list of peptides - sequence, optional modified sequence, charge state.
#' @param fasta - object as created by readPeptideFasta
#' @param peptide - name of column containing peptide sequences default "peptideSeq"
#' @param prefix - default "(([RK])|(^)|(^M))"
#' @param suffix - default ""
#' @return data.frame with columns "peptideSeq", "proteinID","Offset","proteinSequence","matched", "lengthPeptide","proteinlength"
#' @export
#' @examples
#'
#' library(dplyr)
#'
#' file = system.file("extdata/IDResults.txt.gz" , package = "prozor")
#' specMeta <- utils::read.delim(file)
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
    pepinfo[[peptide]] <- as.character(pepinfo[[peptide]])
    pepinfo$lengthPeptide <- nchar(pepinfo[[peptide]])

    pepseq  = unique(as.character(pepinfo[, peptide]))
    restab <- annotateAHO(pepseq, fasta, peptide = peptide)
    if (!is.null(restab)) {

        before <- substr(restab$proteinSequence, 1, restab$Offset - 1)
        after <- substring(restab$proteinSequence, restab$Offset + nchar(restab[[peptide]]))
        restab$matched <- grepl(paste0("(", prefix, ")$"), before) & grepl(paste0("^(", suffix, ")"), after)
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
#' @param peptide - name of the result column holding the peptide sequences, default "peptideSeq"
#' @import AhoCorasickTrie
#' @return A data.frame with proteinID, peptideSeq, Offset and proteinSequence
#' @examples
#'
#' library(dplyr)
#'
#' file = system.file("extdata/IDResults.txt.gz" , package = "prozor")
#' specMeta <- utils::read.delim(file)
#' upeptide <- unique(specMeta$peptideSeq)
#' resCan <-
#'    prozor::readPeptideFasta(
#'        system.file("p1000_db1_example/Annotation_canSeq.fasta.gz" , package = "prozor"))
#' resCanU <- resCan[!duplicated(unlist(resCan))]
#' annotAll = annotateAHO(upeptide[seq_len(20)], resCanU)
#' dim(annotAll)
#'
#' @export
annotateAHO <- function(pepseq, fasta, peptide = "peptideSeq") {
    #100_000 peptides
    #40_000 Proteine

    pepseq <- trimws(unique(pepseq))
    proteinIDS <- names(unlist(fasta))
    fasta <- trimws(unlist(fasta))
    names(fasta) <- proteinIDS

    res <-
        AhoCorasickSearch(unique(pepseq) , unlist(fasta), alphabet = "aminoacid")
    if (length(res) == 0)
    {
        return(NULL)
    }
    xx <- dplyr::bind_rows(Map(
        function(name, x) {
            data.frame(
                proteinID = name,
                dplyr::bind_rows(x),
                stringsAsFactors = FALSE
            )
        },
        names(res),
        res
    ))
    colnames(xx)[colnames(xx) == "Keyword"] <- peptide
    dbframe <- data.frame(
        proteinID = names(fasta),
        proteinSequence = as.character(unlist(fasta)),
        stringsAsFactors = FALSE
    )
    matches <- dplyr::inner_join(xx, dbframe, by = "proteinID")
    return(matches)
}
