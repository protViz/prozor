#' create fasta db from one or more fasta files
#' @export
#' @param databasedirectory directory with fasta files
#' @param useContaminants contaminants to add
#' @param revLab reverse label
#' @param outputdir output directory
#'
create_fgcz_fasta_db <- function(databasedirectory ,
                                 useContaminants = loadContaminantsFasta2021(),
                                 revLab = "REV_",
                                 outputdir = NULL){
    mcall <- match.call()

    dir.exists(databasedirectory)
    dbname <- basename(databasedirectory)

    fasta <- grep("fasta", dir(databasedirectory), value = TRUE)
    files1 <- file.path(databasedirectory, fasta)
    annot <- grep("annotation.txt",dir(databasedirectory), value = TRUE)
    if (length(annot) == 0) { stop("NO annotation file found.") }
    annotation <- readLines(file.path(databasedirectory, annot))[1]

    resDB <- createDecoyDB(files1,useContaminants = useContaminants, annot = annotation, revLab = revLab)
    resDB <- resDB[!duplicated(names(resDB))]

    if (is.null(outputdir)) {
        outputdir <- dirname(databasedirectory)
    }

    if (is.null(revLab)) {
        filepath <- file.path(outputdir, paste(dbname,"_",format(Sys.time(), "%Y%m%d"),".fasta" ,sep = ""))
    } else {
        dbname <- paste0(dbname,"_d")
        filepath <- file.path(outputdir, paste(dbname,"_",format(Sys.time(), "%Y%m%d"),".fasta" ,sep = ""))
    }

    message("writing db to : ", filepath)
    writeFasta(resDB, file = filepath)

    {
        bigstr <- paste(resDB, collapse = "")
        vec <- strsplit(bigstr,split = "")[[1]]

        aafreq <- table(vec)
        aafreq <- paste(utils::capture.output(as.matrix(aafreq)),"\n", sep = "")
        length_s <- summary(sapply(resDB, seqinr::getLength))
        length_s <- paste(utils::capture.output(length_s),"\n", sep = "")


        summary <- paste0(
            "Database created with prozor: ", paste(utils::capture.output(mcall), collapse = "\n"), "\n",
            "where databasedirectory was prepared according to https://fgcz-intranet.uzh.ch/tiki-index.php?page=SOPrequestFASTA \n\n",
            "\n      FASTA name : ", dbname,
            "\n FASTA file name : ", basename(filepath),
            "\n  written to dir : ", outputdir,
            "\n       annotation:\n",
            annotation, "\n",
            "\n      nr sequences: " , length(resDB), "\n")
        summary <- c(summary, "length summary:\n", length_s, "AA frequencies:\n", aafreq)

    }
    return( list(resDB = resDB, filepath = filepath, summary = summary, mcall = mcall, dbname = dbname))
}
