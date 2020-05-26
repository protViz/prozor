#' create fasta db from one or more fasta files
#' @export
#'
create_fgcz_fasta_db <- function(databasedirectory , useContaminants = loadContaminantsFasta(), revLab = "REV_"){
    mcall <- match.call()

    dir.exists(databasedirectory)
    dbname <- basename(databasedirectory)

    fasta <- grep("fasta", dir(databasedirectory), value = T)
    files1 <- file.path(databasedirectory, fasta)
    annot <- grep("annotation",dir(databasedirectory), value = T)

    annotation <- readLines(file.path(databasedirectory,annot))

    resDB <- createDecoyDB(files1,useContaminants = loadContaminantsFasta(), annot = annotation, revLab = revLab)
    resDB <- resDB[!duplicated(names(resDB))]

    if (is.null(revLab)) {
        filepath <- file.path(dirname(databasedirectory), paste(dbname,"_",format(Sys.time(), "%Y%m%d"),".fasta" ,sep = ""))
    } else {
        dbname_decoy <- paste0(dbname,"_d")
        filepath <- file.path(dirname(databasedirectory), paste(dbname_decoy,"_",format(Sys.time(), "%Y%m%d"),".fasta" ,sep = ""))
    }
    writeFasta(resDB, file = filepath)

    {
        bigstr <- paste(resDB, collapse = "")
        vec <- strsplit(bigstr,split = "")[[1]]

        aafreq <- table(vec)
        aafreq <- paste(capture.output(as.matrix(aafreq)),"\n", sep = "")
        length_s <- summary(sapply(resDB, getLength))
        length_s <- paste(capture.output(length_s),"\n", sep = "")


        summary <- paste0(
            "Database created with prozor: ", capture.output(resDB$mcall), "\n",
            "where databasedirectory was prepared according to https://fgcz-intranet.uzh.ch/tiki-index.php?page=SOPrequestFASTA \n\n",
            "FASTA name: ",basename(filepath),
            "\nannotation:\n",
            annotation, "\n",
            "nr sequences: " , length(resDB), "\n")
        summary <- c(summary, "length summary:\n", length_s, "AA frequencies:\n", aafreq)

    }
    return( list(resDB = resDB, filepath = filepath, summary = summary))
}
