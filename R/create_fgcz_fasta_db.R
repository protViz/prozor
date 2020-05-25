#' create fasta db from one or more fasta files
#' @export
#'
create_fgcz_fasta_db <- function(databasedirectory , revLab = "REV_"){
    dir.exists(databasedirectory)
    dbname <- basename(databasedirectory)

    fasta <- grep("fasta", dir(databasedirectory), value = T)
    files1 <- file.path(databasedirectory, fasta)
    annot <- grep("annotation",dir(databasedirectory), value = T)

    annotation <- readLines(file.path(databasedirectory,annot))

    resDB <- createDecoyDB(files1,useContaminants = loadContaminantsFasta(), annot = annotation, revLab = revLab)
    resDB <- resDB[!duplicated(names(resDB))]

    if (is.null(revLab)) {
        xx <- file.path(dirname(databasedirectory), paste(dbname,"_",format(Sys.time(), "%Y%m%d"),".fasta" ,sep = ""))
    } else {
        dbname_decoy <- paste0(dbname,"_d")
        xx <- file.path(dirname(databasedirectory), paste(dbname_decoy,"_",format(Sys.time(), "%Y%m%d"),".fasta" ,sep = ""))
    }
    writeFasta(resDB, file = xx)

    bigstr <- paste(resDB, collapse = "")
    vec <- strsplit(bigstr,split = "")[[1]]

    aafreq <- table(vec)
    aafreq <- paste(capture.output(as.matrix(aafreq)),"\n", sep = "")

    length_s <- summary(sapply(resDB, getLength))
    length_s <- paste(capture.output(length_s),"\n", sep = "")


    summary <- paste0("name: ",basename(xx), "\nannotation:\n", annotation, "\n", "nr sequences: " , length(resDB), "\n")
    summary <- c(summary, "length summary:\n", length_s, "AA frequencies:\n", aafreq)

    return( list(resDB = resDB, summary = summary))
}
