#' make db summary which includes number of sequences, and amino acid statistcis
#' @export
#' @return array of strings which should be passed to the cat function.
#' @examples
#' file = system.file("extdata/fgcz_contaminants2022_20220405.fasta.gz",package="prozor")
#' xx <- prozor::readPeptideFasta(file)
#' cat(make_fasta_summary(xx))
#' rbenchmark::benchmark(make_fasta_summary(xx),replications = 10)
#' rbenchmark::benchmark(make_fasta_summary(xx, old = TRUE), replications = 10)
#' make_fasta_summary(xx, as_string = FALSE)
#' make_fasta_summary(xx, old =TRUE, as_string = FALSE)
make_fasta_summary <- function(resDB, old = FALSE, as_string = TRUE){
    if (old) {
        bigstr <- paste(resDB, collapse = "")
        vec <- strsplit(bigstr,split = "")[[1]]
        aafreq <- table(vec)
    } else {
        res <- list()
        tmp <- lapply(resDB, function(x){table(strsplit(x,split = "")[[1]])})
        for (i in seq_along(tmp)) {
            x <- tmp[[i]]
            for (j in names(x) ) {
                if (is.null(res[[j]])) {
                    res[[j]] <- as.numeric(x[j])
                } else {
                    res[[j]] <- res[[j]] + as.numeric(x[j])
                }
            }
        }
        aafreq <- unlist(res)
        aafreq <- aafreq[order(names(aafreq))]
    }
    length_s <- summary(vapply(resDB, seqinr::getLength, numeric(1)))

    if (as_string) {
        aafreq <- paste(utils::capture.output(as.matrix(aafreq)),"\n", sep = "")
        length_s <- paste(utils::capture.output(length_s),"\n", sep = "")
        summaryRes <- c("nr sequences:\n", length(resDB) , "\n length summary:\n", length_s, "AA frequencies:\n", aafreq)

    } else {
        summaryRes <- list()
        summaryRes$nrSequences <- length(resDB)
        summaryRes$lengthSummary <- length_s
        summaryRes$aafreq <- aafreq
    }
    return(summaryRes)
}


#' create fasta db from one or more fasta files
#' @export
#' @param databasedirectory directory with fasta files
#' @param useContaminants contaminants to add
#' @param revLab reverse label
#' @param outputdir output directory
#' @return list list(resDB, filepath , summary, mcall, dbname)
#' @examples
#' print("NO exmple since function also writes the fasta files")
#'
create_fgcz_fasta_db <- function(databasedirectory ,
                                 useContaminants = loadContaminantsFGCZ2022(),
                                 revLab = "REV_",
                                 outputdir = NULL,
                                 make_summary = TRUE){
    mcall <- match.call()

    dir.exists(databasedirectory)
    dbname <- basename(databasedirectory)

    fasta <- grep("fasta", dir(databasedirectory), value = TRUE)
    files1 <- file.path(databasedirectory, fasta)
    annot <- grep("annotation.txt",dir(databasedirectory), value = TRUE)
    if (length(annot) == 0) { stop("NO annotation file found.") }
    annotation <- readLines(file.path(databasedirectory, annot))[1]

    resDB <- createDecoyDB(files1,
                           useContaminants = useContaminants,
                           annot = annotation,
                           revLab = revLab)
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
    summaryStr <- paste0(
        "Database created with prozor: ", paste(utils::capture.output(mcall), collapse = "\n"), "\n",
        "where databasedirectory was prepared according to https://fgcz-intranet.uzh.ch/tiki-index.php?page=SOPrequestFASTA \n\n",
        "\n      FASTA name : ", dbname,
        "\n FASTA file name : ", basename(filepath),
        "\n  written to dir : ", outputdir,
        "\n       annotation:\n",
        annotation, "\n",
        "\n      nr sequences: " , length(resDB), "\n")

    if ( make_summary ) {
        addSR <- make_fasta_summary(resDB, old = FALSE, as_string = TRUE)
        summaryStr <- c(summaryStr, addSR)
    }
    return( list(resDB = resDB, filepath = filepath, summary = summaryStr, mcall = mcall, dbname = dbname))
}
