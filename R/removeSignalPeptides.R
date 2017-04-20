#' make id for chain compatible with uniprot
#' @export
#' @param sequence - aa sequence as string
#' @param id uniprot id id: sp|P30443|1A01_HUMAN
#' @param sp start position of chain numeric
#' @examples
#' seq <- "MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGR\
#' GEPRFIAVGYVDDTQFVRFDSDAASQKMEPRAPWIEQEGPEYWDQETRN\
#' MKAHSQTDRANLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQ\
#' DAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAVHAAEQRRVYLEGRC\
#' VDGLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEI\
#' TLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQ\
#' HEGLPKPLTLRWELSSQPTIPIVGIIAGLVLLGAVITGAVVAAVMWRRK\
#' SSDRKGGSYTQAASSDSAQGSDVSLTACKV"
#' nam <-"sp|P30443|1A01_HUMAN"
#' sp <- 24
#' makeIDUnip(seq, nam, sp)
makeIDUnip <- function(sequence, id, sp){
  tpos <- unlist(strsplit(id, split="\\|")) [1:2]
  tpos <- paste(tpos, "|", collapse="", sep="")
  sp <- paste( (sp) , "-", nchar(sequence), sep="")
  res <- paste(tpos,  sp, sep="")
  return(res)
}
#' make id for chain in format sp|P30443|1A01_HUMANs25
#' @export
#' @param sequence - aa sequence as string
#' @param id uniprot id id: sp|P30443|1A01_HUMAN
#' @param sp start position of chain numeric
#' @examples
#' seq <- "MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGR\
#' GEPRFIAVGYVDDTQFVRFDSDAASQKMEPRAPWIEQEGPEYWDQETRN\
#' MKAHSQTDRANLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQ\
#' DAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAVHAAEQRRVYLEGRC\
#' VDGLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEI\
#' TLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQ\
#' HEGLPKPLTLRWELSSQPTIPIVGIIAGLVLLGAVITGAVVAAVMWRRK\
#' SSDRKGGSYTQAASSDSAQGSDVSLTACKV"
#' nam <-"sp|P30443|1A01_HUMAN"
#' sp <- 24
#' makeID(seq, nam, sp)
makeID <- function(sequence, id, sp){
  res <- paste(id, "s", sp, sep="")
  return(res)
}

.makesig <- function( sequence, chainstart,  idfun=makeID ){
  attrlist <- attributes(sequence)
  newID <- idfun(sequence , attrlist$name, chainstart)
  attrlist$Annot <- newID
  attrlist$name <- newID
  seq <- substring(sequence , first=chainstart  )
  attributes(seq) <- attrlist
  return(seq)
}

#' remove signal peptides from main chain
#' @export
#' @param db uniprot fasta database as list
#' @param signal tab delimited file with signals
#' @param idfun function to generate id's
#' @examples
#' \dontrun{
#' library(prozor)
#'
#' hsfasta <- loadHomoSapiensFasta()
#' hssignal <- loadHomoSapiensSignalPeptides()
#' xx <- removeSignalPeptide(hsfasta, hssignal)
#' db <- hsfasta
#' signal <- hssignal
#' }
removeSignalPeptide <- function(db, signal, idfun=makeID){
  res <- db
  tmp2 <- as.numeric(gsub("^SIGNAL [0-9]+ ([0-9]+).*"  , "\\1", signal$Signal.peptide))

  cnamessplit <- strsplit(as.character(names(db)),split="\\|")
  protnam <-do.call("rbind",cnamessplit)

  stopifnot(sum(protnam[,2] == signal$Entry) == length(signal$Entry))

  for( i in 1:length(tmp2)){
    if(!is.na(tmp2[i])){
      res[[i]] <- .makesig(db[[i]], tmp2[i] + 1,idfun=idfun  )
    }
  }
  return(res)
}


