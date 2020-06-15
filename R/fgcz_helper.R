#' copy to server
#' @keywords internal
#' @export
#' @param filepath path of file
#' @param to destination folder
#' @param username username
#' @param host host
scp_copy <- function(
    filepath,
    to = "/srv/www/htdocs/FASTA",
    username = "wolski",
    host="fgcz-r-021.uzh.ch")
{
    session <- ssh::ssh_connect(paste0(username, "@", host ))
    if (file.exists(filepath)) {
        ssh::scp_upload(session, filepath, to = to)
    }else{
        warning("file :" , file , " does not exist!")
    }
    ssh::ssh_disconnect(session)

}

#' copy to server
#' @keywords internal
#' @export
#' @param project_id project_id
#' @param filepath filepath
#' @param summary summary (string with info about DB)
#' @param python_path python_path
bfabric_save_fasta <- function(project_id,
                               filepath,
                               summary,
                               python_path = "c:/Program Files/Python38"){
    library(reticulate)
    reticulate::use_python(python_path)
    bfabricpy <- reticulate::import("bfabric.scripts.bfabric_save_fasta")
    bfabricpy$save_fasta(projectid = project_id,
                         fasta_file = filepath,
                         description_resource =  paste(summary, collapse = ""),
                         description_workunit = paste(summary, collapse = "")
    )

}
