.onAttach <- function(lib, pkg) {
  if (interactive()) {
    packageStartupMessage("Package 'prozor' version ", utils::packageVersion("prozor"))
  }
}
