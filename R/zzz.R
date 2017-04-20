#R

.onAttach <- function(lib, pkg){
	if(interactive()){
		version <- utils::packageVersion('prozor')
		packageStartupMessage("Package 'prozor' version ", version)
		# packageStartupMessage("Type 'citation(\"prozor\")' for citing this R package in publications.")
	  invisible()
	}
}
