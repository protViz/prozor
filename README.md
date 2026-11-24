[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/prozor)](https://cran.r-project.org/package=prozor)
[![](https://cranlogs.r-pkg.org/badges/prozor)](https://cran.r-project.org/package=prozor)
[![](https://cranlogs.r-pkg.org/badges/grand-total/prozor)](https://cran.r-project.org/package=prozor)


# prozor

- Determine the minimal protein set explaining peptide spectrum matches (PSMs).
- Create FASTA amino acid databases with decoys and contaminants.
- Estimate the false discovery rate for target-decoy search results on PSM, precursor, peptide and protein level.
- Compute dynamic SWATH (DIA) window sizes from the distribution of MS1 masses.
 
An HTML version of the package documentation can be found here:

https://protviz.github.io/prozor

## How to install:
For CRAN version (not the newest). Please use the GitHub version.

```r
install.packages("prozor")
```

This is how you install the GitHub version.

```r
install.packages("remotes")
remotes::install_github("protviz/prozor")
```

### for Developers

Download the git repo. Document new functions with roxygen2, then run these 2 commands to update NAMESPACE and the Rd files:

```r
library(devtools)
document()
```

## Creating a FASTA database from the command line

prozor installs the script `fgcz_create_fasta.R` (it needs the docopt package). It reads all FASTA files and the `annotation.txt` in a directory, adds contaminants and reversed decoy sequences, and writes one FASTA file:

```bash
script=$(Rscript -e 'cat(system.file("script", "fgcz_create_fasta.R", package = "prozor"))')
Rscript "$script" -h
Rscript "$script" path/to/database_dir -o path/to/output_dir
Rscript "$script" path/to/database_dir nodecoy -c none -o path/to/output_dir
```
