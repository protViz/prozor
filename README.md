[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/prozor)](https://cran.r-project.org/package=prozor)
[![](https://cranlogs.r-pkg.org/badges/grand-total/prozor)](https://cran.r-project.org/package=prozor)
[![R-CMD-check](https://github.com/prolfqua/prozor/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/prolfqua/prozor/actions/workflows/R-CMD-check.yaml)
[![lint](https://github.com/prolfqua/prozor/actions/workflows/lint.yaml/badge.svg)](https://github.com/prolfqua/prozor/actions/workflows/lint.yaml)
[![altdoc](https://github.com/prolfqua/prozor/actions/workflows/altdoc.yaml/badge.svg)](https://prolfqua.github.io/prozor/)


# prozor

- Determine the minimal protein set explaining peptide spectrum matches (PSMs).
- Create FASTA amino acid databases with decoys and contaminants.
- Estimate the false discovery rate for target-decoy search results on PSM, precursor, peptide and protein level.
- Compute dynamic SWATH (DIA) window sizes from the distribution of MS1 masses.

## Documentation

The website https://prolfqua.github.io/prozor/ has the function reference, the articles and the [changelog](https://prolfqua.github.io/prozor/NEWS.html).

Articles:

- [Creating Target Decoy Databases](https://prolfqua.github.io/prozor/vignettes/CreatingDecoyDB.html)
- [Peptide Annotation and Protein Inference](https://prolfqua.github.io/prozor/vignettes/PeptideAnnotationwithProzor.html)
- [Protein Sequence Coverage](https://prolfqua.github.io/prozor/vignettes/ProteinCoverage.html)
- [Target Decoy FDR](https://prolfqua.github.io/prozor/vignettes/TargetDecoyFDR_Example.html)
- [Target Decoy FDR from COMET search results](https://prolfqua.github.io/prozor/vignettes/TargetDecoyComet.html)
- [Computing Dynamic SWATH Windows](https://prolfqua.github.io/prozor/vignettes/cdsw.html)

## Installation

The CRAN version (0.3.1) is older than the GitHub version:

```r
install.packages("prozor")
```

The current version from GitHub:

```r
install.packages("remotes")
remotes::install_github("prolfqua/prozor")
```

### For developers

Clone the repository and use the Makefile targets: `make document` (roxygen2 docs and NAMESPACE), `make test`, `make check`, `make lint` and `make format` (air).

## Creating a FASTA database from the command line

prozor installs the script `fgcz_create_fasta.R` (it needs the docopt package). It reads all FASTA files and the `annotation.txt` in a directory, adds contaminants and reversed decoy sequences, and writes one FASTA file:

```bash
script=$(Rscript -e 'cat(system.file("script", "fgcz_create_fasta.R", package = "prozor"))')
Rscript "$script" -h
Rscript "$script" path/to/database_dir -o path/to/output_dir
Rscript "$script" path/to/database_dir nodecoy -c none -o path/to/output_dir
```
