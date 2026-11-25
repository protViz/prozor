[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/prozor)](https://cran.r-project.org/package=prozor)
[![](https://cranlogs.r-pkg.org/badges/grand-total/prozor)](https://cran.r-project.org/package=prozor)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24535408.svg)](https://doi.org/10.5281/zenodo.24535408)
[![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/prolfqua/prozor/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/prolfqua/prozor/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/prolfqua/prozor/branch/master/graph/badge.svg)](https://app.codecov.io/gh/prolfqua/prozor)
[![lint](https://github.com/prolfqua/prozor/actions/workflows/lint.yaml/badge.svg)](https://github.com/prolfqua/prozor/actions/workflows/lint.yaml)
[![altdoc](https://github.com/prolfqua/prozor/actions/workflows/altdoc.yaml/badge.svg)](https://prolfqua.github.io/prozor/)


# prozor

- Determine the minimal protein set explaining peptide spectrum matches (PSMs).
- Create FASTA amino acid databases with decoys and contaminants.
- Estimate the false discovery rate for target-decoy search results on PSM, precursor, peptide and protein level.
- Compute dynamic SWATH (DIA) window sizes from the distribution of MS1 masses.

Terms:

- PSM (peptide-spectrum match): a fragment spectrum that a search engine assigned to a peptide.
- Decoy: a reversed protein sequence added to the database; matches to decoys estimate the number of false matches.
- DIA (data-independent acquisition, also SWATH): the instrument fragments all precursors within a series of m/z windows.
- MS1: the survey scan that measures the m/z of the intact peptide ions (precursors).

Version 0.5.0 breaks the API of CRAN 0.3.1: functions, arguments and result columns are renamed, and R6 classes compute the DIA windows and the FDR. See [`NEWS.md`](https://github.com/prolfqua/prozor/blob/master/NEWS.md).

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

## Example

```r
library(prozor)

# protein sequences and the peptides identified in a human sample, bundled with prozor
fasta <- read_fasta(system.file("p1000_db1_example/Annotation_canSeq.fasta.gz", package = "prozor"))
peptides <- unique(read.delim(system.file("extdata/IDResults.txt.gz", package = "prozor"))$peptideSeq)

# the proteins each peptide occurs in, then the minimal protein set that explains all peptides
annotation <- annotate_peptides(peptides, fasta)
groups <- greedy_parsimony(annotation)
c(peptides = length(peptides), proteins = length(unique(annotation$protein_id)), groups = length(unique(groups)))
#> peptides proteins   groups
#>     1520      382      244

# FDR of target-decoy search results: the score at a q-value of 1 % (smaller scores are better)
fdr <- TargetDecoyFDR$new(fdr_sample$score, decoy_hit = grepl("^REV_", fdr_sample$protein_id), larger_better = FALSE)
fdr$score_at(q_value = 0.01)
#> [1] 0.03394749
```

The 1,520 peptides match 382 proteins; 244 proteins or protein groups explain all of them.

## Creating a FASTA database from the command line

prozor installs the script `fgcz_create_fasta.R` (it needs the docopt package). It reads all FASTA files and the `annotation.txt` in a directory, adds contaminants and reversed decoy sequences, and writes one FASTA file:

```bash
script=$(Rscript -e 'cat(system.file("script", "fgcz_create_fasta.R", package = "prozor"))')
Rscript "$script" -h
Rscript "$script" path/to/database_dir -o path/to/output_dir
Rscript "$script" path/to/database_dir nodecoy -c none -o path/to/output_dir
```

## Citation

Cite prozor with its DOI [10.32614/CRAN.package.prozor](https://doi.org/10.32614/CRAN.package.prozor); `citation("prozor")` gives the full reference.

Each GitHub release is archived on Zenodo under [10.5281/zenodo.24535408](https://doi.org/10.5281/zenodo.24535408), which resolves to the latest release.
