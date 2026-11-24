# prozor 0.3.4

* `load_special_prot_2024()` works again: the bundled `fgcz_special_prot_20241112.fasta.gz` is a plain gzipped FASTA file again instead of a gzipped tar archive (#7).
* Removed the undocumented debug dataset `mm`.
* The vignette "Peptide Annotation and Protein Inference" builds again.
* Documented the arguments `peptide` of `annotateAHO()`, `resDB`, `old` and `as_string` of `make_fasta_summary()`, `make_summary` of `create_fgcz_fasta_db()` and `use.last.ij` of `prepareMatrix()`.
* `create_fgcz_fasta_db()` now requires `outputdir` (it used to write next to `databasedirectory` by default) and stops if `databasedirectory` does not exist. `inst/script/fgcz_create_fasta.R` still writes next to the FASTA directory when `-o` is not given.
* `plotFDR()` restores the graphics parameters it changes.
* prozor no longer depends on protViz.
* The `Cdsw` help page no longer documents method arguments without a usage section (r-devel NOTE).
* New or runnable examples for `create_fgcz_fasta_db()`, `greedyRes2Matrix()`, `removeSignalPeptide()` and `writeFasta()`; removed debugging code from examples and vignettes.
* `annotatePeptides()` no longer triggers dplyr's `funs()` deprecation warning; results are unchanged.
* Fewer dependencies: purrr, rlang and stringr are no longer needed; readr and docopt moved to Suggests (docopt is only needed by `inst/script/fgcz_create_fasta.R`).

# prolfqua 0.2.15

## Prepare submission to Bioconductor

* address warnings and notes from BiocCheck tool.
