# prozor 0.3.4

* `greedy_parsimony()` assigns every peptide again. When several proteins tie, only those explaining exactly the same peptides are reported together (`"A;B"`); the others compete in later rounds. Since 2023 (GitHub only, not on CRAN), unrelated tied proteins were merged into one group and their remaining peptides were left without a protein.
* `greedyRes2Matrix()` works with modified peptides: the charge is taken after the last dot, so `"AC[+57.0]K.2"` gives peptide `AC[+57.0]K` and charge 2. Such names used to fail with "length of 'dimnames' [2] not equal to array extent".
* `Cdsw$sampling_breaks()` stops with an error when even equal-width windows violate `minwindow` or `maxwindow`; it used to warn and then never return.
* `removeSignalPeptide()` works on a database with a single entry.
* `annotateAHO()` and `annotatePeptides()` no longer print dplyr's "Joining with `by = join_by(proteinID)`" message.
* Breaking: `annotatePeptides()` sets `matched` for each occurrence of a peptide in a protein. It used to give every occurrence in a protein the same value, so a second, non-tryptic occurrence was flagged as tryptic.
* Breaking: `computeFDR()` and `computeFDRwithID()` rank decoys before targets within tied scores, so q-values no longer depend on the input order.
* Breaking: the default `decoy` pattern of `computeFDRwithID()` is `"^REV_"`; only IDs starting with `REV_` count as decoys (before: any ID containing `REV_`).
* Breaking: `predictScoreFDR()` stops with an error when no score reaches the requested FDR, instead of returning `Inf` with a warning.
* `prepareMatrix()` documents only the weightings it implements (`NULL`, `"inverse"` or one numeric weight per row); other values give an error naming them instead of "FALSE is not TRUE".
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
