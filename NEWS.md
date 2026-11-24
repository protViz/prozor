# prozor 0.4.0

## Breaking changes

* All exported functions and their arguments use snake_case. The old names are gone, without aliases:

| Old | New |
|---|---|
| `annotateAHO()`, `annotatePeptides()` | `annotate_aho()`, `annotate_peptides()` |
| `computeFDR()`, `computeFDRwithID()` | `compute_fdr()`, `compute_fdr_with_id()` |
| `plotFDR()`, `predictScoreFDR()` | `plot_fdr()`, `predict_score_fdr()` |
| `createDecoyDB()`, `reverseSeq()` | `create_decoy_db()`, `reverse_seq()` |
| `readPeptideFasta()`, `writeFasta()` | `read_peptide_fasta()`, `write_fasta()` |
| `loadContaminantsFasta2019()`, `loadContaminantsFasta2021()`, `loadContaminantsFGCZ2022()`, `load_universal_contaminants_2024()`, `load_special_prot_2024()` | `load_contaminants("fgcz2019")`, `"fgcz2021"`, `"fgcz2022"`, `"universal2024"`, `"special2024"` |
| `prepareMatrix()`, `greedyRes2Matrix()` | `prepare_matrix()`, `greedy_to_table()` |
| `removeSignalPeptide()`, `makeID()`, `makeIDUnip()` | `remove_signal_peptide()`, `make_chain_id()`, `make_chain_id_uniprot()` |
| `readjustWindows()` | `readjust_windows()` |
| `Cdsw` methods `asTable()`, `showCycle()`, `optimizeWindows()` | `as_table()`, `show_cycle()`, `optimize_windows()` |
| arguments `useContaminants`, `revLab`, `noHuman`, `databasedirectory`, `outputdir`, `resDB`, `fdrObj`, `qValue`, `proteinID`, `peptideID`, `idfun`, `ID` | `contaminants`, `rev_lab`, `no_human`, `database_dir`, `output_dir`, `db`, `fdr`, `q_value`, `protein_id`, `peptide_id`, `id_fun`, `id` |
| data set `fdrSample` | `fdr_sample` |

* Result columns and list elements are snake_case too:
  * `annotate_aho()`, `annotate_peptides()`: `protein_id`, `offset`, `protein_sequence`, `peptide_length`, and `peptide_seq` for vector input (was `proteinID`, `Offset`, `proteinSequence`, `lengthPeptide`, `peptideSeq`)
  * `compute_fdr()`: `fpr`, `simple_fdr`, `q_value_fpr`, `q_value_simple_fdr`; `compute_fdr_with_id()` adds `id` (was `FPR`, `SimpleFDR`, `qValue_FPR`, `qValue_SimpleFDR`, `ID`)
  * `greedy_to_table()`: `peptide`, `precursor_charge`, `protein` (was `Peptide`, `precursorCharge`, `Protein`)
  * `make_fasta_summary(as_string = FALSE)`: `nr_sequences`, `length_summary`, `aa_freq`; `create_fgcz_fasta_db()`: `db`, `file_path`, `summary`, `call`, `db_name`
* `predict_score_fdr(method = )` takes `"simple_fdr"` (default) or `"fpr"` (was `"SimpleFDR"`, `"FPR"`).
* Default column names follow the results: `annotate_peptides(peptide = "peptide_seq")`, `prepare_matrix(protein_id = "protein_id", peptide_id = "peptide_seq")` (was `"peptideSeq"`, `"proteinID"`, `"strippedSequence"`).
* The data sets `fdr_sample`, `protpepmetashort` and `pepprot` have snake_case columns (`protein_id`, `peptide_seq`, `peptide_mod_seq`, `precursor_charge`, `protein_length`, `peptide_length`); `bcd` keeps COMET's column names.
* `greedy()`, last on CRAN in 0.3.1, has no alias; use `greedy_parsimony()`.
* `write_fasta(fasta, file)` takes one list of sequences, data first (was `writeFasta(file, ...)`); combine several lists with `c()`.
* Plots use ggplot2, which moved from Suggests to Imports. `plot_fdr()`, `Cdsw$plot()` and `Cdsw$show_cycle()` return ggplot objects and draw nothing themselves. The `plot` argument of `Cdsw$sampling_breaks()`, `Cdsw$optimize_windows()` and `readjust_windows()` is gone; use `Cdsw$plot_breaks()` and `plot_windows()`.
* `annotate_peptides()` sets `matched` for each occurrence of a peptide in a protein. It used to give every occurrence in a protein the same value, so a second, non-tryptic occurrence was flagged as tryptic.
* `annotate_aho()` and `annotate_peptides()` stop for peptides with characters other than upper-case letters (AhoCorasickTrie 0.1.3 crashed R on them) and for duplicated protein names in `fasta` (they multiplied the result rows).
* `compute_fdr()` and `compute_fdr_with_id()` rank decoys before targets within tied scores, so q-values no longer depend on the input order.
* The default `decoy` pattern of `compute_fdr_with_id()` is `"^REV_"`: only IDs starting with `REV_` count as decoys. Pass `decoy = "REV_"` for the old matching anywhere in the ID.
* `predict_score_fdr()` stops with an error when no score reaches the requested FDR, instead of returning `Inf` with a warning.
* `create_fgcz_fasta_db()` requires `output_dir` (it used to write next to the database directory by default). `inst/script/fgcz_create_fasta.R` still writes next to the FASTA directory when `-o` is not given.
* `greedy_parsimony()` no longer prints its run time.
* prozor requires R 4.4.0 or later.

## New features

* `load_contaminants(set, no_human)` loads any of the bundled contaminant sets.
* `plot_windows()` shows the precursor masses with the DIA window boundaries; `Cdsw$plot_breaks()` compares the current breaks with equal-width and quantile breaks.

## Minor improvements and fixes

* `greedy_parsimony()` assigns every peptide again. When several proteins tie, only those explaining exactly the same peptides are reported together (`"A;B"`); the others compete in later rounds. Since 2023 (GitHub only, not on CRAN), unrelated tied proteins were merged into one group and their remaining peptides were left without a protein.
* `greedy_to_table()` works with modified peptides: the charge is taken after the last dot, so `"AC[+57.0]K.2"` gives peptide `AC[+57.0]K` and charge 2. Such names used to fail with "length of 'dimnames' [2] not equal to array extent".
* `annotate_peptides()` and `prepare_matrix()` accept tibbles; a tibble used to crash R or fail with "cannot xtfrm data frames".
* `write_fasta()` keeps a `>` inside a header (`ADP>ATP` used to become `ADPATP`).
* `remove_signal_peptide()` reads UniProt's current signal peptide format `SIGNAL 1..20`, and works on a database with a single entry.
* `Cdsw$sampling_breaks()` stops with an error when even equal-width windows violate `minwindow` or `maxwindow`; it used to warn and then never return.
* `annotate_aho()` and `annotate_peptides()` no longer print dplyr's "Joining with `by = join_by(proteinID)`" message, nor dplyr's `funs()` deprecation warning.
* Error messages name the argument and the offending value, e.g. for an unknown `method` of `predict_score_fdr()` or `weighting` of `prepare_matrix()`.
* `prepare_matrix()` documents only the weightings it implements (`NULL`, `"inverse"` or one numeric weight per row).
* `load_contaminants("special2024")` works: the bundled `fgcz_special_prot_20241112.fasta.gz` is a plain gzipped FASTA file again instead of a gzipped tar archive (#7).
* `inst/script/fgcz_create_fasta.R` stops early with a clear message for an unknown `--contamin` value.
* The summary of `create_fgcz_fasta_db()` no longer contains a link to the FGCZ intranet.
* `load_contaminants()` cites the source of the `"universal2024"` set; `inst/COPYRIGHTS` lists its origin.
* Help pages use markdown, link related functions in families, and describe the actual return values (e.g. the elements of `compute_fdr()`'s result).
* New or runnable examples for `create_fgcz_fasta_db()`, `greedy_to_table()`, `remove_signal_peptide()` and `write_fasta()`; removed debugging code from examples and vignettes.
* The vignettes use the new names and ggplot2; "Peptide Annotation and Protein Inference" builds again.
* The `Cdsw` help page no longer documents method arguments without a usage section (r-devel NOTE).
* Removed the undocumented debug dataset `mm`.
* Fewer dependencies: protViz, purrr, rlang and stringr are no longer needed; readr and docopt moved to Suggests (docopt is only needed by `inst/script/fgcz_create_fasta.R`).

# prozor 0.2.15

## Prepare submission to Bioconductor

* address warnings and notes from BiocCheck tool.
