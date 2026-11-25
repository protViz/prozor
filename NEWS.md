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
  * `compute_fdr()`: `fdr_concatenated`, `simple_fdr`, `q_value_fdr_concatenated`, `q_value_simple_fdr`; `compute_fdr_with_id()` adds `id` (was `FPR`, `SimpleFDR`, `qValue_FPR`, `qValue_SimpleFDR`, `ID`). `FPR` became `fdr_concatenated` because 2 * decoys / all hits estimates the FDR of a concatenated target-decoy search, not a false positive rate; the legend of `plot_fdr()` uses the new names too
  * `greedy_to_table()`: `peptide`, `precursor_charge`, `protein` (was `Peptide`, `precursorCharge`, `Protein`); `precursor_charge` is an integer instead of a character
  * `make_fasta_summary(as_string = FALSE)`: `nr_sequences`, `length_summary`, `aa_freq`; `create_fgcz_fasta_db()`: `db`, `file_path`, `summary`, `call`, `db_name`
* `predict_score_fdr(method = )` takes `"simple_fdr"` (default) or `"fdr_concatenated"` (was `"SimpleFDR"`, `"FPR"`).
* `predict_score_fdr(q_value = )` is a fraction like the q-values it reads, `0.01` for 1 % (the default); `qValue` was in percent. Values outside 0 to 1 stop with an error.
* Default column names follow the results: `annotate_peptides(peptide = "peptide_seq")`, `prepare_matrix(protein_id = "protein_id", peptide_id = "peptide_seq")` (was `"peptideSeq"`, `"proteinID"`, `"strippedSequence"`).
* The data sets `fdr_sample`, `protpepmetashort` and `pepprot` have snake_case columns (`protein_id`, `peptide_seq`, `peptide_mod_seq`, `precursor_charge`, `protein_length`, `peptide_length`); `bcd` keeps COMET's column names.
* `greedy()`, last on CRAN in 0.3.1, has no alias; use `greedy_parsimony()`.
* `prepare_matrix()` has no `sep` argument (it was not used), and `make_fasta_summary()` no `old` argument (it selected a slower implementation with the same counts).
* `Cdsw()` has no `digits` argument and the object no `digits` field: they were stored and never read. `constant_breaks()`, `quantile_breaks()`, `sampling_breaks()` and `optimize_windows()` take their own `digits`.
* `Cdsw$sampling_breaks()` is deterministic: it mixes the masses with uniformly spread data and takes the quantile breaks of the mixture, with the smallest share of uniform data for which all windows lie within `minwindow` and `maxwindow` (the quantile breaks themselves when they fit). It used to add random uniform masses until the quantile windows fitted, so the breaks changed from call to call unless a seed was set, and it took seconds to minutes, or never returned, when `maxwindow` was close to the width of equal windows (more than 120 s at 1.02 times that width; now 0.05 s). The breaks differ from those of earlier versions; the counts per window are as even as before: on `masses` with `maxwindow = 100` they range from 3,120 to 6,391 (coefficient of variation 0.14), with the random mixing from about 3,050 to 6,420 (0.14).
* Precursors per window are counted one way: the masses `m` with `from < m <= to`, the intervals of `hist()`. `Cdsw$as_table()` left out masses on a window boundary (`from < m < to`), and `readjust_windows()` and `Cdsw$optimize_windows()` counted them in both neighboring windows (`from <= m <= to`). `Cdsw$plot()` and `Cdsw$error()` count as before for breaks that enclose the masses, as all `*_breaks()` methods return them. Counts change only for masses exactly on a boundary, e.g. in one of 25 windows of `quantile_breaks()` on `masses` with `as_table(overlap = 0)`.
* `write_fasta(fasta, file)` takes one list of sequences, data first (was `writeFasta(file, ...)`); combine several lists with `c()`.
* Plots use ggplot2, now in Imports. `plot_fdr()`, `Cdsw$plot()` and `Cdsw$show_cycle()` return ggplot objects and draw nothing themselves. The `plot` argument of `Cdsw$sampling_breaks()`, `Cdsw$optimize_windows()` and `readjust_windows()` is gone; use `Cdsw$plot_breaks()` and `plot_windows()`.
* `annotate_peptides()` sets `matched` for each occurrence of a peptide in a protein. It used to give every occurrence in a protein the same value, so a second, non-tryptic occurrence was flagged as tryptic.
* `annotate_aho()` and `annotate_peptides()` stop for peptides with characters other than upper-case letters (AhoCorasickTrie 0.1.3 crashed R on them) and for duplicated protein names in `fasta` (they multiplied the result rows).
* `annotate_aho()` no longer removes whitespace around the peptides: it stops and shows the peptide escaped, e.g. `"PEPTIDEK\r"`. `annotate_peptides()` still accepts such peptides: it removes the whitespace first.
* `annotate_peptides()` stops, naming the columns, when `pepinfo` has a column named like a result column (`protein_id`, `offset`, `protein_sequence`, `matched`, `peptide_length`). A `protein_id` or `offset` column, e.g. of search engine output, gave `protein_id.x` and `protein_id.y`, on which `prepare_matrix()` failed; a `peptide_length` column was overwritten.
* `compute_fdr()` and `compute_fdr_with_id()` rank decoys before targets within tied scores, so q-values no longer depend on the input order.
* The default `decoy` pattern of `compute_fdr_with_id()` is `"^REV_"`: only IDs starting with `REV_` count as decoys. Pass `decoy = "REV_"` for the old matching anywhere in the ID.
* `predict_score_fdr()` stops with an error when no score reaches the requested FDR, instead of returning `Inf` with a warning.
* `create_fgcz_fasta_db()` requires `output_dir` (it used to write next to the database directory by default). `inst/script/fgcz_create_fasta.R` still writes next to the FASTA directory when `-o` is not given.
* `greedy_parsimony()` no longer prints its run time.
* `greedy_parsimony()` ranks proteins by the number of peptides they explain, as documented, also for a weighted matrix such as `prepare_matrix(weighting = "inverse")`: it gives the same result as the matrix with weights 1. It ranked by the sum of the weights, so a protein with fewer but less shared peptides could win the shared ones, and ties were decided by comparing sums of fractions.
* prozor requires R 4.4.0 or later.
* `prepare_matrix()` applies `use.last.ij` to numeric weights too: by default a repeated peptide-protein pair keeps the weight of its last row, as it keeps a single entry for the other weightings. It summed the weights; pass `use.last.ij = FALSE` for the sum.
* prozor no longer imports seqinr: `read_peptide_fasta()` and `write_fasta()` use base R. The entries are still `SeqFastaAA` objects with the attributes `name` and `Annot`, so seqinr's `getAnnot()`, `getName()` and `getLength()` work on them as before. `library(prozor)` no longer loads seqinr, so `summary()` of an entry gives R's default summary of a string unless seqinr is loaded (it used to dispatch to seqinr's method, which returns the length, the composition and the amino-acid properties).
* `create_decoy_db()` and `create_fgcz_fasta_db()` add the contaminant sets `"special2024"` and `"universal2024"` by default, like `inst/script/fgcz_create_fasta.R` (was `"fgcz2022"`).

## New features

* `load_contaminants(set, no_human)` loads any of the bundled contaminant sets.
* `plot_windows()` shows the precursor masses with the DIA window boundaries; `Cdsw$plot_breaks()` compares the current breaks with equal-width and quantile breaks.
* `residue_coverage()`, `sequence_coverage()` and `plot_coverage()` compute and plot the peptide coverage of protein sequences from the result of `annotate_peptides()`; the new vignette "Protein Sequence Coverage" shows them, including coverage of the minimal protein set and of isoforms.

## Minor improvements and fixes

* `greedy_parsimony()` assigns every peptide that matches a protein again. When several proteins tie, only those explaining exactly the same peptides are reported together (`"A;B"`); the others compete in later rounds. Since 2023 (GitHub only, not on CRAN), unrelated tied proteins were merged into one group and their remaining peptides were left without a protein.
* `greedy_parsimony()` is 15 to 60 times faster on thousands of peptides, e.g. 0.16 s instead of 5.6 s for 6,000 and 0.44 s instead of 21 s for 12,000 peptides. It compared tied proteins one at a time, which made it quadratic in the number of peptides.
* `greedy_parsimony()` no longer fails with "NA subscripts in x[i,j] not supported" when tied proteins leave only peptides without a protein. Peptides without a protein are left out of the result, as the help page now says.
* `greedy_parsimony()` stops for a matrix without peptide row names or protein column names; it returned an empty list or empty protein names.
* `greedy_to_table()` works with modified peptides: the charge is taken after the last dot, so `"AC[+57.0]K.2"` gives peptide `AC[+57.0]K` and charge 2. Such names used to fail with "length of 'dimnames' [2] not equal to array extent".
* `greedy_to_table()` stops, naming the offending IDs, when a peptide ID does not end in a dot and a charge, e.g. for a `greedy_parsimony()` result on plain peptide sequences. It used to copy the peptide into `precursor_charge`, or split a modified peptide such as `"AC[+57.0]K"` into `AC[+57` and `0]K`.
* `compute_fdr()` and `compute_fdr_with_id()` stop for `NA` scores or decoy flags, a `decoy_hit` that is not logical, and a `score` and `decoy_hit` (or `id`) of different lengths. One `NA` decoy flag used to make all q-values `NA`, and `predict_score_fdr()` then returned `NA` without an error; `"TRUE"`/`"FALSE"` strings gave a warning and `NA` q-values.
* `plot_fdr()` uses a logarithmic score axis only when all scores are positive. Scores of 0 (e.g. PEP or rounded E-values) used to vanish from the histogram and the q-value curves with a warning.
* `annotate_peptides()` and `prepare_matrix()` accept tibbles; a tibble used to crash R or fail with "cannot xtfrm data frames".
* `write_fasta()` keeps a `>` inside a header (`ADP>ATP` used to become `ADPATP`).
* `remove_signal_peptide()` reads UniProt's current signal peptide format `SIGNAL 1..20`, and works on a database with a single entry.
* `remove_signal_peptide()` stops when a name of `db` has no `|` and so no UniProt accession, or `signal$Entry` is `NA`; such entries skipped the check that `signal` lists the entries in the order of `db`.
* `remove_signal_peptide()` gives a shortened sequence its new ID as list name too, so `annotate_aho()` and the written FASTA file report the same ID, and keeps the description in its header (it used to keep the old ID as list name and drop the description).
* `Cdsw$sampling_breaks()` stops with an error when even equal-width windows violate `minwindow` or `maxwindow`; it used to warn and then never return.
* `Cdsw$constant_breaks()`, `readjust_windows()` and `Cdsw$optimize_windows()` work for any `digits`. The mass range was rounded to 2 decimals in `constant_breaks()` and to 1 in `readjust_windows()`, whatever `digits` was, so the outer breaks could lie inside the range of the masses: `constant_breaks(digits = 3)` followed by `$plot()` or `$error()`, and `readjust_windows(digits = 2)`, failed with "some 'x' not counted; maybe 'breaks' do not span range of 'x'" for about 38 % of data sets. The range is now rounded outward at `digits` and padded by `1 / 10^digits`. With the default `digits` the breaks are the same; other values give slightly different breaks (on `masses`, by at most one unit of the last digit for `digits` 0 and 1, and 0.004 m/z for 3 and 4), and the bins of `readjust_windows(digits = 0)` start at whole m/z values.
* `readjust_windows()` and `Cdsw$optimize_windows()` also move a window start down to the first empty bin and a window end up to the last one; such boundaries stayed in place when that bin was the nearest. On `masses`, the first window of `optimize_windows(maxbin = 10)` starts at 350.05 instead of 350.13. When no bin is empty, e.g. with `maxbin = 0`, the windows are returned unchanged; this failed with "argument is of length zero".
* `Cdsw$quantile_breaks()` and `Cdsw$sampling_breaks()` return unnamed breaks, like `constant_breaks()`, and `Cdsw$as_table()` then has row numbers instead of the quantile names `0%`, `4%`, ..., which did not describe mixed breaks.
* `Cdsw$copy()` works; it failed with 'argument "masses" is missing, with no default'.
* `annotate_peptides()` removes whitespace around the peptides, e.g. the `"\r"` of a file with Windows line endings, before it counts, searches and merges them. Such peptides were found but then dropped from the result without a warning; when all peptides had whitespace, the result was an empty data frame.
* `annotate_aho()` and `annotate_peptides()` search proteins with lower-case residues or characters other than letters, such as the stop codon `*` of translated databases: the sequences are searched in upper case, with `X` for each other character, so the offsets stay the same. One such protein used to stop the whole search with "character '*' is not in the trie's alphabet".
* The help of `annotate_peptides()` says what `matched` checks: by default, that the peptide is preceded by K or R or the protein N-terminus (or by the Met at position 1). The end of the peptide is not checked.
* `annotate_aho()` is faster: 0.02 s instead of 0.12 s for the 1,520 bundled peptides against 1,072 proteins.
* `annotate_aho()` and `annotate_peptides()` no longer print dplyr's "Joining with `by = join_by(proteinID)`" message, nor dplyr's `funs()` deprecation warning.
* Error messages name the argument and the offending value, e.g. for an unknown `method` of `predict_score_fdr()` or `weighting` of `prepare_matrix()`.
* `prepare_matrix()` documents only the weightings it implements (`NULL`, `"inverse"` or one numeric weight per row).
* `prepare_matrix()` drops unused factor levels of the ID columns, so it works on subsets of the bundled data sets. It failed with "length of Dimnames[[2]] (250) is not equal to Dim[2] (239)", or gave a peptide row or protein column without entries.
* `write_fasta()`, `make_fasta_summary()` and `reverse_seq()` are faster and give the same output: for 82,000 proteins and their decoys, writing takes 0.9 s instead of 11 s, the summary 0.6 s instead of 10 s and the reversal 1.2 s instead of 3.8 s, and `create_fgcz_fasta_db()` 4.7 s instead of 30 s.
* `read_peptide_fasta()` gives an entry without sequence lines (a header followed directly by the next header, or a header on the last line) the empty string. It used to return the next header pasted to this one, e.g. `">b>a"`, or `"NA>b"` for the last entry, which `write_fasta()` then wrote as an extra header line.
* `make_fasta_summary()` gives an empty sequence the length 0 in the length summary (it counted 1).
* `prepare_matrix(weighting = "inverse")` no longer builds a dense peptide x protein array: for 12,000 peptides and 5,759 proteins it peaks at 0.28 GB instead of 2.5 GB, the same as without weighting, and takes 0.01 s instead of 0.7 s.
* `load_contaminants("special2024")` works: the bundled `fgcz_special_prot_20241112.fasta.gz` is a plain gzipped FASTA file again instead of a gzipped tar archive (#7).
* `create_fgcz_fasta_db()` reads only FASTA files (`*.fasta`, `*.fa` or `*.faa`, optionally gzipped) and the file `annotation.txt`, and stops when `database_dir` holds no FASTA file or when `output_dir` is `database_dir`. It used to read every file with `fasta` in its name: a `*.fa.gz` database was skipped and the output held only the contaminants, and a second run read the database written into `database_dir` as targets, so target sequences were labeled as decoys (`REV_REV_`).
* `create_decoy_db()` keeps only the first entry of each name, as `create_fgcz_fasta_db()` did. It used to return duplicated names, e.g. 1010 when a target file holds the default contaminants, and `annotate_aho()` refused the result.
* `load_contaminants(no_human = TRUE)` also drops entries whose header says `Homo sapiens`; it kept 7 human keratins of the `"fgcz2019"`, `"fgcz2021"` and `"fgcz2022"` sets.
* `inst/script/fgcz_create_fasta.R --contamin` takes every set of `load_contaminants()`, also `universal2024` and `special2024`, besides `special`, `none` and a FASTA file, and stops early naming the sets for an unknown value. The parameter log of a `nodecoy` run says `revLab: NULL` instead of `REV_`.
* The summary of `create_fgcz_fasta_db()` and the parameter log of `fgcz_create_fasta.R` record the argument values, e.g. `contaminants : 505 sequences`, instead of the call, which showed `contaminants = contaminants` when run from the script. The summary no longer contains a link to the FGCZ intranet.
* The dummy entry `CRAPCRAPCRAP` that `create_decoy_db()` puts first has the attributes it gets when the written file is read back: `Annot` starts with `>` and `name` is the first word of the header (it was the whole annotation line). The written files are unchanged.
* The FGCZ shell wrappers (`inst/script/*.sh`) and the notes in `inst/creatingDBs/` are no longer part of the package; `inst/script/fgcz_create_fasta.R` stays.
* `load_contaminants()` cites the source of the `"universal2024"` set; `inst/COPYRIGHTS` lists its origin.
* Help pages use markdown, link related functions in families, and describe the actual return values (e.g. the elements of `compute_fdr()`'s result).
* `?create_fgcz_fasta_db` gives the memory a build needs (about 4 KB per written sequence), and `?greedy_parsimony` says that protein groups are formed on the peptides left in each round. `pepprot` is titled "Precursors of identified peptides", as it holds no proteins.
* New or runnable examples for `create_fgcz_fasta_db()`, `greedy_to_table()`, `remove_signal_peptide()` and `write_fasta()`; removed debugging code from examples and vignettes.
* The vignettes use the new names and ggplot2; "Peptide Annotation and Protein Inference" builds again. "Creating Target Decoy Databases" finds the FASTA files as `create_fgcz_fasta_db()` does, writes to `tempdir()` and ends with the one-call `create_fgcz_fasta_db()`.
* The `Cdsw` help page no longer documents method arguments without a usage section (r-devel NOTE).
* Removed the undocumented debug data set `mm`.
* Fewer dependencies: protViz, purrr, rlang and seqinr are no longer needed, so a default install also skips the five packages only seqinr needed, RcppArmadillo among them; dplyr and docopt moved to Suggests (dplyr is only used by the vignettes and tests, docopt only by `inst/script/fgcz_create_fasta.R`).

# prozor 0.2.15

## Prepare submission to Bioconductor

* address warnings and notes from BiocCheck tool.
