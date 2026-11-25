# prozor 0.5.0

prozor 0.4.0 was not released to CRAN; this section lists all changes since CRAN 0.3.1.

## Breaking changes

* All exported functions and their arguments use snake_case. The old names are gone, without aliases:

| Old | New |
|---|---|
| `annotateAHO()`, `annotatePeptides()` | `annotate_aho()`, `annotate_peptides()` |
| `computeFDR()`, `computeFDRwithID()`, `predictScoreFDR()`, `plotFDR()` | `TargetDecoyFDR$new()`, `$score_at()`, `$plot()` (see below) |
| `createDecoyDB()`, `reverseSeq()` | `create_decoy_db()`, `reverse_seq()` |
| `readPeptideFasta()`, `writeFasta()` | `read_fasta()`, `write_fasta()` |
| `loadContaminantsFasta2019()`, `loadContaminantsFasta2021()`, `loadContaminantsFGCZ2022()`, `load_universal_contaminants_2024()`, `load_special_prot_2024()` | `load_contaminants("fgcz2019")`, `"fgcz2021"`, `"fgcz2022"`, `"universal2024"`, `"special2024"` |
| `greedyRes2Matrix(res)` | `greedy_to_table(groups)` |
| `removeSignalPeptide()`, `makeID()`, `makeIDUnip()` | `remove_signal_peptide()`, `make_chain_id()`, `make_chain_id_uniprot()` |
| `readjustWindows()` | `readjust_windows()` |
| `Cdsw` methods `asTable()`, `showCycle()`, `optimizeWindows()` | `as_table()`, `plot_cycle()`, `optimize_windows()` |
| `Cdsw` methods `constant_breaks()`, `quantile_breaks()`, `sampling_breaks()` | `Cdsw$new(method = "constant")`, `"quantile"`, `"mixture"` |
| `Cdsw$error()` elements `score1`, `score2`, `counts`, `optimumN` | `manhattan`, `euclidean`, `normalized_counts`, `equal_counts` |
| arguments `useContaminants`, `revLab`, `noHuman`, `databasedirectory`, `outputdir`, `resDB`, `qValue`, `idfun`, `ID` | `contaminants`, `rev_lab`, `no_human`, `database_dir`, `output_dir`, `fasta`, `q_value`, `id_fun`, `id` |
| arguments `pepinfo` of `annotatePeptides()`, `pepseq` of `annotateAHO()` | `peptides` |
| arguments `dbs`, `annot` of `createDecoyDB()`, `db` of `removeSignalPeptide()` | `files`, `header`, `fasta` |
| arguments `wind`, `ms1data`, `maxbin` of `readjustWindows()` and `Cdsw$optimizeWindows()` | `windows`, `masses`, `min_count` |
| column `counts` of `Cdsw$asTable()`, `Cdsw$optimizeWindows()` and `readjustWindows()` | `n_precursors` |
| `Cdsw` arguments and field `nbins`, `maxwindow`, `minwindow` | `n_bins`, `max_window`, `min_window` |
| argument `sp` of `makeID()`, `makeIDUnip()` | `start` |
| data set `fdrSample`, its column `score2` | `fdr_sample`, column `neg_log10_score` |
| column `m` of the data sets `pepprot` and `protpepmetashort` | `precursor_mz` |

* Renamed since 0.4.0 (GitHub only), without aliases:
  * `read_peptide_fasta()` → `read_fasta()`
  * `annotate_peptides(pepinfo)`, `annotate_aho(pepseq)` → `annotate_peptides(peptides)`, `annotate_aho(peptides)`
  * `create_decoy_db(dbs, annot)` → `create_decoy_db(files, header)`
  * `make_fasta_summary(db)`, `remove_signal_peptide(db)`, element `db` of `create_fgcz_fasta_db()` → `fasta`
  * `make_fasta_summary()` element and summary label `nr_sequences`, `nr sequences` → `n_sequences`
  * `readjust_windows(wind, ms1data)`, `plot_windows(windows, ms1data)` → `readjust_windows(windows, masses)`, `plot_windows(windows, masses)`
  * `maxbin` of `readjust_windows()` and `Cdsw$optimize_windows()` → `min_count`
  * column `counts` of `Cdsw$as_table()`, `Cdsw$optimize_windows()` and `readjust_windows()` → `n_precursors`
  * `Cdsw$new(nbins, maxwindow, minwindow)`, field `nbins` → `n_bins`, `max_window`, `min_window`, field `n_bins`
  * `Cdsw$new(method = "sampling")` → `method = "mixture"`
  * `Cdsw$show_cycle()` → `Cdsw$plot_cycle()`
  * `Cdsw$error()` elements `score1`, `score2`, `counts`, `optimum_n` → `manhattan`, `euclidean`, `normalized_counts`, `equal_counts`
  * `greedy_to_table(res)`, column `protein` → `greedy_to_table(groups)`, column `protein_group`
  * `TargetDecoyFDR`: `simple_fdr` → `fdr_simple` (method value, columns `fdr_simple` and `q_value_fdr_simple`)
  * `make_chain_id(sp)`, `make_chain_id_uniprot(sp)` → `start`
  * column `m` of `pepprot` and `protpepmetashort` → `precursor_mz`
  * column `score2` of `fdr_sample` → `neg_log10_score`

* Result columns and list elements are snake_case too:
  * `annotate_aho()`, `annotate_peptides()`: `protein_id`, `offset`, `protein_sequence`, `peptide_length`, and `peptide_seq` for vector input (was `proteinID`, `Offset`, `proteinSequence`, `lengthPeptide`, `peptideSeq`)
  * `greedy_to_table()`: `peptide`, `precursor_charge`, `protein_group` (was `Peptide`, `precursorCharge`, `Protein`); `precursor_charge` is an integer instead of a character
  * `make_fasta_summary()`: `n_sequences`, `length_summary`, `aa_freq`; `create_fgcz_fasta_db()`: `fasta`, `file_path`, `summary`, `call`, `db_name`
* The data sets `fdr_sample`, `protpepmetashort` and `pepprot` have snake_case columns (`protein_id`, `peptide_seq`, `peptide_mod_seq`, `precursor_charge`, `precursor_mz`, `protein_length`, `peptide_length`, `neg_log10_score`); `bcd` keeps COMET's column names.
* Default column names follow the results: `annotate_peptides(peptide = "peptide_seq")`, and `greedy_parsimony(peptide = "peptide_seq")` reads the proteins from `protein_id` (was `annotatePeptides(peptide = "peptideSeq")` and `prepareMatrix(proteinID = "proteinID", peptideID = "strippedSequence")`).
* The R6 class `TargetDecoyFDR` replaces `computeFDR()`, `computeFDRwithID()`, `predictScoreFDR()` and `plotFDR()`, without aliases (see `?TargetDecoyFDR`):
  * `computeFDR(score, decoy_hit)` → `fdr <- TargetDecoyFDR$new(score, decoy_hit)`; the data frame `fdr$table` holds the hits in score order with both FDR estimates and their q-values (was a list of vectors).
  * `FPR`, `SimpleFDR`, `qValue_FPR`, `qValue_SimpleFDR` → `fdr_concatenated`, `fdr_simple`, `q_value_fdr_concatenated`, `q_value_fdr_simple`; 2 * decoys / all hits estimates the FDR of a concatenated search, not a false positive rate.
  * `predictScoreFDR(fdr, qValue = 1)` → `fdr$score_at(q_value = 0.01)`: `q_value` is a fraction, not a percentage. `method = "SimpleFDR"` (default) and `"FPR"` → `"fdr_simple"` and `"fdr_concatenated"`.
  * `plotFDR(fdr)` → `fdr$plot()`.
  * `computeFDRwithID(score, ID)` → `TargetDecoyFDR$new(score, decoy_hit = grepl("^REV_", ID), id = ID)`; only IDs starting with `REV_` count as decoys (`decoy = "REV_"` matched anywhere in the ID).
* `TargetDecoyFDR` computes the FDR per score, counting all tied hits: they share one FDR and q-value, whatever the input order, and `score_at()` no longer returns a cut-off whose FDR exceeds `q_value` (1.07 % at 0.01 on `bcd` with rounded scores).
* `TargetDecoyFDR$score_at()` stops when no score reaches the requested FDR; `predictScoreFDR()` returned `Inf` with a warning.
* `greedy_parsimony(annotation, peptide = "peptide_seq")` replaces `greedy()` and `prepareMatrix()` (`prepare_matrix()` on GitHub), without aliases. It takes the peptide-protein matches as a data frame, e.g. from `annotate_peptides()`, instead of a peptide-by-protein matrix (see `?greedy_parsimony`).
  * For precursors, add a column such as `paste(peptide_seq, precursor_charge, sep = ".")` and pass its name as `peptide`.
  * The weightings of `prepareMatrix()` are gone: they no longer change the result, since proteins are ranked by their number of peptides.
* `greedy_parsimony()` returns a named character vector instead of a named list of single strings; drop `unlist()` on its result.
* `greedy_parsimony()` ranks proteins by the number of peptides they explain, as documented. `greedy()` ranked them by the sum of their matrix entries, so with `prepareMatrix(weighting = "inverse")` a protein with fewer but less shared peptides could win the shared ones.
* `greedy_parsimony()` no longer prints its run time.
* `annotate_aho()` and `annotate_peptides()` return a data frame with 0 rows and the result columns when no peptide occurs in `fasta`, instead of `NULL`; `annotate_peptides()` still warns. `residue_coverage()` and `sequence_coverage()` return 0 rows for such a result.
* `annotate_aho()` and `annotate_peptides()` stop for peptides with characters other than upper-case letters (AhoCorasickTrie 0.1.3 crashed R on them) and for duplicated protein names in `fasta` (they multiplied the result rows).
* `annotate_aho()` stops for peptides with whitespace around them and shows the peptide escaped, e.g. `"PEPTIDEK\r"`; it used to remove the whitespace. `annotate_peptides()` still removes it.
* `annotate_peptides()` stops, naming the columns, when `peptides` has a column named like a result column (`protein_id`, `offset`, `protein_sequence`, `matched`, `peptide_length`).
  * Such columns, e.g. of search engine output, gave `protein_id.x` and `protein_id.y` instead of `protein_id`, or were overwritten.
* `annotate_peptides()` sets `matched` for each occurrence of a peptide in a protein; a second, non-tryptic occurrence was flagged as tryptic.
* `make_fasta_summary(fasta)` always returns the list; the arguments `as_string` and `old` are gone. The summary of `create_fgcz_fasta_db(make_summary = TRUE)` lists the sequence lengths and amino-acid counts one per line.
* `create_fgcz_fasta_db(database_dir, output_dir, ...)` requires `output_dir`, now its second argument; it used to write next to `database_dir` by default. `inst/script/fgcz_create_fasta.R` still does when `-o` is not given.
* `create_decoy_db()` and `create_fgcz_fasta_db()` add the contaminant sets `"special2024"` and `"universal2024"` by default, like `inst/script/fgcz_create_fasta.R` (was `"fgcz2021"` in CRAN 0.3.1, `"fgcz2022"` in 0.4.0).
* `write_fasta(fasta, file)` takes one list of sequences, data first, and returns `file` invisibly (was `writeFasta(file, ...)`); combine several lists with `c()`.
* `Cdsw` is an R6 class whose window breaks are computed once, by `Cdsw$new()` (see `?Cdsw`). It was a reference class whose `*_breaks()` methods overwrote the breaks, so the results of the other methods depended on which of them ran last.
  * `constant_breaks()`, `quantile_breaks()`, `sampling_breaks()` → `Cdsw$new(masses, method = "constant")`, `"quantile"`, `"mixture"` (the default); create one object per method. `Cdsw(masses)` started with `"constant"`.
  * `copy()` → `clone()`. The fields `masses`, `n_bins`, `method` and `breaks` are read-only.
* `Cdsw$new(method = "mixture")` is deterministic and fast (see `?Cdsw`). It added random masses until the windows fitted, so the breaks changed from call to call, and it could take minutes or never return.
  * The breaks differ from those of earlier versions; the counts per window are as even as before.
* `readjust_windows()`, `plot_windows()` and `Cdsw$optimize_windows()` take the bin width of the masses in m/z, `bin_width = 0.1`, instead of `digits = 1`; `digits` of `Cdsw$new()` sets only the decimal places of the breaks.
  * Replace `digits = d` with `bin_width = 10^-d`.
* `readjust_windows()` returns a data frame instead of a list, like `Cdsw$as_table()` and `Cdsw$optimize_windows()`.
* `Cdsw$optimize_windows()` widens the windows by `overlap = 1` m/z by default, as `as_table()` and `plot_cycle()` do (was 0).
* Precursors per window are counted as the masses `m` with `from < m <= to`, the intervals of `hist()`. `Cdsw$as_table()` left out masses on a window boundary, and `readjust_windows()` and `Cdsw$optimize_windows()` counted them in both windows. Counts change only for masses exactly on a boundary.
* Plots use ggplot2: `TargetDecoyFDR$plot()`, `Cdsw$plot()` and `Cdsw$plot_cycle()` return ggplot objects and draw nothing themselves.
  * The `plot` argument of `readjust_windows()`, `Cdsw$optimize_windows()` and the old `sampling_breaks()` is gone; use `plot_windows()` and `Cdsw$plot_breaks()`.
* `read_fasta()` and `write_fasta()` use base R instead of seqinr. The entries are still `SeqFastaAA` objects with the attributes `name` and `Annot`, so seqinr's `getAnnot()`, `getName()` and `getLength()` work on them.
  * `library(prozor)` no longer loads seqinr; load it to get seqinr's `summary()` of an entry.
* prozor requires R 4.4.0 or later (was R 3.1.0).

## New features

* `load_contaminants(set, no_human)` loads any of the bundled contaminant sets.
* `load_contaminants()` loads the new set `"fgcz2024"` by default: `"special2024"` and `"universal2024"` together, the default contaminants of `create_decoy_db()`, `create_fgcz_fasta_db()` and `fgcz_create_fasta.R --contamin=special`.
* `plot_windows()` shows the precursor masses with the DIA window boundaries.
* `Cdsw$plot_breaks()` compares the breaks of the object with the constant and quantile breaks of the same masses.
* `residue_coverage()`, `sequence_coverage()` and `plot_coverage()` compute and plot the peptide coverage of protein sequences from the result of `annotate_peptides()`; the new vignette "Protein Sequence Coverage" shows them.
* `citation("prozor")` cites the package with its CRAN DOI.

## Minor improvements and fixes

* `greedy_parsimony()` assigns every peptide that matches a protein again: only proteins that explain exactly the same peptides are reported together (`"A;B"`). Since 2023 (GitHub only), unrelated tied proteins were merged into one group and their other peptides were left without a protein.
* `greedy_parsimony()` is 250 to 400 times faster on thousands of peptides, e.g. 0.05 s instead of 21 s for 12,000 peptides.
* `greedy_parsimony()` no longer scans all proteins in every round: 2.2 s instead of 36 s for 400,000 peptides that match many proteins.
* `greedy_parsimony()` breaks ties between character IDs in C-locale order (upper case first) in every locale; it followed the session's collation, so a tie between proteins `"a"` and `"B"` depended on the locale.
* `greedy_parsimony()` ignores unused factor levels of the ID columns, so it works on subsets of the bundled data sets; `prepareMatrix()` failed on them or gave rows or columns without entries.
* `annotate_peptides()` and `greedy_parsimony()` accept tibbles; `annotatePeptides()` and `prepareMatrix()` crashed R or failed on them.
* `greedy_to_table()` takes the charge after the last dot, so modified peptides such as `"AC[+57.0]K.2"` work; they failed with "length of 'dimnames' [2] not equal to array extent".
* `greedy_to_table()` stops, naming the IDs, when a peptide ID does not end in a dot and a charge, e.g. for a `greedy_parsimony()` result on plain peptide sequences; it used to return a wrong `precursor_charge`.
* `TargetDecoyFDR$new()` checks its input: it stops for a `score` that is not numeric, a `decoy_hit` that is not logical, `NA` in either, a `decoy_hit` or `id` of another length than `score`, and a `larger_better` other than `TRUE` or `FALSE`.
  * Character scores were sorted as text, and one `NA` decoy flag made all q-values `NA`.
* `TargetDecoyFDR$new()` stops for zero hits; `$plot()` failed on them with "arguments imply differing number of rows: 0, 1".
* `TargetDecoyFDR` caps both FDR estimates and their q-values at 1; when the best hit was a decoy, `fdr_simple` was `Inf` and q-values exceeded 1.
* `TargetDecoyFDR$plot()` uses a logarithmic score axis only when all scores are positive; scores of 0 (e.g. PEP or rounded E-values) vanished from the plot with a warning.
* `annotate_peptides()` removes whitespace around the peptides, e.g. the `"\r"` of a file with Windows line endings, before it searches them; such peptides were dropped from the result without a warning.
* `annotate_aho()` and `annotate_peptides()` search proteins with lower-case residues or characters other than letters, such as the stop codon `*` of translated databases, with unchanged offsets; one such protein stopped the search with "character '*' is not in the trie's alphabet".
* `annotate_aho()` and `annotate_peptides()` give the offsets in the sequences of `fasta` also for a sequence that starts with blanks; they removed such blanks and gave offsets shifted by their number. Blanks are searched as `X`, like `*`.
* `annotate_peptides()` stops when `peptides` lacks the column named by `peptide`, and `annotate_aho()` when `fasta` has no names; they failed with "replacement has 0 rows, data has 1" and "arguments imply differing number of rows: 0, 1".
* `annotate_aho()` is about 6 times faster on the bundled example.
* `annotate_aho()` and `annotate_peptides()` no longer print dplyr's "Joining with `by = join_by(proteinID)`" message, nor dplyr's `funs()` deprecation warning.
* `readjust_windows()` and `plot_windows()` stop for `masses` with `NA` or non-numeric values, as `Cdsw$new()` does, and for `windows` without the column `from` or `to`.
* `plot_coverage()` stops unless `width` is a positive whole number.
* `residue_coverage()` and `sequence_coverage()` take 0.5 s and 0.2 s instead of 15 s for 190,000 peptides, and 0.6 and 0.3 GB of memory instead of 6 GB.
* `Cdsw$new(method = "mixture")` stops when even equal-width windows violate `min_window` or `max_window`; `sampling_breaks()` warned and then never returned.
* `Cdsw$new(method = "mixture")` keeps the rounded and padded breaks within `min_window` and `max_window`; the last window was up to `2 / 10^digits` wider than `max_window`.
* `Cdsw$new()` stops for an infinite `n_bins`, `max_window`, `min_window` or `digits`, `readjust_windows()` and `Cdsw$optimize_windows()` for a `min_count` that is not a non-negative number, and `Cdsw$as_table()`, `$optimize_windows()` and `$plot_cycle()` for such an `overlap`; `digits = Inf` gave `NaN` breaks, and `min_count = "a"` compared the counts as text.
* Constant breaks enclose the masses for any `digits`, and `readjust_windows()` and `Cdsw$optimize_windows()` work for any bin width. `constant_breaks(digits = 3)` followed by `$plot()` or `$error()`, and `readjust_windows(digits = 2)`, failed with "some 'x' not counted" for about 38 % of data sets.
  * With the default `digits` the breaks are unchanged; other values give slightly different breaks.
* `readjust_windows()` and `Cdsw$optimize_windows()` also move a window start down to the first empty bin and a window end up to the last one; such boundaries stayed in place when that bin was the nearest.
* `readjust_windows()` and `Cdsw$optimize_windows()` return the windows unchanged when no bin is empty, e.g. with `min_count = 0`; they failed with "argument is of length zero".
* `readjust_windows()` leaves a boundary that is already on the mid of an empty bin in place, so it returns its own result unchanged; it moved such a boundary to the next empty bin, and a second pass widened the windows by up to 0.7 m/z.
* The breaks of `Cdsw` are unnamed for all methods, so `Cdsw$as_table()` has row numbers instead of the quantile names `0%`, `4%`, ...
* `write_fasta()` keeps a `>` inside a header (`ADP>ATP` used to become `ADPATP`).
* `write_fasta()` stops, naming the first such entry, unless every entry has one `Annot` string, before it creates the file; it failed with "values must be length 1" and left an empty or partial file.
* `write_fasta()`, `make_fasta_summary()` and `reverse_seq()` are faster and give the same output; `create_fgcz_fasta_db()` takes 4.7 s instead of 30 s for 82,000 proteins.
* `read_fasta()` gives an entry without sequence lines the empty string; it returned the next header pasted to this one, e.g. `">b>a"`, which `write_fasta()` then wrote as an extra header line.
* `make_fasta_summary()` counts an empty sequence as length 0 (was 1).
* `read_fasta()` ends an entry's name at the first whitespace, so a tab-separated header `>sp|P1|A_HUMAN<TAB>...` gives `sp|P1|A_HUMAN`; only a space ended it. `reverse_seq()` takes the decoy names from the decoy headers in the same way.
* `read_fasta()` and `create_fgcz_fasta_db()` read a FASTA file or `annotation.txt` without a final newline without the warning "incomplete final line found".
* `remove_signal_peptide()` reads UniProt's current signal peptide format `SIGNAL 1..20`, and works on a database with a single entry.
* `remove_signal_peptide()` stops when a name of `fasta` has no `|` and so no UniProt accession, or `signal$Entry` is `NA`; such entries skipped the check that `signal` lists the entries in the order of `fasta`.
* `remove_signal_peptide()` gives a shortened sequence its new ID as list name too, so `annotate_aho()` and the written FASTA file report the same ID, and keeps the description in its header.
* `remove_signal_peptide()` stops, naming the missing columns, when `signal` lacks `Entry` or `Signal.peptide`; with UniProt's `Signal peptide` (as `read.delim(check.names = FALSE)` names it) it returned `fasta` unchanged, and without `Entry` it skipped the order check.
* `remove_signal_peptide()` takes linear time: 0.5 s instead of 5 s for 82,000 entries.
* `load_contaminants("special2024")` works: the bundled `fgcz_special_prot_20241112.fasta.gz` is a plain gzipped FASTA file again instead of a gzipped tar archive (#7).
* `load_contaminants(no_human = TRUE)` also drops entries whose header says `Homo sapiens`; it kept 7 human keratins of the `"fgcz2019"`, `"fgcz2021"` and `"fgcz2022"` sets.
* `create_fgcz_fasta_db()` reads only FASTA files (`*.fasta`, `*.fa` or `*.faa`, optionally gzipped) and `annotation.txt`, and stops when `database_dir` holds no FASTA file or is `output_dir`.
  * It read every file with `fasta` in its name, so a `*.fa.gz` database was skipped, and a second run read the written database as targets (`REV_REV_` decoys).
* `create_fgcz_fasta_db()` stops before it reads any file when `output_dir` does not exist; it failed with "cannot open file" after building the whole database.
* `create_decoy_db()` keeps only the first entry of each name, as `create_fgcz_fasta_db()` did; duplicated names made `annotate_aho()` refuse the result.
* The summary of `create_fgcz_fasta_db()` and the parameter log of `fgcz_create_fasta.R` record the argument values, e.g. `contaminants : 505 sequences`, instead of the call, and no longer link to the FGCZ intranet.
* The dummy entry `CRAPCRAPCRAP` that `create_decoy_db()` puts first has the attributes it gets when the written file is read back: `Annot` starts with `>` and `name` is the first word of the header. The written files are unchanged.
* `inst/script/fgcz_create_fasta.R --contamin` takes every set of `load_contaminants()` besides `special`, `none` and a FASTA file, and stops early, naming the sets, for an unknown value. The parameter log of a `nodecoy` run says `revLab: NULL` instead of `REV_`.
* `inst/script/fgcz_create_fasta.R` writes its parameter log next to the FASTA file in the output directory, named like it, e.g. `p1000_db1_example_d_20261125.txt`; it wrote `<db_name>.txt` into the working directory.
* The FGCZ shell wrappers (`inst/script/*.sh`) and the notes in `inst/creatingDBs/` are no longer part of the package; `inst/script/fgcz_create_fasta.R` stays.
* The example database `inst/p1000_db1_example` is described as human in its `annotation.txt`; it said mouse.
* Errors raised by exported functions show their call; errors from internal checks and R6 methods show none, so `residue_coverage()` no longer reports `.check_annotation()`.
* Error messages name the argument and the offending value, and put column names in backticks.
* `library(prozor)` no longer prints the package version.
* `library(prozor)` no longer loads ggplot2, only the first plot does: 0.002 s instead of 0.33 s.
* Removed the undocumented debug data set `mm`.

## Documentation

* Help pages use markdown, link related functions in families, and describe the actual return values, e.g. the columns of `TargetDecoyFDR$table`. The data sets are listed with the functions of their family, e.g. `fdr_sample` and `bcd` with `TargetDecoyFDR`.
* `?prozor` gives an overview of the functions and vignettes by topic.
* `?TargetDecoyFDR` cites Elias and Gygi (2007) and Käll et al. (2008), with DOIs. `?load_contaminants` cites the source of the `"universal2024"` set, and `inst/COPYRIGHTS` lists its origin.
* The `Cdsw` help page documents each method with its usage and arguments; it documented method arguments without a usage section (r-devel NOTE).
* `?annotate_peptides` says what `matched` checks: by default, that the peptide is preceded by K or R or the protein N-terminus (or by the Met at position 1). The end of the peptide is not checked.
* `?create_fgcz_fasta_db` gives the memory a build needs, and `?greedy_parsimony` says that protein groups are formed on the peptides left in each round. `pepprot` is titled "Precursors of identified peptides", as it holds no proteins.
* New or runnable examples for `create_fgcz_fasta_db()`, `greedy_to_table()`, `remove_signal_peptide()` and `write_fasta()`; removed debugging code from examples and vignettes.
* The vignettes use the new names, snake_case code and ggplot2; "Peptide Annotation and Protein Inference" builds again. "Creating Target Decoy Databases" finds the FASTA files as `create_fgcz_fasta_db()` does, writes to `tempdir()` and ends with the one-call `create_fgcz_fasta_db()`.
* The README has an R example and a citation section.

## Dependencies

Compared with CRAN 0.3.1:

* Removed from Imports: Matrix, methods, purrr, readr, rlang, seqinr and stringr.
* Moved to Suggests: docopt, used only by `inst/script/fgcz_create_fasta.R`, and dplyr, used only by the vignettes and tests.
* Added to Imports: ggplot2 and R6; the base packages graphics, stats and utils are now declared too.
* Added to Suggests: altdoc, testthat and withr.
* protViz, imported by the GitHub versions from 2023 on, is no longer needed; no CRAN release depended on it.

# prozor 0.3.1

CRAN release after 0.2.11 (0.2.15 was prepared for Bioconductor and not released to CRAN). Changes since 0.2.15:

* Removed `loadContaminantsFasta()` and `loadContaminantsNoHumanFasta()`, the FGCZ contaminants of 2015. `createDecoyDB()` adds `loadContaminantsFasta2021()` by default, as `create_fgcz_fasta_db()` did already.
* The example FASTA files moved from `inst/extdata` to `inst/p1000_db1_example`; several bundled files were removed and the data files made smaller to reduce the package size.
* Removed the internal helpers that copied FASTA files to an FGCZ server and registered them in B-Fabric; ssh and reticulate are no longer suggested.
* prozor requires R 3.1.0 or later (was R 4.1.0).
* The maintainer address is wewolski@gmail.com.

# prozor 0.2.15

## Prepare submission to Bioconductor

* address warnings and notes from BiocCheck tool.
