# AGENTS.md

Guide for coding agents working in prozor. The workspace guide `../AGENTS.md` applies too; this file wins for this repository.

## What prozor is

- R package: minimal protein set from peptide-spectrum matches (greedy parsimony), FASTA databases with decoys and contaminants, target-decoy FDR, dynamic DIA (SWATH) windows
- On CRAN as 0.3.1 (2021); GitHub `prolfqua/prozor` master (moved from protViz on 25 Nov 2026; site https://prolfqua.github.io/prozor/) is ahead and heading for 0.5.0
- Used at FGCZ through `inst/script/fgcz_create_fasta.R`, the command-line tool that builds FASTA databases
- Current work plan: `../TODO/prozor/TODO_cran_health.md` (review reports next to it)

## Commands

- `make document`: roxygen2 docs and NAMESPACE
- `make test`: testthat suite (3rd edition) in `tests/testthat/`
- `make check`: full `R CMD check` via devtools (`--as-cran`); must stay at 0 errors, 0 warnings, 0 notes
- `make site`: altdoc site into `docs/` (config `altdoc/quarto_website.yml`); CI deploys it to gh-pages
- To test-render the site against the working tree, install into a scratch library (`R CMD INSTALL --library=<scratch>`) and run with `R_LIBS=<scratch>`; never overwrite the user library

## Conventions

- Edit roxygen comments in `R/` (markdown is on), never `NAMESPACE` or `man/` by hand
- Format with `make format` (air, `air.toml`: 120 columns, 2 spaces); `.lintr` configures `make lint`. The air reformat is listed in `.git-blame-ignore-revs`
- Every user-visible change gets a `NEWS.md` bullet under the heading of the current DESCRIPTION version; tests and tooling do not
- Vignettes live in `vignettes/`
- No `debug()`, `browser()`, `:::` to prozor itself, or commented-out code in examples and vignettes
- Examples write only to `tempdir()` and restore any `par()` they change
- Declare every package used in examples, vignettes and tests in DESCRIPTION (Imports or Suggests)
- Bundled FASTA files in `inst/extdata/` are plain gzip (`gzip -9n`), never tar archives; `tests/testthat/test-fasta-data.R` guards this

## Decisions

- Maintainer address stays `wewolski@gmail.com`
- DESCRIPTION is at 0.5.0, the next CRAN version. Don't run `make new-version` before CRAN accepts it (it bumps the patch, tags and pushes)
- API is snake_case: functions, arguments, `Cdsw` methods, result columns and list elements, default column names, method values and the bundled data sets (except `bcd`, which keeps COMET's column names); renamed in 0.4.0 without aliases, `greedy()` included
- `remove_signal_peptide()` reads the UniProt columns `Entry` and `Signal.peptide` of `signal` as `read.delim()` names them: an exception to snake_case, like `bcd`
- Plots use ggplot2 and return ggplot objects; no base graphics, no `plot = TRUE` arguments
- `Depends: R (>= 4.4.0)`
- The CLI options of `inst/script/fgcz_create_fasta.R` (`--revLab`, `--contamin`, ...) stay as they are: the FGCZ wrappers (slurmworker `config/prozor/script/`, moved out of prozor in 0.4.0) and FGCZ jobs call them
- `annotate_aho()` rejects peptides outside `^[A-Z]+$`: modifications, whitespace and lower case (AhoCorasickTrie 0.1.3 crashes R on them, 0.1.4 gives a less specific error), and `""` and `NA` (wrong hits in both versions)
- Dev tools other than altdoc stay out of Suggests, unlike the workspace rule: the Makefile calls devtools, covr and lintr from the local library, and CI installs its tools through `extra-packages`

## Git

- Commit on `master`, one commit per batch of work, then push to `prolfqua/prozor`
- A push to `master` runs the altdoc workflow, which deploys the site
