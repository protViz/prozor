# AGENTS.md

Guide for coding agents working in prozor. The workspace guide `../AGENTS.md` applies too; this file wins for this repository.

## What prozor is

- R package: minimal protein set from peptide-spectrum matches (greedy parsimony), FASTA databases with decoys and contaminants, target-decoy FDR, dynamic DIA (SWATH) windows
- On CRAN as 0.3.1 (2021); GitHub `protviz/prozor` master is ahead and heading for 0.4.0
- Used at FGCZ through `inst/script/fgcz_create_fasta.R`, the command-line tool that builds FASTA databases
- Current work plan: `../TODO/prozor/TODO_cran_health.md` (review reports next to it)

## Commands

- `make document`: roxygen2 docs and NAMESPACE
- `make test`: testthat suite (3rd edition) in `tests/testthat/`
- `make check`: full `R CMD check` via devtools (`--as-cran`); must stay at 0 errors, 0 warnings, 0 notes
- `make site`: altdoc site into `docs/` (config `altdoc/quarto_website.yml`); CI deploys it to gh-pages
- To test-render the site against the working tree, install into a scratch library (`R CMD INSTALL --library=<scratch>`) and run with `R_LIBS=<scratch>`; never overwrite the user library

## Conventions

- Edit roxygen comments in `R/`, never `NAMESPACE` or `man/` by hand
- Every user-visible change gets a `NEWS.md` bullet under the heading of the current DESCRIPTION version; tests and tooling do not
- Vignettes live in `vignettes/`; `vignettes/PlottingProteinCoverage.qmd` is a draft under the maintainer's review: do not edit its content
- No `debug()`, `browser()`, `:::` to prozor itself, or commented-out code in examples and vignettes
- Examples write only to `tempdir()` and restore any `par()` they change
- Declare every package used in examples, vignettes and tests in DESCRIPTION (Imports or Suggests)
- Bundled FASTA files in `inst/extdata/` are plain gzip (`gzip -9n`), never tar archives; `tests/testthat/test-fasta-data.R` guards this

## Decisions

- Maintainer address stays `wewolski@gmail.com`
- Next CRAN version is 0.4.0: bump with `usethis::use_version("minor")`, not `make new-version` (which bumps the patch, tags and pushes before CRAN acceptance)
- `greedy()` returns as a deprecated alias of `greedy_parsimony()` in 0.4.0 and becomes defunct in 0.5.0

## Git

- Commit on `master`, one commit per batch of work, then push to `protviz/prozor`
- A push to `master` runs the altdoc workflow, which deploys the site
