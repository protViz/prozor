# Contributing to prozor

Bug reports, questions and pull requests are welcome.

## Reporting a problem

- Open an issue at https://github.com/prolfqua/prozor/issues
- Include a minimal example that reproduces the problem, ideally with the bundled data (`fdr_sample`, `masses`, the files in `inst/p1000_db1_example/`)
- Include the output of `sessionInfo()`

## Proposing a change

- For anything larger than a typo, open an issue first so the change can be discussed
- Fork the repository, branch from `master`, and open a pull request against `master`
- One topic per pull request

## Development

The `Makefile` wraps the usual steps:

- `make document`: roxygen2 documentation and `NAMESPACE`; edit the roxygen comments in `R/`, never `man/` or `NAMESPACE`
- `make test`: testthat suite in `tests/testthat/`
- `make format`: format with air (`air.toml`: 120 columns, 2 spaces)
- `make lint`: lintr (`.lintr`)
- `make check`: `R CMD check --as-cran`; it must stay at 0 errors, 0 warnings, 0 notes

Conventions:

- snake_case for functions, arguments, result columns and list elements
- New code comes with tests; bug fixes come with a test that fails without the fix
- Every user-visible change gets a bullet in `NEWS.md` under the current version
- Plots use ggplot2 and return the ggplot object
- Examples and tests write only to temporary directories

CI runs R CMD check on Linux, macOS and Windows, lint, and the website build for every pull request.
