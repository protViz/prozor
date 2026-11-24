test_that("all bundled FASTA files are present", {
  # 5 contaminant files in inst/extdata, 2 in inst/p1000_db1_example
  expect_length(bundled_fasta_files(), 7)
})

test_that("every bundled FASTA file is a plain gzipped FASTA of amino-acid sequences", {
  for (file in bundled_fasta_files()) {
    # a tar archive starts with a tar header instead of ">" and makes read_peptide_fasta() warn
    expect_match(readLines(file, n = 1), "^>", label = paste(basename(file), "first line"))
    expect_no_warning(fasta <- read_peptide_fasta(file))
    seqs <- unlist(fasta)
    expect_gt(length(fasta), 0)
    expect_false(anyDuplicated(names(fasta)) > 0, label = paste(basename(file), "has duplicated names"))
    expect_true(all(nchar(seqs) > 0), label = paste(basename(file), "has only non-empty sequences"))
    expect_true(all(grepl("^[A-Za-z*]+$", seqs)), label = paste(basename(file), "has only amino-acid letters"))
  }
})

test_that("load_contaminants() returns each set and no_human removes human entries", {
  for (set in c("fgcz2019", "fgcz2021", "fgcz2022", "universal2024", "special2024")) {
    all_entries <- load_contaminants(set)
    no_human <- load_contaminants(set, no_human = TRUE)
    annot <- vapply(no_human, seqinr::getAnnot, character(1))
    expect_gt(length(all_entries), 0)
    expect_lt(length(no_human), length(all_entries))
    expect_false(any(grepl("HUMAN", annot)), label = paste(set, "(no_human = TRUE) has HUMAN entries"))
  }
})

test_that("load_contaminants() rejects unknown sets", {
  expect_error(load_contaminants("fgcz2023"), "`set` must be one of .* not \"fgcz2023\"")
})
