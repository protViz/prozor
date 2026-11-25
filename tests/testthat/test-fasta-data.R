test_that("all bundled FASTA files are present", {
  # 5 contaminant files in inst/extdata, 2 in inst/p1000_db1_example
  expect_length(bundled_fasta_files(), 7)
})

test_that("every bundled FASTA file is a plain gzipped FASTA of amino-acid sequences", {
  for (file in bundled_fasta_files()) {
    # a tar archive starts with a tar header instead of ">" and makes read_fasta() warn
    expect_match(readLines(file, n = 1), "^>", label = paste(basename(file), "first line"))
    expect_no_warning(fasta <- read_fasta(file))
    seqs <- unlist(fasta)
    expect_gt(length(fasta), 0)
    expect_false(anyDuplicated(names(fasta)) > 0, label = paste(basename(file), "has duplicated names"))
    expect_true(all(nchar(seqs) > 0), label = paste(basename(file), "has only non-empty sequences"))
    expect_true(all(grepl("^[A-Za-z*]+$", seqs)), label = paste(basename(file), "has only amino-acid letters"))
  }
})
